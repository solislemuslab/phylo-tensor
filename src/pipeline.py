import pandas as pd
import numpy as np
from pathlib import Path
import sys
from .sparsification import *

def run_pipeline(network_dir, num_seeds, method, num_samples):

    # 1. Load Data and Build Tensor
    network_dir = Path(network_dir)
    df = pd.read_csv(network_dir / "CFs.csv")
    T, taxa = build_tensor(df)

    # 2. Compute sampling probabilities
    W = collapse_tensor(T, method)
    L = compute_laplacian(W)
    R = compute_effective_resistance(L)
    edge_list = compute_edge_sampling_probabilities(W, R, taxa)

    # 3. Run sparsification experiment
    weights_by_edge = {
        (i, j): []
        for i, j, w, r, s, p in edge_list
    }

    for seed in range(num_seeds):
        W_sparse = sample_edges(
            edge_list,
            num_samples=num_samples,
            n=len(taxa),
            taxa=taxa,
            seed=seed
        )

        for i, j in weights_by_edge:
            # Includes zero when the edge wasn't retained
            weights_by_edge[(i, j)].append(W_sparse[i, j])

    # 4. Summarize results
    rows = []

    W_mean = np.zeros_like(W, dtype=float)

    for (i, j), weights in weights_by_edge.items():
        weights = np.array(weights)

        count = np.count_nonzero(weights > 0)
        mean_weight = np.mean(weights)

        W_mean[i, j] = mean_weight
        W_mean[j, i] = mean_weight

        rows.append({
            "edge": f"{taxa[i]}-{taxa[j]}",
            "runs_retained": count,
            "retention_percent": 100 * count / num_seeds,
            "original_weight": W[i, j],
            "total_weight": np.sum(weights),
            "mean_weight": mean_weight
        })

    weight_df = (
        pd.DataFrame(rows)
        .sort_values("mean_weight", ascending=False)
        .reset_index(drop=True)
    )

    # 5. Display results
    print(weight_df.round(3).to_string(index=False))
    return weight_df