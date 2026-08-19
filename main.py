import pandas as pd

from sparsification import (
    build_tensor,
    collapse_tensor,
    compute_laplacian,
    compute_effective_resistance,
    compute_edge_sampling_probabilities,
    sample_edges,
    visualize_graph
)

df = pd.read_csv("cfs.csv")

T, taxa = build_tensor(df)

W = collapse_tensor(T, method="sum")

L = compute_laplacian(W)

R = compute_effective_resistance(L)

edge_list = compute_edge_sampling_probabilities(
    W,
    R,
    taxa
)

W_sparse = sample_edges(
    edge_list,
    num_samples=20,
    n=len(taxa),
    taxa=taxa,
    seed=42
)

graph = visualize_graph(W_sparse, taxa,20)