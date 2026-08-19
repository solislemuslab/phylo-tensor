import pandas as pd
import numpy as np
from scipy.linalg import pinv
import networkx as nx
import matplotlib.pyplot as plt

def build_tensor(df):
    taxa = set()

    # iterates through the 4 columns and adds them into taxa set
    for col in ["t1", "t2", "t3", "t4"]:
        for name in df[col]:
            taxa.add(name)

    taxa = sorted(taxa, key=lambda x: int(x[1:]))
    n = len(taxa)
    q = len(df)
    #print(taxa) # Current taxas present in the csv
    #print("n =", len(taxa))

    # Creating dictionary and assigning each taxa its index
    taxon_to_i = {}
    for i, taxon in enumerate(taxa):
        taxon_to_i[taxon] = i
    #print(taxon_to_i)

    # Creating the tensor
    T = np.zeros((n, n, q), dtype=float) #rows, columns, quartets
    #print(T.shape)

    # Filling the tensor with quartet CFs
    for k,row in df.iterrows():
        a = row["t1"]
        b = row["t2"]
        c = row["t3"]
        d = row["t4"]

        #print(a,b,c,d)
        ia = taxon_to_i[a]
        ib = taxon_to_i[b]
        ic = taxon_to_i[c]
        id = taxon_to_i[d]
        #print(ia, ib, ic, id)

        cf12_34 = row["CF12_34"]
        cf13_24 = row["CF13_24"]
        cf14_23 = row["CF14_23"]
        #print(cf12_34)

        # Adds CF at index ia,ib,k... respectively

        T[ia,ib,k] += cf12_34
        T[ib,ia,k] += cf12_34
        T[ic,id,k] += cf12_34
        T[id,ic,k] += cf12_34

        T[ia,ic,k] += cf13_24
        T[ic,ia,k] += cf13_24
        T[ib,id,k] += cf13_24
        T[id,ib,k] += cf13_24

        T[ia,id,k] += cf14_23
        T[id,ia,k] += cf14_23
        T[ib,ic,k] += cf14_23
        T[ic,ib,k] += cf14_23

        #print("Quartet:", a, b, c, d)
        #print(T[:, :, k])
    print("Tensor shape:", T.shape)
    return T, taxa

# Step 1 - Collapsing the tensor across the quartet dimension to get a single adjacency matrix W
# Sum across all quartet layers to get total pairwise support
def collapse_tensor(T, method):
    if method == "sum":
        W = np.sum(T, axis=2)   
    elif method == "mean":
        W = np.mean(T, axis=2)  
    elif method == "max":
        W = np.max(T, axis=2)

    elif method == "l2":
        W = np.sqrt(np.sum(T**2, axis=2))

    else:
        raise ValueError("Invalid method. Use 'sum', 'mean', 'max', or 'l2'.")
    print("Step 1 - Collapsed adjacency matrix W (sum across all quartet layers)")
    print(W)
    return W


# Step 2 - Compute the Laplacian matrix L from W and D 
# D is a diagonal matrix where D[i,i] = sum of all edge weights connected to node i
def compute_laplacian(W):
    D = np.diag(np.sum(W, axis=1))  # Degree matrix 
    L = D - W  # Unnormalized Laplacian
    print("Step 2 - Graph Laplacian L = D - W")
    print("Degree matrix diagonal:", np.diag(D))
    print("Laplacian: ")
    print(np.round(L, 4))
    
    return L


# Step 3 - Computing the effective resistance matrix R from L
# The effective resistance between nodes i and j is:
#   R[i,j] = L_pinv[i,i] + L_pinv[j,j] - 2 * L_pinv[i,j]
# where L_pinv is the pseudoinverse of the Laplacian
def compute_effective_resistance(L):
    L_pinv = pinv(L) 
    n = L.shape[0]
    R = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            R[i,j] = L_pinv[i,i] + L_pinv[j,j] - 2 * L_pinv[i,j]
            R[j,i] = R[i,j]  # symmetric
     
    print("Step 3 — Effective resistance matrix R")
    print(np.round(R, 4))
    
    return R


# Step 4 - Compute edge sampling probabilities
# p[i,j] is proportional to W[i,j] * R[i,j]
# Edges that are both heavy AND structurally important get high probability
def compute_edge_sampling_probabilities(W, R, taxa):
    scores = W * R 
    n = W.shape[0]
    # Normalize so probabilities sum to 1 (only upper triangle to avoid double-counting)
    total_score = 0
    edge_list = []
    for i in range(n):
        for j in range(i+1, n):
            if W[i,j] > 0: # Redundant
                total_score += scores[i,j]
                edge_list.append((i, j, W[i,j], R[i,j], scores[i,j]))

    # Add normalized probability to each edge
    for idx in range(len(edge_list)):
        i, j, w, r, s = edge_list[idx]
        prob = s / total_score
        edge_list[idx] = (i, j, w, r, s, prob)

    print("Step 4 - Edge sampling probabilities")
    print(f"{'Edge':<12} {'Weight':<10} {'Eff.Res.':<10} {'Score':<10} {'Prob':<10}")
    for i, j, w, r, s, p in sorted(edge_list, key=lambda x: -x[5]):
        name_i = taxa[i]
        name_j = taxa[j]
        print(f"{name_i}-{name_j:<8} {w:<10.4f} {r:<10.4f} {s:<10.4f} {p:<10.4f}")
    return edge_list

# Step 5 - Sample edges to create sparsified graph
# Keep a subset of edges by sampling based on probabilities
# num_samples controls how sparse the result is
def sample_edges(edge_list, num_samples, n, taxa,seed):
    np.random.seed(seed) 
    probs = np.array([e[5] for e in edge_list])
    indices = np.arange(len(edge_list))
    # Sample edges (with replacement, then reweight)
    sampled_indices = np.random.choice(indices, size=num_samples, p=probs, replace=True)
    # Build the sparsified adjacency matrix
    W_sparse = np.zeros((n, n))
    for s_idx in sampled_indices:
        i, j, w, r, s, p = edge_list[s_idx]   
        # Reweight: divide by (num_samples * probability) to preserve expected value
        reweight = w / (num_samples * p)
        W_sparse[i,j] += reweight
        W_sparse[j,i] += reweight
    print("Step 5 - Sparsified adjacency matrix")
    print(np.round(W_sparse, 4))
    # Count how many unique edges survived
    surviving_edges = set()
    for s_idx in sampled_indices:
        i, j = edge_list[s_idx][0], edge_list[s_idx][1]
        surviving_edges.add((i,j))
     
    print(f"Original edges: {len(edge_list)}")
    print(f"Surviving unique edges: {len(surviving_edges)}")
    print(f"Edges removed: {len(edge_list) - len(surviving_edges)}")
    # Print surviving edges with species names
    print("Surviving edges:")
    for i, j in sorted(surviving_edges):
        print(f"  {taxa[i]} — {taxa[j]}: {W_sparse[i,j]:.4f}")
    
    return W_sparse



def visualize_graph(W_sparse, taxa, num_samples):
    G = nx.from_numpy_array(W_sparse)
     
    # Relabel nodes from numeric indices (0,1,2,...) to species names (t1, t2, t3,...)
    mapping = {i: taxa[i] for i in range(len(taxa))}
    G = nx.relabel_nodes(G, mapping)
     
    # Compute layout positions for nodes
    pos = nx.spring_layout(G, seed=42)  
     
    # Extract edge weights for visual encoding
    edges = G.edges(data=True)
    weights = [d['weight'] for (u, v, d) in edges]
     
    # Scale edge widths so the plot looks readable
    max_weight = max(weights) if weights else 1
    edge_widths = [3 * w / max_weight for w in weights]
     
    # Draw the graph
    plt.figure(figsize=(10, 8))
    nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=800)
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')
    nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color='gray', alpha=0.7)
     
    # Add edge weight labels (rounded to 2 decimals for readability)
    edge_labels = {(u, v): f"{d['weight']:.2f}" for (u, v, d) in edges}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)
     
    plt.title(f"Sparsified Quartet Graph (num_samples = {num_samples})")
    plt.axis('off')
    plt.tight_layout()
    plt.show()
