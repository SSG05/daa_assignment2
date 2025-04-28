# Exact and Core-Exact Algorithms

---

## Exact Algorithm

The Exact algorithm is designed to find the optimal solution for the Clique Densest Subgraph (CDS) problem. It systematically explores all possible subsets of vertices to identify the subgraph that maximizes the ratio of the number of cliques to the number of vertices. 

While it guarantees the best possible result, it can be computationally expensive for large graphs due to its exhaustive search.

---

## Core-Exact Algorithm

The Core-Exact algorithm improves the efficiency of the Exact approach by first reducing the graph using a core decomposition. It filters out vertices that cannot be part of the densest subgraph, significantly shrinking the search space. 

Then, the Exact algorithm is applied only to this reduced graph, achieving the same optimal result with much faster performance on real-world datasets.

---
