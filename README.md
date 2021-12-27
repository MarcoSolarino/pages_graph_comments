# Graph of the pages using comments: Maximal Cliques

- **Academic year:** 2020-2021
- **Course:** Advanced Algorithms and Graph Mining
- **Students:** [Marco Solarino](https://github.com/MarcoSolarino)
- **CFUs:** 6

# Overview
Given an undirecred graph G, a Maximal Clique [1] is an induced subgraph of G that is complete (clique) and to which no other adjacent vertex can be added maintaining the property of completeness. They have an important role in graph theory and computer science and have many applications. The problem of enumerating all the maximal cliques of a graph is NP-complete but many algorithms to solve this task exist. In this project I choose Bron-Kerbosch[3] (and two variants) that runs in exponential time.

In particular in this project are implemented:
- **Algorithm that builds an undirected graph given a subset of OEIS®[2] in JSON format**
- **Algorithm that given a graph G and a vertex v, finds a non trivial maximal clique**
- **Bron-Kerosch without pivoting**
- **Bron-Kerbosch with pivoting**
- **Bron-Kerbosch using a degeneracy order of the nodes of the graph[4]**
- **Algorithm that finds a maximal clique of max size**

# Bibliography
- [1] Maximal Clique: https://mathworld.wolfram.com/MaximalClique.html
- [2] Online Encyclopedia of Integer Sequences: https://oeis.org/
- [3] Bron-Kerbosch: https://dl.acm.org/doi/10.1145/362342.362367
- [4] Listing All Maximal Cliques in Sparse Graphs in Near-optimal Time: https://arxiv.org/abs/1006.5440
