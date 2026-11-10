# Social Network Analysis Report — Friendship Network (60 People)

**Author:** Yoga Vignesh V · **Course:** Network Analysis (KTU) Lab-8· **Tool:** R + igraph

---

## 1. Abstract

An anonymized ego-network of 60 people (Facebook friends, family, and university classmates) using R and igraph. The analysis covers visualization, centrality measures, clique calculation, and community detection, revealing a star-shaped ego-network with one dense classmate clique of size 25 and three communities (modularity Q = 0.440).

---

## 2. Introduction

Social network analysis (SNA) studies relationships between people (vertices) and their connections (edges) using graph theory. In this project, my personal friendship network was collected from three sources — Facebook, family members, and my university (Moodle) classmates — and modelled as an undirected weighted graph G = (V, E).

The goal is to answer four questions about the network:

1. How is the network structured visually? (visualization)
2. Who are the most central people? (centrality)
3. Which groups are fully interconnected? (cliques)
4. What natural communities exist? (community detection)

The analysis was performed in R with the `igraph` package (full code in `network_analysis.R`).

### 2.1 Data Description (Analytical Methods of Dataset)

**Sources and size.** The dataset contains **60 vertices** (people) and **336 edges** (relationships):

| Group | Count | Source |
|---|---|---|
| ME (ego) | 1 | self |
| Family | 2 | Facebook "Family members" |
| Classmates (incl. 8 friends) | 24 | Moodle course roster (MGDMDiM-6) |
| Facebook friends | 33 | Facebook friends list (top by mutual-friend count) |

**Privacy.** Every person is anonymized using the first letters of their first name and family name (e.g., A.B, J.T). Numeric suffixes resolve collisions (M.M, M.M2, M.M3).

**Edge weights.** Each edge carries a relationship type (family, friend, classmate, facebook_friend) and a weight: Facebook edges are weighted by the number of mutual friends, other ties by tie strength (family = 5, friend = 3, classmate = 1).

**Mathematical methods used.** Let G = (V, E) be an undirected graph with n = |V| vertices, m = |E| edges, and adjacency matrix A (a_ij = 1 if {i,j} ∈ E).

- **Degree:** k_i = Σ_j a_ij; normalized degree centrality C_D(i) = k_i / (n−1)
- **Betweenness:** C_B(v) = Σ_{s≠v≠t} σ_st(v) / σ_st, where σ_st = number of shortest paths from s to t and σ_st(v) = how many pass through v
- **Closeness:** C_C(v) = (n−1) / Σ_{t≠v} d(v,t), where d(v,t) is the shortest-path distance
- **Eigenvector centrality:** the principal eigenvector of A, A x = λ_max x
- **PageRank:** PR(i) = (1−d)/n + d Σ_{j→i} PR(j)/k_j with damping d ≈ 0.85
- **Clique:** a subset S ⊆ V where every pair of vertices is adjacent; ω(G) is the maximum clique size
- **Modularity:** Q = (1/2m) Σ_{i,j} [ a_ij − k_i·k_j/(2m) ] · δ(c_i, c_j); communities were detected with Louvain (modularity maximization)

Because edge weights represent tie **strength** (not distance), betweenness and closeness were computed on the unweighted version of the graph.

### 2.2 Visualizations (Results and Findings)

**Network overview.** The static plot shows a classic ego-network: the red ME node (YV) sits in the centre connected to everyone; a dense cluster of blue/green classmates and friends occupies one region; magenta Facebook friends radiate outward as pendant vertices; the two yellow family nodes form a small triangle with ME. An interactive force-directed version (networkD3, zoomable, colored by role) is saved as `networkInteractive.html`.

<img width="1180" height="946" alt="Rplot1" src="https://github.com/user-attachments/assets/152fd06b-4be4-4d94-9831-f0cee084db32" />
                                                        Figure 1 — Network overview.

**Table 2 — Centrality comparison (top 15 by degree).**

| id | role | degree | degree_norm | betweenness | closeness | eigenvector | pagerank |
|---|---|---|---|---|---|---|---|
| YV | ME | 59 | 1.0000 | 1434 | 0.01695 | 0.7848 | 0.27649 |
| B.B | friend | 31 | 0.5254 | 0 | 0.01064 | 1.0000 | 0.02254 |
| A.B | friend | 31 | 0.5254 | 0 | 0.01064 | 1.0000 | 0.02254 |
| Y.D | friend | 31 | 0.5254 | 0 | 0.01064 | 1.0000 | 0.02254 |
| G.J | friend | 31 | 0.5254 | 0 | 0.01064 | 1.0000 | 0.02254 |
| M.M2 | friend | 31 | 0.5254 | 0 | 0.01064 | 1.0000 | 0.02254 |
| S.N | friend | 31 | 0.5254 | 0 | 0.01064 | 1.0000 | 0.02254 |
| E.O | friend | 31 | 0.5254 | 0 | 0.01064 | 1.0000 | 0.02254 |
| J.T | friend | 31 | 0.5254 | 0 | 0.01064 | 1.0000 | 0.02254 |
| K.B | classmate | 24 | 0.4068 | 0 | 0.01064 | 0.7477 | 0.01469 |

<img width="1180" height="946" alt="plot3" src="https://github.com/user-attachments/assets/c6b9c31f-35b7-4f38-9f4e-b2fe7b3899ef" />
Figure 2 - Centrality comparison plot


**Communities (Louvain, Q = 0.440).** Three communities were detected: (1) ME together with the Facebook friends, (2) the family pair, and (3) the dense classmate/friend cluster.

<img width="1180" height="946" alt="Rplot-2" src="https://github.com/user-attachments/assets/77071b01-1b93-421c-80c1-03b84390f8d2" />

Figure 3 — Communities (Louvain, Q = 0.440)

**Cliques.**

- Clique number ω(G) = **25** (a single maximum clique consisting of ME + all 24 classmates — everyone in the course is connected to everyone else)
- Number of maximal cliques = 2: the size-25 classmate clique and one size-3 family triangle (S., V.R, YV)
- Clique size distribution shows a sharp peak at sizes 12–14 (over 5.2 million cliques each), all contained inside the classmate cluster

---

## 3. Discussion

The results describe a typical **ego-network with a star core and one dense community**:

- **ME dominates every centrality measure.** Degree 59 (connected to everyone), betweenness 1434 (nearly all shortest paths between Facebook friends and classmates pass through ME), and the highest PageRank (0.276). This is expected: the network was collected from ME's point of view, so ME is a structural bottleneck.
- **The 8 university friends are the second tier.** They combine the classmate clique (degree 31 = 24 classmates + ME + family + friend-circle edges) with perfect eigenvector centrality (1.0000) because they sit in the densest part of the graph — being connected to the well-connected classmate cluster boosts their score beyond even ME (0.7848).
- **Facebook friends are pendant vertices.** They connect only to ME, giving degree 1 and zero betweenness. They widen the network but add no internal structure.
- **Cliques mirror the data design.** The single clique of size 25 (ME + all classmates) shows the class group is a complete subgraph — every classmate knows every other. The only other maximal clique is the family triangle.
- **Community detection separates life domains.** Louvain (Q = 0.440, a moderately strong partition) found three communities matching real social circles: university (classmates/friends), family, and online friends (grouped around ME).

**Limitations.** The network is ego-centric, so alter–alter ties are only known inside the classmate group and family; Facebook friends' mutual connections were approximated by mutual-friend counts (edge weights) rather than observed edges. The classmate clique is complete by construction (same course), which inflates clique size and eigenvector scores.

---

## 4. Conclusions

The friendship network of 60 people is a **star-shaped ego-network with a dense university cluster**: ME is the single most central actor (degree 59, betweenness 1434, PageRank 0.276), the university classmates form one giant clique of size 25, and the network splits into three natural communities (university, family, online friends) with modularity Q = 0.440. Eigenvector centrality highlights the university friend group as the network's structural core, while Facebook friends form a wide but shallow periphery. Overall, the analysis demonstrates how different centrality measures and community detection reveal complementary views of the same social structure.

---

## 5. Appendix (Links)

**Dataset and code (project files):** `friendship_dataset.csv`, `network_analysis.R`, `networkInteractive.html`, `centrality_results.csv`, `community_results.csv`

**References:**

1. Network analysis with R — https://bookdown.org/jdholster1/idsr/network-analysis.html
2. Network Analysis and Visualization with R — Kateto.net — https://kateto.net/netscix2016.html
3. Network Analysis in R — Schochastics — https://schochastics.github.io/netAnaR/
4. igraph package vignette — https://cran.r-project.org/web/packages/igraph/vignettes/igraph.html
