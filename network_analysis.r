# =====================================================================
# Social Network Analysis of a Facebook / Classmate Friendship Network
# Dataset: friendship_dataset.csv  (60 people, 364 relationships)
# Tasks: a) visualization  b) centrality  c) cliques  d) communities
# =====================================================================

# ---- 0. Setup -------------------------------------------------------
# install.packages(c("igraph", "networkD3", "htmlwidgets"))  # run once if needed
library(igraph)

set.seed(42)  # reproducible layouts

# ---- 1. Load dataset ------------------------------------------------
# CSV columns: from, from_name, to, to_name, relationship, weight
df <- read.csv("friendship_dataset.csv", stringsAsFactors = FALSE)

# keep only the anonymized columns (privacy: names are for the author only)
edges <- df[, c("from", "to", "relationship", "weight")]

# ---- 2. Build the graph ---------------------------------------------
# Undirected weighted graph  G = (V, E), w: E -> R+
g <- graph_from_data_frame(
  d  = edges[, c("from", "to", "relationship", "weight")],
  directed = FALSE
)
E(g)$relationship <- as.character(edges$relationship)
E(g)$weight       <- as.numeric(edges$weight)

# Vertex role: how each person is tied to the ego (YV)
role <- rep("classmate", vcount(g))
names(role) <- V(g)$name
ego_edges  <- incident_edges(g, "YV")[[1]]
for (e in seq_along(ego_edges)) {
  v <- ends(g, ego_edges[e])[, 2]
  if (v != "YV") role[v] <- as.character(E(g)$relationship[ego_edges[e]])
}
role["YV"] <- "ego"
V(g)$role <- role

cat("Graph summary\n")
print(g)
cat("Vertices:", vcount(g), " Edges:", ecount(g), "\n")

# =====================================================================
# 3. METHODS (mathematical description)
# =====================================================================
# Let G = (V, E) be an undirected graph with n = |V| vertices and
# m = |E| edges, and A its adjacency matrix (a_ij = 1 if {i,j} in E).
#
# (1) DEGREE of vertex i:
#         k_i = SUM_j a_ij
#     Normalized degree centrality:
#         C_D(i) = k_i / (n - 1)
#
# (2) BETWEENNESS centrality of vertex v:
#         C_B(v) = SUM_{s != v != t} sigma_st(v) / sigma_st
#     where sigma_st = number of shortest paths s -> t and
#     sigma_st(v) = number of those passing through v.
#
# (3) CLOSENESS centrality of vertex v:
#         C_C(v) = (n - 1) / SUM_{t != v} d(v, t)
#     where d(v,t) is the graph-theoretic distance.
#
# (4) EIGENVECTOR centrality: principal eigenvector of A:
#         A x = lambda_max x ,  x_i = (1/lambda) SUM_j a_ij x_j
#
# (5) PAGERANK (with damping d ~ 0.85):
#         PR(i) = (1-d)/n + d * SUM_{j -> i} PR(j)/k_j
#
# (6) CLIQUE: a subset S of V such that every pair in S is adjacent.
#     Maximum clique size = omega(G).
#
# (7) MODULARITY of a partition C = {c_1,...,c_k}:
#         Q = (1/2m) SUM_{i,j} [ a_ij - k_i k_j / (2m) ] * delta(c_i, c_j)
#     Communities are found by maximizing Q (Louvain, greedy) or by
#     iteratively removing high-betweenness edges (edge betweenness).

# =====================================================================
# 4. TASK (a): DATA VISUALIZATION
# =====================================================================
cols <- c(ego = "red", family = "gold", friend = "forestgreen",
          classmate = "steelblue", facebook_friend = "orchid")

# ---------- (a1) static plot (goes into the PDF report) --------------
plot(g,
     layout = layout_with_fr(g),          # Fruchterman-Reingold
     vertex.color  = cols[V(g)$role],
     vertex.size   = 3 + 2 * sqrt(degree(g)),
     vertex.label  = V(g)$name,
     vertex.label.cex = 0.55,
     edge.width    = 0.4 + 0.25 * E(g)$weight,
     edge.color    = "grey65",
     main = "Friendship network (60 people)")

legend("topleft", legend = names(cols), col = cols, pch = 16, cex = 0.8)

# ---------- (a2) INTERACTIVE graph (networkD3) -----------------------
library(networkD3)
library(htmlwidgets)

# convert the igraph object to networkD3 format
# (nodes: name + role group; links: 0-based source/target indices)
d3 <- igraph_to_networkD3(g, group = V(g)$role)
d3$nodes$degree <- degree(g)              # used for node size
d3$links$value  <- E(g)$weight            # used for link width

# note: simpleNetwork() supports only ONE colour for all nodes,
# so we use forceNetwork() to keep the per-role colours + zoom
p <- forceNetwork(
  Links  = d3$links,                      # edge list
  Nodes  = d3$nodes,                      # node list
  Source = "source",                      # 0-based index of source node
  Target = "target",                      # 0-based index of target node
  Value  = "value",                       # edge weight -> link width
  NodeID = "name",                        # node label (anonymized initials)
  Nodesize = "degree",                    # node size by degree k_i
  Group  = "group",                       # colour by role
  colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"),
  charge = -150,                          # node repulsion (negative)
  linkDistance = 60,                      # space between nodes
  fontSize = 14,
  fontFamily = "serif",
  opacity = 0.9,
  zoom = TRUE,                            # zoom + drag enabled
  legend = TRUE,                          # role legend
  bounded = FALSE,
  width = 900, height = 700
)
p                                        # view in RStudio Viewer

# save as a standalone HTML file (open in any browser)
saveWidget(p, file = "networkInteractive.html")

# =====================================================================
# 5. TASK (b): CENTRALITY CALCULATIONS
# =====================================================================
# NOTE: weights are tie STRENGTHS, not distances. igraph's weighted
# betweenness/closeness treat weights as distances, so we compute the
# unweighted versions (weights = NA) for a correct interpretation.
deg   <- degree(g)                          # k_i
deg_c <- degree(g, normalized = TRUE)       # C_D(i) = k_i/(n-1)
bet   <- betweenness(g, weights = NA)       # C_B(v)
clo   <- closeness(g, weights = NA)         # C_C(v)
eig   <- eigen_centrality(g, weights = NA)$vector  # A x = lambda x
pr    <- page_rank(g, directed = FALSE)$vector    # PageRank

centrality <- data.frame(
  id = V(g)$name, role = V(g)$role,
  degree = round(deg, 0), degree_norm = round(deg_c, 4),
  betweenness = round(bet, 2), closeness = round(clo, 5),
  eigenvector = round(eig, 4), pagerank = round(pr, 5)
)
centrality <- centrality[order(-centrality$degree), ]
print(head(centrality, 10))     # top-10 most central people
write.csv(centrality, "centrality_results.csv", row.names = FALSE)

# quick comparison plot of the three main measures (top 15 by degree)
top <- head(centrality$id, 15)
op <- par(mfrow = c(1, 3), mar = c(7, 4, 3, 1))
barplot(deg[top],  names.arg = top, las = 2, col = "steelblue",
        main = "Degree", border = NA)
barplot(bet[top],  names.arg = top, las = 2, col = "forestgreen",
        main = "Betweenness", border = NA)
barplot(clo[top],  names.arg = top, las = 2, col = "orchid",
        main = "Closeness", border = NA)
par(op)

# =====================================================================
# 6. TASK (c): CLIQUES CALCULATION
# =====================================================================
clique_sizes <- clique_size_counts(g)        # how many cliques per size
cat("Clique number omega(G) =", clique_num(g), "\n")
cat("Number of maximal cliques =", count_max_cliques(g), "\n")
cat("Clique size distribution:\n"); print(clique_sizes)

# all cliques of maximum size
max_cl <- largest_cliques(g)
cat("Largest clique size:", length(max_cl[[1]]),
    " count:", length(max_cl), "\n")
for (cl in max_cl) cat("  ", paste(sort(names(cl)), collapse = ", "), "\n")

# top-10 largest maximal cliques by size
mc <- max_cliques(g, min = 3)                # cliques of size >= 3
mc <- mc[order(-lengths(mc))]
cat("\nTen largest maximal cliques (size >= 3):\n")
for (i in seq_len(min(10, length(mc))))
  cat(sprintf("  size %d: %s\n", lengths(mc)[i],
              paste(sort(names(mc[[i]])), collapse = ", ")))

# =====================================================================
# 7. TASK (d): COMMUNITY DETECTION
# =====================================================================
# Louvain modularity optimization
com_louv  <- cluster_louvain(g, weights = E(g)$weight)
# Greedy modularity (Clauset-Newman-Moore)
com_greedy<- cluster_fast_greedy(g, weights = E(g)$weight)
# Edge betweenness (Girvan-Newman)
com_eb    <- cluster_edge_betweenness(g, weights = NA, directed = FALSE)
# Infomap
com_info  <- cluster_infomap(g, e.weights = E(g)$weight)

cat("Modularity Q (Louvain):        ", round(modularity(com_louv), 4), "\n")
cat("Modularity Q (Greedy):         ", round(modularity(com_greedy), 4), "\n")
cat("Modularity Q (Edge betw.):     ", round(modularity(com_eb), 4), "\n")
cat("Modularity Q (Infomap):        ", round(modularity(com_info), 4), "\n")

for (nm in c("Louvain", "Greedy", "EdgeBetweenness", "Infomap")) {
  com <- switch(nm, Louvain = com_louv, Greedy = com_greedy,
                EdgeBetweenness = com_eb, Infomap = com_info)
  cat(sprintf("\n--- %s: %d communities ---\n", nm, length(com)))
  for (i in seq_len(no.clusters(com)))
    cat(sprintf("  community %2d (%2d members): %s\n", i, sizes(com)[i],
                paste(sort(names(membership(com))[membership(com) == i]),
                      collapse = ", ")))
}

# visualize the best partition (Louvain)
plot(com_louv, g,
     layout = layout_with_fr(g),
     vertex.label = V(g)$name, vertex.label.cex = 0.55,
     vertex.size = 3 + 2 * sqrt(degree(g)),
     main = sprintf("Communities (Louvain, Q = %.3f)", modularity(com_louv)))

# export memberships
communities_out <- data.frame(
  id = V(g)$name, role = V(g)$role,
  louvain = as.integer(membership(com_louv)),
  greedy  = as.integer(membership(com_greedy)),
  edge_betweenness = as.integer(membership(com_eb)),
  infomap = as.integer(membership(com_info))
)
write.csv(communities_out, "community_results.csv", row.names = FALSE)

# =====================================================================
# 8. SESSION INFO (for the report)
# =====================================================================
sessionInfo()
