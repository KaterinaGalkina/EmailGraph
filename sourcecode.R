
# ---------------------------------------------------
#                 Given graph Email-EU
# ---------------------------------------------------
# Part a - Results for 2a and 2b
# ---------------------------------------------------

# To load the package
library(igraph)

path <- file.choose()
opinions <- read.table(path)

optab <-as.matrix(opinions)
n <- nrow(optab)
v1 <- optab[1:n,1]
v2 <- optab[1:n,2]

relations<- data.frame(from=v1,to=v2)
relations

g <- graph_from_data_frame(relations,directed=TRUE) 

# Attempts to plot the whole graph :
plot(g) # not visible at all
# plot.igraph(g) # not visible at all also
# str(mygraph) # too much information, not very useful 

plot(g, edge.arrow.size=0.3, edge.color="gray", vertex.size=5, vertex.label=NA) # slightly better 

# The whole graph's representation but with communities 
# We have to make sure that all weights are positif
E(g)$weight <- rep(1, ecount(g))
V(g)$weight <- rep(1, vcount(g))
wc <- cluster_walktrap(g)
plot(wc, g, vertex.size = 1, layout = layout.fruchterman.reingold, vertex.label = NA, edge.arrow.size = 0.01)

# First way to simplify the graph, based on random normal variable to eliminate some vertices
E(g)$weight <- rnorm(ecount(g))
V(g)$weight <- rnorm(vcount(g))

sg <- igraph::induced_subgraph(g, which(igraph::V(g)$weight > 1.55))

# Remove edges with weight <= 0
sg <- delete_edges(sg, E(sg)[weight <= 0])

# Remove isolated vertices after edge deletion
sg <- delete_vertices(sg, igraph::degree(sg) == 0)
# plotting the simplified graph 
plot(sg, vertex.size=12, edge.arrow.size=0.5) 

# Identifying the communities in the simplified graph
wc2 <- cluster_walktrap(sg)
plot(wc2, sg, vertex.size=12, layout = layout.fruchterman.reingold, edge.arrow.size = 0.1)

# Second way to simplify : based on each node's number of neighbors
g_simplified <- delete_vertices(g, igraph::degree(g)<=250)
E(g_simplified)$weight <- sample(1, ecount(g_simplified), replace = T)
V(g_simplified)$weight <- sample(1, vcount(g_simplified), replace = T)
plot(g_simplified, vertex.color = "pink", vertex.label.color = "black", 
     edge.color = "darkslateblue", vertex.size = 14, edge.arrow.size=0.1)

# Third way to simplify 
# compute clusters ----
cl <- walktrap.community(g, steps = 5)

#contract vertices----
E(g)$weight <- 1
V(g)$weight <- 1

mygraph_simple <- contract(g, cl$membership, vertex.attr.comb = list(weight = "sum", 
                  name = function(x)x[1], "ignore"))

mygraph_simplified <- induced_subgraph(mygraph_simple, V(mygraph_simple)$weight > 30)
V(mygraph_simplified)$degree <- unname(igraph::degree(mygraph_simplified))
mygraph_simplified <- simplify(mygraph_simplified)
plot(mygraph_simplified, vertex.size = 14, edge.arrow.size=0.1)

hist(igraph::degree(mygraph_simplified))

# ---------------------------------------------------
# Functions from "An introduction to Graph Analytics"
# ---------------------------------------------------
# Part b - Results for 3
# ---------------------------------------------------

library(sna)

# Number of nodes
V(mygraph_simplified)

# Number of edges
E(mygraph_simplified)

# Adjacency matrix
mygraph_simplified.adj = igraph::as_adjacency_matrix(mygraph_simplified)
mygraph_simplified.adj

# Edge density
igraph::edge_density(mygraph_simplified)
igraph::edge_density(mygraph_simplified, loops = T)

# Graph density
mygraph_simplified_matrix <- as.matrix(mygraph_simplified.adj)
sna::gden(mygraph_simplified_matrix)

# Degree of the graph's nodes
degree <- igraph::degree(mygraph_simplified)
degree

# Histogram
hist(igraph::degree(mygraph_simplified))

#Betweenness Centrality
centr_betw(mygraph_simplified)

# Closeness Centrality 
centr_clo(mygraph_simplified)

# Check if the graph is simple
is_simple(mygraph_simplified)
# And to simplify it : mygraph_simplified <- simplify(mygraph_simplified)

# Geodesic
mygraph_simplified.edge_list <- as_edgelist(mygraph_simplified)
mygraph_simplified.geo <- geodist(mygraph_simplified.edge_list)
mygraph_simplified.geo


# ---------------------------------------------------
#               Part c - Answers for 4
# ---------------------------------------------------

# -------------------------------
# a - Central node in the graph 
# -------------------------------

# Degree centrality counting all neighbors
degree_number = igraph::degree(g, mode = "total")
sort(degree_number)[length(degree_number)]

# Degree centrality counting only successors
degree_number = igraph::degree(g, mode = "out")
sort(degree_number)[length(degree_number)]

# Eigenvector centrality 
eigen_centrality <- igraph::eigen_centrality(g)
highest_node <- which.max(eigen_centrality$vector)
highest_node

# Betweenneess centrality 
betweenness_scores <- igraph::betweenness(g)
highest_betweenness_node <- which.max(betweenness_scores)
highest_betweenness_node

# Closeness centrality
closeness_scores <- igraph::closeness(g)
highest_closeness_node <- which.max(closeness_scores)
highest_closeness_node


# ---------------------
# b - Longest path(s)
# ---------------------

# Calculating the graph's diameter
g_simplified.d <- igraph::diameter(g_simplified)

# Calculating the matrix of shortest paths
g_simplified.sp <- igraph::distances(g_simplified)

# Finding each pair of nodes between which the shortest distance is the diameter
longest_distances_nodes <- list()
for (i in 1:nrow(g_simplified.sp)) {
  for (j in 1:ncol(g_simplified.sp)) {
    if (g_simplified.sp[i, j] == g_simplified.d - 1) {
      longest_distances_nodes <- append(longest_distances_nodes, 
                    list(c(V(g_simplified)[i]$name, V(g_simplified)[j]$name)))
    }
  }
}

# The actual path between some pairs of vertices found
for (i in c(1:2)){
  pair = longest_distances_nodes[[i]]
  end_node = V(g_simplified)[name == pair[1]] 
  start_node = V(g_simplified)[name == pair[2]] 
  g_simplified.spvalues <- igraph::get.shortest.paths(graph = g_simplified, 
                    from = start_node, to = end_node)
  g_simplified.spvalues <- g_simplified.spvalues$vpath[[1]]
  cat("Path number", i, ":\n")
  print(g_simplified.spvalues)
  cat("\n")
}

# -----------------------
# c - Largest clique(s)
# -----------------------

# Finding the size of largest clique(s)
g.largest_clique_size = igraph::clique_num(g) 
g.largest_clique_size # maximum clique's size is 12 

# Finding all largest cliques
g.largest_cliques <- igraph::largest_cliques(g)

length(g.largest_cliques) # 24 largest cliques
g.largest_cliques # their list

# plotting some of them 
largest_clique_nodes <- g.largest_cliques[[1]]
largest_clique_subgraph <- igraph::induced_subgraph(g, largest_clique_nodes)
plot(largest_clique_subgraph, vertex.color = "lightgreen", edge.color = "black", 
     vertex.size = 30, edge.arrow.size = 0.5)

largest_clique_nodes <- g.largest_cliques[[5]]
largest_clique_subgraph <- igraph::induced_subgraph(g, largest_clique_nodes)
plot(largest_clique_subgraph, vertex.color = "pink", edge.color = "black", 
     vertex.size = 30, edge.arrow.size = 0.5)

largest_clique_nodes <- g.largest_cliques[[10]]
largest_clique_subgraph <- igraph::induced_subgraph(g, largest_clique_nodes)
plot(largest_clique_subgraph, vertex.color = "darkslategray1", edge.color = "black", 
     vertex.size = 30, edge.arrow.size = 0.5)


# ------------
# d - Ego(s)
# ------------

# Finding the ego(s) of the given graph
egos <- igraph::ego(g, mode = "out") 

# Finding the largest ego(s)' length
max_len <- max(sapply(egos, length))

# Identifying "root" nodes of the largest ego(s)
ego_nodes <- which(sapply(egos, length) == max_len)

# Getting their name(s) in the graph
ego_nodes_names <- V(g)$name[ego_nodes]
ego_nodes_names # In the given graph it is only the node "14648"

# Extracting the first (and here unique) largest ego
largest_ego_subgraph <- induced_subgraph(g, egos[[ego_nodes[1]]])

V(largest_ego_subgraph)$color <- "lightgreen"
E(largest_ego_subgraph)$color <- "gray"

# Emphasizing the egocentric node (the "root" of the subgraph) and its ongoing edges
V(largest_ego_subgraph)[name == ego_nodes_names[1]]$color <- "magenta"
E(largest_ego_subgraph)[.from(ego_nodes_names[1])]$color <- "blue"

plot(largest_ego_subgraph, vertex.size = 6, vertex.label = NA, edge.arrow.size = 0.1)

# How many nodes and edges are in this subgraph 
ecount(largest_ego_subgraph)
vcount(largest_ego_subgraph)


# ----------------------
# e - Power centrality 
# ----------------------

power_centrality_scores <- power_centrality(g_simplified)
power_centrality_scores <- sort(power_centrality_scores)
power_centrality_scores[length(power_centrality_scores)]






