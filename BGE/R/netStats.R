# Load required libraries
library(WGCNA,quietly = TRUE)
library(igraph,quietly = TRUE)
library(dplyr,quietly = TRUE)

load("expressionData.Rdata")

# Set global options
options(stringsAsFactors = FALSE)

# Enable multi-threading within WGCNA
enableWGCNAThreads()

# Load metadata and batch-corrected counts
metadata <- read.csv("metadata.csv", stringsAsFactors = TRUE)
metadata$Time <- factor(metadata$Time) # Ensure Time is treated as a factor
stillmeta <- metadata %>% filter(Condition == "still")
normmeta <- metadata %>% filter(Condition == "norm")
pupmeta <- metadata %>% filter(Stage == "pupae")
adultmeta <- metadata %>% filter(Stage == "adult")
stillexp <- expressionData[row.names(expressionData) %in% stillmeta$Sample_ID,]
normexp <- expressionData[row.names(expressionData) %in% normmeta$Sample_ID,]
pupexp <- expressionData[row.names(expressionData) %in% pupmeta$Sample_ID,]
adultexp <- expressionData[row.names(expressionData) %in% adultmeta$Sample_ID,]

# Create adjacency matrix and signed, undirected graph
adjacencyMatrix<- adjacency(pupexp,power = 4,type = "signed")
graph <- graph_from_adjacency_matrix(adjacencyMatrix, mode = "undirected", weighted = TRUE)

# Extract edge weights
edge_weights <- E(graph)$weight

# Identify positive and negative edges
positive_edges <- edge_weights[edge_weights > 0]
negative_edges <- edge_weights[edge_weights < 0]

# Calculate strength (sum of weights) for positive and negative interactions separately
# For each node, sum weights of positive and negative edges connected to it
positive_strength <- rep(0, vcount(graph))
negative_strength <- rep(0, vcount(graph))

# Get all edge weights
edge_weights <- E(graph)$weight

# Get indices of positive and negative weights
positive_indices <- which(edge_weights > 0)
negative_indices <- which(edge_weights < 0)

# For positive weights
if (length(positive_indices) > 0) {
  edges <- get.edgelist(graph)[positive_indices, ]
  for (i in 1:nrow(edges)) {
    positive_strength[edges[i, ]] <- positive_strength[edges[i, ]] + edge_weights[positive_indices[i]]
  }
}

# For negative weights (store absolute values)
if (length(negative_indices) > 0) {
  edges <- get.edgelist(graph)[negative_indices, ]
  for (i in 1:nrow(edges)) {
    negative_strength[edges[i, ]] <- negative_strength[edges[i, ]] + abs(edge_weights[negative_indices[i]])
  }
}


V(graph)$positive_strength <- positive_strength
V(graph)$negative_strength <- negative_strength


# Create sub-networks
positive_subgraph <- subgraph.edges(graph, which(E(graph)$weight > 0), delete.vertices = FALSE)
negative_subgraph <- subgraph.edges(graph, which(E(graph)$weight < 0), delete.vertices = FALSE)

# Calculate the clustering coefficient for each sub-network
positive_clustering <- transitivity(positive_subgraph, type = "localaverage")
negative_clustering <- transitivity(negative_subgraph, type = "localaverage")

# Local clustering coefficients for each node
positive_local_clustering <- transitivity(positive_subgraph, type = "local")
negative_local_clustering <- transitivity(negative_subgraph, type = "local")

# Betweenness centrality
betweenness_centrality <- betweenness(graph, directed = FALSE, weights = E(graph)$weight)

# Closeness centrality
closeness_centrality <- closeness(graph, mode = "all", weights = E(graph)$weight)

# Modularity
community <- cluster_louvain(graph, weights = E(graph)$weight)
modularity_score <- modularity(community)

# Eigenvector centrality
eig_centrality <- eigen_centrality(graph, weights = E(graph)$weight)$vector
#
pupal_graph <- graph
pupal_positive_strength <- positive_strength
pupal_negative_strength <- negative_strength

save(eig_centrality,
     modularity_score,
     closeness_centrality,
     betweenness_centrality,
     positive_clustering,
     negative_clustering,
     positive_local_clustering,
     negative_local_clustering,
     pupal_positive_strength,
     pupal_negative_strength,
     file = "pupae.stats.Rdata")
