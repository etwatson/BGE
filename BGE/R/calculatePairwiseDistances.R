calculatePairwiseDistances <- function(chromSeparatedList) {
  # Initialize an empty list to store the distance matrices for each chromosome
  distanceMatrices <- list()
  # Loop through each chromosome in the chromSeparatedList
  for (chrom in names(chromSeparatedList)) {
    grangesList <- chromSeparatedList[[chrom]]
    n <- length(grangesList)

    # Initialize a matrix to store distances for the current chromosome
    distanceMatrix <- matrix(nrow = n, ncol = n, data = NA)
    names_c <-vector()
    for (k in 1:length(grangesList)){ names_c <- c(names_c,names(grangesList[[k]]))}
    rownames(distanceMatrix) <- names_c
    colnames(distanceMatrix) <- names_c

    # Calculate pairwise distances
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i == j) {
          distanceMatrix[i, j] <- 0  # Distance to self is 0
        } else {
          # Calculate the minimum distance between genes i and j
          dist_ij <- min(abs(start(grangesList[[i]]) - end(grangesList[[j]])),
                         abs(start(grangesList[[j]]) - end(grangesList[[i]])))/1000000

          # If genes overlap, set distance to 0
          if (any(overlapsAny(grangesList[[i]], grangesList[[j]]))) {
            dist_ij <- 0
          }


          distanceMatrix[i, j] <- dist_ij
        }
      }
    }

    # Add the computed distance matrix to the distanceMatrices list
    distanceMatrices[[chrom]] <- distanceMatrix
  }

  # Set the names of the list elements to be the chromosome IDs
  names(distanceMatrices) <- names(chromSeparatedList)

  return(distanceMatrices)
}
