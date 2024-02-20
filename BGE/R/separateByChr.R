separateByChr <- function(grangesList) {
  chromLists <- list()

  for (i in seq_along(grangesList)) {
    currentElement <- grangesList[[i]]
    geneName <- names(grangesList)[i]  # Extract the name (GENEID) of the current element

    uniqueChroms <- unique(seqnames(currentElement))

    for (chrom in uniqueChroms) {
      chromElement <- currentElement[seqnames(currentElement) == chrom]
      names(chromElement) <- geneName  # Set the name (GENEID) for the current chromElement

      if (is.null(chromLists[[as.character(chrom)]])) {
        chromLists[[as.character(chrom)]] <- GRangesList(chromElement)
      } else {
        chromLists[[as.character(chrom)]] <- c(chromLists[[as.character(chrom)]], GRangesList(chromElement))
      }
    }
  }

  names(chromLists) <- unique(seqnames(unlist(grangesList)))

  return(chromLists)
}
