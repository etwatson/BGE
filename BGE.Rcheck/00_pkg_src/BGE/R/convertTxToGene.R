convertTxToGene <- function(transcriptIds, tx2geneFilePath, header = FALSE, sep = "\t") {
  # Load the tx2gene mapping file
  tx2gene <- read.csv(tx2geneFilePath, header = header, sep = sep, stringsAsFactors = FALSE)

  # Ensure the tx2gene object has appropriate column names
  if (!header) {
    colnames(tx2gene) <- c("TXNAME", "GENEID")
  }

  # Create a named vector for the conversion
  tx2geneVector <- setNames(tx2gene$gene, tx2gene$transcript)

  # Convert transcript IDs to gene IDs
  geneIds <- sapply(transcriptIds, function(txId) tx2geneVector[txId])

  return(geneIds)
}
