convertTxToGene <- function(txName, tx2gene) {
  # Check if txName is a character vector
  if (!is.character(txName)) {
    stop("txNames must be a character vector.")
  }

  # Check if tx2gene has the required columns
  if (!("TXNAME" %in% names(tx2gene)) || !("GENEID" %in% names(tx2gene))) {
    stop("tx2gene must have columns 'TXNAME' and 'GENEID'.")
  }
  # Map TXNAMEs to GENEIDs
  geneId <- tx2gene[tx2gene$TXNAME %in% txName,]$GENEID

  # Handle missing values
  missing <- is.na(geneId)
  if (any(missing)) {
    warning("Some TXNAMEs have no matching GENEID: ", paste(txName[missing], collapse = ", "))
  }

  return(geneId)
}
