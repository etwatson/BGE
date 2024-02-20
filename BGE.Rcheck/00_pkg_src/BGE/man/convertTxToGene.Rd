\name{convertTxToGene}
\alias{convertTxToGene}
\title{Convert a TXNAME to a GENEID}
\usage{
convertTxToGene(transcriptIds, tx2geneFilePath, header = FALSE, sep = "\t")
}
\description{
This script requires a tx2gene file to already be created, with a GENEID column
and a TXNAME column. When TXNAMEs are provided, they are converted to GENEIDs,

If a tx2gene file is not available, one can be made in this manner:
library(GenomicFeatures)
gffFile <- "file.gff"
txdb <- makeTxDbFromGFF(file = gffFile, format = "gff3")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
}
\examples{
convertTxToGene(transcriptIds, tx2geneFilePath, header = FALSE, sep = "\t")
}
