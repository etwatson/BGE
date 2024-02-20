\name{convertTxToGene}
\alias{convertTxToGene}
\title{Convert a TXNAME to a GENEID}
\usage{
convertTxToGene(transcriptIds, tx2gene)
}
\description{
This script requires a tx2gene object to already be created, with a GENEID column
and a TXNAME column. When TXNAMEs are provided, they are converted to GENEIDs.


}
\examples{
# Create a tx2gene object if you haven't already
library(GenomicFeatures)
gffFile <- "file.gff"
txdb <- makeTxDbFromGFF(file = gffFile, format = "gff3")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Now run the function
convertTxToGene(transcriptIds, tx2gene)
}
