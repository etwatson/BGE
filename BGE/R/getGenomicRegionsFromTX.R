getGenomicRegionsFromTX <- function(transcriptId, txdb) {
  # Ensure transcriptIds is a character vector
  if (!is.character(transcriptId)) {
    stop("transcriptIds must be a character vector.")
  }
  # Create txdb for project
  #gffFile <- "/media/ark/tri/expression/Tcas5.2_rna_mt.gff"
  #txdb <- makeTxDbFromGFF(file = gffFile, format = "gff3")

  # Extract transcripts information from the TxDb object
  transcripts <- transcriptsBy(txdb, by = "gene")

  # Get GeneID for inquiry
  GeneID <- vector()
  for(i in 1:length(transcriptId)){
    GeneID <- c(GeneID,convertTxToGene(transcriptId[i],tx2gene))
  }


  # Subset the transcripts GRangesList by the provided transcript IDs
  selectedTranscripts <- transcripts[GeneID]

  # Reduce the GRangesList to a single GRanges object (if needed)

  reducedTranscripts <- GenomicRanges::reduce(selectedTranscripts)

  # Combine all reduced GRanges objects into a single GRanges object
  genomicRegions <- do.call(c, reducedTranscripts)
  return(genomicRegions)
}


