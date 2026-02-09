
# Libraries ---------------------------------------------------------------

library(ChIPseeker)

# Parameters --------------------------------------------------------------
dar_res_path <- as.character(commandArgs(TRUE)[1])
comparison <- as.character(commandArgs(TRUE)[2])
direction <- as.character(commandArgs(TRUE)[3])
log_thr <- as.numeric(commandArgs(TRUE)[4])
output_path <- as.character(commandArgs(TRUE)[5])

dars <- read.table(dar_res_path, sep = "\t", header = T)

if(nrow(dars) != 0){
  peak_gr <- Signac::StringToGRanges(dars$peak, sep = c(":", "-"))
  if(grepl("KPlung", dar_res_path)){
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    annot <- "org.Mm.eg.db"
  }else if(grepl("human", dar_res_path)){ 
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    annot <- "org.Hs.eg.db"
  }else if(grepl("NfatKO", dar_res_path)){ 
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    annot <- "org.Mm.eg.db"
  }
  
  anno <- try({
    ChIPseeker::annotatePeak(
      peak_gr,
      TxDb = txdb,
      addFlankGeneInfo = TRUE,
      overlap = "all", #any nearest gene is reported regardless of the overlap with the TSS
      annoDb = annot
    )
  }, silent = TRUE)
  
  # Check if an error occurred
  if (inherits(anno, "try-error")) {
    message("An error occurred, saving dars without annotation.")
    merged_data <- dars
  }else{
    
    # Save annotation
    anno_df = as.data.frame(anno)
    head(anno_df)
    
    anno_df$annotation_type = anno_df$annotation
    anno_df$annotation_type[grep("Intron",anno_df$annotation_type)]="Intron"
    anno_df$annotation_type[grep("Exon",anno_df$annotation_type)]="Exon"
    anno_df$annotation_type[grep("Promoter",anno_df$annotation_type)]="Promoter"
    anno_df$annotation_type[grep("Downstream",anno_df$annotation_type)]="Downstream"
    
    anno_df$peak <- paste0(anno_df$seqnames, ":",anno_df$start, "-", anno_df$end)
    table(anno_df$annotation_type)
    
    merged_data <- merge(dars, anno_df, by = "peak")
  }
  write.table(merged_data, file = paste0(output_path, "annotated_",comparison, "_", direction, "_logFC", log_thr, ".txt"), row.names = F, col.names = T, sep = "\t", quote = F)
  
}else{
  write.table(dars, file = paste0(output_path, "annotated_",comparison, "_", direction, "_logFC", log_thr, ".txt"), row.names = F, col.names = T, sep = "\t", quote = F)
  
}

