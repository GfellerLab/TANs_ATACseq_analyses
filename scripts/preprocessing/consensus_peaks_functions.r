
library(PEPATACr)
library(dplyr)
library(yaml)
library(rtracklayer)
library(Biostrings)
library(SummarizedExperiment)
library(GenomeInfoDb)
library(GenomicRanges)
library(readr)
library(edgeR)
library(BSgenome)


extendedPeakSet <- function(df, BSgenome = NULL, blacklist = NULL, extend = 250, scorePerMillion = 2, selectionRules = "(n+1)/2", final_converge = T){
  
  #a lot faster than import.bed from rtracklayer
  readSummits <- function(file){

    df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr","start","end","name","score1", "strand", "signifValue", "score", "qval","peak"))))
    df <- df[,c(1,2,3,8)] #do not keep name column it can make the size really large

    return(GenomicRanges::makeGRangesFromDataFrame(df=df,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE))
  }


  #Error-------
  stopifnot(extend > 0)
  stopifnot("samples" %in% colnames(df))
  stopifnot("groups" %in% colnames(df))
  stopifnot("summits" %in% colnames(df))
  stopifnot(!is.null(BSgenome))
  # stopifnot(all(apply(df,1,function(x){file.exists(paste0(x[3]))})))
  #--------------
  #Deal with blacklist
  if(is.null(blacklist)){
    blacklist <- GRanges()
  }else if(is.character(blacklist)){
    blacklist <- rtracklayer::import.bed(blacklist)
  }
  stopifnot(inherits(blacklist,"GenomicRanges"))
  #--------------
  
  chromSizes <- getChromSizes(BSgenome = BSgenome)
  groups <- unique(df$groups)
  
  groupSEList <- lapply(seq_along(groups), function(i){
    
    message(groups[i])
    
    df_group = df[which(df$groups==groups[i]),]
    
    #--- 1. read in summit files
    grList <- GenomicRanges::GRangesList(lapply(paste0(df_group$summits),readSummits))
    
    #--- 2. resize
    #--- 3. within chromsizes
    #--- 4. not in blacklist
    #--- 5. non-overlapping
    #--- 6. score-per-million
    grList <- lapply(seq_along(grList),function(x){
      extended_summits <- grList[[x]] %>%
        resize(., width = 2 * extend + 1, fix = "center") %>%
        subsetByOverlaps(.,chromSizes,type="within") %>%
        subsetByOverlaps(.,blacklist,invert=TRUE) %>%
        convergeClusterGRanges(., by = "score", decreasing = T)
      mcols(extended_summits)$score <- edgeR::cpm(mcols(extended_summits)$score)
      return(extended_summits)
    })
    
    #--- 7. non-overlapping overlaps determined score-per-million
    grNonOverlapping <- GenomicRanges::GRangesList(grList) %>%
      unlist %>%
      convergeClusterGRanges(., by = "score", decreasing = T) %>%
      sortSeqlevels %>%
      sort
    
    #--- 8. select those meeting score-per-million threshold
    grNonOverlapping <- grNonOverlapping[which(mcols(grNonOverlapping)$score > scorePerMillion),]
    
    #--- 9. get score-per-million for each sample
    scores <- lapply(seq_along(grList), function(x){
      columnOverlaps(query = grNonOverlapping,
                     subject = grList[[x]],
                     colname = "score",
                     decreasing = TRUE)
    }) %>% data.frame()
    colnames(scores) <- df_group$samples
    
    #--- 10. create summarized experiment
    se <- SummarizedExperiment(assays = SimpleList(scores = as.matrix(scores)), rowRanges = grNonOverlapping)
    
    #free up some memory
    remove(grList,grNonOverlapping)
    
    #--- 11. apply selection rules
    n = ncol(se)
    minSamples <- eval(parse(text=selectionRules))
    if(n >= minSamples){
      message("Applying selection criteria: minimum of ", minSamples," samples with a score-per-million of ", scorePerMillion)
      keep <- which(rowSums(assay(se) > scorePerMillion) >= minSamples)
      message("Final ", groups[i], " Peak Set ", length(keep))
      se <- se[keep,]
    }else{
      warning("group with less than min samples required by selectionRules (most likely 1 summit for group)")
    }
    
    #--- 12. re-normalize score-per-million
    mcols(se)$score <- edgeR::cpm(mcols(se)$score) #Re-scale for final comparison #10^6*mcols(se)$score/sum(mcols(se)$score)
    mcols(se)$name <- paste0(groups[i],"_",seq_len(nrow(se))) #add short name
    return(se)
    
  }) #END groupSEList lapply
  
  names(groupSEList) <- groups
  
  #Create GRanges List for all Groups
  groupGRlist <- lapply(groupSEList, rowRanges)
  #Get Non Overlapping Peaks for all groups
  message("Now Getting Non Overlapping Peak Set Between Groups...")
  
  if(final_converge){
    #there may be some R version-specific errors from this command
    finalGR <- Reduce("c",groupGRlist) %>%
      convergeClusterGRanges(., by = "score", decreasing = T) %>%
      sortSeqlevels %>%
      sort
    
    message("Determined ", length(finalGR), " Peaks")
    
    return(list(groupSE = groupSEList, finalGR = finalGR))
  }else{
    return(list(groupSE = groupSEList, groupGRlist = groupGRlist))
  }
  
  
}

collapse_peaks_from_GRList = function(grList, group_name, scorePerMillion = 2, selectionRules = "(n+1)/2"){
  
  #--- non-overlapping overlaps determined score-per-million
  grNonOverlapping <- GenomicRanges::GRangesList(grList) %>%
    unlist %>%
    convergeClusterGRanges(., by = "score", decreasing = T) %>%
    sortSeqlevels %>%
    sort
  
  #--- select those meeting score-per-million threshold
  grNonOverlapping <- grNonOverlapping[which(mcols(grNonOverlapping)$score > scorePerMillion),]
  
  #--- get score-per-million for each sample
  scores <- lapply(seq_along(grList), function(x){
    columnOverlaps(query = grNonOverlapping,
                   subject = grList[[x]],
                   colname = "score",
                   decreasing = TRUE)
  }) %>% data.frame()
  colnames(scores) <- names(grList)
  
  #--- create summarized experiment
  se <- SummarizedExperiment(assays = SimpleList(scores = as.matrix(scores)), rowRanges = grNonOverlapping)
  
  #free up some memory
  remove(grList,grNonOverlapping)
  
  #--- apply selection rules
  n = ncol(se)
  minSamples <- eval(parse(text=selectionRules))
  if(n >= minSamples){
    message("Applying selection criteria: minimum of ", minSamples," samples with a score-per-million of ", scorePerMillion)
    keep <- which(rowSums(assay(se) > scorePerMillion) >= minSamples)
    message("Final Peak Set ", length(keep))
    se <- se[keep,]
  }else{
    warning("group with less than min samples required by selectionRules (most likely 1 summit for group)")
  }
  
  #--- re-normalize score-per-million
  mcols(se)$score <- edgeR::cpm(mcols(se)$score) #Re-scale for final comparison #10^6*mcols(se)$score/sum(mcols(se)$score)
  mcols(se)$name <- paste0(group_name,"_",seq_len(nrow(se))) #add short name
  
  
  #Create GRanges List for all Groups
  groupGRlist <- rowRanges(se)
  return(groupGRlist)
}


columnOverlaps <- function(query, subject, colname = "score", decreasing = TRUE){
  #First get overlaps
  o <- findOverlaps(query, subject) %>% data.frame()
  #Then append information
  o$col <- mcols(subject)[[colname]][o[,2]]
  #Order it by the factor to rank
  o <- o[order(o$col, decreasing = decreasing),]
  #Deduplicate
  o <- o[!duplicated(o$queryHits),]
  #Initialize
  val <- rep(0, length(query))
  #Fill Values
  val[o[,1]] <- o$col
  return(val)
}

convergeClusterGRanges <- function(gr, by = "score", decreasing = TRUE, verbose = FALSE){
  stopifnot(by %in% colnames(mcols(gr)))
  i = 0
  gr_initial <- gr #initial gr i
  if(verbose){
    message("Converging", appendLF = FALSE)
  }
  while(length(gr_initial) > 0){
    if(verbose){
      message(".", appendLF = FALSE)
    }
    i = i + 1
    gr_clustered <- clusterGRanges(gr = gr_initial, filter = TRUE, by = by, decreasing = decreasing, verbose = verbose) #initial gr i
    gr_initial <- subsetByOverlaps(gr_initial ,gr_clustered, invert=TRUE) #blacklist called cluster
    if(i == 1){ #if i=1 then set gr_all to clustered
      gr_all <- gr_clustered
    }else{
      gr_all <- c(gr_all, gr_clustered)
    }
  }
  return(gr_all)
}

clusterGRanges <- function(gr, filter = TRUE, by = "score", decreasing = TRUE, verbose = TRUE){
  #reduce first
  gr <- sort(sortSeqlevels(gr))
  r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
  o <- findOverlaps(gr,r)
  mcols(gr)$cluster <- subjectHits(o)
  if(verbose){
    message(sprintf("found %s overlaps...",length(gr) - max(subjectHits(o))))
  }
  #filter by
  if(filter){
    if(any(toupper(colnames(mcols(gr))) %in% toupper(by))){
      if(verbose){
        message(sprintf("filtering overlaps by %s...",by))
      }
      gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
      grn <- gr[!duplicated(mcols(gr)$cluster),]
      gr <- sort(sortSeqlevels(grn))
    }else{
      if(verbose){
        message(sprintf("filtering by order..."))
      }
      gr <- gr[!duplicated(mcols(gr)$cluster),]
    }
    mcols(gr)$cluster <- NULL
  }
  
  return(gr)
}

getChromSizes <- function (BSgenome, filter = TRUE) {
  stopifnot(!is.null(BSgenome))
  if (inherits(BSgenome, "BSgenome")) {
    BSgenome <- BSgenome
  }
  else if (is.character(BSgenome)) {
    BSgenome <- BSgenome::getBSgenome(BSgenome, masked = FALSE)
  }
  else {
    stop("Error cannot validate BSgenome options are a valid BSgenome or character for getBSgenome")
  }
  grCS <- GRanges(seqlevels(BSgenome), IRanges(1, seqlengths(BSgenome)))
  if (filter) {
    grCS <- keepFilteredChromosomes(grCS)
  }
  seqlengths(grCS) <- end(grCS)
  return(grCS)
}

keepFilteredChromosomes <- function (x, remove = c("chrM"), underscore = TRUE, standard = TRUE, pruning.mode = "coarse") {
  if (standard) {
    x <- GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = pruning.mode)
  }
  seq_names <- seqlevels(x)
  chr_remove <- c()
  if (underscore) {
    chr_remove <- c(chr_remove, which(grepl("_", seq_names)))
  }
  chr_remove <- c(chr_remove, which(seq_names %in% remove))
  if (length(chr_remove) > 0) {
    chr_keep <- seq_names[-chr_remove]
  }
  else {
    chr_keep <- seq_names
  }
  seqlevels(x, pruning.mode = pruning.mode) <- chr_keep
  return(x)
}

removeCliffedPeaks <- function(gr, BSgenome) {
  chromSizes <- getChromSizes(BSgenome = BSgenome)
  return(subsetByOverlaps(gr, chromSizes, type = "within"))
}

removeSpecificChromosome <- function(gr, remove = "chrY"){
  seq_names <- seqlevels(gr)
  if(!(remove %in% seq_names)) {
    print(paste0("Warning - ",remove," does not exist in provided GRange"))
    return(gr)
  }
  else {
    gr <- gr[which(!(as.vector(seqnames(gr)) == remove))]
    return(gr)
  }
}

removeNbasePeaks <- function(gr, BSgenome, percent = 0){
  BSgenome <- BSgenome::getBSgenome(BSgenome, masked = FALSE)
  seqs <- getSeq(BSgenome, gr)
  nuc_freqs <- as.data.frame(letterFrequency(seqs, c("A","C","G","T","N")))
  return(gr[which(nuc_freqs$N <= percent)])
}

gather_stats <- function(summary_stats_file , annotation_path , TSS_thr = 10){
  
  stats_files = list.files(
    summary_stats_file,
    pattern = "stats.tsv",
    recursive = T,
    full.names = T
  )
  samples_annot = read.table(annotation_path, header = T, sep="\t")
  
  summary_df = data.frame(sample=vector(), NRF=vector(), PBC1=vector(), PBC2=vector(), TSS_score=vector(), FRiP=vector(), Peak_count=vector(), Alignment_rate=vector(), Alignment_rate_chrM2x= vector())
  for(i in stats_files){
    sample = unlist(strsplit(i, "/"))
    sample = sample[length(sample)-1]
    stats_df = read.table(i, header = F, sep = "\t")
    summary_df = rbind(summary_df,data.frame(sample=sample,
                                             NRF=as.numeric(stats_df$V2[which(stats_df$V1 == "NRF")]),
                                             PBC1=as.numeric(stats_df$V2[which(stats_df$V1 == "PBC1")]),
                                             PBC2=as.numeric(stats_df$V2[which(stats_df$V1 == "PBC2")]),
                                             TSS_score=as.numeric(stats_df$V2[which(stats_df$V1 == "TSS_score")]),
                                             FRiP=as.numeric(stats_df$V2[which(stats_df$V1 == "FRiP")]),
                                             Peak_count=as.numeric(stats_df$V2[which(stats_df$V1 == "Peak_count")]),
                                             Alignment_rate=as.numeric(stats_df$V2[which(stats_df$V1 == "Alignment_rate")]),
                                             Alignment_rate_chrM2x=as.numeric(stats_df$V2[which(stats_df$V1 == "Alignment_rate_mouse_chrM2x")])))
  }
  
  merged_data = merge(summary_df, samples_annot, by.x = "sample", by.y = "SRA_ID")
  
  merged_data$to_keep = F
  merged_data$to_keep[which(merged_data$TSS_score >= TSS_thr)] = T
  
  return(merged_data)
  
} 
get_stats_figure <- function(stats_df, group_var = "protocol"){
  
  pNRF <- ggplot2::ggplot(data = stats_df, aes(x = factor(sample, level=stats_df$sample), y = NRF, fill = get(group_var))) +
    geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") +
    theme(legend.title=element_blank())
  
  pPBC1 <- ggplot2::ggplot(data = stats_df, aes(x = factor(sample, level=stats_df$sample), y = PBC1, fill = get(group_var))) +
    geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")+
    theme(legend.title=element_blank())
  
  PBC2 <- ggplot2::ggplot(data = stats_df, aes(x = factor(sample, level=stats_df$sample), y = PBC2, fill = get(group_var))) +
    geom_bar(stat="identity") + theme_bw()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")+
    theme(legend.title=element_blank())
  
  pTSS_score <- ggplot2::ggplot(data = stats_df,  aes(x = factor(sample, level=stats_df$sample), y = TSS_score, fill = get(group_var))) +
    geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")+
    theme(legend.title=element_blank())
  
  pFRiP <- ggplot2::ggplot(data = stats_df, aes(x = factor(sample, level=stats_df$sample), y = FRiP, fill = get(group_var))) +
    geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")+
    theme(legend.title=element_blank())
  
  pPeak_count <- ggplot2::ggplot(data = stats_df, aes(x = factor(sample, level=stats_df$sample), y = Peak_count, fill = get(group_var))) +
    geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")+
    theme(legend.title=element_blank())
  
  pAlignment_rate <- ggplot2::ggplot(data = stats_df, aes(x = factor(sample, level=stats_df$sample), y = Alignment_rate, fill = get(group_var))) +
    geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")+
    theme(legend.title=element_blank())
  
  pAlignment_rate_chrM2x <- ggplot2::ggplot(data = stats_df, aes(x = factor(sample, level=stats_df$sample), y = Alignment_rate_chrM2x, fill = get(group_var))) +
    geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")+
    theme(legend.title=element_blank())
  
  p <- gridExtra::grid.arrange(pNRF, pPBC1, PBC2, pFRiP, pPeak_count, pTSS_score, pAlignment_rate, pAlignment_rate_chrM2x, nrow = 2)
  return(p)
} 
