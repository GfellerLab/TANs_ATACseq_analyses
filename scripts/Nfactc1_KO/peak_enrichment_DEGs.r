library(biomaRt)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicFeatures)
library(dplyr)
library(ggplot2)
library(tidyr)
set.seed(2025)

# Functions ---------------------------------------------------------------

get_tss_info <- function(genes, mart_obj, window_size){
  tss_info <- getBM(attributes = c('external_gene_name', 'chromosome_name', 'transcription_start_site'),
                    filters = 'external_gene_name',
                    values = genes,
                    mart = mart_obj,useCache = F)
  
  tss_info <- tss_info[!duplicated(tss_info$external_gene_name), ]
  tss_info <- tss_info[which(tss_info$chromosome_name != "MT"), ]
  tss_info$start <- tss_info$transcription_start_site - window_size
  tss_info$start <- ifelse(tss_info$start < 0, 0, tss_info$start)
  tss_info$end <- tss_info$transcription_start_site + window_size
  tss_info$chromosome_name <- paste0("chr", tss_info$chromosome_name)
  gr_ranges <- makeGRangesFromDataFrame(tss_info, keep.extra.columns = T)
  
  return(gr_ranges)
}


# Parameters --------------------------------------------------------------

options(scipen = 999)
rna_deg <- as.character(commandArgs(TRUE)[1])
output_path <- as.character(commandArgs(TRUE)[2])
window_size <- as.numeric(commandArgs(TRUE)[3])
DAR_path <- as.character(commandArgs(TRUE)[4])
cluster_peaks_path <- as.character(commandArgs(TRUE)[5])

dir.create(output_path, recursive = T)
# Load mouse genome data --------------------------------------------------
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# retrieve all genes info
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
all_genes <- genes(txdb)
gene_ids <- mcols(all_genes)$gene_id
gene_symbols <- mapIds(org.Mm.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
mcols(all_genes)$symbol <- gene_symbols

print(window_size)
tss_gr <- get_tss_info(all_genes$symbol, mart_obj = mouse, window_size = window_size)


# Run test using the DARs  ------------------------------------------------

summary_tests <- data.frame(DEG_direction = vector(),
                            DAR_direction = vector(),
                            nb_genes = vector(),
                            window_size = vector(),
                            RNA_genes = vector(),
                            Random_genes = vector(),
                            pvalue = vector())

for(DEG_direction in c("up", "down")){
  print(paste0("DEG", DEG_direction))
  # Load RNAseq diff genes  -------------------------------
  diff_genes <- as.data.frame(readxl::read_excel(rna_deg))
  if(DEG_direction == "up"){
    diff_genes <- diff_genes[which(diff_genes$log2FoldChange>0 & diff_genes$padj < 0.05), ]
  }else{
    diff_genes <- diff_genes[which(diff_genes$log2FoldChange<0 & diff_genes$padj < 0.05), ]
  }
  
  diff_genes <- diff_genes[order(abs(diff_genes$log2FoldChange), decreasing = T), ]
  
  for(DAR_direction in c("up", "down")){
    print(paste0("DAR", DAR_direction))
    
    peaks <- read.table(DAR_path, header = T, sep = "\t",quote = "")
    peaks <- peaks[peaks$peak_cluster == paste0("C", DAR_direction),]
    peaks$peak <- paste0(peaks$seqnames, ":", peaks$start, "-", peaks$end)
   
    
    peaks_df <- as.data.frame(Signac::StringToGRanges(peaks$peak, sep =c(":", "-")))
    peaks_df$peak_id <- peaks$peak 
    colnames(peaks_df) <- c("chr", "start", "end", "length","strand" ,"peak_id")
    peaks <- makeGRangesFromDataFrame(peaks_df, keep.extra.columns = T)
    
    # Match genes to peaks -------------------------------
    nb_genes <- nrow(diff_genes)
    tss_subset_gr <- tss_gr[tss_gr$external_gene_name %in% diff_genes$gene, ]
    
    overlapping_peaks <- findOverlaps(peaks,  tss_subset_gr)
    overlapping_peaks
    
    nb_overlapping_peaks_random = vector()
    for(iter in 1:1000){
      # random genes:
      random_genes <- sample(all_genes$symbol, nb_genes, replace = F) 
      tss_random_gr <- tss_gr[tss_gr$external_gene_name %in% random_genes, ]
      
      overlapping_peaks_random <- findOverlaps(peaks,  tss_random_gr)
      overlapping_peaks_random
      nb_overlapping_peaks_random <- c(nb_overlapping_peaks_random, length(overlapping_peaks_random))
    }
    
    observed_count <- length(overlapping_peaks)
    rand_counts <- nb_overlapping_peaks_random
    
    z_score <- (observed_count - mean(rand_counts)) / sd(rand_counts)
    p_value_z <- pnorm(z_score, lower.tail = FALSE)
    print(p_value_z)
    
    summary_tests <- rbind(summary_tests, data.frame(DEG_direction = DEG_direction,
                                                     DAR_direction = DAR_direction,
                                                     nb_genes = nb_genes,
                                                     window_size = window_size,
                                                     RNA_genes = length(overlapping_peaks),
                                                     Random_genes = mean(nb_overlapping_peaks_random),
                                                     pvalue = p_value_z))
    
    
    
    
  }
  
}


pdf(paste0(output_path, "/NfatKODARs_DEGs_overlap_window", window_size, ".pdf"), w = 11, h = 5)

summary_subset <- summary_tests[which(summary_tests$window_size == window_size ), ]
summary_subset$nb_genes <- as.character(summary_subset$nb_genes)
summary_subset <- pivot_longer(summary_subset,
                               cols = c(RNA_genes, Random_genes),
                               names_to = "Analysis",
                               values_to =  "Nb_overlaping_peaks")
summary_subset$Analysis <- paste0(summary_subset$Analysis, "_", summary_subset$DEG_direction)

write.table(summary_subset, file = paste0(output_path, "/NfatKODARs_DEGs_overlap_window", window_size, ".txt"), col.names = T, row.names = F, quote = F, sep = "\t")

diff_genes <- as.data.frame(readxl::read_excel(rna_deg))
diff_genes <- diff_genes[which(diff_genes$log2FoldChange>0 & diff_genes$padj < 0.05), ]
nb_induced <- nrow(diff_genes)
diff_genes <- as.data.frame(readxl::read_excel(rna_deg))
diff_genes <- diff_genes[which(diff_genes$log2FoldChange<0 & diff_genes$padj < 0.05), ]
nb_repressed <- nrow(diff_genes)


summary_subset$DEG_direction <- plyr::revalue(summary_subset$DEG_direction, c("up" = paste0("Induced genes \n(n = ",nb_induced,")") , "down" = paste0("Repressed genes \n(n = ",nb_repressed,")")))
summary_subset$DEG_direction <- factor(summary_subset$DEG_direction , levels = c(paste0("Induced genes \n(n = ",nb_induced,")"), paste0("Repressed genes \n(n = ",nb_repressed,")")))
summary_subset$DAR_direction <- plyr::revalue(summary_subset$DAR_direction, c("up" = "Opening" , "down" = "Closing"))
summary_subset$DAR_direction <- factor(summary_subset$DAR_direction , levels = c("Opening", "Closing" ))


barplot_colors <- c("RNA_genes_up" = "#FF6666", "RNA_genes_down" = "#FF6666", "Rnadom_genes_up" = "gray", "Rnadom_genes_down" = "gray")
ggplot(summary_subset, aes(x = DAR_direction, y = Nb_overlaping_peaks, fill = Analysis)) +
  geom_bar(position = position_dodge(width = 0.9), stat = "identity", show.legend = F) +
  geom_label(inherit.aes = FALSE, data = . %>%
               dplyr::select(DEG_direction, DAR_direction, pvalue) %>%
               distinct(), aes(label = paste0( "p = ", sprintf("%.2e", pvalue)),
                               x = DAR_direction, 
                               y = summary_subset$Nb_overlaping_peaks[summary_subset$Analysis %in% c("RNA_genes_up", "RNA_genes_down")]+50), size = 4) +
  facet_wrap(~ DEG_direction) +
  scale_fill_manual(values = barplot_colors) + 
  theme_bw() + xlab("") + ylab("Number of peaks") +
  theme(text = element_text(size=17), legend.title = element_blank()) 
dev.off()


# Run test using the cluster peaks ----------------------------------------

summary_tests <- data.frame(DEG_direction = vector(),
                            Cluster = vector(),
                            nb_genes = vector(),
                            window_size = vector(),
                            RNA_genes = vector(),
                            Random_genes = vector(),
                            pvalue = vector())

for(DEG_direction in c("up", "down")){
  
  # Load RNAseq diff genes  -------------------------------
  diff_genes <- as.data.frame(readxl::read_excel(rna_deg))
  if(DEG_direction == "up"){
    diff_genes <- diff_genes[which(diff_genes$log2FoldChange>0 & diff_genes$padj < 0.05), ]
  }else{
    diff_genes <- diff_genes[which(diff_genes$log2FoldChange<0 & diff_genes$padj < 0.05), ]
  }
  
  diff_genes <- diff_genes[order(abs(diff_genes$log2FoldChange), decreasing = T), ]
  cluster_peaks <- as.data.frame(data.table::fread(cluster_peaks_path))
  for(cluster in paste0("C", 1:length(unique(cluster_peaks$peak_cluster)))){
    
    peaks <- cluster_peaks
    table(peaks$peak_cluster)
    
    peaks$peak <- paste0(peaks$seqnames, ":", peaks$start, "-", peaks$end)
    peaks <- peaks[peaks$peak_cluster %in% cluster, ]
    
    peaks_df <- as.data.frame(Signac::StringToGRanges(peaks$peak, sep =c(":", "-")))
    peaks_df$peak_id <- peaks$peak 
    colnames(peaks_df) <- c("chr", "start", "end", "length","strand" ,"peak_id")
    peaks <- makeGRangesFromDataFrame(peaks_df, keep.extra.columns = T)
    
    # Match genes to peaks -------------------------------
    nb_genes_vector <- nrow(diff_genes)
    
    for(nb_genes in nb_genes_vector){
      print(paste0("nb_genes: ", nb_genes))
      print(paste0("window_size: ", window_size))
      
      genes <- diff_genes$gene[1:nb_genes]
      tss_subset_gr <- tss_gr[tss_gr$external_gene_name %in% genes, ]
      
      overlapping_peaks <- findOverlaps(peaks,  tss_subset_gr)
      
      nb_overlapping_peaks_random = vector()
      for(iter in 1:1000){
        # random genes:
        random_genes <- sample(all_genes$symbol, nb_genes, replace = F) 
        tss_random_gr <- tss_gr[tss_gr$external_gene_name %in% random_genes, ]
        
        overlapping_peaks_random <- findOverlaps(peaks,  tss_random_gr)
        nb_overlapping_peaks_random <- c(nb_overlapping_peaks_random, length(overlapping_peaks_random))
      }
      
      observed_count <- length(overlapping_peaks)
      rand_counts <- nb_overlapping_peaks_random
      z_score <- (observed_count - mean(rand_counts)) / sd(rand_counts)
      p_value_z <- pnorm(z_score, lower.tail = FALSE)
      print(p_value_z)
      
      summary_tests <- rbind(summary_tests, data.frame(DEG_direction = DEG_direction,
                                                       Cluster = cluster,
                                                       nb_genes = nb_genes,
                                                       window_size = window_size,
                                                       RNA_genes = length(overlapping_peaks),
                                                       Random_genes = mean(nb_overlapping_peaks_random),
                                                       pvalue = p_value_z))
      
      
      
    }
    
    
  }
  
}


pdf(paste0(output_path, "/KPclusterPeaks_DEGs_overlap_window", window_size, ".pdf"), w = 11, h = 5)

summary_subset <- summary_tests[which(summary_tests$window_size == window_size ), ]
summary_subset <- pivot_longer(summary_subset,
                               cols = c(RNA_genes, Random_genes),
                               names_to = "Analysis",
                               values_to =  "Nb_overlaping_peaks")
summary_subset$Analysis <- paste0(summary_subset$Analysis, "_", summary_subset$DEG_direction)

write.table(summary_subset, file = paste0(output_path, "/KPclusterPeaks_DEGs_overlap_window", window_size, ".txt"), col.names = T, row.names = F, quote = F, sep = "\t")


diff_genes <- as.data.frame(readxl::read_excel(rna_deg))
diff_genes <- diff_genes[which(diff_genes$log2FoldChange>0 & diff_genes$padj < 0.05), ]
nb_induced <- nrow(diff_genes)
diff_genes <- as.data.frame(readxl::read_excel(rna_deg))
diff_genes <- diff_genes[which(diff_genes$log2FoldChange<0 & diff_genes$padj < 0.05), ]
nb_repressed <- nrow(diff_genes)

summary_subset$DEG_direction <- plyr::revalue(summary_subset$DEG_direction, c("up" = paste0("Induced genes \n(n = ",nb_induced,")") , "down" = paste0("Repressed genes \n(n = ",nb_repressed,")")))
summary_subset$DEG_direction <- factor(summary_subset$DEG_direction , levels = c(paste0("Induced genes \n(n = ",nb_induced,")"), paste0("Repressed genes \n(n = ",nb_repressed,")")))

summary_subset$pvalue_star <- ifelse(summary_subset$pvalue < 0.001, "***",
                                     ifelse(summary_subset$pvalue < 0.01, "**",
                                            ifelse(summary_subset$pvalue < 0.05, "*", "ns")))

barplot_colors <- c("RNA_genes_up" = "#FF6666", "RNA_genes_down" = "#FF6666", "Rnadom_genes_up" = "gray", "Rnadom_genes_down" = "gray")
ggplot(summary_subset, aes(x = Cluster, y = Nb_overlaping_peaks, fill = Analysis)) +
  geom_bar(position = position_dodge(width = 0.9), stat = "identity", show.legend = F) +
  geom_label(inherit.aes = FALSE, data = . %>%
               dplyr::select(DEG_direction, Cluster, pvalue) %>%
               distinct(), 
             aes(label = paste0( "p=", sprintf("%.1e", pvalue)),
                 x = Cluster, 
                 y = summary_subset$Nb_overlaping_peaks[summary_subset$Analysis %in% c("RNA_genes_up", "RNA_genes_down")]+50), 
             size = 4) +
  facet_wrap(~ DEG_direction) +
  scale_fill_manual(values = barplot_colors) + 
  theme_bw() + xlab("") + ylab("Number of peaks") +
  theme(text = element_text(size=17), legend.title = element_blank()) 
dev.off()
