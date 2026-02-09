
# Libraries ---------------------------------------------------------------

library(GenomicRanges)
library(ggplot2)
library(ComplexHeatmap)


# Parameters --------------------------------------------------------------
counts_path <- as.character(commandArgs(TRUE)[1])
metadata_path <- as.character(commandArgs(TRUE)[2])
dar_res_path <- as.character(commandArgs(TRUE)[3])
fasta_file <- as.character(commandArgs(TRUE)[4])
output_path <- as.character(commandArgs(TRUE)[5])


dir.create(output_path, recursive = T)

descriminative_colors = c("WT"="darkgray",
                          "Nfatc1KO"="#EE6677")

# Read the raw counts -----------------------------------------------------

counts_data <- read.table(counts_path, header = T, sep = "\t", row.names = 1)

# Load samples metadata ---------------------------------------------------

metadata <- read.table(metadata_path, header = T, sep = "\t")
rownames(metadata) <- metadata$sample_name
metadata <- metadata[colnames(counts_data),] 

condition <- metadata$groups


# Get diff_peaks ----------------------------------------------------------

diff_res <- read.table(dar_res_path, header = T, sep = "\t", quote = "")
colnames(diff_res)[1] = "IDs"
summary(diff_res$IDs %in% rownames(counts_data))

pval <- diff_res[, grepl("adj..p.value", colnames(diff_res)), drop = F]
min_pval <- apply(pval, 1, min)
logFC <- abs(diff_res[, grepl("Log2.Fold.Change", colnames(diff_res)), drop = F])

max_logFC <- apply(logFC, 1, max)

diff_peaks <- which(min_pval < 0.05 )

length(unique(diff_res$IDs[diff_peaks]))

gc_norm <- function(input_matrix, grouping_var, fa_file, by_group = F){
  ff <- Rsamtools::FaFile(fa_file)
  gr <- Signac::StringToGRanges(rownames(input_matrix), sep = c(":", "-"))
  gr = GRanges(seqnames=seqnames(gr), ranges=IRanges(start(gr), end(gr)), strand="*", mcols=data.frame(peakID=rownames(input_matrix)))
  
  peakSeqs <- Biostrings::getSeq(x=ff, gr)
  gcContentPeaks = Biostrings::letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]
  
  if(by_group){
    dataNorm <- vector()
    for(group in unique(grouping_var)){
      print(group)
      input_matrix_subset <- input_matrix[, which(grouping_var == group)] 
      dataWithin <- EDASeq::withinLaneNormalization(input_matrix_subset, y = gcContentPeaks, num.bins = 20, which = "full")#adjusting for gene length and gc content  
      dataNorm_subset <- EDASeq::betweenLaneNormalization(dataWithin, which = "full")# adjusting for sequencing depth 
      
      dataNorm <- cbind(dataNorm, dataNorm_subset)
    } 
  } else{
    dataWithin = EDASeq::withinLaneNormalization(input_matrix, y = gcContentPeaks, num.bins = 20, which = "full")#adjusting for gene length and gc content
    dataNorm = EDASeq::betweenLaneNormalization(dataWithin, which = "full")# adjusting for sequencing depth
  } 
  
  return(dataNorm)
} 
dataNorm <- gc_norm(input_matrix = as.matrix(counts_data), 
                    grouping_var = metadata$groups, 
                    fa_file = fasta_file, 
                    by_group = F)
norm_data <- log(dataNorm+1)
norm_data_subset <- norm_data[diff_res$IDs[diff_peaks] , ]



# Run PCA and UMAP --------------------------------------------------------
cluster_assignment <- ifelse(diff_res[diff_peaks,]$Nfatc1KO.vs..WT.Log2.Fold.Change > 0 , "down", "up")


table(cluster_assignment)
cluster_assignment <- as.character(cluster_assignment)

# Plot clustering ---------------------------------------------------------



cl_colors <- c("#4DAF4A", "#984EA3")
names(cl_colors) <- paste0("C", c("up", "down"))

groups_factor <- factor(metadata$groups, levels = c("Nfatc1KO", "WT"))
sorted_indices <- order(groups_factor)
col_an = HeatmapAnnotation(Condition = metadata$groups[sorted_indices], 
                           annotation_name_rot = 45,
                           col = list(Condition = descriminative_colors
                           )
)
ha = rowAnnotation(Peak_cluster = paste0("C", cluster_assignment[order(cluster_assignment)] ),
                   col = list(Peak_cluster = cl_colors))


pdf(paste0(output_path,"peak_clustering_heatmap.pdf"), w= 10, h=15)
Heatmap(norm_data_subset[order(cluster_assignment), sorted_indices] , show_row_names = F,
        left_annotation = ha, top_annotation = col_an,
        show_row_dend = F, cluster_rows = F, cluster_columns = F)

Heatmap(t(scale(t(norm_data_subset[order(cluster_assignment), sorted_indices]))) , show_row_names = F,
        left_annotation = ha, top_annotation = col_an,
        show_row_dend = F, cluster_rows = F, cluster_columns = F)
dev.off()

# Write list of cell-type specific peaks with annotations
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

for(i in unique(cluster_assignment)){
  print(paste0("C", i))
  cluster_peaks <- rownames(norm_data_subset)[which(cluster_assignment == i)]
  peaks_bed <- as.data.frame( Signac::StringToGRanges(cluster_peaks, sep = c(":", "-")) )
  peaks_bed$peak_id <- cluster_peaks
  print(dim(peaks_bed))
  dir.create(paste0(output_path, paste0("C", i)), recursive = T)
  write.table(peaks_bed[,c(1:3,6,4,5)], file = paste0(output_path, paste0("C", i), "/cluster_peaks.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
  
}

anno_all <- vector()
for(i in unique(cluster_assignment)){
  print(paste0("C", i))
  cluster_peaks <- rownames(norm_data_subset)[which(cluster_assignment == i)]
  peak_gr <- Signac::StringToGRanges(cluster_peaks, sep = c(":", "-")) 
  
  anno = ChIPseeker::annotatePeak(
    peak_gr,
    TxDb = txdb,
    # tssRegion=c(-1000, 1000),
    addFlankGeneInfo = TRUE,
    overlap = "all", #any nearest gene is reported regardless of the overlap with the TSS
    annoDb = "org.Mm.eg.db"
  )
  
  # Save annotation
  anno_df = as.data.frame(anno)
  head(anno_df)
  
  anno_df$annotation_type = anno_df$annotation
  anno_df$annotation_type[grep("Intron",anno_df$annotation_type)]="Intron"
  anno_df$annotation_type[grep("Exon",anno_df$annotation_type)]="Exon"
  anno_df$annotation_type[grep("Promoter",anno_df$annotation_type)]="Promoter"
  anno_df$annotation_type[grep("Downstream",anno_df$annotation_type)]="Downstream"
  
  table(anno_df$annotation_type)
  anno_df$peak_cluster <- paste0("C", i)
  anno_all <- rbind(anno_all, anno_df)
}

write.table(anno_all, file = paste0(output_path, "cluster_peaks_annotation.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

library(dplyr)
proportions_table <- anno_all %>%
  group_by(peak_cluster, annotation_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  # Calculate proportions within each cluster
  group_by(peak_cluster) %>%
  mutate(proportion = count / sum(count))

# View the table to ensure it looks correct
print(proportions_table)

pdf(paste0(output_path,"peak_clustering_annotations.pdf"))

ggplot(proportions_table, aes(x = peak_cluster, y = proportion, fill = annotation_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "Proportion of Annotation Types by Peak Cluster",
       x = "Peak Cluster",
       y = "Proportion (%)",
       fill = "Annotation Type") + scale_fill_manual(values = c("Promoter"="indianred2","Downstream"="palegoldenrod", "3' UTR"="deepskyblue3", "Exon" = "green4", 
                                                                "Intron"="lightblue", "Distal Intergenic"="darkorange1", "5' UTR"="plum3", "thistle", "saddlebrown"))
dev.off()
