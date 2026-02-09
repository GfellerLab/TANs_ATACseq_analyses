
library(ComplexHeatmap)
library(circlize)
library(GenomicRanges)
library(ggplot2)
library(dplyr)

peak_clustering <- as.character(commandArgs(TRUE)[1])
peak_clustering_pairwise <- as.character(commandArgs(TRUE)[2])
DAR_res <- as.character(commandArgs(TRUE)[3])
metadata_path <- as.character(commandArgs(TRUE)[4])
counts_path <- as.character(commandArgs(TRUE)[5])
fasta_file <- as.character(commandArgs(TRUE)[6])
homer_path <- as.character(commandArgs(TRUE)[7])
motif_annotation_path <- as.character(commandArgs(TRUE)[8])
output_path <- as.character(commandArgs(TRUE)[9])


# Peak annotations --------------------------------------------------------
anno_peak_clustering <- read.table(peak_clustering, header = T, sep = "\t",quote = "")

diff_res <- read.table(DAR_res, header = T, sep = "\t", quote = "")
colnames(diff_res)[1] = "IDs"

pval <- diff_res[, grepl("adj..p.value", colnames(diff_res)), drop = F]
min_pval <- apply(pval, 1, min)
nonDARs <- which(min_pval > 0.05 )

peak_gr <- Signac::StringToGRanges(diff_res$IDs[nonDARs], sep = c(":", "-")) 

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

anno = ChIPseeker::annotatePeak(
  peak_gr,
  TxDb = txdb,
  addFlankGeneInfo = TRUE,
  overlap = "all", 
  annoDb = "org.Hs.eg.db"
)

anno_nonDARs = as.data.frame(anno)
head(anno_nonDARs)

anno_nonDARs$annotation_type = anno_nonDARs$annotation
anno_nonDARs$annotation_type[grep("Intron",anno_nonDARs$annotation_type)]="Intron"
anno_nonDARs$annotation_type[grep("Exon",anno_nonDARs$annotation_type)]="Exon"
anno_nonDARs$annotation_type[grep("Promoter",anno_nonDARs$annotation_type)]="Promoter"
anno_nonDARs$annotation_type[grep("Downstream",anno_nonDARs$annotation_type)]="Downstream"

table(anno_nonDARs$annotation_type)
anno_nonDARs$peak_cluster <- "non-DARs"

anno_peak_clustering$peak_cluster <- ifelse(anno_peak_clustering$peak_cluster == "Cup", "Opening", "Closing")
anno_all <- rbind(anno_peak_clustering, anno_nonDARs)

table(anno_all$peak_cluster)

proportions_table <- anno_all %>%
  group_by(peak_cluster, annotation_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  # Calculate proportions within each cluster
  group_by(peak_cluster) %>%
  mutate(proportion = count / sum(count))


# View the table to ensure it looks correct
print(proportions_table)
proportions_table$peak_cluster <- factor(proportions_table$peak_cluster, levels = c("Opening", "Closing", "non-DARs"))
pdf(paste0(output_path,"peak_clustering_annotations_withNonDARs.pdf"))

ggplot(proportions_table, aes(x = peak_cluster, y = proportion, fill = annotation_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "Proportion of Annotation Types by Peak Cluster",
       x = "Peak Cluster",
       y = "Proportion (%)",
       fill = "Annotation Type") + scale_fill_manual(values = c("Promoter"="indianred2","Downstream"="palegoldenrod", "3' UTR"="deepskyblue3", "Exon" = "green4", 
                                                                "Intron"="lightblue", "Distal Intergenic"="darkorange1", "5' UTR"="plum3", "thistle", "saddlebrown")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


# Peak cluster heatmap + summarized heatmaps ------------------------------


descriminative_colors = c("Blood"="#3288bdff",
                          "Tumor"="#d90017ff",
                          "Adjacent_lung"="#ffd92fff",
                          "P1" = "#f39800ff", "P2" = "#009944ff", "P3" = "#0068b7ff", "P4" = "#8f76d6ff")

# Read the raw counts 

counts_data <- read.table(counts_path, header = T, sep = "\t", row.names = 1)
colnames(counts_data) <- gsub("[.]","-", colnames(counts_data))

# Load samples metadata 

metadata <- read.table(metadata_path, header = T, sep = "\t")
rownames(metadata) <- metadata$replicate
metadata <- metadata[colnames(counts_data),]  
metadata$patient2 <- plyr::revalue(metadata$patient, replace = c("G-LC-01" = "P1", "G-LC-03" = "P2","G-LC-05" = "P3", "G-LC-08" = "P4"))

batch <- metadata$batch
condition <- metadata$groups

# normalize data
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

anno_peak_clustering <- read.table(peak_clustering_pairwise, header = T, sep = "\t",quote = "")
anno_peak_clustering$ID <- paste0(anno_peak_clustering$seqnames, ":", anno_peak_clustering$start, "-", anno_peak_clustering$end)
norm_data_subset <- norm_data[anno_peak_clustering$ID, ]

cluster_assignment <- plyr::revalue(anno_peak_clustering$peak_cluster, c("C1" = "C2", "C2" = "C4", "C3" = "C5", "C4" = "C1", "C5" = "C3" ))
cl_colors <- c("#4DAF4A", "#984EA3", "#1B9E77","#6A3D9A",
               "#D95F02", "#7570B3", "#E7298A", "#66A61E",
               "#E6AB02", "#A6761D", "black", "gray")
names(cl_colors) <- paste0("C", 1:12)

groups_factor <- factor(metadata$groups, levels = c("Blood", "Adjacent_lung", "Tumor"))
sorted_indices <- order(groups_factor)
col_an = HeatmapAnnotation(Condition = metadata$groups[sorted_indices], 
                           Sample = metadata$patient2[sorted_indices], 
                           annotation_name_rot = 45,
                           col = list(Condition = descriminative_colors, Sample = descriminative_colors
                           )
)
ha = rowAnnotation(Peak_cluster = cluster_assignment[order(cluster_assignment)] ,
                   col = list(Peak_cluster = cl_colors))

gene_names_to_add <- c("RHOG", "MBLN1", "RHOB", "FCGR2A", "CXCL8", "IL1RAP", "OSM", "IL1B", "PLAU", "EGR1", "OLR1",
                       "CCL3L3", "CXCR4", "CAMK1G", "KCNH4", "IRAK2", "NOD2", "NFKBIE", "SEMA6B")
df <- norm_data_subset[order(cluster_assignment), sorted_indices]
rownames(anno_all) <- paste0(anno_all$seqnames, ":", anno_all$start, "-", anno_all$end)


panel_fun = function(index, nm) {
  grid.rect()
  grid.text(paste(gene_names_to_add[gene_names_to_add %in% anno_all[rownames(df)[index], "SYMBOL"]], 
                  collapse = "\n"), 0.5, 0.5,gp = gpar(fontface = "italic", fontsize = 10))
}
align_to = lapply(sort(unique(cluster_assignment)), function(x) which(cluster_assignment[order(cluster_assignment)] == x))
names(align_to) <- sort(unique(cluster_assignment))

right_annotation <- rowAnnotation(foo = anno_block(
  align_to = align_to,
  panel_fun = panel_fun, 
  width = unit(3, "cm")
))


pdf(paste0(output_path,"peak_clustering_heatmap_reordered.pdf"), w=8, h=10)

Heatmap(df, show_row_names = F, show_column_names = F, #col = col_fun1,
        left_annotation = ha, right_annotation = right_annotation, top_annotation = col_an, 
        show_row_dend = F, cluster_rows = F, cluster_columns = F)


dev.off()

col_an = HeatmapAnnotation(Condition = metadata$groups[sorted_indices], 
                           annotation_name_rot = 45,
                           col = list(Condition = descriminative_colors
                           )
)
data <- as.data.frame(norm_data_subset[order(cluster_assignment), sorted_indices])
df_avg <- data %>%
  mutate(Group = factor(cluster_assignment[order(cluster_assignment)])) %>%     
  group_by(Group) %>%              
  summarise(across(everything(), mean)) %>%  
  as.data.frame()

df_avg
rownames(df_avg) <- df_avg$Group


pdf(paste0(gsub("peaks_clustering","peaks_clustering_pairwise",output_path),"peak_clustering_heatmap_averaged.pdf"), w=5, h=4)

col_fun = circlize::colorRamp2(quantile(unlist(df_avg)), c("#0000ffff", "#996efbff","white","#ff7a5aff","red"))

avg_heatmap <- Heatmap(as.matrix(df_avg[,2:ncol(df_avg)]), show_row_names = T, show_column_names = F, 
        top_annotation = col_an,row_names_side = "left", 
        row_names_gp = gpar(fontsize = 15),
        heatmap_legend_param = list(direction = "horizontal"), col = col_fun,
        show_row_dend = F, cluster_rows = F, cluster_columns = F, width = unit(3, "cm"))

print(avg_heatmap)

dev.off()


df_avg <- t(as.data.frame(t(df_avg[,2:ncol(df_avg)])) %>%
              mutate(Group = factor(metadata$groups, levels = c("Blood", "Adjacent_lung", "Tumor"))) %>%     
              group_by(Group) %>%              
              summarise(across(everything(), mean)) %>%  
              as.data.frame())

colnames(df_avg) <- df_avg[1,]
df_avg <- as.data.frame(df_avg[-1,])
df_avg <- mutate_all(df_avg, function(x) as.numeric(as.character(x)))
col_an = HeatmapAnnotation(Condition = colnames(df_avg), 
                           annotation_name_rot = 45,
                           col = list(Condition = descriminative_colors)
)

pdf(paste0(gsub("peaks_clustering","peaks_clustering_pairwise",output_path),"peak_clustering_heatmap_averaged.pdf"), w=5, h=4)

col_fun = circlize::colorRamp2(quantile(unlist(df_avg)), c("#0000ffff", "#996efbff","white","#ff7a5aff","red"))
avg_heatmap2 <- Heatmap(as.matrix(df_avg), show_row_names = T, show_column_names = F, #df_avg[,2:ncol(df_avg)]
                       top_annotation = col_an,row_names_side = "left", 
                       row_names_gp = gpar(fontsize = 15),
                       heatmap_legend_param = list(direction = "horizontal"), col = col_fun,
                       show_row_dend = F, cluster_rows = F, cluster_columns = F, width = unit(3, "cm"))

print(avg_heatmap2)

dev.off()

# Retrieve Homer results --------------------------------------------
homer_res <- read.table(homer_path, header = T)
homer_res$significance <- -log10(homer_res$pval)

homer_res$cluster <- plyr::revalue(homer_res$cluster, c("C1" = "C2", "C2" = "C4", "C3" = "C5", "C4" = "C1", "C5" = "C3" ))

motif_annotation <- read.table(motif_annotation_path, header = T, sep = "\t")
homer_res$motifID <- sapply(homer_res$motif_name, function(i) strsplit(i,"\\(")[[1]][1])
if(grepl("mouse", motif_annotation_path)){
  homer_res$gene_name <- sapply(homer_res$motifID, function(i) motif_annotation$mouse_gene_symbol[motif_annotation$name == i])
}else{
  homer_res$gene_name <- sapply(homer_res$motifID, function(i) motif_annotation$human_gene_symbol[motif_annotation$name == i])
}

signif_tfs <- unique(homer_res[homer_res$significance >=12, "gene_name"])
length(signif_tfs)
signif_tfs <- vector()
clusters_names <- unique(homer_res$cluster)
for(c in clusters_names){
  signif_tfs <- c(signif_tfs, homer_res[homer_res$significance >=12 & homer_res$cluster == c, "gene_name"])
}
signif_tfs <- unique(signif_tfs)

signif_matrix <- matrix(0, nrow = length(signif_tfs), ncol = length(clusters_names), dimnames = list(signif_tfs, clusters_names))

for(i in 1:length(signif_tfs)){
  for(j in clusters_names){
    signif_matrix[signif_tfs[i], j] <- ifelse(j %in% homer_res[homer_res$gene_name == signif_tfs[i],"cluster"],
                                              homer_res[homer_res$gene_name == signif_tfs[i] & homer_res$cluster == j,"significance"],
                                              0)
  }
}
signif_matrix <- t(signif_matrix)
signif_matrix <- signif_matrix[sort(rownames(signif_matrix)),]
# Get plot ----------------------------------------------------------------

pdf(paste0(gsub("peaks_clustering","peaks_clustering_pairwise",output_path),"homer_res_heatmap.pdf"), 
    w=8, h=3)

col_fun = colorRamp2(c(min(signif_matrix), 12, 30, max(signif_matrix)), c("white","#F4BA58","orange","#A03124"))

ht <-Heatmap(signif_matrix, name = "-log10(pval)", 
             cluster_rows = F, col = col_fun,
             rect_gp = gpar(col = "black", lwd = 1),
             column_names_rot = 45, row_names_side = "left",
             heatmap_legend_param = list(direction = "horizontal"), 
             row_names_gp = gpar(fontface = "bold", fontsize = 15),
             column_names_gp = gpar(fontsize = 15),
             show_column_names = T, column_title = NULL, cluster_columns = F,
             show_row_names = T, row_dend_gp = gpar(lwd = 0.5))

draw(avg_heatmap + ht,  merge_legend = TRUE, heatmap_legend_side = "right")

dev.off()



pdf(paste0(gsub("peaks_clustering","peaks_clustering_pairwise",output_path),"homer_res_heatmap2.pdf"), 
    w=8, h=3)

draw(avg_heatmap2 + ht,  merge_legend = TRUE, heatmap_legend_side = "right")

dev.off()

