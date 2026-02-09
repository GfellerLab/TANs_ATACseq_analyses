# Libraries ---------------------------------------------------------------
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(rstatix)
library(ggpubr)
# Functions ---------------------------------------------------------------

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

# Parameters --------------------------------------------------------------
metadata_path <- as.character(commandArgs(TRUE)[1])
counts_path <- as.character(commandArgs(TRUE)[2])
fasta_file <- as.character(commandArgs(TRUE)[3])
output_path <- as.character(commandArgs(TRUE)[4])
cluster_peaks_path <- as.character(commandArgs(TRUE)[5])

descriminative_colors <- c("Blood"="#3288bdff",
                           "Tumor" = "#d90017ff",
                           "Bone marrow" = "#ffd92fff",
                           "Spleen" = "#70c4b4ff",
                           "KP_lung_SiglecF_low" = "#ff7f00ff",
                           "IMM" = "darkgreen", "MAT" = "purple","Bone marrow_IMM" = "#4DAF4A", "Spleen_IMM" = "#984EA3", "Blood_IMM"="#E6AB02",
                           "Bone marrow_MAT" = "#1B9E77", "Spleen_MAT" = "#6A3D9A", "Blood_MAT"="#D95F02",
                           "Tumor_IMM" = "pink", "Tumor_MAT" = "red")


# Load the count matrix ------------------------------
counts <- read.table(counts_path, header = T, sep = "\t", row.names = 1)

# Load metadata ------------------------------

# Load samples metadata
metadata = read.table(metadata_path, header = T, sep = ",")
head(metadata)
summary(colnames(counts) %in% metadata$Run)


rownames(metadata) = metadata$Run
metadata = metadata[colnames(counts),] 
table(metadata$tissue)
table(metadata$cell_type)
metadata$group <- paste0(metadata$tissue, "_", metadata$cell_type)


# Run PCA -----------------------------------------------------------------
counts <- counts[which(rowSums(counts) != 0),]
NormFactor <- edgeR::calcNormFactors(object = counts, method="TMM") 
d0 <- edgeR::DGEList(counts, norm.factors = NormFactor)
drop <- which(apply(edgeR::cpm(d0), 1, max) < 5)
if(length(drop)!=0){
  d <- d0[-drop,] 
}else{
  d <- d0
}

peaks_gr = as(rownames(counts),"GRanges")
unique(peaks_gr@ranges@width) 

# Normalize data (considering GC content) ---------------------------------

dataNorm <- gc_norm(input_matrix = as.matrix(counts), 
                    grouping_var = metadata$cell_type, 
                    fa_file = paste0(fasta_file), 
                    by_group = F)
norm_data <- log(dataNorm+1)

# Keep most variable peaks 
rvdm = apply(norm_data, 1, var)
norm_data_subset = t(norm_data[order(rvdm, decreasing = T)[1:10000] ,])

# Run PCA and UMAP 
pca = prcomp(norm_data_subset,scale. = T)
print(summary(pca))
pcaData = as.data.frame(pca$x)
rownames(pcaData) = rownames(norm_data_subset)

metadata$annotation = metadata$cell_type


pdf(paste0(output_path, "PCA_PC1-3.pdf"), w=4, h=4)
variance_explained <- summary(pca)
to_plot=merge(pcaData, metadata[,c("Run","annotation","tissue")],by.x="row.names",by.y="row.names")
base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=annotation, shape= tissue)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  scale_shape_manual(values=c(16,17,18,8,15,7,9,10,12,13,14,1,2,3)) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot + theme(legend.position = "right") #c(0.7, 0.5)
scatter_only_plot 

base_plot <- ggplot(to_plot, aes(x=PC1, y=PC3, color=annotation, shape= tissue)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  scale_shape_manual(values=c(16,17,18,8,15,7,9,10,12,13,14,1,2,3)) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + ylab(paste0("PC3 (", round(variance_explained$importance[2,"PC3"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot + theme(legend.position = "right") #c(0.7, 0.5)
scatter_only_plot 

base_plot <- ggplot(to_plot, aes(x=PC2, y=PC3, color=annotation, shape= tissue)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  scale_shape_manual(values=c(16,17,18,8,15,7,9,10,12,13,14,1,2,3)) +
  theme_bw() + xlab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) + ylab(paste0("PC3 (", round(variance_explained$importance[2,"PC3"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot + theme(legend.position = "right") #c(0.7, 0.5)
scatter_only_plot 


dev.off()



# Plot clustering ---------------------------------------------------------
peak_clusters <- read.table(cluster_peaks_path, 
                            header = T, sep = "\t", quote = "")
rownames(peak_clusters) <- paste0(peak_clusters$seqnames, ":", peak_clusters$start, "-", peak_clusters$end)
summary(rownames(peak_clusters) %in% rownames(counts))
# head(rownames(peak_clusters)[!rownames(peak_clusters) %in% rownames(counts)])

cl_colors <- c("#4DAF4A", "#984EA3",  "#1B9E77","#6A3D9A",
               "#D95F02", "#7570B3", "#E7298A", "#66A61E",
               "#E6AB02", "#A6761D", "black", "gray")
names(cl_colors) <- paste0("C", 1:12)

groups_factor <- factor(metadata$group, levels = c("Bone marrow_IMM", "Bone marrow_MAT", 
                                                   "Spleen_IMM", "Spleen_MAT",
                                                   "Blood_IMM", "Blood_MAT",
                                                   "Tumor_IMM", "Tumor_MAT"))
sorted_indices <- order(groups_factor)
library(ComplexHeatmap)
col_an = HeatmapAnnotation(Condition = metadata$group[sorted_indices], 
                           Cell_type = metadata$cell_type[sorted_indices],
                           annotation_name_rot = 45,
                           col = list(Condition = descriminative_colors,
                                      Tissue = descriminative_colors,
                                      Cell_type = descriminative_colors
                           )
)
# C4_peaks <- row.names(peak_clusters)[peak_clusters$peak_cluster == "C4"]
norm_data_subset <- norm_data[row.names(peak_clusters),]

cluster_assignment = peak_clusters[rownames(norm_data_subset), "peak_cluster"]
ha = rowAnnotation(Peak_cluster = cluster_assignment[order(cluster_assignment)] ,
                   col = list(Peak_cluster = cl_colors))

pdf(paste0(output_path, "pancreas_heatmap.pdf"), w=10,h=10)

Heatmap(t(scale(t(norm_data_subset[order(cluster_assignment), sorted_indices]))) , show_row_names = F,
        left_annotation = ha, top_annotation = col_an,
        show_row_dend = F, cluster_rows = F, cluster_columns = F)

Heatmap(norm_data_subset[order(cluster_assignment), sorted_indices] , show_row_names = F,
        left_annotation = ha, top_annotation = col_an,
        show_row_dend = F, cluster_rows = F, cluster_columns = F)
dev.off()

# Generate boxplots -------------------------------------------------------

avg_data <- as.data.frame(sapply(unique(metadata$tissue), function(g) rowMeans(norm_data_subset[, metadata$tissue == g, drop = FALSE])))
avg_data$peak_cluster <- cluster_assignment

avg_df <- as.data.frame(norm_data_subset) %>%
  t() %>%
  as.data.frame() %>%
  mutate(tissue = metadata$tissue, annotation = metadata$annotation) %>%
  group_by(tissue, annotation) %>%
  summarise(across(everything(), mean)) 
avg_df$tissue <- factor(avg_df$tissue, levels = c("Bone marrow", "Spleen", "Blood", "Tumor"))

pdf(paste0(output_path, "pancreas_boxplots.pdf"), w=7, h=4)
for(i in 1:length(unique(cluster_assignment))){
  
  print(paste0("C",i))
  avg_subset <- avg_df[,c(1:2, which(colnames(avg_df) %in% rownames(norm_data_subset)[cluster_assignment == paste0("C",i)]))]

  
  df_long <- reshape2::melt(
    avg_subset,
    id.vars = c("tissue", "annotation"),
    variable.name = "Peaks",
    value.name = "Accessibility"
  )

  df_long$tissue <- as.factor(df_long$tissue)
  df_long$annotation <- as.factor(df_long$annotation)
  stats <- df_long %>%
    group_by(annotation) %>%
    pairwise_wilcox_test(Accessibility ~ tissue, p.adjust.method = "BH",paired=T) %>% #comparisons = list(c("Sell+", "Ccl3-"), c("Ccl3-", "Ccl3+"), c("Sell+", "Ccl3+"))
    ungroup() %>%
    add_xy_position(x = "tissue") 
  
  stats$subtype = stats$group1 
  stats <- stats[stats$group1 == "Tumor" | stats$group2 == "Tumor", ]
  
  stats$y.position = rep(c(max(df_long[,"Accessibility"])+0.2,
                           max(df_long[,"Accessibility"])+0.9, 
                           max(df_long[,"Accessibility"])+1.6
                           ),
                         length(unique(df_long$annotation)))
  
  stats$effsize <- sapply(1:nrow(stats), function(i){
    df_tmp <- df_long[df_long$annotation == stats$annotation[i] & df_long$tissue %in% c(stats$group1[i], stats$group2[i]),]
    df_tmp$tissue <- factor(df_tmp$tissue, levels = unique(df_tmp$tissue))
    test = rstatix::wilcox_effsize(data.frame(tissue = df_tmp[, "tissue"],
                                       Accessibility = df_tmp[, "Accessibility"],
                                       Peaks = df_tmp[, "Peaks"]),
                            formula = Accessibility ~ tissue)

    return(test$effsize)
  })
  
  stats$label_star_r <- paste0(stats$p.adj.signif, ", r=", round(stats$effsize,digits = 1))
 
  p <- ggplot(df_long, aes(x = tissue, y = Accessibility)) +
    geom_boxplot(position = position_dodge(width = 0.8), size = 0.5) +
    ylim(c(min(df_long[,"Accessibility"])-0.3, max(df_long[,"Accessibility"])+2))+
    facet_grid(~ annotation, scales = "free_x", space = "free_x")+ 
    stat_pvalue_manual(stats, label = "p.adj.signif", hide.ns = T, tip.length = 0.01, size = 4)+
    theme_bw(base_size = 14) + ylab("Gene Signature") + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
    theme(axis.text.x = element_text( size = 10)) + #angle = 25, vjust = 1, hjust=1,
    labs(title = paste0("C", i, " accessibility"), y = "Signature", x = "Study")
  print(p)

}
dev.off()

