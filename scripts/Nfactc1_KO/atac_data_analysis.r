# Libraries ---------------------------------------------------------------
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)

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
pfm_file <- as.character(commandArgs(TRUE)[5])
cluster_peaks_path <- as.character(commandArgs(TRUE)[6])
counts_KP_peaks_path <- as.character(commandArgs(TRUE)[7])


descriminative_colors <- c("WT"="darkgray",
                           "NfatKO" = "#EE6677","Nfatc1KO" = "#EE6677")

dir.create(output_path, recursive = T)

# Load the count matrix ------------------------------
counts <- read.table(counts_path, header = T, sep = "\t", row.names = 1)
dim(counts)

# Load metadata ------------------------------

# Load samples metadata
metadata = read.table(metadata_path, header = T, sep = "\t")
head(metadata)
summary(colnames(counts) %in% metadata$sample_name)

rownames(metadata) = metadata$sample_name
metadata = metadata[colnames(counts),] 


# Run PCA -----------------------------------------------------------------
counts <- counts[which(rowSums(counts) != 0),]

peaks_gr = as(rownames(counts),"GRanges")
unique(peaks_gr@ranges@width) 

# Normalize data (considering GC content) ---------------------------------

dataNorm <- gc_norm(input_matrix = as.matrix(counts), 
                    grouping_var = metadata$groups, 
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

metadata$annotation = metadata$groups

dir.create(paste0(output_path, "/PCA/"))
pdf(paste0(output_path, "/PCA/PCA_PC1-3.pdf"), w=5.5,h=4)
variance_explained <- summary(pca)
to_plot=merge(pcaData, metadata[,c("sample_name","annotation")],by.x="row.names",by.y="row.names")

base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=annotation)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  stat_ellipse(geom = "polygon",type = "norm",level = 0.7, alpha = 0.3) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot + theme(legend.position = "right") #c(0.7, 0.5)
scatter_only_plot 

base_plot <- ggplot(to_plot, aes(x=PC1, y=PC3, color=annotation)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  stat_ellipse(geom = "polygon",type = "norm",level = 0.7, alpha = 0.3) +

  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + ylab(paste0("PC3 (", round(variance_explained$importance[2,"PC3"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot + theme(legend.position = "right") #c(0.7, 0.5)
scatter_only_plot 

dev.off()


# Distances and correlation plots -------------------------------------------------------
dir.create(paste0(output_path, "/corrPlots/"))
sampleDists <- dist(t(norm_data))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- sapply(metadata$sample_name, function(i) unlist(strsplit(i,"_"))[[1]] )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

ha = rowAnnotation(Condition = metadata$groups, 
                   annotation_name_rot = 45,
                   col = list(Condition = descriminative_colors)
)
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- colnames(norm_data)
dist_plot_correction <- Heatmap(sampleDistMatrix, name = "Distance", 
                                       left_annotation = ha, col = colors, 
                                       rect_gp = gpar(col = "black", lwd = 0.5))


pdf(paste0(output_path, "/corrPlots/samples_distances.pdf"), w = 10, h = 9)
draw(dist_plot_correction,
     column_title="",
     column_title_gp=grid::gpar(fontsize=16))
dev.off()


# Compute correlations instead of distances 
sampleCor <- cor(norm_data)
sampleCorMatrix <- as.matrix(sampleCor)
rownames(sampleCorMatrix) <- sapply(metadata$sample_name, function(i) unlist(strsplit(i,"_"))[[1]] )
colnames(sampleCorMatrix) <- NULL
colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)

ha = rowAnnotation(Condition = metadata$groups,
                   annotation_name_rot = 45,
                   col = list(Condition = descriminative_colors
                   )
)
rownames(sampleCorMatrix) <- colnames(sampleCorMatrix) <- colnames(norm_data)
cor_plot_correction <- Heatmap(sampleCorMatrix, name = "Correlation", 
                                      left_annotation = ha, col = colors, show_column_names = T,
                                      rect_gp = gpar(col = "black", lwd = 0.5))


pdf(paste0(output_path, "/corrPlots/samples_correlation.pdf"), w = 10, h = 9)
draw(cor_plot_correction,
     column_title="",
     column_title_gp=grid::gpar(fontsize=16))
dev.off()



# Remove N12 --------------------------------------------------------------
counts <- read.table(counts_path, header = T, sep = "\t", row.names = 1)
dim(counts)
metadata = read.table(metadata_path, header = T, sep = "\t")
head(metadata)
write.table(cbind(data.frame(IDs = rownames(counts)) ,counts[,-which(metadata$sample_name == "N12")]), 
            file = paste0(output_path, "raw_counts_N12out.txt"), 
            row.names = F, col.names = T, quote = F, sep = "\t")

write.table(metadata[-which(metadata$sample_name == "N12"),], 
            file = paste0(output_path, "metadata_N12out.txt"), 
            row.names = F, col.names = T, quote = F, sep = "\t")

counts_N12out <- counts[,-which(metadata$sample_name == "N12")]
metadata_N12out <- metadata[-which(metadata$sample_name == "N12"),]
rownames(metadata_N12out) <- metadata_N12out$sample_name

dataNorm <- gc_norm(input_matrix = as.matrix(counts_N12out), 
                    grouping_var = metadata_N12out$groups, 
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

metadata_N12out$annotation = metadata_N12out$groups

dir.create(paste0(output_path, "/PCA/"))
pdf(paste0(output_path, "/PCA/PCA_PC1-3_N12out.pdf"), w=5.5, h=4)
variance_explained <- summary(pca)
to_plot=merge(pcaData, metadata_N12out[,c("sample_name","annotation")],by.x="row.names",by.y="row.names")


base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, fill=annotation)) + 
  geom_point(size=5, color = "black", shape = 21, stroke = 1.2) + 
  scale_fill_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata_N12out$annotation))], name = "Group" ) +
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata_N12out$annotation))], name = "Group" ) +
  # stat_ellipse(geom = "polygon",type = "norm",level = 0.7, alpha = 0.3, ) +
  ggforce::geom_mark_ellipse(aes(fill = annotation, color = annotation), expand = unit(5, "mm")) +
  xlim(c(-110,100)) + ylim(c(-70,70))+
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=17,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=17,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
base_plot + theme(legend.position = "right") #c(0.5, 0.5)


base_plot <- ggplot(to_plot, aes(x=PC1, y=PC3, fill=annotation)) + 
  geom_point(size=5, color = "black", shape = 21, stroke = 1.2) + 
  scale_fill_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata_N12out$annotation))], name = "Group" ) +
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata_N12out$annotation))], name = "Group" ) +
  # stat_ellipse(geom = "polygon",type = "norm",level = 0.7, alpha = 0.3, ) +
  ggforce::geom_mark_ellipse(aes(fill = annotation, color = annotation), expand = unit(5, "mm")) +
  xlim(c(-110,110)) + ylim(c(-75,90))+
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=17,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=17,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
base_plot + theme(legend.position = "right") #c(0.5, 0.5)


dev.off()


# Run chromvar ------------------------------------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2025)

dir.create(paste0(output_path, "/chromVAR/"))
rowRanges <- Signac::StringToGRanges(rownames(counts_N12out), sep = c(":", "-"))
counts_save=counts_N12out
counts <- SummarizedExperiment(assays=list(counts=as.matrix(counts_N12out)),
                               rowRanges=rowRanges, colData=metadata_N12out)
# add GC bias
counts <- addGCBias(counts, 
                    genome = BSgenome.Mmusculus.UCSC.mm10)
bg <- getBackgroundPeaks(object = counts, niterations=2000)

head(rowData(counts))

# get PFM data and match motifs
pfm <- readRDS(pfm_file)
motif_ix <- matchMotifs(pfm, counts, 
                        genome = BSgenome.Mmusculus.UCSC.mm10)

dev <- chromVAR::computeDeviations(object = counts, annotations = motif_ix, background_peaks = bg)
dev_matrix <- assays(dev)$deviations

motifs_annot <- data.frame(motifs = names(pfm), 
                           TF =  sapply(rownames(dev_matrix), function(i) pfm@listData[[i]]@tags$tf), 
                           family =  sapply(rownames(dev_matrix), function(i) pfm@listData[[i]]@tags$family), 
                           subfamily = sapply(rownames(dev_matrix), function(i) pfm@listData[[i]]@tags$sub.family))
rownames(motifs_annot) <- motifs_annot$TF
motifs_annot$subfamily[motifs_annot$subfamily %in% c('', "Unclassified", "Other")] <- motifs_annot$family[motifs_annot$subfamily %in% c('', "Unclassified", "Other")]

rownames(dev_matrix) <- sapply(rownames(dev_matrix), function(i) pfm@listData[[i]]@tags$tf)

z_matrix <- assays(dev)$z
rownames(z_matrix) <- sapply(rownames(z_matrix), function(i) pfm@listData[[i]]@tags$tf)

saveRDS(dev_matrix, file = paste0(output_path, "/chromVAR/dev_scores.rds"))
saveRDS(z_matrix, file = paste0(output_path, "/chromVAR/z_scores.rds"))

z_matrix_toSave <- rbind(t(data.frame(metadata_N12out$group)),z_matrix )
write.csv2(cbind(data.frame(row_names = rownames(z_matrix_toSave)),z_matrix_toSave), row.names = F, file = paste0(output_path, "/chromVAR/chromvar_scores.csv"))
write.table(cbind(data.frame(row_names = rownames(z_matrix_toSave)),z_matrix_toSave), row.names = F, file = paste0(output_path, "/chromVAR/chromvar_scores.txt"), col.names = F, sep = "\t")

# chromvar_mat <- scale(z_matrix) 
chromvar_mat <- z_matrix

rvdm = apply(chromvar_mat, 1, var)

colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
            '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
            '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000", '#E4C755', '#F7F398',
            '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
            '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
            '#968175', '#FF6347', '#4682B4', '#D2B48C', '#008080',
            '#D8BFD8', '#FF4500', '#DA70D6', '#EEE8AA', '#98FB98',
            '#AFEEEE', '#DB7093', '#FFEFD5', '#FFDAB9', '#CD853F')



rvdm = apply(chromvar_mat, 1, var)
Top_chromVAR <- rownames(chromvar_mat[order(rvdm, decreasing = T)[1:100] ,])

mat <- chromvar_mat[Top_chromVAR,]
col_fun = circlize::colorRamp2(c(-2, -1, 0 ,1, 2), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124"))
split.vector <- factor(metadata_N12out[colnames(mat), "groups"], levels = c("WT","NfatKO"))
colours <- list('Condition' = descriminative_colors)
colAnn <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(Condition =  split.vector),
                                            which = 'col',
                                            col = colours, 
                                            show_annotation_name = c(Condition = FALSE) )
indices <- which(rownames(mat) %in% c("Cebpbe", "Cebpd", "Cebpa", "Cebpb", "Ddit3", "Atf4",
                                              "Fos", "Jun", "Atf3", "Nfkb1", "Rel",
                                              "Mafg", "Maff", "Mafk",
                                              "Crem", "Creb1", "Xbp1", "Rara", "Rarb",
                                              "Foxo6", "Foxj2", "Foxo4", "Rfx1", "Irf1",
                                              "Elf1", "Elf3", "Nfyc", "Nfya", "Elk1", "E2f1", "E2f3",
                                              "Myc", "E2f4", "Nfatc4","Nfatc1", "Mef2a",
                                              "Runx3","Stat5a","Stat5b", "Nfatc3", "Bcl6", "Smad3"))
ha = rowAnnotation(foo = anno_mark(at = indices, 
                                   labels = rownames(mat)[indices],link_gp = gpar(lwd = 0.5), 
                                   labels_gp = gpar(fontsize = 6) ))
pdf(paste0(output_path, "/chromVAR/top100TF_chromVar.pdf"), w=10, h=10)

ht <-Heatmap(t(scale(t(mat))), name = "ChromVar \n z-score", right_annotation = ha,
             clustering_method_rows = "ward.D2", show_row_dend = T,
             column_split = split.vector, top_annotation = colAnn, #col = col_fun,
             show_column_names = T, column_title = NULL, cluster_columns = F,
             show_row_names = T,row_names_gp = gpar(fontsize = 5), row_dend_gp = gpar(lwd = 0.5))
draw(ht)

dev.off()


# Plot clustering ---------------------------------------------------------
counts <- read.table(counts_KP_peaks_path, header = T, sep = "\t", row.names = 1)
dim(counts)
counts <- counts[, colnames(counts) %in% metadata_N12out$sample_name]
summary(colnames(counts) %in% metadata$sample)

metadata_N12out = metadata_N12out[colnames(counts),]
dataNorm <- gc_norm(input_matrix = as.matrix(counts),
                    grouping_var = metadata_N12out$group,
                    fa_file = paste0(fasta_file),
                    by_group = F)
norm_data <- log(dataNorm+1)

peak_clusters <- read.table(cluster_peaks_path,
                            header = T, sep = "\t", quote = "")
rownames(peak_clusters) <- paste0(peak_clusters$seqnames, ":", peak_clusters$start, "-", peak_clusters$end)
summary(rownames(peak_clusters) %in% rownames(counts))

norm_data_subset <- norm_data[row.names(peak_clusters),]
cluster_assignment = peak_clusters[rownames(norm_data_subset), "peak_cluster"]


# Transpose to apply across columns
averaged <- sapply(levels(as.factor(metadata_N12out$group)), function(g) {
  rowMeans(norm_data_subset[, metadata_N12out$group == g, drop = FALSE])
})

# Convert to data.frame
averaged_df <- as.data.frame(averaged)
# Add row group as a column
averaged_df$row_group <- cluster_assignment


# Reshape to long format
plot_data <- tidyr::pivot_longer(averaged_df, cols = c("Nfatc1KO", "WT"),
                                 names_to = "Condition",
                                 values_to = "Value")

GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin", ggplot2::GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x'])
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", ggplot2::GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

pdf(paste0(output_path, "KP_peaks_V4_boxplots.pdf"), w=9,h=5)
plot_data$Condition <- factor(plot_data$Condition, levels = c("WT", "Nfatc1KO"))
ggplot(plot_data, aes(x = row_group, y = Value, fill = Condition)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  ggpubr::stat_compare_means(aes(group = Condition),
                             method = "t.test",
                             label = "p.signif") +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  #             alpha = 0.6) +
  labs(x = "KP clusters", y = "Averaged Accessibility", title = "Accessibility signal of KP cluster peaks in NfatKO samples") +
  theme_minimal()+
  scale_fill_manual(values=descriminative_colors)

ggplot(plot_data, aes(x = row_group, y = Value, fill = Condition)) + 
  geom_split_violin(trim = TRUE) + 
  geom_boxplot(width = 0.25, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0) +
  labs(x=NULL,y="GM Attitude Score") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  ggpubr::stat_compare_means(aes(group = Condition),
                             method = "t.test",
                             label = "p.signif") +
  scale_fill_manual(values=descriminative_colors)

dev.off()


