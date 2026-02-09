
# Library -----------------------------------------------------------------
library(GenomicRanges)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)
library(limma)
library(umap)
library(edgeR)
library(grid)
library(Rsamtools)
library(Biostrings)
library(EDASeq)
library(Signac)

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


# metadata_path <- "/mnt/curnagl//work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/results/PEPATAC_processing/KP_lung/consensus_peaks/metadata.txt"
# counts_path <- paste0("/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/results/PEPATAC_processing/KP_lung/consensus_peaks/raw_counts.txt")
# fasta_file <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/results/PEPATAC_processing/genome_folder/data/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/fasta/default/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1.fa"

dir.create(output_path, recursive = T)

descriminative_colors <- c("Healthy_blood"="#3288bdff",
                           "KP_lung_SiglecF_high" = "#d90017ff",
                           "Healthy_lung" = "#ffd92fff",
                           "KP_blood" = "#70c4b4ff",
                           "KP_lung_SiglecF_low" = "#ff7f00ff",
                           "batch1" = "#a0451fff", "batch2" = "black", "batch3" = "gray", "batch4" = "#984ea3ff", "batch5" = "yellow4")


# Load counts -------------------------------------------------------------
counts <- read.table(counts_path, header = T, sep = "\t", row.names = 1)

# Load samples metadata ---------------------------------------------------
metadata <- read.table(metadata_path, header = T, sep = "\t")
rownames(metadata) <- metadata$sample_name

colnames(counts) <- gsub("[.]", "-",colnames(counts))
metadata <- metadata[colnames(counts),] 

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

umap_res = umap::umap((pcaData[, 1:3]), verbose = T, n_neighbors = min(nrow(pcaData)-1, 15))
umap_dr = umap_res$layout
dim(umap_dr)
colnames(umap_dr) = c("UMAP1", "UMAP2")

metadata$annotation = metadata$groups


pdf(paste0(output_path, "UMAP.pdf"), w = 8, h = 5)

to_plot=merge(umap_dr,metadata[,c("sample_name","annotation","batch")],by.x="row.names",by.y="sample_name")
ggplot(to_plot, aes(x=UMAP1, y=UMAP2, color=annotation, shape=batch)) +
  geom_point(size=5) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))] ) +
  scale_shape_manual(values=c(16,17,18,8,15,7,9,10,12,13,14,1,2,3)) +
  theme_bw() +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))+
  guides(col=guide_legend("Experiment",override.aes = list(size=4)),shape=guide_legend("Cell types",override.aes = list(size=4)))

# no shape use
ggplot(to_plot, aes(x=UMAP1, y=UMAP2, color=annotation)) +
  geom_point(size=5) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))] ) +
  theme_bw() +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))+
  guides(col=guide_legend("Experiment",override.aes = list(size=4)),shape=guide_legend("Cell types",override.aes = list(size=4))) 

dev.off()


pdf(paste0(output_path, "PCA.pdf"), w = 7, h = 7)

barplot(summary(pca)$importance[2,][1:10], ylab = "Proportion of variance")

to_plot=merge(pcaData, metadata[,c("sample_name","annotation","batch")], by.x="row.names", by.y="sample_name")
shapes = c(16,17,18,8,15)
names(shapes) = sort(unique(metadata$batch))
pairs(pcaData[,1:5] , col = sapply(metadata$annotation, function(i) descriminative_colors[i]), upper.panel = NULL, pch = sapply(metadata$batch, function(i) shapes[i]), cex = 1)

dev.off()


pdf(paste0(output_path, "PCA_PC1-3.pdf"), w=4, h=4)
variance_explained <- summary(pca)
to_plot=merge(pcaData, metadata[,c("sample_name","annotation","batch")],by.x="row.names",by.y="sample_name")
base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, fill=annotation)) + #, shape=batch
  geom_point(size=4, color = "black", shape = 21, stroke = 1.2) + 
  scale_fill_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  # scale_shape_manual(values=c(16,17,18,8,15,7,9,10,12,13,14,1,2,3)) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot + theme(legend.position = "none") #c(0.7, 0.5)
scatter_only_plot 

ggplot(to_plot, aes(x=PC1, y=PC2, color=annotation, shape=batch)) +
  geom_point(size=5) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))] ) +
  scale_shape_manual(values=c(16,17,18,8,15,7,9,10,12,13,14,1,2,3)) +
  theme_bw() +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))+
  guides(col=guide_legend("Experiment",override.aes = list(size=4)),shape=guide_legend("Cell types",override.aes = list(size=4)))


base_plot <- ggplot(to_plot, aes(x=PC1, y=PC3, fill=annotation)) + #, shape=batch
  geom_point(size=4, color = "black", shape = 21, stroke = 1.2) + 
  scale_fill_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  # scale_shape_manual(values=c(16,17,18,8,15,7,9,10,12,13,14,1,2,3)) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + ylab(paste0("PC3 (", round(variance_explained$importance[2,"PC3"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot + theme(legend.position = "none") #c(0.5, 0.5)
scatter_only_plot 

base_plot <- ggplot(to_plot, aes(x=PC2, y=PC3, fill=annotation)) + #, shape=batch
  geom_point(size=4, color = "black", shape = 21, stroke = 1.2) + 
  scale_fill_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  # scale_shape_manual(values=c(16,17,18,8,15,7,9,10,12,13,14,1,2,3)) +
  theme_bw() + xlab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) + ylab(paste0("PC3 (", round(variance_explained$importance[2,"PC3"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot + theme(legend.position = "none") #c(0.5, 0.5)
scatter_only_plot 


dev.off()


# Distances and correlation plots -------------------------------------------------------

sampleDists <- dist(t(norm_data))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- sapply(metadata$sample_name, function(i) unlist(strsplit(i,"_"))[[1]] )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

ha = rowAnnotation(Condition = metadata$groups, 
                   Batch = metadata$batch, 
                   annotation_name_rot = 45,
                   col = list(Condition = descriminative_colors,
                              Batch = descriminative_colors
                   )
)
dist_plot_before_correction <- Heatmap(sampleDistMatrix, name = "Distance", left_annotation = ha, col = colors, rect_gp = gpar(col = "black", lwd = 0.5))


# Remove batch effect using limma::removeBatchEffect 
mm <- model.matrix(~metadata$groups)
mat <- norm_data
mat <- limma::removeBatchEffect(mat, batch = metadata$batch, design = mm)  #http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-after-vst-are-there-still-batches-in-the-pca-plot

sampleDists <- dist(t(mat))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- sapply(metadata$sample_name, function(i) unlist(strsplit(i,"_"))[[1]] )
colnames(sampleDistMatrix) <- NULL

dist_plot_after_correction <- Heatmap(sampleDistMatrix, name = "Distance", left_annotation = ha, col = colors, rect_gp = gpar(col = "black", lwd = 0.5))

pdf(paste0(output_path, "samples_distances.pdf"), w = 10, h = 9)
draw(dist_plot_before_correction,
     column_title="Before batch effect correction",
     column_title_gp=grid::gpar(fontsize=16))
draw(dist_plot_after_correction,
     column_title="After batch effect correction",
     column_title_gp=grid::gpar(fontsize=16))
dev.off()



# Compute correlations instead of distances --------------------

sampleCor <- cor(norm_data)
sampleCorMatrix <- as.matrix(sampleCor)
rownames(sampleCorMatrix) <- sapply(metadata$sample_name, function(i) unlist(strsplit(i,"_"))[[1]] )
colnames(sampleCorMatrix) <- NULL
colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)

ha = rowAnnotation(Condition = metadata$groups, 
                   Batch = metadata$batch, 
                   annotation_name_rot = 45,
                   col = list(Condition = descriminative_colors,
                              Batch = descriminative_colors
                   )
)
cor_plot_before_correction <- Heatmap(sampleCorMatrix, name = "Correlation", left_annotation = ha, col = colors, rect_gp = gpar(col = "black", lwd = 0.5))


# Remove batch effect using limma::removeBatchEffect 

sampleCor <- cor(mat)
sampleCorMatrix <- as.matrix(sampleCor)
rownames(sampleCorMatrix) <- sapply(metadata$sample_name, function(i) unlist(strsplit(i,"_"))[[1]] )
colnames(sampleCorMatrix) <- NULL

cor_plot_after_correction <- Heatmap(sampleCorMatrix, name = "Correlation", left_annotation = ha, col = colors, rect_gp = gpar(col = "black", lwd = 0.5))

pdf(paste0(output_path, "samples_correlations.pdf"), w = 10, h = 9)
draw(
  cor_plot_before_correction,
  column_title = "Before batch effect correction",
  column_title_gp = grid::gpar(fontsize = 16)
)
draw(
  cor_plot_after_correction,
  column_title = "After batch effect correction",
  column_title_gp = grid::gpar(fontsize = 16)
)
dev.off()

# Save PCA coordinates ----------------------------------------------------
saveRDS(pcaData, file = paste0(output_path, "pca_coordinates.rds"))

