
# Libraries ---------------------------------------------------------------

library(GenomicRanges)
library(ggplot2)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(dplyr)
library(SummarizedExperiment)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(BiocParallel)
set.seed(12345)

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
get_TF_heatmap <- function(gene_list, genes_matrix, metadata, mat.name = "Deviation score",
                           split.vector = NULL, colors.vector = descriminative_colors, row_clustering = T,clustering_method_rows = "ward.D",
                           target_annotation = NULL, colors_annotation=NULL, text_size=8, show_row_names = T, show_col_names = T, scale = F){
  gene_list <- gene_list[which(gene_list %in% rownames(genes_matrix))]
  genes_matrix <- genes_matrix[, order(metadata$groups)]
  split.vector <- factor(metadata[colnames(genes_matrix), "groups"], levels = c("Blood", "Adjacent_lung", "Tumor"))
  colours <- list('Condition' = colors.vector,
                  'Batch' = colors.vector)
  colAnn <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(Condition =  split.vector),
                                              which = 'col',
                                              col = colours, show_annotation_name = c(Condition = FALSE) )
  if(!is.null(target_annotation)){
    rowAnn <- ComplexHeatmap::HeatmapAnnotation(df = target_annotation,
                                                which = 'row',
                                                col = colors_annotation)
    if(!scale){
      ht <-ComplexHeatmap::Heatmap(genes_matrix[gene_list , ],
                                   name = mat.name, column_split = split.vector,
                                   top_annotation = colAnn, right_annotation = rowAnn,
                                   column_title = NULL, cluster_columns = F, cluster_rows = row_clustering, clustering_method_rows = clustering_method_rows,
                                   show_column_names = show_col_names, show_row_names = show_row_names, rect_gp = grid::gpar(col = "black", lwd = 0.5),
                                   row_names_gp = grid::gpar(fontsize = text_size))
    }else{
      ht <-ComplexHeatmap::Heatmap(t(scale(t(genes_matrix[gene_list , ]))),
                                   name = mat.name, column_split = split.vector,
                                   top_annotation = colAnn, right_annotation = rowAnn,
                                   column_title = NULL, cluster_columns = F, cluster_rows = row_clustering, clustering_method_rows = clustering_method_rows,
                                   show_column_names = show_col_names, show_row_names = show_row_names, rect_gp = grid::gpar(col = "black", lwd = 0.5),
                                   row_names_gp = grid::gpar(fontsize = text_size))
    }
    
  } else{
    
    if(!scale){
      ht <-ComplexHeatmap::Heatmap(genes_matrix[gene_list , ],
                                   name = mat.name, column_split = split.vector,
                                   top_annotation = colAnn,
                                   column_title = NULL, cluster_columns = F,cluster_rows = row_clustering, clustering_method_rows = clustering_method_rows,
                                   show_column_names = show_col_names, show_row_names = show_row_names, rect_gp = grid::gpar(col = "black", lwd = 0.5),
                                   row_names_gp = grid::gpar(fontsize = text_size))
    }else{
      ht <-ComplexHeatmap::Heatmap(t(scale(t(genes_matrix[gene_list , ]))),
                                   name = mat.name, column_split = split.vector,
                                   top_annotation = colAnn,
                                   column_title = NULL, cluster_columns = F,cluster_rows = row_clustering, clustering_method_rows = clustering_method_rows,
                                   show_column_names = show_col_names, show_row_names = show_row_names, rect_gp = grid::gpar(col = "black", lwd = 0.5),
                                   row_names_gp = grid::gpar(fontsize = text_size))
    }
    
  }
  
  
  return(ht)
}

# Parameters --------------------------------------------------------------
counts_path <- as.character(commandArgs(TRUE)[1])
metadata_path <- as.character(commandArgs(TRUE)[2])
fasta_file <- as.character(commandArgs(TRUE)[3])
output_path <- as.character(commandArgs(TRUE)[4])
pfm_file <- as.character(commandArgs(TRUE)[5])
motif_annotation_path <- as.character(commandArgs(TRUE)[6])
remove_donor9 <- as.logical(commandArgs(TRUE)[7])

# metadata_path <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/results/PEPATAC_processing/Human_data/consensus_peaks/metadata.txt"
# counts_path <- paste0("/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/results/PEPATAC_processing/Human_data/consensus_peaks/raw_counts.txt")
# fasta_file <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/results/PEPATAC_processing/genome_folder/data/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/fasta/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.fa"
# output_path <- paste0("/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/results_human_data/")
# # ga_path <- paste0("/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/results/PEPATAC_processing/Human_data/gene_activity/gene_activity.rds")
# pfm_file <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/SuperCellMultiomicsAnalyses/H12CORE_human_pfm.rds"
# motif_annotation_path <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/input_data/HOCOMOCO_v12_annotation_human.txt"


# counts_path <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/results/PEPATAC_processing/Human_data/consensus_peaks_pval/raw_counts.txt"
# metadata_path <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/results/PEPATAC_processing/Human_data/consensus_peaks_pval/metadata.txt"
# fasta_file <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/results/PEPATAC_processing/genome_folder/data/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/fasta/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.fa"
# output_path <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/results_human_data/"
# pfm_file <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/SuperCellMultiomicsAnalyses/H12CORE_human_pfm.rds"
# motif_annotation_path <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/backup/KP_results_v1/final_experiment/hocomoco_version/motifs_data/H12CORE_motifs.tsv"

# output_path <- paste0("/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/results/PEPATAC_processing/Human_data/figures/")
dir.create(output_path, recursive = T)

descriminative_colors = c("Blood"="#3288bdff",
                          "Tumor"="#d90017ff",
                          "Adjacent_lung"="#ffd92fff",
                          "batch_1" = "#a0451fff", "batch_2" = "black", "batch3" = "gray", "batch4" = "#984ea3ff", "batch5" = "yellow4")
# register(MulticoreParam(3))

# Read the raw counts -----------------------------------------------------

counts_data <- read.table(counts_path, header = T, sep = "\t", row.names = 1)

# Load samples metadata ---------------------------------------------------

metadata <- read.table(metadata_path, header = T, sep = "\t")
rownames(metadata) <- metadata$sample_name

colnames(counts_data) <- sapply(colnames(counts_data), function(i) unlist(strsplit(i, "_"))[[1]])
counts_data <- counts_data[, gsub("[.]", "-",colnames(counts_data)) %in% metadata$sample_name]
colnames(counts_data) <- gsub("[.]", "-",colnames(counts_data))
metadata <- metadata[colnames(counts_data),]
metadata$batch <- paste0("batch_",c(rep(1,6), rep(2,6),1,1,2,1,3,3,3))
metadata$patient <- gsub("-1|-2","", metadata$patient)
batch <- metadata$batch
condition <- metadata$groups

if(remove_donor9){
  counts_data <- counts_data[, metadata$patient != "G-LC-09"]
  metadata <- metadata[metadata$patient != "G-LC-09", ]
}

# Normalize data ----------------------------------------------------------
dataNorm <- gc_norm(input_matrix = as.matrix(counts_data), 
                    grouping_var = metadata$groups, 
                    fa_file = fasta_file, 
                    by_group = F)
norm_data <- log(dataNorm+1)

norm_data_uncorrected <- norm_data

# Remove batch effect using limma::removeBatchEffect 
mm <- model.matrix(~metadata$groups)
mat <- norm_data
mat <- limma::removeBatchEffect(mat, batch = metadata$batch, design = mm)  #http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-after-vst-are-there-still-batches-in-the-pca-plot
norm_data_corrected <- mat
norm_data <- norm_data_corrected

# PCA on top 5000 features ------------------------------------------------
dir.create(paste0(output_path,"PCA"), recursive = T)
rvdm = apply(norm_data_uncorrected, 1, var)
norm_data_subset = t(norm_data_uncorrected[order(rvdm, decreasing = T)[1:5000] ,])

pca = prcomp(norm_data_subset,scale. = T)
print(summary(pca))
pcaData = as.data.frame(pca$x)
rownames(pcaData) = rownames(norm_data_subset)

metadata$annotation = metadata$groups

pdf(paste0(output_path, "/PCA/PCA_pairwise_before_batch.pdf"), w = 7, h = 7)

barplot(summary(pca)$importance[2,][1:10], ylab = "Proportion of variance")

to_plot=merge(pcaData, metadata[,c("sample_name","annotation","batch","patient")], by.x="row.names", by.y="sample_name")
shapes = c(16,17,18,8,7)
names(shapes) = sort(unique(metadata$patient))
print(pairs(pcaData[,1:5] , col = sapply(metadata$annotation, function(i) descriminative_colors[i]), upper.panel = NULL, 
            pch = sapply(metadata$patient, function(i) shapes[i]), cex = 1))
dev.off()


pdf(paste0(output_path, "/PCA/PCA_plots_all_samples_before_batch.pdf"), w=6, h=4)
variance_explained <- summary(pca)
to_plot <- merge(pcaData, metadata[,c("sample_name","annotation","batch","patient")],by.x="row.names",by.y="sample_name")
to_plot$annotation <- as.factor(to_plot$annotation)
to_plot$sample_name <- as.factor(to_plot$Row.names)

base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=annotation, fill=annotation, shape=sample_name)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  scale_shape_manual(values=c(16,17,18,8,15,7,9,10,12,13,14,1,2,3,4,5,6,11,19)) +
  scale_fill_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))] ) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + 
  ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot 
print(scatter_only_plot )

base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=annotation,shape=patient)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  scale_shape_manual(values=c(16,17,18,15, 8,7,9,10,12,13,14,1,2,3,4,5,6,11,19)) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + 
  ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot 
print(scatter_only_plot )


base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=annotation, fill=annotation, shape=patient)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  scale_shape_manual(values=c(16,17,18,15,8,7,9,10,12,13,14,1,2,3,4,5,6,11,19)) +
  scale_fill_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))] ) +
  ggforce::geom_mark_ellipse(aes(fill = annotation, color = annotation, group = annotation), expand = unit(5, "mm")) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + 
  ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  # xlim(c(-40,50)) + ylim(c(-60,40)) +
  xlim(c(min(to_plot$PC1) - 10,max(to_plot$PC1) + 10)) + ylim(c(min(to_plot$PC2) - 10,max(to_plot$PC2) + 10)) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot #+ theme(legend.position = "none") #c(0.5, 0.5)
print(scatter_only_plot )

dev.off()

# PCA on top features after correction ------------------------------------
dir.create(paste0(output_path,"PCA"), recursive = T)
rvdm = apply(norm_data_corrected, 1, var)
norm_data_subset = t(norm_data_corrected[order(rvdm, decreasing = T)[1:5000] ,])

pca = prcomp(norm_data_subset,scale. = T)
print(summary(pca))
pcaData = as.data.frame(pca$x)
rownames(pcaData) = rownames(norm_data_subset)

metadata$annotation = metadata$groups

pdf(paste0(output_path, "/PCA/PCA_pairwise_after_batch.pdf"), w = 7, h = 7)

barplot(summary(pca)$importance[2,][1:10], ylab = "Proportion of variance")

to_plot=merge(pcaData, metadata[,c("sample_name","annotation","batch","patient")], by.x="row.names", by.y="sample_name")
shapes = c(16,17,18,8,7)
names(shapes) = sort(unique(metadata$patient))
print(pairs(pcaData[,1:5] , col = sapply(metadata$annotation, function(i) descriminative_colors[i]), upper.panel = NULL, 
            pch = sapply(metadata$patient, function(i) shapes[i]), cex = 1))
dev.off()


pdf(paste0(output_path, "/PCA/PCA_plots_all_samples_after_batch.pdf"), w=6, h=4)
variance_explained <- summary(pca)
to_plot <- merge(pcaData, metadata[,c("sample_name","annotation","batch","patient")],by.x="row.names",by.y="sample_name")
to_plot$annotation <- as.factor(to_plot$annotation)
to_plot$sample_name <- as.factor(to_plot$Row.names)

base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=annotation, fill=annotation, shape=sample_name)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  scale_shape_manual(values=c(16,17,18,8,15,7,9,10,12,13,14,1,2,3,4,5,6,11,19)) +
  scale_fill_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))] ) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + 
  ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot 
print(scatter_only_plot )

base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=annotation,shape=patient)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  scale_shape_manual(values=c(16,17,18,15, 8,7,9,10,12,13,14,1,2,3,4,5,6,11,19)) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + 
  ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot 
print(scatter_only_plot )


base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=annotation, fill=annotation, shape=patient)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  scale_shape_manual(values=c(16,17,18,15,8,7,9,10,12,13,14,1,2,3,4,5,6,11,19)) +
  scale_fill_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))] ) +
  ggforce::geom_mark_ellipse(aes(fill = annotation, color = annotation, group = annotation), expand = unit(5, "mm")) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + 
  ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  # xlim(c(-40,50)) + ylim(c(-60,40)) +
  xlim(c(min(to_plot$PC1) - 10,max(to_plot$PC1) + 10)) + ylim(c(min(to_plot$PC2) - 10,max(to_plot$PC2) + 10)) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot #+ theme(legend.position = "none") #c(0.5, 0.5)
print(scatter_only_plot )

dev.off()

############################################################################
############################################################################
###                                                                      ###
###                      MERGE TECHNICAL REPLICATES                      ###
###                                                                      ###
############################################################################
############################################################################

metadata$replicate <- plyr::revalue(metadata$sample, c("G-LC-05-B-1"="G-LC-05-B", "G-LC-05-B-2"="G-LC-05-B", 
                                                       "G-LC-05-T-1"="G-LC-05-T", "G-LC-05-T-2"="G-LC-05-T",
                                                       "G-LC-08-T-2" = "G-LC-08-T", "G-LC-08-B-2" = "G-LC-08-B"))

counts_data <- t(rowsum(t(counts_data), group = metadata$replicate))
norm_data_corrected  <- t(rowsum(t(norm_data_corrected), group = metadata$replicate) / as.vector(table(metadata$replicate)))

metadata <- metadata[-which(duplicated(metadata$replicate)),]
rownames(metadata) <- metadata$replicate
metadata <- metadata[colnames(counts_data),]

df_output <- cbind(data.frame(IDs = rownames(counts_data)), counts_data)
write.table(df_output, file=paste0(output_path, "/merged_data/raw_counts_merged.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# df_output <- cbind(data.frame(IDs = rownames(counts_data)), counts_data[, metadata$patient != "G-LC-09"])
# write.table(df_output, file=paste0(output_path, "/merged_data/raw_counts_merged2.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


df_output <- cbind(data.frame(IDs = rownames(norm_data_corrected)), norm_data_corrected)
write.table(df_output, file=paste0(output_path, "/merged_data/corrected_counts_merged.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

write.table(metadata, file=paste0(output_path, "/merged_data/metadata_merged.txt"), quote = F, row.names = F, col.names = T, sep = "\t")



if(remove_donor9){
  norm_data <- norm_data_corrected
}else{
  dataNorm <- gc_norm(input_matrix = as.matrix(counts_data),
                      grouping_var = metadata$groups,
                      fa_file = fasta_file,
                      by_group = F)
  norm_data <- log(dataNorm+1)
}

score <- sapply(1:nrow(norm_data), function(i){
  fit <- lm(norm_data[i,] ~ metadata$groups)
  anova <- anova(fit)
  return(anova$`F value`[1])
}  )
norm_data_subset = t(norm_data[order(score, decreasing = T)[1:10000] ,])

# Run PCA and UMAP 
pca = prcomp(norm_data_subset,scale. = T)
print(summary(pca))
pcaData = as.data.frame(pca$x)
rownames(pcaData) = rownames(norm_data_subset)

metadata$annotation = metadata$groups

pdf(paste0(output_path, "/PCA/PCA_pairwise_merged_samples.pdf"), w = 7, h = 7)

barplot(summary(pca)$importance[2,][1:10], ylab = "Proportion of variance")

to_plot=merge(pcaData, metadata[,c("sample_name","annotation","patient")], by.x="row.names", by.y="sample_name")
shapes = c(16,17,18,15,8)
names(shapes) = sort(unique(metadata$patient))
print(pairs(pcaData[,1:5] , col = sapply(metadata$annotation, function(i) descriminative_colors[i]), upper.panel = NULL, pch = sapply(metadata$patient, function(i) shapes[i]), cex = 1))
dev.off()


pdf(paste0(output_path, "/PCA/PCA_plots_merged_samples.pdf"), w=6, h=4)
variance_explained <- summary(pca)
to_plot=merge(pcaData, metadata[,c("replicate","annotation","batch", "patient")],by.x="row.names",by.y="replicate")
to_plot$annotation <- as.factor(to_plot$annotation)
to_plot$sample_name <- as.factor(to_plot$Row.names)

base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=annotation, fill=annotation, shape=sample_name)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  scale_shape_manual(values=c(16,17,18,15,8,7,9,10,12,13,14,1,2,3,4,5)) +
  scale_fill_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))] ) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + 
  ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot 
print(scatter_only_plot )

base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=annotation, fill=annotation, shape=patient)) + #, shape=batch
  geom_point(size=4) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))], name = "Group" ) +
  scale_shape_manual(values=c(16,17,18,15,8,7,9,10,12,13,14,1,2,3)) +
  scale_fill_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(metadata$annotation))] ) +
  ggforce::geom_mark_ellipse(aes(fill = annotation, color = annotation, group = annotation), expand = unit(5, "mm")) +
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + 
  ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  # xlim(c(-80,100)) + ylim(c(-80,50)) +
  xlim(c(min(to_plot$PC1) - 10,max(to_plot$PC1) + 10)) + ylim(c(min(to_plot$PC2) - 10,max(to_plot$PC2) + 10)) +
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
scatter_only_plot <- base_plot #+ theme(legend.position = "none") #c(0.5, 0.5)
print(scatter_only_plot )

dev.off()



# total_var <- matrixStats::rowVars(norm_data)
# within_var <- sapply(metadata$groups, function(g) {
#   matrixStats::rowVars(norm_data[, metadata$groups == g, drop = FALSE])
# }) %>% rowMeans()
# between_var <- total_var - within_var
# f_score <- between_var / (within_var + 1e-8)  

score <- sapply(1:nrow(norm_data), function(i){
  fit <- lm(norm_data[i,] ~ metadata$groups)
  anova <- anova(fit)
  return(anova$`F value`[1])
}  )
norm_data_subset = t(norm_data[order(score, decreasing = T)[1:10000] ,])
write.table(data.frame(IDs = rownames(norm_data[order(score, decreasing = T)[1:10000] ,])), file=paste0(output_path, "10000_variablePeaks.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

pca = prcomp(norm_data_subset,scale. = T)
top_pc1_up_merged <- rownames(pca$rotation)[order((pca$rotation[,1]), decreasing = T)][1:500]
top_pc1_down_merged <- rownames(pca$rotation)[order((pca$rotation[,1]), decreasing = F)][1:500]

top_pc2_up_merged <- rownames(pca$rotation)[order((pca$rotation[,2]), decreasing = T)][1:500]
top_pc2_down_merged <- rownames(pca$rotation)[order((pca$rotation[,2]), decreasing = F)][1:500]

write.table(data.frame(IDs = top_pc1_up_merged), file=paste0(output_path, "/PCA/top_pc1_up.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
write.table(data.frame(IDs = top_pc1_down_merged), file=paste0(output_path, "/PCA/top_pc1_down.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
write.table(data.frame(IDs = top_pc2_up_merged), file=paste0(output_path, "/PCA/top_pc2_up.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
write.table(data.frame(IDs = top_pc2_down_merged), file=paste0(output_path, "/PCA/top_pc2_down.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


# Plot top PCs peaks ---------------------------------------------------------
groups_factor <- factor(metadata$groups, levels = c("Blood", "Adjacent_lung", "Tumor"))
sorted_indices <- order(groups_factor)
col_an = HeatmapAnnotation(Group = metadata$groups[sorted_indices], 
                           annotation_name_rot = 45,
                           col = list(Group = descriminative_colors
                           )
)


norm_data_subset = norm_data[c(top_pc1_up_merged, top_pc1_down_merged) ,]

pdf(paste0(output_path, "/PCA/top_peaks_PC1_mergedData.pdf"), w=8, h=8)
Heatmap(t(scale(t(norm_data_subset[, sorted_indices]))) , show_row_names = F,
        top_annotation = col_an, 
        show_row_dend = F, cluster_rows =T, cluster_columns = F)
dev.off()


norm_data_subset = norm_data[c(top_pc2_up_merged, top_pc2_down_merged) ,]

pdf(paste0(output_path, "/PCA/top_peaks_PC2_mergedData.pdf"), w=8, h=8)
Heatmap(t(scale(t(norm_data_subset[, sorted_indices]))) , show_row_names = F,
        # left_annotation = ha,
        top_annotation = col_an, 
        show_row_dend = F, cluster_rows =T, cluster_columns = F)
dev.off()


# Run chromVAr ------------------------------------------------------------

# Get chromVAR deviation scores
# total_var <- matrixStats::rowVars(norm_data)
# within_var <- sapply(metadata$groups, function(g) {
#   matrixStats::rowVars(norm_data[, metadata$groups == g, drop = FALSE])
# }) %>% rowMeans()
# between_var <- total_var - within_var
# f_score <- between_var / (within_var + 1e-8)  
score <- sapply(1:nrow(norm_data), function(i){
  fit <- lm(norm_data[i,] ~ metadata$groups)
  anova <- anova(fit)
  return(anova$`F value`[1])
}  )
counts_subset <- counts_data[rownames(norm_data)[order(score, decreasing = T)[1:10000]],]

rowRanges <- Signac::StringToGRanges(rownames(counts_subset), sep = c(":", "-"))
counts <- SummarizedExperiment(assays=list(counts=as.matrix(counts_subset)),
                               rowRanges=rowRanges, colData=metadata)
# add GC bias
counts <- addGCBias(counts, 
                    genome = BSgenome.Hsapiens.UCSC.hg38)
bg <- getBackgroundPeaks(object = counts) #, niterations=2000

head(rowData(counts))

# get PFM data and match motifs
pfm <- readRDS(pfm_file)
motif_ix <- matchMotifs(pfm, counts, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)

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

z_matrix_toSave <- rbind(t(data.frame(metadata$groups)),z_matrix )
write.csv2(cbind(data.frame(row_names = rownames(z_matrix_toSave)),z_matrix_toSave), row.names = F, file = paste0(output_path, "/chromVAR/chromvar_scores.csv"))


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


# Heatmap with avreaged signal per TF family 
dir.create(paste0(output_path, "/chromVAR"))
pdf(paste0(output_path, "/chromVAR/chromVar_TFfamilies_heatmap.pdf"), w=10, h=9)
annot_df <- data.frame(row.names = rownames(chromvar_mat), 
                       TF_cluster = motifs_annot[rownames(chromvar_mat), "family"])
annot_colors <- as.list(colors[1:length(unique(motifs_annot$subfamily))])
names(annot_colors) <- unique(motifs_annot$subfamily)

agg_chromvar_mat <- aggregate(. ~ factor(motifs_annot[rownames(chromvar_mat), "family"]), data=as.data.frame(chromvar_mat), FUN=mean)
write.csv2(agg_chromvar_mat, row.names = F, file = paste0(output_path, "/chromVAR/chromvar_FamilyScores.csv"))

rownames(agg_chromvar_mat) <- agg_chromvar_mat[,1]
agg_chromvar_mat <- agg_chromvar_mat[,-1]
print(get_TF_heatmap(gene_list = rownames(agg_chromvar_mat), 
                     genes_matrix = t(scale(t(agg_chromvar_mat))), metadata = metadata,
                     target_annotation = NULL, 
                     row_clustering = T, clustering_method_rows = "ward.D2",
                     show_col_names = T, text_size = 6))
dev.off()

pdf(paste0(output_path, "/chromVAR/chromVar_TFsubfamilies_heatmap.pdf"), w=10, h=25)
annot_df <- data.frame(row.names = rownames(chromvar_mat), 
                       TF_cluster = motifs_annot[rownames(chromvar_mat), "subfamily"])
annot_colors <- as.list(colors[1:length(unique(motifs_annot$subfamily))])
names(annot_colors) <- unique(motifs_annot$subfamily)

agg_chromvar_mat <- aggregate(. ~ factor(motifs_annot[rownames(chromvar_mat), "subfamily"]), data=as.data.frame(chromvar_mat), FUN=mean)
write.csv2(agg_chromvar_mat, row.names = F, file = paste0(output_path, "/chromVAR/chromvar_subFamilyScores.csv"))

rownames(agg_chromvar_mat) <- agg_chromvar_mat[,1]
agg_chromvar_mat <- agg_chromvar_mat[,-1]
print(get_TF_heatmap(gene_list = rownames(agg_chromvar_mat), 
                     genes_matrix = t(scale(t(agg_chromvar_mat))), metadata = metadata,
                     target_annotation = NULL, 
                     row_clustering = T, clustering_method_rows = "ward.D2",
                     show_col_names = T, text_size = 6))
dev.off()

# total_var <- matrixStats::rowVars(chromvar_mat)
# within_var <- sapply(metadata$groups, function(g) {
#   matrixStats::rowVars(chromvar_mat[, metadata$groups == g, drop = FALSE])
# }) %>% rowMeans()
# between_var <- total_var - within_var
# f_score <- between_var / (within_var + 1e-8)  # Add small constant to avoid divide-by-zero
score <- sapply(1:nrow(chromvar_mat), function(i){
  fit <- lm(chromvar_mat[i,] ~ metadata$groups)
  anova <- anova(fit)
  return(anova$`F value`[1])
}  )
Top_chromVAR <- rownames(chromvar_mat[order(score, decreasing = T)[1:300] ,])

rvdm = apply(chromvar_mat, 1, var)
Top_chromVAR <- rownames(chromvar_mat[order(rvdm, decreasing = T)[1:300] ,])

mat <- chromvar_mat[Top_chromVAR,]
col_fun = circlize::colorRamp2(c(-2, -1, 0 ,1, 2), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124"))
split.vector <- factor(metadata[colnames(mat), "groups"], levels = c("Blood","Adjacent_lung","Tumor"))
colours <- list('Condition' = descriminative_colors)
colAnn <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(Condition =  split.vector),
                                            which = 'col',
                                            col = colours, 
                                            show_annotation_name = c(Condition = FALSE) )
indices <- which(rownames(mat) %in% toupper(c("Cebpbe", "Cebpd", "Cebpa", "Cebpb", "Ddit3", "Atf4",
                                              "Fos", "Jun", "Atf3", "Nfkb1", "Rel",
                                              "Mafg", "Maff", "Mafk",
                                              "Crem", "Creb1", "Xbp1", "Rara", "Rarb",
                                              "Foxo6", "Foxj2", "Foxo4", "Rfx1", "Irf1",
                                              "Elf1", "Elf3", "Nfyc", "Nfya", "Elk1", "E2f1", "E2f3",
                                              "Myc", "E2f4", "Nfac4","Nfac1", "Mef2a",
                                              "Runx3","Sta5a","Sta5b", "Nfac3", "Bcl6", "Smad3")))
ha = rowAnnotation(foo = anno_mark(at = indices, 
                                   labels = rownames(mat)[indices],link_gp = gpar(lwd = 0.5), 
                                   labels_gp = gpar(fontsize = 6) ))
pdf(paste0(output_path, "/chromVAR/top300TF_chromVar.pdf"), w=10, h=10)

ht <-Heatmap(t(scale(t(mat))), name = "ChromVar \n z-score", right_annotation = ha,
             clustering_method_rows = "ward.D2", show_row_dend = T,
             column_split = split.vector, top_annotation = colAnn, #col = col_fun,
             show_column_names = T, column_title = NULL, cluster_columns = F,
             show_row_names = F,row_names_gp = gpar(fontsize = 4), row_dend_gp = gpar(lwd = 0.5))
draw(ht)


# ht <-Heatmap(mat, name = "ChromVar \n z-score", right_annotation = ha,
#              clustering_method_rows = "ward.D2", show_row_dend = T,
#              column_split = split.vector, top_annotation = colAnn, #col = col_fun,
#              show_column_names = T, column_title = NULL, cluster_columns = F,
#              show_row_names = F,row_names_gp = gpar(fontsize = 4), row_dend_gp = gpar(lwd = 0.5))
# draw(ht)
dev.off()

