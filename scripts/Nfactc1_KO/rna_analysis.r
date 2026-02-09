
# Libraries ----------------------------------------------------------------

library(readxl)
library(SummarizedExperiment)
library(edgeR)
library(dplyr)
library(DESeq2)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)

# Parameters --------------------------------------------------------------

RNA_seq_data_path <- as.character(commandArgs(TRUE)[1])
output_path <- as.character(commandArgs(TRUE)[2])
dir.create(output_path)

# Metadata ----------------------------------------------------------------

metadata <- data.frame(samples = c("KO5SHI","KO7SHI","KO8SHI","KO9SHI","WT2SHI","WT5SHI","WT7SHI","WT8SHI"))
metadata$Sample_ID <- metadata$samples
metadata$Group <- c(rep("KO", 4), rep("WT", 4))
rownames(metadata) <- metadata$Sample_ID

# Load counts  -------------------------------------------------------------

kallisto_res <- readRDS(paste0(RNA_seq_data_path, "kallisto/kallisto.merged.gene_counts.rds"))

counts <- assay(kallisto_res)
colnames(counts) <- sapply(colnames(counts), function(i) unlist(strsplit(i, "_"))[[1]])
counts <- counts[, rownames(metadata)]

kallisto_res <- read.table(paste0(RNA_seq_data_path, "kallisto/kallisto.merged.gene_tpm.tsv"))
colnames(kallisto_res) <- kallisto_res[1,]
kallisto_res <- kallisto_res[-1,]
tpm <- kallisto_res[, 3:ncol(kallisto_res)]
tpm <- apply(tpm, 2, as.numeric)
rownames(tpm) <- kallisto_res$gene_id
colnames(tpm) <- sapply(colnames(tpm), function(i) unlist(strsplit(i, "_"))[[1]])
tpm <- tpm[, rownames(metadata)]

gene_matching <- data.frame(row.names = kallisto_res$gene_id, ensembl_id = kallisto_res$gene_id, symbol = kallisto_res$gene_name)
gene_matching$nbcounts <- rowSums(counts)[rownames(gene_matching)]

# Get vst data from DESeq2 ------------------------------------------------------------
load(paste0(RNA_seq_data_path, "kallisto/deseq2_qc/deseq2.dds.RData"))
dds
vst <- assay(dds, "vst")
colnames(vst) <- sapply(colnames(vst), function(i) unlist(strsplit(i, "_"))[[1]])
vst <- vst[, metadata$samples]

# Get samples distance heatmap --------------------------------------------

descriminative_colors = c("KO" = "#ff7f00ff",
                          "WT" = "#d90017ff")
sampleDists <- stats::dist(t(vst))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- grDevices::colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

ha = rowAnnotation(Group = metadata$Group, 
                   annotation_name_rot = 45,
                   col = list(Group = descriminative_colors)
)

ha2 = columnAnnotation(Group = metadata$Group, 
                       annotation_name_rot = 45,
                       col = list(Group = descriminative_colors)
)

dist_plot <- Heatmap(sampleDistMatrix, name = "Distance", left_annotation = ha, col = colors, rect_gp = gpar(col = "black", lwd = 0.5))
dist_plot

# Remove WT5 --------------------------------------------------------------
counts_subset <- counts[-which(rowSums(counts) == 0),metadata$samples != "WT5SHI"]
metadata_subset <- metadata[metadata$samples != "WT5SHI",]


# Remove lowly expressed genes --------------------------------------------

NormFactor <- edgeR::calcNormFactors(object = counts_subset, method="TMM")
d0 <- DGEList(counts_subset, norm.factors = NormFactor)
keep_genes <- rowSums(edgeR::cpm(d0) >= 5) >= 2
counts_subset <- counts_subset[keep_genes, ]
dim(counts_subset)

to_save <- cbind(data.frame(gene_id  = rownames(counts_subset), symbol =  gene_matching[rownames(counts_subset), "symbol"]),tpm[rownames(counts_subset),])
to_save <- to_save %>%
  group_by(symbol) %>%
  summarise(across(where(is.numeric), sum))
write.table(to_save, file = paste0(output_path, "tpm_values.txt"), row.names = F, col.names = T, quote = F)



norm_data <- vst[rownames(counts_subset), colnames(vst) != "WT5SHI"]
rvdm = apply(norm_data,1,var)
selectdm = order(rvdm, decreasing = TRUE)[1:2000]
norm_data_subset = t(norm_data[selectdm,])
norm_data_subset[1:5,1:5]
dim(norm_data_subset)

# Run PCA and UMAP 
pca = prcomp(norm_data_subset, scale. = T)
print(summary(pca))
plot(pca)
pcaData = as.data.frame(pca$x)
loadings = pca$rotation
rownames(pcaData)=rownames(norm_data_subset)

to_plot=merge(pcaData, metadata, by.x="row.names",by.y="Sample_ID")
variance_explained <- summary(pca)
pdf(paste0(output_path,"/pca_plot.pdf"), w=5, h=5)
base_plot <- ggplot(to_plot, aes(x=PC1, y=PC2, fill=Group)) + 
  geom_point(size=5, color = "black", shape = 21, stroke = 1.2) + 
  scale_fill_manual(values = c("KO" = "#EE6677", "WT" = "#BBBBBB"), name = "Group" ) +
  scale_color_manual(values = c("KO" = "#EE6677", "WT" = "#BBBBBB"), name = "Group" ) +
  ggforce::geom_mark_ellipse(aes(fill = Group, color = Group), expand = unit(5, "mm")) +
  xlim(c(-50,50)) + ylim(c(-40,40))+
  theme_bw() + xlab(paste0("PC1 (", round(variance_explained$importance[2,"PC1"], 2)*100, "%)")) + ylab(paste0("PC2 (", round(variance_explained$importance[2,"PC2"], 2)*100, "%)")) +
  theme(axis.text=element_text(size=17,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=17,face="bold"),legend.background = element_rect(fill = "white", color = "black"))
base_plot + theme(legend.position = "none") #c(0.5, 0.5)
dev.off()

pdf(paste0(output_path,"/pca_plots_noWT5.pdf"), w=7, h=5)
ggplot(to_plot, aes(x=PC1, y=PC2, color=Group, shape = samples)) +
  geom_point(size=5) + scale_shape_manual(values = c(16,17,18,1,2,3,7,9)) +
  scale_color_manual(values = c("KO" = "#EE6677", "WT" = "darkgray")) + 
  theme_bw() +
  theme(axis.text=element_text(size=15,face="bold"),
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))+
  guides(col=guide_legend("Group",override.aes = list(size=4)))

ggplot(to_plot, aes(x=PC1, y=PC3, color=Group,shape = samples)) +
  geom_point(size=5) + scale_shape_manual(values = c(16,17,18,1,2,3,7,9)) +
  scale_color_manual(values = c("KO" = "#EE6677", "WT" = "darkgray")) + 
  theme_bw() +
  theme(axis.text=element_text(size=15,face="bold"),
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))+
  guides(col=guide_legend("Group",override.aes = list(size=4)))

ggplot(to_plot, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=5) + 
  scale_color_manual(values = c("KO" = "#EE6677", "WT" = "darkgray")) + 
  theme_bw() +
  theme(axis.text=element_text(size=15,face="bold"),
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))+
  guides(col=guide_legend("Group",override.aes = list(size=4)))

ggplot(to_plot, aes(x=PC1, y=PC3, color=Group)) +
  geom_point(size=5) + 
  scale_color_manual(values = c("KO" = "#EE6677", "WT" = "darkgray")) + 
  theme_bw() +
  theme(axis.text=element_text(size=15,face="bold"),
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold"),legend.background = element_rect(fill = "white", color = "black"))+
  guides(col=guide_legend("Group",override.aes = list(size=4)))
dev.off()


# Run DESeq2 --------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_subset),
  colData = metadata_subset,
  design = ~ Group 
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Group", "KO", "WT"))
head(res)
res$gene <- gene_matching[rownames(res), "symbol"]
sig_genes <- as.data.frame(res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 0 ), ])
head(sig_genes)
dim(sig_genes)
grep("Igl", sig_genes$gene, value = T)

openxlsx::write.xlsx(as.data.frame(res), file = paste0(output_path, "noWT5_diff_genes_all.xlsx"))


df_to_plot <- vst[rownames(sig_genes),-6]
rownames(df_to_plot) <- gene_matching[rownames(df_to_plot), "symbol"]
split.vector1 <- factor(metadata[colnames(vst[rownames(sig_genes),-6]), "Group"], levels = c("WT","KO"))


descriminative_colors = c("KO" = "#EE6677",
                          "WT" = "darkgray")
colours <- list('Group' = descriminative_colors)
colAnn <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(Group =  split.vector1),
                                            which = 'col',
                                            col = colours)
ht <-ComplexHeatmap::Heatmap(t(scale(t(df_to_plot))),
                             name = "vst", column_split = split.vector1,
                             top_annotation = colAnn,
                             column_title = NULL, cluster_columns = F, cluster_rows = T, clustering_method_rows = "ward.D",
                             show_column_names = T, show_row_names = T, #rect_gp = grid::gpar(col = "black", lwd = 0.5),
                             row_names_gp = grid::gpar(fontsize = 8))
pdf(paste0(output_path, "heatmap_DEGs_WT5removed.pdf"), w=7,h=25)
draw(ht)
dev.off()

openxlsx::write.xlsx(sig_genes, file = paste0(output_path, "noWT5_sig_DEGs.xlsx"))

