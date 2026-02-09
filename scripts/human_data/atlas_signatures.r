# library(BPCells)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(rstatix)
library(ggpubr)
library(cowplot)
set.seed(2025)


# Functions ---------------------------------------------------------------

get_facet_plot <- function(df, yvar, title, comparison_categ = "subtype"){
  pos_map <- df %>%
    distinct(disease, study) %>%
    group_by(disease) %>%
    arrange(study, .by_group = TRUE) %>%   # ensure alphabetical order
    mutate(xnum = row_number()) %>%
    ungroup()
  
  # tests: per disease AND per study, comparing subtypes
  stats <- df %>%
    group_by(disease, study) %>%
    pairwise_wilcox_test(reformulate(comparison_categ, response = yvar), p.adjust.method = "BH") %>% #comparisons = list(c("Sell+", "Ccl3-"), c("Ccl3-", "Ccl3+"), c("Sell+", "Ccl3+"))
    ungroup() %>%
    add_xy_position(x = "study") %>% 
    left_join(pos_map, by = c("disease","study")) %>%
    mutate(xmin = xnum + xmin - x, xmax = xnum + xmax - x) 
  
  stats$subtype = stats$group1 
  if(comparison_categ == "subtype"){
    stats$y.position = rep(c(max(df[,yvar])+0.2,
                             max(df[,yvar])+0.9, 
                             max(df[,yvar])+1.6),
                           length(unique(df$study)))
    extra_y = 2
    colors <- c("#1E88E5", "#FFC107", "#DA3C08")
  }else{
    stats$y.position = rep(c(max(df[,yvar])+0.2),
                           length(unique(df$study)))
    df$tissue <- factor(df$tissue, levels = c("blood", "tumor"))
    extra_y = 0.5
    colors <- c("#3288bdff","#d90017ff")
  }
  
  p <- ggplot(df, aes(x = study, y = get(yvar), color = get(comparison_categ))) +
    geom_boxplot(position = position_dodge(width = 0.8), size = 0.5) +
    ylim(c(min(df[,yvar])-0.3, max(df[,yvar])+extra_y))+
    scale_color_manual(values = colors) +
    facet_grid(~ disease, scales = "free_x", space = "free_x")+
    # facet_wrap(~ disease, scales = "free_x",ncol = length(unique(df$disease))) +
    stat_pvalue_manual(stats, label = "p.adj.signif", hide.ns = T, tip.length = 0.01, size = 4)+
    theme_bw(base_size = 14) + ylab("Gene Signature") + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
    theme(axis.text.x = element_text( size = 10)) + #angle = 25, vjust = 1, hjust=1,
    labs(title = title, y = "Signature", x = "Study", color = "Category")
  p
  
  
  
  
  # one comparison per disease
  stats <- df %>%
    group_by(disease) %>%
    pairwise_wilcox_test(reformulate(comparison_categ, response = yvar), p.adjust.method = "BH") %>% #comparisons = list(c("Sell+", "Ccl3-"), c("Ccl3-", "Ccl3+"), c("Sell+", "Ccl3+"))
    add_xy_position(x = comparison_categ)
  
  stats$subtype = stats$group1 
  if(comparison_categ == "subtype"){
    stats$y.position = rep(c(max(df[,yvar])+0.2,
                             max(df[,yvar])+0.9, 
                             max(df[,yvar])+1.6),
                           length(unique(df$disease)))
    extra_y = 2
    colors <- c("#1E88E5", "#FFC107", "#DA3C08")
  }else{
    stats$y.position = rep(c(max(df[,yvar])+0.2),
                           length(unique(df$disease)))
    df$tissue <- factor(df$tissue, levels = c("blood", "tumor"))
    extra_y = 0.5
    colors <- c("#3288bdff","#d90017ff")
  }
  
  p2 <- ggplot(df, aes(x = get(comparison_categ), y = get(yvar), color = get(comparison_categ))) +
    geom_boxplot(position = position_dodge(width = 0.8), size = 0.5) +
    ylim(c(min(df[,yvar])-0.3, max(df[,yvar])+extra_y))+
    scale_color_manual(values = colors) +
    facet_grid(~ disease, scales = "free_x", space = "free_x")+
    stat_pvalue_manual(stats, label = "p.adj.signif", hide.ns = T, tip.length = 0.01, size = 4)+
    theme_bw(base_size = 14) + ylab("Gene Signature") + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
    theme(axis.text.x = element_text( size = 10)) + #angle = 25, vjust = 1, hjust=1,
    labs(title = title, y = "Signature", x = "", color = "Category")
  p2
  return(list(p,p2))
}
###########################################################################
###########################################################################
###                                                                     ###
###                             MOUSE ATLAS                             ###
###                                                                     ###
###########################################################################
###########################################################################


data_path <- as.character(commandArgs(TRUE)[1])
cluster_peaks_path <- as.character(commandArgs(TRUE)[2])
output_folder <- as.character(commandArgs(TRUE)[3])
dar_res_path <- as.character(commandArgs(TRUE)[4])

dir.create(paste0(output_folder, "/mouse_atlas/"), recursive = T)
dir.create(paste0(output_folder, "/human_atlas/"), recursive = T)
# convert mouse to human genes 
genes_info <- read.table(paste0(data_path,"HOM_MouseHumanSequence.rpt"), header = T, sep = "\t")
df_mouse <- genes_info %>%
  filter(Common.Organism.Name == "mouse, laboratory") %>%
  dplyr::select(DB.Class.Key, mouse_symbol = Symbol)
df_human <- genes_info %>%
  filter(Common.Organism.Name == "human") %>%
  dplyr::select(DB.Class.Key, human_symbol = Symbol)
df_map <- df_mouse %>%
  dplyr::inner_join(df_human, by = "DB.Class.Key")
df_map

data <- readRDS(paste0(data_path, "/mouse_expression.rds"))
metadata <- readRDS(paste0(data_path, "/mouse_annotation.rds"))


seurat.obj <- CreateSeuratObject(counts = data, data = data, meta.data = metadata)
seurat.obj$subtype <- seurat.obj$Ntype
seurat.obj$study <- sapply(colnames(seurat.obj), function(i) unlist(strsplit(i, "@"))[1])
table(seurat.obj$study)

head(seurat.obj@meta.data)

seurat.obj$group <- paste0(seurat.obj$study, "_", seurat.obj$tissue)
averaged_data_perStudy_perTissue <- Seurat::AverageExpression(object = seurat.obj, assays = "RNA", layer = "data", 
                                                              group.by = "group")$RNA

seurat.obj$group <- paste0(seurat.obj$study, "_", seurat.obj$subtype)
averaged_data_perStudy_perSubtype <- Seurat::AverageExpression(object = seurat.obj, assays = "RNA", layer = "data", 
                                                               group.by = "group")$RNA

# Select cluster nearest genes
cluster_peaks <- read.table(cluster_peaks_path, header = F)

# get nearest genes
peak_gr <- Signac::StringToGRanges(cluster_peaks$V4, sep = c(":", "-")) 
txdb <-  TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
anno = ChIPseeker::annotatePeak(
  peak_gr,
  TxDb = txdb,
  # tssRegion=c(-1000, 1000),
  addFlankGeneInfo = TRUE,
  overlap = "all", #any nearest gene is reported regardless of the overlap with the TSS
  annoDb = "org.Hs.eg.db"
)

anno_df = as.data.frame(anno)
head(anno_df)

anno_df$annotation_type = anno_df$annotation
anno_df$annotation_type[grep("Intron",anno_df$annotation_type)]="Intron"
anno_df$annotation_type[grep("Exon",anno_df$annotation_type)]="Exon"
anno_df$annotation_type[grep("Promoter",anno_df$annotation_type)]="Promoter"
anno_df$annotation_type[grep("Downstream",anno_df$annotation_type)]="Downstream"
anno_df$ID <- paste0(anno_df$seqnames, ":", anno_df$start, "-", anno_df$end)

diff_res <- read.table(dar_res_path, header = T, sep = "\t", quote = "")
colnames(diff_res)[1] = "IDs"
rownames(diff_res) <- diff_res$IDs
summary(cluster_peaks$V4 %in% diff_res$IDs)

if(grepl("Cup", cluster_peaks_path)){
  diff_res <- diff_res[diff_res$Normal.vs..Tumor.Log2.Fold.Change > 0 & diff_res$Normal.vs..Tumor.adj..p.value < 0.05,] 
  diff_res <- diff_res[order(abs(diff_res$Normal.vs..Tumor.Log2.Fold.Change), decreasing = T),]
  top_peaks <- diff_res$IDs[1:500]
}else{
  diff_res <- diff_res[diff_res$Normal.vs..Tumor.Log2.Fold.Change < 0 & diff_res$Normal.vs..Tumor.adj..p.value < 0.05,] 
  diff_res <- diff_res[order(abs(diff_res$Normal.vs..Tumor.Log2.Fold.Change), decreasing = T),]
  top_peaks <- diff_res$IDs[1:500]
}
anno_df <- anno_df[anno_df$ID %in% top_peaks,]
# Compute signature from the cluster peaks -------------------------------------
my_theme <- theme(
  axis.title.x = element_text(size = 16),  # x-axis label size
  axis.title.y = element_text(size = 16),  # y-axis label size
  axis.text.x  = element_text(size = 14),  # x tick labels
  axis.text.y  = element_text(size = 14),   # y tick labels
  strip.text = element_text(size = 16)  # increase facet title size
)

for(w in c(2000,10000,50000,100000)){
  
  cluster_genes <- anno_df$SYMBOL[abs(anno_df$distanceToTSS) < w]
  cluster_genes <- unique(df_map[df_map$human_symbol %in% cluster_genes, "mouse_symbol"])
  cluster_genes <- cluster_genes[ cluster_genes %in% rownames(seurat.obj)]
  length(cluster_genes)
  
  seurat.obj <- AddModuleScore(object = seurat.obj, features = list("cluster_closest_genes" = cluster_genes), name = "cluster", ctrl = 200)
  seurat.obj$cluster_closest_genes <- scale(seurat.obj$cluster1)
  

  seurat.obj <- AddModuleScore(seurat.obj, features = list("Random_genes" = sample(rownames(seurat.obj), length(cluster_genes), replace=F)), 
                               name = "Random", ctrl = 200)
  seurat.obj$Random_genes <- scale(seurat.obj$Random1)

  feature_data <- FetchData(seurat.obj, vars = c("cluster_closest_genes","Random_genes", "subtype", "tissue","disease", "study"))
  feature_data$subtype <- factor(feature_data$subtype, levels = c("Sell+","Ccl3-", "Ccl3+"))
  feature_data$tissue <- factor(feature_data$tissue, levels = unique(feature_data$tissue))
  feature_data$disease <- factor(feature_data$disease, levels = unique(feature_data$disease))
  feature_data$cell <- rownames(feature_data)
  feature_data <- feature_data[order(feature_data$disease),]
  write.table(feature_data[,c(7,1:6)], file = paste0(output_folder, "/mouse_atlas/TSSdist", w,"_signatures_perCell.txt"), quote = F, row.names = F, col.names = T)
  

  study_infos <- as.data.frame(feature_data %>% group_by(study, disease) %>% count())
  rownames(study_infos) <- study_infos$study
  as.data.frame(feature_data %>% group_by(study, disease, subtype) %>% count())
  
  feature_data$study <- factor(feature_data$study, levels = rownames(study_infos)[order(study_infos$disease)])

  studies_with_blood <- unique(seurat.obj$study[seurat.obj$tissue == "blood"])
  pdf(paste0(output_folder, "/mouse_atlas/TSSdist", w,"_signatures_byTissue.pdf"), w=4, h=5)
  tissue_plot <- get_facet_plot(df = feature_data[which(feature_data$study %in% studies_with_blood),], 
                                yvar = "cluster_closest_genes", title = "DAR-associated genes signature",
                                comparison_categ = "tissue")[[1]]
  print(tissue_plot)
  print(get_facet_plot(df = feature_data[which(feature_data$study %in% studies_with_blood),], 
                       yvar = "Random_genes", title = "Random signature",
                       comparison_categ = "tissue")[[1]])
  dev.off()

  
  pdf(paste0(output_folder, "/mouse_atlas/TSSdist", w,"_signatures_bySubtype.pdf"), w=10, h=5)
  print(get_facet_plot(df = feature_data, yvar = "cluster_closest_genes", title = "DAR-associated genes signature")[[1]])
  print(get_facet_plot(df = feature_data, yvar = "Random_genes", title = "Random signature")[[1]])
  dev.off()
  
  pdf(paste0(output_folder, "/mouse_atlas/TSSdist", w,"_signatures_bySubtype2.pdf"), w=10, h=5)
  subtype_plot <- get_facet_plot(df = feature_data, yvar = "cluster_closest_genes", title = "DAR-associated genes signature")[[2]]
  print(subtype_plot)
  print(get_facet_plot(df = feature_data, yvar = "Random_genes", title = "Random signature")[[2]])
  dev.off()
  
  pdf(paste0(output_folder, "/mouse_atlas/TSSdist", w,"_signatures_combined.pdf"), w=14, h=5)
  print(plot_grid(tissue_plot + my_theme, subtype_plot + my_theme, ncol = 2, rel_widths = c(0.25, 0.7)))
  dev.off()
}


###########################################################################
###########################################################################
###                                                                     ###
###                             HUMAN ATLAS                             ###
###                                                                     ###
###########################################################################
###########################################################################


data <- readRDS(paste0(data_path, "/human_expression.rds"))
metadata <- readRDS(paste0(data_path, "/human_annotation.rds"))

seurat.obj <- CreateSeuratObject(counts = data, data = data, meta.data = metadata)
seurat.obj$subtype <- seurat.obj$Ntype
seurat.obj$study <- sapply(colnames(seurat.obj), function(i) unlist(strsplit(i, "@"))[1])
table(seurat.obj$study)

seurat.obj <- seurat.obj[, !is.na(seurat.obj$subtype)]
head(seurat.obj@meta.data)

seurat.obj$group <- paste0(seurat.obj$study, "_", seurat.obj$tissue)
averaged_data_perStudy_perTissue <- Seurat::AverageExpression(object = seurat.obj, assays = "RNA", layer = "data", 
                                                              group.by = "group")$RNA

seurat.obj$group <- paste0(seurat.obj$study, "_", seurat.obj$subtype)
averaged_data_perStudy_perSubtype <- Seurat::AverageExpression(object = seurat.obj, assays = "RNA", layer = "data", 
                                                               group.by = "group")$RNA

# Select cluster nearest genes
# get nearest genes
peak_gr <- Signac::StringToGRanges(cluster_peaks$V4, sep = c(":", "-")) 
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
anno = ChIPseeker::annotatePeak(
  peak_gr,
  TxDb = txdb,
  addFlankGeneInfo = TRUE,
  overlap = "all", 
  annoDb = "org.Hs.eg.db"
)

anno_df = as.data.frame(anno)
head(anno_df)

anno_df$annotation_type = anno_df$annotation
anno_df$annotation_type[grep("Intron",anno_df$annotation_type)]="Intron"
anno_df$annotation_type[grep("Exon",anno_df$annotation_type)]="Exon"
anno_df$annotation_type[grep("Promoter",anno_df$annotation_type)]="Promoter"
anno_df$annotation_type[grep("Downstream",anno_df$annotation_type)]="Downstream"
anno_df$ID <- paste0(anno_df$seqnames, ":", anno_df$start, "-", anno_df$end)

diff_res <- read.table(dar_res_path, header = T, sep = "\t", quote = "")
colnames(diff_res)[1] = "IDs"
rownames(diff_res) <- diff_res$IDs
summary(cluster_peaks$V4 %in% diff_res$IDs)

if(grepl("Cup", cluster_peaks_path)){
  diff_res <- diff_res[diff_res$Normal.vs..Tumor.Log2.Fold.Change > 0 & diff_res$Normal.vs..Tumor.adj..p.value < 0.05,] 
  diff_res <- diff_res[order(abs(diff_res$Normal.vs..Tumor.Log2.Fold.Change), decreasing = T),]
  top_peaks <- diff_res$IDs[1:500]
}else{
  diff_res <- diff_res[diff_res$Normal.vs..Tumor.Log2.Fold.Change < 0 & diff_res$Normal.vs..Tumor.adj..p.value < 0.05,] 
  diff_res <- diff_res[order(abs(diff_res$Normal.vs..Tumor.Log2.Fold.Change), decreasing = T),]
  top_peaks <- diff_res$IDs[1:500]
}
anno_df <- anno_df[anno_df$ID %in% top_peaks,]
# Compute signature from the cluster peaks -------------------------------------


for(w in c(2000,10000,50000,100000)){
  cluster_genes <- anno_df$SYMBOL[abs(anno_df$distanceToTSS) < w]
  cluster_genes <- cluster_genes[ cluster_genes %in% rownames(seurat.obj)]
  length(cluster_genes)
  
  seurat.obj <- AddModuleScore(seurat.obj, features = list("cluster_closest_genes" = cluster_genes), name = "cluster")
  seurat.obj$cluster_closest_genes <- scale(seurat.obj$cluster1)
  
  seurat.obj <- AddModuleScore(seurat.obj, features = list("Random_genes" = sample(rownames(seurat.obj), length(cluster_genes), replace=F)), name = "Random")
  seurat.obj$Random_genes <- scale(seurat.obj$Random1)
  
  
  feature_data <- FetchData(seurat.obj, vars = c("cluster_closest_genes","Random_genes", "subtype", "tissue","disease", "study"))
  feature_data$subtype <- factor(feature_data$subtype, levels = c("SELL+","CCL3-", "CCL3+"))
  feature_data$tissue <- factor(feature_data$tissue, levels = unique(feature_data$tissue))
  feature_data$disease <- factor(feature_data$disease, levels = unique(feature_data$disease))
  feature_data$cell <- rownames(feature_data)
  feature_data <- feature_data[order(feature_data$disease),]
  
  write.table(feature_data[,c(7,1:6)], file = paste0(output_folder, "/human_atlas/TSSdist", w,"_signatures_perCell.txt"), quote = F, row.names = F, col.names = T)
  
  
  study_infos <- as.data.frame(feature_data %>% group_by(study, disease) %>% count())
  rownames(study_infos) <- study_infos$study
  as.data.frame(feature_data %>% group_by(study, disease, subtype) %>% count())
  
  feature_data$study <- factor(feature_data$study, levels = rownames(study_infos)[order(study_infos$disease)])
  
  studies_with_blood <- unique(seurat.obj$study[seurat.obj$tissue == "blood"])
  pdf(paste0(output_folder, "/human_atlas/TSSdist", w,"_signatures_byTissue.pdf"), w=6, h=5)
  tissue_plot <- get_facet_plot(df = feature_data[which(feature_data$study %in% studies_with_blood),], 
                                yvar = "cluster_closest_genes", title = "DAR-associated genes signature",
                                comparison_categ = "tissue")[[1]]
  print(tissue_plot)
  print(get_facet_plot(df = feature_data[which(feature_data$study %in% studies_with_blood),], 
                       yvar = "Random_genes", title = "Random signature",
                       comparison_categ = "tissue")[[1]])
  dev.off()

  pdf(paste0(output_folder, "/human_atlas/TSSdist", w,"_signatures_bySubtype.pdf"), w=10, h=5)
  subtype_plot <- get_facet_plot(df = feature_data, yvar = "cluster_closest_genes", title = "DAR-associated genes signature")[[1]]
  print(subtype_plot)
  print(get_facet_plot(df = feature_data, yvar = "Random_genes", title = "Random signature")[[1]])
  dev.off()
  
  pdf(paste0(output_folder, "/human_atlas/TSSdist", w,"_signatures_bySubtype2.pdf"), w=10, h=5)
  print(get_facet_plot(df = feature_data, yvar = "cluster_closest_genes", title = "DAR-associated genes signature")[[2]])
  print(get_facet_plot(df = feature_data, yvar = "Random_genes", title = "Random signature")[[2]])
  dev.off()
  
  pdf(paste0(output_folder, "/human_atlas/TSSdist", w,"_signatures_combined.pdf"), w=15, h=5)
  print(plot_grid(tissue_plot + my_theme, subtype_plot + my_theme, ncol = 2, rel_widths = c(0.35, 0.65)))
  dev.off()
}
