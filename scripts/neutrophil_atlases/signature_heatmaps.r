# library(BPCells)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
set.seed(2025)

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
cluster_peaks <- read.table(cluster_peaks_path, header = T, sep = "\t", quote = "")

# get nearest genes
peak_gr <- Signac::StringToGRanges(paste0(cluster_peaks$seqnames, ":", cluster_peaks$start, "-", cluster_peaks$end), sep = c(":", "-")) 
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
anno = ChIPseeker::annotatePeak(
  peak_gr,
  TxDb = txdb,
  addFlankGeneInfo = TRUE,
  overlap = "all", 
  annoDb = "org.Mm.eg.db"
)

anno_df = as.data.frame(anno)
head(anno_df)

anno_df$annotation_type = anno_df$annotation
anno_df$annotation_type[grep("Intron",anno_df$annotation_type)]="Intron"
anno_df$annotation_type[grep("Exon",anno_df$annotation_type)]="Exon"
anno_df$annotation_type[grep("Promoter",anno_df$annotation_type)]="Promoter"
anno_df$annotation_type[grep("Downstream",anno_df$annotation_type)]="Downstream"
anno_df$cluster <- cluster_peaks$peak_cluster
  
# Compute signature from the cluster peaks -------------------------------------

for(w in c(2000,10000,50000,100000)){
  
  cluster_list <- list()
  for(i in 1:length(unique(anno_df$cluster))){
    cluster_genes <- anno_df$SYMBOL[abs(anno_df$distanceToTSS) < w & anno_df$cluster == paste0("C", i)]
    cluster_genes <- cluster_genes[ cluster_genes %in% rownames(seurat.obj)]
    length(cluster_genes)
    cluster_list[[paste0("C",i)]] <- cluster_genes
  }
  seurat.obj <- AddModuleScore(object = seurat.obj, features = cluster_list, name = "cluster", ctrl = 200)

  feature_data <- FetchData(seurat.obj, vars = c(paste0("cluster", 1:length(unique(anno_df$cluster))), "subtype", "tissue","disease", "study"))
  feature_data$subtype <- factor(feature_data$subtype, levels = c("Sell+","Ccl3-", "Ccl3+"))
  feature_data$tissue <- factor(feature_data$tissue, levels = unique(feature_data$tissue))
  feature_data$disease <- factor(feature_data$disease, levels = unique(feature_data$disease))
  feature_data$cell <- rownames(feature_data)

  library(ComplexHeatmap)
  library(circlize)
  library(tidyr)
  library(dplyr)
  library(tibble)
  
  cluster_means <- feature_data %>%
    group_by(subtype, study) %>%
    summarise(across(starts_with("cluster"), mean, .names = "mean_{.col}")) %>%
    ungroup()
  annotations <- cluster_means %>%
    dplyr::select(subtype, study)
  cluster_means <- cluster_means %>%
    unite("group", subtype, study, sep = "_")
  rownames(annotations) <- cluster_means$group
  mat <- cluster_means %>%
    column_to_rownames("group") %>%
    as.matrix()
  row_ha <- rowAnnotation(
    Subtype = annotations$subtype,
    Study = annotations$study,
    # col = list(
    #   Subtype = subtype_colors,
    #   Study = study_colors
    # ),
    show_annotation_name = TRUE
  )
  pdf(paste0(output_folder, "/mouse_atlas/TSSdist", w,"_heatmap.pdf"))
  print(Heatmap(mat,
          name = "Mean\nValue",
          col = colorRamp2(c(min(mat), 0.2, max(mat)), c("blue", "white", "red")),
          cluster_rows = FALSE,
          cluster_columns = TRUE,
          column_title = "Clusters",
          show_row_names = TRUE,
          left_annotation = row_ha))
  
  col_ha <- columnAnnotation(
    Subtype = annotations$subtype,
    Study = annotations$study,
    # col = list(
    #   Subtype = subtype_colors,
    #   Study = study_colors
    # ),
    show_annotation_name = TRUE
  )
  print(Heatmap(t(scale(t(t(mat)))),
                name = "Mean\nValue",
                # col = colorRamp2(c(min(mat), 0.1, max(mat)), c("blue", "white", "red")),
                cluster_rows = T,
                cluster_columns = F,
                column_title = "Clusters",
                show_row_names = TRUE, top_annotation = col_ha))
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

# get nearest genes
peak_gr <- Signac::StringToGRanges(paste0(cluster_peaks$seqnames, ":", cluster_peaks$start, "-", cluster_peaks$end), sep = c(":", "-")) 
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
anno = ChIPseeker::annotatePeak(
  peak_gr,
  TxDb = txdb,
  # tssRegion=c(-1000, 1000),
  addFlankGeneInfo = TRUE,
  overlap = "all", #any nearest gene is reported regardless of the overlap with the TSS
  annoDb = "org.Mm.eg.db"
)

anno_df = as.data.frame(anno)
head(anno_df)

anno_df$annotation_type = anno_df$annotation
anno_df$annotation_type[grep("Intron",anno_df$annotation_type)]="Intron"
anno_df$annotation_type[grep("Exon",anno_df$annotation_type)]="Exon"
anno_df$annotation_type[grep("Promoter",anno_df$annotation_type)]="Promoter"
anno_df$annotation_type[grep("Downstream",anno_df$annotation_type)]="Downstream"
anno_df$cluster <- cluster_peaks$peak_cluster

# Compute signature from the cluster peaks -------------------------------------


for(w in c(2000,10000,50000,100000)){

  cluster_list <- list()
  for(i in 1:length(unique(anno_df$cluster))){
    cluster_genes <- anno_df$SYMBOL[abs(anno_df$distanceToTSS) < w & anno_df$cluster == paste0("C", i)]
    cluster_genes <- unique(df_map[df_map$mouse_symbol %in% cluster_genes, "human_symbol"])
    cluster_genes <- cluster_genes[ cluster_genes %in% rownames(seurat.obj)]
    length(cluster_genes)
    cluster_list[[paste0("C",i)]] <- cluster_genes
  }
  seurat.obj <- AddModuleScore(object = seurat.obj, features = cluster_list, name = "cluster", ctrl = 200)
  
  feature_data <- FetchData(seurat.obj, vars = c(paste0("cluster", 1:length(unique(anno_df$cluster))), "subtype", "tissue","disease", "study"))
  feature_data$subtype <- factor(feature_data$subtype, levels = c("SELL+","CCL3-", "CCL3+"))
  feature_data$tissue <- factor(feature_data$tissue, levels = unique(feature_data$tissue))
  feature_data$disease <- factor(feature_data$disease, levels = unique(feature_data$disease))
  feature_data$cell <- rownames(feature_data)

  library(ComplexHeatmap)
  library(circlize)
  library(tidyr)
  library(dplyr)
  library(tibble)
  
  cluster_means <- feature_data %>%
    group_by(subtype, study) %>%
    summarise(across(starts_with("cluster"), mean, .names = "mean_{.col}")) %>%
    ungroup()
  annotations <- cluster_means %>%
    dplyr::select(subtype, study)
  cluster_means <- cluster_means %>%
    unite("group", subtype, study, sep = "_")
  rownames(annotations) <- cluster_means$group
  mat <- cluster_means %>%
    column_to_rownames("group") %>%
    as.matrix()
  row_ha <- rowAnnotation(
    Subtype = annotations$subtype,
    Study = annotations$study,
    show_annotation_name = TRUE
  )
  pdf(paste0(output_folder, "/human_atlas/TSSdist", w,"_heatmap.pdf"))
  print(Heatmap(mat,
          name = "Mean\nValue",
          col = colorRamp2(c(min(mat), 0.1, max(mat)), c("blue", "white", "red")),
          cluster_rows = FALSE,
          cluster_columns = TRUE,
          column_title = "Clusters",
          show_row_names = TRUE,
          left_annotation = row_ha))
  
  col_ha <- columnAnnotation(
    Subtype = annotations$subtype,
    Study = annotations$study,
    show_annotation_name = TRUE
  )
  print(Heatmap(t(scale(t(t(mat)))),
          name = "Mean\nValue",
          cluster_rows = T,
          cluster_columns = F,
          column_title = "Clusters",
          show_row_names = TRUE, top_annotation = col_ha))
  dev.off()
  

}
