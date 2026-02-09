
# Libraries ---------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(dplyr)
library(circlize)


# Functions ---------------------------------------------------------------


get_indiv_heatmap <- function(df, col, levels, binary, row_order, title, colors = NULL){
  stopifnot(col %in% colnames(df))
  vals <- as.character(df[[col]])
  vals[is.na(vals)] <- "NA"
  lvls <- unique(vals)
  rids <- if (!is.null(rownames(df))) rownames(df) else as.character(seq_len(nrow(df)))
  m <- sapply(lvls, function(l) as.integer(vals == l))
  dimnames(m) <- list(rids, lvls)
  
  m <- m[row_order, levels, drop = FALSE]
  
  if(binary){
    ht <- Heatmap(
      t(m),
      name = title,
      col = colorRamp2(c(0, 1), c("white", "black")),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = T, show_column_names = FALSE,
      column_title = title,border = T, use_raster = T, raster_by_magick = F, raster_quality = 20
    )
  }else{
    for(i in 1:ncol(m)){
      m[,i] <- ifelse(m[,i] == 1, i, 0)
    }
    cols <- c('white',colors)
    names(cols) <- 0:ncol(m)
    ht <- Heatmap(
      t(m),
      name = title,
      col = cols,#colorRamp2(0:ncol(m), c("white", colors[levels])),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = T, show_column_names = FALSE,
      column_title = title,border = T,
      show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
    )
  }
  
  return(ht)
}

###########################################################################
###########################################################################
###                                                                     ###
###                             MOUSE ATLAS                             ###
###                                                                     ###
###########################################################################
###########################################################################

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


# Load atlas data ---------------------------------------------------------------

data <- readRDS(paste0(data_path, "/mouse_expression.rds"))
metadata <- readRDS(paste0(data_path, "/mouse_annotation.rds"))

seurat.obj <- CreateSeuratObject(counts = data, data = data, meta.data = metadata)
seurat.obj$subtype <- seurat.obj$Ntype
seurat.obj$study <- sapply(colnames(seurat.obj), function(i) unlist(strsplit(i, "@"))[1])
table(seurat.obj$study)

head(seurat.obj@meta.data)
metadata <- seurat.obj@meta.data

colors <- c("#1E88E5", "#FFC107", "#DA3C08")
names(colors) <- c("Sell+", "Ccl3-", "Ccl3+")


# Load signatures ---------------------------------------------------------

for(cluster in c("C2", "C5")){
  signature <- read.csv2(paste0(signature_path_KP,cluster,"_associated_genes/mouse_atlas/TSSdist1e+05_signatures_perCell.txt"), 
                         header = T, sep = " ")
  rownames(signature) <- signature$cell
  
  metadata[[cluster]] <- as.numeric(signature[rownames(metadata),"cluster_closest_genes"])
}

for(cluster in c("Cup", "Cdown")){
  signature <- read.csv2(paste0(signature_path_human,cluster,"_associated_genes/mouse_atlas/TSSdist1e+05_signatures_perCell.txt"), 
                         header = T, sep = " ")
  rownames(signature) <- signature$cell
  
  metadata[[cluster]] <- as.numeric(signature[rownames(metadata),"cluster_closest_genes"])
}

# Build heatmaps ----------------------------------------------------------

ht_list <- list()
levels <- list("subtype" = c("Sell+", "Ccl3-", "Ccl3+"), 
               "tissue" = c("blood", "tumor"), 
               "disease" = c("lung", "breast", "glioma", "HCC", "pancreas"))
titles <- list("subtype" = "Neutrophil state", 
               "tissue" = "Sample site", 
               "disease" = "Tumor type")
for(colname in c("subtype", "tissue", "disease")){
  if(colname == "subtype"){
    ht_list[[colname]] <- get_indiv_heatmap(df = metadata, col = colname, 
                                            levels = levels[[colname]],
                                            binary = F, 
                                            row_order = rownames(metadata)[order(metadata$tree_order)],
                                            title = titles[[colname]], colors = colors)
  }else{
    ht_list[[colname]] <- get_indiv_heatmap(df = metadata, col = colname, 
                                            levels = levels[[colname]],
                                            binary = T, 
                                            row_order = rownames(metadata)[order(metadata$tree_order)],
                                            title = titles[[colname]])
  }

}

ht_list[["subtype_1row"]] <- Heatmap(
  t(metadata[order(metadata$tree_order), c("subtype"), drop = F]),
  name = "Neutrophil state",
  col = colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T, show_column_names = FALSE,
  column_title = "Signatures", border = T,
  show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
)

ht_list[["Signatures_KP"]] <- Heatmap(
  t(metadata[order(metadata$tree_order), c("C2", "C5")]),
  name = "Signatures",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T, show_column_names = FALSE,
  column_title = "Signatures", border = T,
  show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
)

ht_list[["Signatures_human"]] <- Heatmap(
  t(metadata[order(metadata$tree_order), c("Cup", "Cdown")]),
  name = "Signatures",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T, show_column_names = FALSE,
  column_title = "Signatures", border = T,
  show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
)

# Combine heatmaps --------------------------------------------------------

pdf(paste0(output_folder_KPpeaks, "mouse_metadata_heatmap.pdf"), h=5, w=5)
ht_combined = ht_list$subtype_1row %v% ht_list$disease %v% ht_list$tissue %v% ht_list$Signatures_KP
draw(ht_combined)
dev.off()

pdf(paste0(output_folder_HUMANpeaks, "mouse_metadata_heatmap.pdf"), h=5, w=5)
ht_combined = ht_list$subtype_1row %v% ht_list$disease %v% ht_list$tissue %v% ht_list$Signatures_human
draw(ht_combined)
dev.off()

###########################################################################
###########################################################################
###                                                                     ###
###                             HUMAN ATLAS                             ###
###                                                                     ###
###########################################################################
###########################################################################


# Load atlas data ---------------------------------------------------------------

data <- readRDS(paste0(data_path, "/human_expression.rds"))
metadata <- readRDS(paste0(data_path, "/human_annotation.rds"))

seurat.obj <- CreateSeuratObject(counts = data, data = data, meta.data = metadata)
seurat.obj$subtype <- seurat.obj$Ntype
seurat.obj$study <- sapply(colnames(seurat.obj), function(i) unlist(strsplit(i, "@"))[1])
table(seurat.obj$study)

head(seurat.obj@meta.data)
metadata <- seurat.obj@meta.data

colors <- c("#1E88E5", "#FFC107", "#DA3C08")
names(colors) <- c("SELL+","CCL3-", "CCL3+")


# Load signatures ---------------------------------------------------------

for(cluster in c("C2", "C5")){
  signature <- read.csv2(paste0(signature_path_KP,cluster,"_associated_genes/human_atlas/TSSdist1e+05_signatures_perCell.txt"), 
                         header = T, sep = " ")
  rownames(signature) <- signature$cell
  
  metadata[[cluster]] <- as.numeric(signature[rownames(metadata),"cluster_closest_genes"])
}

for(cluster in c("Cup", "Cdown")){
  signature <- read.csv2(paste0(signature_path_human,cluster,"_associated_genes/human_atlas/TSSdist1e+05_signatures_perCell.txt"), 
                         header = T, sep = " ")
  rownames(signature) <- signature$cell
  
  metadata[[cluster]] <- as.numeric(signature[rownames(metadata),"cluster_closest_genes"])
}


# Build heatmaps ----------------------------------------------------------

ht_list <- list()
levels <- list("subtype" = c("SELL+","CCL3-", "CCL3+"), 
               "tissue" = c("blood", "tumor"), 
               "disease" = c("HNSCC", "NSCLC", "PDAC", "ESCC"),
               "study" = c("Salcher2022",  "Zilionis2019", "Bill2023", "Kurten2021",  "Wang2023","Zhang2021"))
titles <- list("subtype" = "Neutrophil state", 
               "tissue" = "Sample site", 
               "disease" = "Tumor type",
               "study" = "Study")
for(colname in c("subtype", "tissue", "disease","study")){
  if(colname == "subtype"){
    ht_list[[colname]] <- get_indiv_heatmap(df = metadata, col = colname, 
                                            levels = levels[[colname]],
                                            binary = F, 
                                            row_order = rownames(metadata)[order(metadata$tree_order)],
                                            title = titles[[colname]], colors = colors)
  }else{
    ht_list[[colname]] <- get_indiv_heatmap(df = metadata, col = colname, 
                                            levels = levels[[colname]],
                                            binary = T, 
                                            row_order = rownames(metadata)[order(metadata$tree_order)],
                                            title = titles[[colname]])
  }
  
}

ht_list[["subtype_1row"]] <- Heatmap(
  t(metadata[order(metadata$tree_order), c("subtype"), drop = F]),
  name = "Neutrophil state",
  col = colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T, show_column_names = FALSE,
  column_title = "Signatures", border = T,
  show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
)

ht_list[["Signatures_KP"]] <- Heatmap(
  t(metadata[order(metadata$tree_order), c("C2", "C5")]),
  name = "Signatures",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T, show_column_names = FALSE,
  column_title = "Signatures", border = T,
  show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
)

ht_list[["Signatures_human"]] <- Heatmap(
  t(metadata[order(metadata$tree_order), c("Cup", "Cdown")]),
  name = "Signatures",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T, show_column_names = FALSE,
  column_title = "Signatures", border = T,
  show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
)

# Combine heatmaps --------------------------------------------------------

pdf(paste0(output_folder_KPpeaks, "human_metadata_heatmap.pdf"), h=5, w=5)
ht_combined = ht_list$subtype_1row %v% ht_list$study %v% ht_list$tissue %v% ht_list$Signatures_KP
draw(ht_combined)
dev.off()

pdf(paste0(output_folder_HUMANpeaks, "human_metadata_heatmap.pdf"), h=5, w=5)
ht_combined = ht_list$subtype_1row %v% ht_list$study %v% ht_list$tissue %v% ht_list$Signatures_human
draw(ht_combined)
dev.off()
