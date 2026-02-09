library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(cowplot)
set.seed(2025)


library(decoupleR)
library(OmnipathR)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(limma)
library(rstatix)

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
    scale_color_manual(values = c("#1E88E5", "#FFC107", "#DA3C08")) +
    facet_grid(~ disease, scales = "free_x", space = "free_x")+
    # facet_wrap(~ disease, scales = "free_x",ncol = length(unique(df$disease))) +
    stat_pvalue_manual(stats, label = "p.adj.signif", hide.ns = T, tip.length = 0.01, size = 4)+
    theme_bw(base_size = 14) + ylab("Gene Signature") + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
    theme(axis.text.x = element_text( size = 10)) + #angle = 25, vjust = 1, hjust=1,
    labs(title = title, y = "TF activity", x = "Study", color = "Category")
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
    scale_color_manual(values = c("#1E88E5", "#FFC107", "#DA3C08")) +
    facet_grid(~ disease, scales = "free_x", space = "free_x")+
    stat_pvalue_manual(stats, label = "p.adj.signif", hide.ns = T, tip.length = 0.01, size = 4)+
    theme_bw(base_size = 14) + ylab("Gene Signature") + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
    theme(axis.text.x = element_text( size = 10)) + #angle = 25, vjust = 1, hjust=1,
    labs(title = title, y = "TF activity", x = "", color = "Category")
  p2
  return(list(p,p2))
}
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

data_path <- as.character(commandArgs(TRUE)[1])
output_folder <- as.character(commandArgs(TRUE)[3])

dir.create(paste0(output_path, "/mouse_atlas/"), recursive = T)
dir.create(paste0(output_path, "/human_atlas/"), recursive = T)


# Load mouse human genes info ---------------------------------------------

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
# df_map


# Load mouse atlas --------------------------------------------------------

data <- readRDS(paste0(data_path, "/mouse_expression.rds"))
metadata <- readRDS(paste0(data_path, "/mouse_annotation.rds"))

seurat.obj <- CreateSeuratObject(counts = data, data = data, meta.data = metadata)
seurat.obj$subtype <- seurat.obj$Ntype
seurat.obj$study <- sapply(colnames(seurat.obj), function(i) unlist(strsplit(i, "@"))[1])
table(seurat.obj$study)

metadata <- seurat.obj@meta.data


# Load omnipath db --------------------------------------------------------

net <- get_collectri(organism='mouse', split_complexes=FALSE)
unique(net[grep("Nfa",net$source),"source"])

map_mm <- uniprot_full_id_mapping_table(
  to       = "genesymbol",
  from     = "uniprot",
  organism = 10090,   
  reviewed = NULL     
)


# Run ulm
sample_acts <- run_ulm(mat=seurat.obj@assays$RNA$data, 
                       net=net, .source='source', .target='target', .mor='mor', minsize = 5)

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()
colnames(sample_acts_mat) <- sapply(colnames(sample_acts_mat), function(i){
  if(i %in% map_mm$From){
    return(unlist(as.vector(map_mm[map_mm$From == i, "To"]))[1])
  }else{
    return(i)
  }
})
to_save <- t(sample_acts_mat)
saveRDS(to_save,
        file = paste0(output_path, "/mouse_atlas/mouse_TF_scores.rds"))


# Get Pvalue matrix
pvalue_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'p_value') %>%
  column_to_rownames('condition') %>%
  as.matrix()
colnames(pvalue_mat) <- sapply(colnames(pvalue_mat), function(i){
  if(i %in% map_mm$From){
    return(unlist(as.vector(map_mm[map_mm$From == i, "To"]))[1])
  }else{
    return(i)
  }
})
to_save <- t(pvalue_mat)
saveRDS(to_save,
        file = paste0(output_path, "/mouse_atlas/mouse_TF_pvalues.rds"))


feature_data <- FetchData(seurat.obj, vars = c("subtype", "tissue","disease", "study"))
feature_data$subtype <- factor(feature_data$subtype, levels = c("Sell+","Ccl3-", "Ccl3+"))
feature_data$tissue <- factor(feature_data$tissue, levels = unique(feature_data$tissue))
feature_data$disease <- factor(feature_data$disease, levels = unique(feature_data$disease))
feature_data$cell <- rownames(feature_data)
for(gene in c("Nfatc1","Nfatc2", "Nfatc3", "Nfatc4", "Nfat5")){
  feature_data[[gene]] <- sample_acts_mat[,gene]
}
feature_data <- feature_data[order(feature_data$disease),]
write.table(feature_data, file = paste0(output_path, "/mouse_atlas/Nfat_TFs_scores.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


# Identify most variable TFs ----------------------------------------------


# Get top most variable TFs
n_tfs <- 100
tfs <- sample_acts %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)
tfs <- sapply(tfs, function(i){
  if(i %in% map_mm$From){
    return(unlist(as.vector(map_mm[map_mm$From == i, "To"]))[1])
  }else{
    return(i)
  }
})


descriminative_colors = c("Sell+" = "#1E88E5", "Ccl3-" = "#FFC107", "Ccl3+" = "#DA3C08")
col_an = HeatmapAnnotation(Group = feature_data$subtype,
                           annotation_name_rot = 45,
                           col = list(Group = descriminative_colors))


# Build heatmaps ----------------------------------------------------------

colors <- c("#1E88E5", "#FFC107", "#DA3C08")
names(colors) <- c("Sell+", "Ccl3-", "Ccl3+")

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


for(gene in c("Nfatc1","Nfatc2", "Nfatc3", "Nfatc4", "Nfat5")){ 
  metadata[[gene]] <- sample_acts_mat[,gene]
}

ht_list[["Nfat_TFs"]] <- Heatmap(
  t(metadata[order(metadata$tree_order), c("Nfatc1","Nfatc2", "Nfatc3", "Nfatc4", "Nfat5")]),
  name = "TFs activity",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T, show_column_names = FALSE,
  column_title = "Nfat_TFs", border = T,
  show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
)



# Combine heatmaps --------------------------------------------------------

pdf(paste0(output_path, "/mouse_atlas/mouse_metadata_heatmap.pdf"), h=5, w=5)
ht_combined = ht_list$subtype_1row %v% ht_list$disease %v% ht_list$tissue %v% ht_list$Nfat_TFs
draw(ht_combined)
dev.off()


# Compute signature from the cluster peaks -------------------------------------
my_theme <- theme(
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16),  
  axis.text.x  = element_text(size = 14),  
  axis.text.y  = element_text(size = 14),   
  strip.text = element_text(size = 16)  
)


study_infos <- as.data.frame(feature_data %>% group_by(study, disease) %>% count())
rownames(study_infos) <- study_infos$study
as.data.frame(feature_data %>% group_by(study, disease, subtype) %>% count())

feature_data$study <- factor(feature_data$study, levels = rownames(study_infos)[order(study_infos$disease)])

studies_with_blood <- unique(seurat.obj$study[seurat.obj$tissue == "blood"])
pdf(paste0(output_path, "/mouse_atlas/boxplots_Nfat_TF_activities_byTissue.pdf"), w=4, h=3)
for(tf in c("Nfatc1","Nfatc2", "Nfatc3", "Nfatc4", "Nfat5")){
  tissue_plots[[tf]] <- get_facet_plot(df = feature_data[which(feature_data$study %in% studies_with_blood),], 
                                       yvar = tf, title = paste0(tf, " TF activity"),
                                       comparison_categ = "tissue")[[1]]
  print(tissue_plots[[tf]])
}
dev.off()

pdf(paste0(output_path, "/mouse_atlas/boxplots_Nfat_TF_activities_bySubtype.pdf"), w=10, h=3)
subtype_plots <- list()
for(tf in c("Nfatc1","Nfatc2", "Nfatc3", "Nfatc4", "Nfat5")){
  subtype_plots[[tf]] <- get_facet_plot(df = feature_data, yvar = tf, title = paste0(tf, " TF activity"))[[2]]
  print(subtype_plots[[tf]])
}
dev.off()

pdf(paste0(output_path, "/mouse_atlas/boxplots_Nfat_TF_activities_combined.pdf"), w=14, h=3)
for(tf in c("Nfatc1","Nfatc2", "Nfatc3", "Nfatc4", "Nfat5")){
  print(plot_grid(tissue_plots[[tf]] + my_theme, subtype_plots[[tf]] + my_theme, ncol = 2, rel_widths = c(0.25, 0.7)))
}
dev.off()


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

metadata <- seurat.obj@meta.data


# Load omnipath db --------------------------------------------------------

net <- get_collectri(organism='human', split_complexes=FALSE)
unique(net[grep("NFA",net$source),"source"])
unique(net[grep("E2F",net$source),"source"])
unique(net[grep("STAT",net$source),"source"])

map_mm <- uniprot_full_id_mapping_table(
  to       = "genesymbol",
  from     = "uniprot",
  organism = 9606,   
  reviewed = NULL     
)


# Run ulm
sample_acts <- run_ulm(mat=seurat.obj@assays$RNA$data, 
                       net=net, .source='source', .target='target', .mor='mor', minsize = 5)

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()
colnames(sample_acts_mat) <- sapply(colnames(sample_acts_mat), function(i){
  if(i %in% map_mm$From){
    return(unlist(as.vector(map_mm[map_mm$From == i, "To"]))[1])
  }else{
    return(i)
  }
})
to_save <- t(sample_acts_mat)
saveRDS(to_save,
        file = paste0(output_path, "/human_atlas/human_TF_scores.rds"))


# Get Pvalue matrix
pvalue_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'p_value') %>%
  column_to_rownames('condition') %>%
  as.matrix()
colnames(pvalue_mat) <- sapply(colnames(pvalue_mat), function(i){
  if(i %in% map_mm$From){
    return(unlist(as.vector(map_mm[map_mm$From == i, "To"]))[1])
  }else{
    return(i)
  }
})
to_save <- t(pvalue_mat)
saveRDS(to_save,
        file = paste0(output_path, "/human_atlas/human_TF_pvalues.rds"))


feature_data <- FetchData(seurat.obj, vars = c("subtype", "tissue","disease", "study"))
feature_data$subtype <- factor(feature_data$subtype, levels = c("SELL+","CCL3-", "CCL3+"))
feature_data$tissue <- factor(feature_data$tissue, levels = unique(feature_data$tissue))
feature_data$disease <- factor(feature_data$disease, levels = unique(feature_data$disease))
feature_data$cell <- rownames(feature_data)
for(gene in c("NFATC1","NFATC2", "NFATC3", "NFATC4", "NFAT5", 
              "STAT1","STAT2","STAT3","STAT4","STAT5A", "STAT5B", 
              "JUN", "JUNB", "JUND", 
              "RUNX1", "RUNX2", "RUNX3",
              "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7", "SMAD9",
              "ETS1", "ETS2", "EGR1", "EGR2", "EGR3", "EGR4",
              "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7")){
  feature_data[[gene]] <- sample_acts_mat[,gene]
}
feature_data <- feature_data[order(feature_data$disease),]
write.table(feature_data, file = paste0(output_path, "/human_atlas/Nfat_TFs_scores.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


# Identify most variable TFs ----------------------------------------------

# Get top most variable TFs
n_tfs <- 100
tfs <- sample_acts %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)
tfs <- sapply(tfs, function(i){
  if(i %in% map_mm$From){
    return(unlist(as.vector(map_mm[map_mm$From == i, "To"]))[1])
  }else{
    return(i)
  }
})


descriminative_colors = c("SELL+" = "#1E88E5", "CCL3-" = "#FFC107", "CCL3+" = "#DA3C08")
col_an = HeatmapAnnotation(Group = feature_data$subtype,
                           annotation_name_rot = 45,
                           col = list(Group = descriminative_colors))


# Build heatmaps ----------------------------------------------------------

colors <- c("#1E88E5", "#FFC107", "#DA3C08")
names(colors) <- c("SELL+","CCL3-", "CCL3+")

ht_list <- list()
levels <- list("subtype" = c("SELL+","CCL3-", "CCL3+"), 
               "tissue" = c("blood", "tumor"), 
               "disease" = c("HNSCC", "NSCLC", "PDAC", "ESCC"),
               "study" = c("Bill2023", "Kurten2021", "Salcher2022", "Zilionis2019", "Wang2023","Zhang2021"))
titles <- list("subtype" = "Neutrophil state", 
               "tissue" = "Sample site", 
               "disease" = "Tumor type",
               "study" = "Study")
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
  column_title = "Nfat_TFs", border = T,
  show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
)


for(gene in c("NFATC1","NFATC2", "NFATC3", "NFATC4", "NFAT5")){ 
  metadata[[gene]] <- sample_acts_mat[,gene]
}

ht_list[["Nfat_TFs"]] <- Heatmap(
  t(metadata[order(metadata$tree_order), c("NFATC1","NFATC2", "NFATC3", "NFATC4", "NFAT5")]),
  name = "TFs activity",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T, show_column_names = FALSE,
  column_title = "Nfat_TFs", border = T,
  show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
)

all_TFs <- c("NFATC1","NFATC2", "NFATC3", "NFATC4", "NFAT5", 
  "STAT1","STAT2","STAT3","STAT4","STAT5A", "STAT5B", 
  "JUN", "JUNB", "JUND", 
  "RUNX1", "RUNX2", "RUNX3",
  "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7", "SMAD9",
  "ETS1", "ETS2", "EGR1", "EGR2", "EGR3", "EGR4",
  "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7")

for(gene in all_TFs){ 
  metadata[[gene]] <- sample_acts_mat[,gene]
}

ht_list[["all_TFs"]] <- Heatmap(
  t(metadata[order(metadata$tree_order), all_TFs]),
  name = "TFs activity",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T, show_column_names = FALSE,
  column_title = "Nfat_TFs", border = T,
  show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
)

ht_list[["Nfatc1"]] <- Heatmap(
  t(metadata[order(metadata$tree_order), "NFATC1"]),
  name = "TFs activity",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T, show_column_names = FALSE,
  column_title = "NFATC1", border = T,
  show_heatmap_legend = T, use_raster = T, raster_by_magick = F, raster_quality = 20
)


# Combine heatmaps --------------------------------------------------------

pdf(paste0(output_path, "/human_atlas/human_metadata_heatmap.pdf"), h=5, w=5)
ht_combined = ht_list$subtype_1row %v% ht_list$disease %v% ht_list$tissue %v% ht_list$Nfat_TFs
draw(ht_combined)
dev.off()


pdf(paste0(output_path, "/human_atlas/human_metadata_heatmap_NFATC1.pdf"), h=5, w=5)
ht_combined = ht_list$subtype_1row %v% ht_list$disease %v% ht_list$tissue %v% ht_list$Nfatc1
draw(ht_combined)
dev.off()

pdf(paste0(output_path, "/human_atlas/human_metadata_heatmap_allTFs.pdf"), h=10, w=5)
ht_combined = ht_list$subtype_1row %v% ht_list$disease %v% ht_list$tissue %v% ht_list$all_TFs
draw(ht_combined)
dev.off()


# Boxplots by study and tissue -------------------------------------
my_theme <- theme(
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16),  
  axis.text.x  = element_text(size = 14),  
  axis.text.y  = element_text(size = 14),  
  strip.text = element_text(size = 16)  
)


study_infos <- as.data.frame(feature_data %>% group_by(study, disease) %>% count())
rownames(study_infos) <- study_infos$study
as.data.frame(feature_data %>% group_by(study, disease, subtype) %>% count())

feature_data$study <- factor(feature_data$study, levels = rownames(study_infos)[order(study_infos$disease)])

studies_with_blood <- unique(seurat.obj$study[seurat.obj$tissue == "blood"])
pdf(paste0(output_path, "/human_atlas/boxplots_Nfat_TF_activities_byTissue.pdf"), w=4, h=5)
tissue_plots <- list()
for(tf in c("NFATC1","NFATC2", "NFATC3", "NFATC4", "NFAT5")){
  tissue_plots[[tf]] <- get_facet_plot(df = feature_data[which(feature_data$study %in% studies_with_blood),], 
                                yvar = tf, title = paste0(tf, " TF activity"),
                                comparison_categ = "tissue")[[1]]
  print(tissue_plots[[tf]])
}

dev.off()


pdf(paste0(output_path, "/human_atlas/boxplots_Nfat_TF_activities_bySubtype.pdf"), w=10, h=5)
subtype_plots <- list()
for(tf in c("NFATC1","NFATC2", "NFATC3", "NFATC4", "NFAT5")){
  subtype_plots[[tf]] <- get_facet_plot(df = feature_data, yvar = tf, title = paste0(tf, " TF activity"))[[1]]
  print(subtype_plots[[tf]])
}
dev.off()

pdf(paste0(output_path, "/human_atlas/boxplots_Nfat_TF_activities_combined.pdf"), w=15, h=5)
for(tf in c("NFATC1","NFATC2", "NFATC3", "NFATC4", "NFAT5")){
  print(plot_grid(tissue_plots[[tf]] + my_theme, subtype_plots[[tf]] + my_theme, ncol = 2, rel_widths = c(0.35, 0.65)))
}
dev.off()


