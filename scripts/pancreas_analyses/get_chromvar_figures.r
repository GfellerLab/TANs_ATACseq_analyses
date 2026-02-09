
# Libraries ---------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)
library(limma)
library(ggplot2)

# Functions ---------------------------------------------------------------


# Get heatmap with ChromVar most variable features 
get_TF_heatmap <- function(gene_list, genes_matrix, metadata, mat.name = "Deviation score", col_fun = NULL,
                           split.vector = NULL, colors.vector = descriminative_colors, row_clustering = T,clustering_method_rows = "ward.D",
                           target_annotation = NULL, colors_annotation=NULL, text_size=8, show_row_names = T, show_col_names = T, scale = F){
  gene_list <- gene_list[which(gene_list %in% rownames(genes_matrix))]
  genes_matrix <- genes_matrix[, order(metadata$groups)]
  split.vector <- factor(metadata[colnames(genes_matrix), "groups"], levels = c("Bone_marrow_IMM", "Bone_marrow_MAT", 
                                                                                "Spleen_IMM", "Spleen_MAT",
                                                                                "Blood_IMM", "Blood_MAT",
                                                                                "Tumor_IMM", "Tumor_MAT"))
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
                                   name = mat.name, column_split = split.vector, col = col_fun,
                                   top_annotation = colAnn, right_annotation = rowAnn,
                                   column_title = NULL, cluster_columns = F, cluster_rows = row_clustering, clustering_method_rows = clustering_method_rows,
                                   show_column_names = show_col_names, show_row_names = show_row_names, rect_gp = grid::gpar(col = "black", lwd = 0.5),
                                   row_names_gp = grid::gpar(fontsize = text_size))
    }else{
      ht <-ComplexHeatmap::Heatmap(t(scale(t(genes_matrix[gene_list , ]))),
                                   name = mat.name, column_split = split.vector,col = col_fun,
                                   top_annotation = colAnn, right_annotation = rowAnn,
                                   column_title = NULL, cluster_columns = F, cluster_rows = row_clustering, clustering_method_rows = clustering_method_rows,
                                   show_column_names = show_col_names, show_row_names = show_row_names, rect_gp = grid::gpar(col = "black", lwd = 0.5),
                                   row_names_gp = grid::gpar(fontsize = text_size))
    }
    
  } else{
    
    if(!scale){
      ht <-ComplexHeatmap::Heatmap(genes_matrix[gene_list , ],
                                   name = mat.name, column_split = split.vector,col = col_fun,
                                   top_annotation = colAnn,
                                   column_title = NULL, cluster_columns = F,cluster_rows = row_clustering, clustering_method_rows = clustering_method_rows,
                                   show_column_names = show_col_names, show_row_names = show_row_names, rect_gp = grid::gpar(col = "black", lwd = 0.5),
                                   row_names_gp = grid::gpar(fontsize = text_size))
    }else{
      ht <-ComplexHeatmap::Heatmap(t(scale(t(genes_matrix[gene_list , ]))),
                                   name = mat.name, column_split = split.vector,col = col_fun,
                                   top_annotation = colAnn,
                                   column_title = NULL, cluster_columns = F,cluster_rows = row_clustering, clustering_method_rows = clustering_method_rows,
                                   show_column_names = show_col_names, show_row_names = show_row_names, rect_gp = grid::gpar(col = "black", lwd = 0.5),
                                   row_names_gp = grid::gpar(fontsize = text_size))
    }
    
  }
  
  
  return(ht)
}



# Parameters --------------------------------------------------------------

dev_matrix_path <- as.character(commandArgs(TRUE)[1])
z_matrix_path <- as.character(commandArgs(TRUE)[2])
metadata_path <- as.character(commandArgs(TRUE)[3])
motif_annotation_path <- as.character(commandArgs(TRUE)[4])
output_path <- as.character(commandArgs(TRUE)[5])

# Load the deviations scores and TF annotations ----------------------------------------------
dev_matrix <- readRDS(dev_matrix_path)
z_matrix <- readRDS(z_matrix_path)

motifs_annot <- read.table(motif_annotation_path, sep = '\t', header = T)
rownames(motifs_annot) <- motifs_annot$mouse_gene_symbol 

# Scale the z_matrix to create chromvar_mat
chromvar_mat <- z_matrix

# Load metadata -----------------------------------------------------------
metadata <- read.table(metadata_path, header = T, sep = "\t")
rownames(metadata) <- metadata$sample_name
metadata <- metadata[colnames(chromvar_mat),] 

# Calculate the row-wise variance of chromvar_mat
rvdm = apply(chromvar_mat, 1, var)


colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
            '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
            '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000", '#E4C755', '#F7F398',
            '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
            '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
            '#968175', '#FF6347', '#4682B4', '#D2B48C', '#008080',
            '#D8BFD8', '#FF4500', '#DA70D6', '#EEE8AA', '#98FB98',
            '#AFEEEE', '#DB7093', '#FFEFD5', '#FFDAB9', '#CD853F')
descriminative_colors <- c("Blood"="#3288bdff",
                           "Tumor" = "#d90017ff",
                           "Bone marrow" = "#ffd92fff",
                           "Spleen" = "#70c4b4ff",
                           "IMM" = "darkgreen", "MAT" = "purple","Bone_marrow_IMM" = "#4DAF4A", "Spleen_IMM" = "#984EA3", "Blood_IMM"="#E6AB02",
                           "Bone_marrow_MAT" = "#1B9E77", "Spleen_MAT" = "#6A3D9A", "Blood_MAT"="#D95F02",
                           "Tumor_IMM" = "pink", "Tumor_MAT" = "red")


annot_colors <- as.list(colors[1:length(unique(motifs_annot$tfclass_family))])
names(annot_colors) <- unique(motifs_annot$tfclass_family)

# Get the top 200 most variable features from chromvar_mat
Top_chromVAR <- rownames(chromvar_mat[order(rvdm, decreasing = T)[1:300] ,])

# Heatmap with avreaged signal per TF family ------------------------------
split.vector <- factor(metadata[colnames(z_matrix), "groups"],
                       levels = c("Bone_marrow_IMM", "Bone_marrow_MAT", 
                                  "Spleen_IMM", "Spleen_MAT",
                                  "Blood_IMM", "Blood_MAT",
                                  "Tumor_IMM", "Tumor_MAT"))
colours <- list('Condition' = descriminative_colors)
colAnn <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(Condition =  split.vector),
                                            which = 'col',
                                            col = colours,
                                            show_annotation_name = c(Condition = FALSE) )


pdf(paste0(output_path, "chromVar_TFfamilies_heatmap.pdf"), w=10, h=7)
agg_chromvar_mat <- aggregate(. ~ factor(motifs_annot[rownames(chromvar_mat), "tfclass_family"]), 
                              data=as.data.frame(chromvar_mat), FUN=mean)
rownames(agg_chromvar_mat) <- agg_chromvar_mat[,1]
agg_chromvar_mat <- agg_chromvar_mat[,-1]

print(get_TF_heatmap(gene_list = rownames(agg_chromvar_mat),
                     genes_matrix = t(scale(t(agg_chromvar_mat))), metadata = metadata,
                     target_annotation = NULL, mat.name = "Stand. \n chromVAR \n score",
                     row_clustering = T, clustering_method_rows = "ward.D2",
                     show_col_names = F, text_size = 6))

col_fun = colorRamp2(c(-2, -1, 0 ,1, 2), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124"))
print(get_TF_heatmap(gene_list = rownames(agg_chromvar_mat),
                     genes_matrix = t(scale(t(agg_chromvar_mat))), metadata = metadata,
                     target_annotation = NULL, col_fun = col_fun, mat.name = "Stand. \n chromVAR \n score",
                     row_clustering = T, clustering_method_rows = "ward.D2",
                     show_col_names = F, text_size = 6))

print(get_TF_heatmap(gene_list = rownames(agg_chromvar_mat),
                     genes_matrix = (agg_chromvar_mat - min(agg_chromvar_mat)) / (max(agg_chromvar_mat) - min(agg_chromvar_mat)), metadata = metadata,
                     target_annotation = NULL,mat.name = "chromVAR \n score",
                     row_clustering = T, clustering_method_rows = "ward.D2",
                     show_col_names = F, text_size = 6))

col_fun = colorRamp2(c(0, 0.25, 0.5 ,0.75, 1), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124"))
print(get_TF_heatmap(gene_list = rownames(agg_chromvar_mat),
                     genes_matrix = (agg_chromvar_mat - min(agg_chromvar_mat)) / (max(agg_chromvar_mat) - min(agg_chromvar_mat)), metadata = metadata,
                     target_annotation = NULL, col_fun = col_fun,mat.name = "chromVAR \n score",
                     row_clustering = T, clustering_method_rows = "ward.D2",
                     show_col_names = F, text_size = 6))
dev.off()

write.csv2(cbind(data.frame(row_names = rownames(agg_chromvar_mat)),agg_chromvar_mat), row.names = F, file = paste0(output_path, "z_scores_perFamily.csv"))

pdf(paste0(output_path, "chromVar_top200.pdf"),w=6, h = 4)

mat <- chromvar_mat[Top_chromVAR,]
col_fun = colorRamp2(c(-2, -1, 0 ,1, 2), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124"))

indices <- which(rownames(mat) %in% c("Cebpbe", "Cebpd", "Cebpa", "Cebpb", "Ddit3", "Atf4",
                                      "Fos", "Jun", "Atf3", "Nfkb1", "Rel",
                                      "Mafg", "Maff", "Mafk",
                                      "Crem", "Creb1", "Xbp1", "Rara", "Rarb",
                                      "Foxo6", "Foxj2", "FOxo4", "Rfx1", "Irf1","Tfeb","Ctcf",
                                      "Elf1", "Elf3", "Nfyc", "Nfya", "Elk1", "E2f1", "E2f3",
                                      "Myc", "E2f4", "Nfatc1","Nfatc4", "Mef2a",
                                      "Runx3","Stat5a","Stat5b", "Nfatc3", "Bcl6", "Smad3"))
ha = rowAnnotation(foo = anno_mark(at = indices,
                                   labels = rownames(mat)[indices],link_gp = gpar(lwd = 0.5),
                                   labels_gp = gpar(fontsize = 6) ))
ht <-Heatmap(t(scale(t(mat))), name = "Stand. \n chromVAR \n score",
             right_annotation = ha, clustering_method_rows = "ward.D2", show_row_dend = T,
             column_split = split.vector, top_annotation = colAnn, col = col_fun,
             show_column_names = F, column_title = NULL, cluster_columns = F,
             show_row_names = F, row_dend_gp = gpar(lwd = 0.5))
draw(ht)

col_fun = colorRamp2(c(-40, -20, 0 ,20, 40), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124"))
ht <-Heatmap(mat, name = "chromVAR \n z-score",
             right_annotation = ha, clustering_method_rows = "ward.D2", show_row_dend = T,
             column_split = split.vector, top_annotation = colAnn,# col = col_fun,
             show_column_names = F, column_title = NULL, cluster_columns = F,
             show_row_names = F, row_dend_gp = gpar(lwd = 0.5))
draw(ht)
dev.off()


# Run differential analysis on chromVAR deviation scores -------------------

comp_list = list(c("KP_lung_SiglecF_high", "KP_lung_SiglecF_low"), c("KP_lung_SiglecF_high", "Healthy_lung"),
                 c("KP_lung_SiglecF_high", "KP_blood"), c("KP_lung_SiglecF_high", "Healthy_blood"),
                 c("KP_lung_SiglecF_low", "Healthy_lung"), c("KP_lung_SiglecF_low", "KP_blood"), c("KP_lung_SiglecF_low", "Healthy_blood"),
                 c("Healthy_lung", "KP_blood"), c("Healthy_lung", "Healthy_blood"),
                 c("Healthy_blood","KP_blood"))

chromVAR_summary <- data.frame(group1 = vector(),
                               group2 = vector(),
                               TF = vector(),
                               method = vector(),
                               significance = vector(),
                               mean_dev_difference = vector(),
                               direction = vector())

chrom_var_scores <- scale(z_matrix)
dimnames(chrom_var_scores) <- dimnames(z_matrix)

conditions <- metadata$groups
batch <- metadata$batch

groups <- unique(conditions)

mm <- model.matrix( ~ 0 + conditions + batch)
limma_fit <- lmFit(chrom_var_scores, mm)

print(head(coef(limma_fit)))
chromVAR_summary <- vector()
for(i in 1:length(comp_list)){
  
  cond1 = comp_list[[i]][1]
  cond2 = comp_list[[i]][2]
  print(cond1); print(cond2)
  
  contr <- makeContrasts(get(paste0("conditions", cond1)) - get(paste0("conditions", cond2)), levels = mm)
  tmp <- contrasts.fit(limma_fit, contr)
  tmp <- eBayes(tmp)
  
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  print(length(which(top.table$adj.P.Val < 0.01)))
  da_tests <- top.table[which(top.table$adj.P.Val < 0.05), ]
  da_tests$group1 <- rep(cond1, nrow(da_tests))
  da_tests$group2 <- rep(cond2, nrow(da_tests))
  da_tests$TF <- rownames(da_tests)
  da_tests$Mean_diff <- sapply(da_tests$TF, function(tf) mean(z_matrix[tf, metadata$groups == cond1]) - mean(z_matrix[tf, metadata$groups == cond2]))
  chromVAR_summary <- rbind(chromVAR_summary, da_tests)
}
head(chromVAR_summary)

# Keep high vs low comparison ---------------------------------------------

comp_df <- chromVAR_summary[chromVAR_summary$group1 == "KP_lung_SiglecF_high" & 
                              chromVAR_summary$group2 == "KP_lung_SiglecF_low", ]

comp_df$absMaxDeviation <- matrixStats::rowMaxs(abs(dev_matrix[rownames(comp_df),]))
comp_df$absMaxZ <- matrixStats::rowMaxs(abs(z_matrix[rownames(comp_df),]))
comp_df <- comp_df[order(comp_df$adj.P.Val),]

comp_df_subset <- comp_df[which(abs(comp_df$Mean_diff)>2 & comp_df$adj.P.Val<0.001 ),]  #& comp_df$absMaxDeviation>0.005
interesting_TFs <- c("Stat5a", "Stat5b", "Nfatc1", "Nfatc3", "Smad3", "Runx3", "Mef2d", "E2f4",
                     "Tfeb", "Ctcf",
                     "Bcl6", "Cebpb", "Jun", "Junb", "Fos", "Nfe2l2", "Bach1", "Maff", "Atf4", "Ddit3")
comp_df_subset <- comp_df[which(comp_df$TF %in% interesting_TFs),]

pdf(paste0(output_path, "chromVar_meanDiff-vs-maxDeviation.pdf"))
ggplot(comp_df, aes(absMaxZ, Mean_diff, colour=-log10(adj.P.Val), label=TF)) +
  geom_point() + scale_color_viridis_c(option="rocket",direction=-1) + theme_bw() +
  ggrepel::geom_text_repel(data=comp_df_subset,
                           min.segment.length=0, #nudge_y = 1,
                           nudge_x = -1,
                           # max.overlaps=30,
                           nudge_y=sign(comp_df_subset$Mean_diff)*2,
                           size=3) +
  labs(x="Maximum absolute chromVAR deviation", y="Mean_diff") +
  theme(legend.position = "bottom")
dev.off()
