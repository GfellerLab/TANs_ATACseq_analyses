

# Libraries ---------------------------------------------------------------

library(circlize)
library(ComplexHeatmap)

homer_path <- as.character(commandArgs(TRUE)[1])
motif_annotation_path <- as.character(commandArgs(TRUE)[2])
output_path  <- homer_path

# Gather homer results ----------------------------------------------------

# read homer de novo results:
homer_enrichment_summary <- data.frame(cluster = vector(),
                                       description = vector(),
                                       motif_name = vector(),
                                       pval = vector(),
                                       target = vector(),
                                       background =vector())

entries <- list.files(homer_path,pattern = "C", full.names = TRUE, recursive = FALSE, all.files = FALSE)
info <- file.info(entries)
n_clusters <- sum(info$isdir)
clusters_names <- paste0("C",sapply(entries, function(i) unlist(strsplit(i,"/C"))[[2]]))
for(cluster in clusters_names){
  motifs_files <- list.files(paste0(homer_path, cluster,"/homerResults/"), pattern = ".motif", recursive = T)
  motifs_files <- motifs_files[grep("RV|similar", motifs_files, invert = T)] 
  for(file in motifs_files){
    motif = read.table(paste0(homer_path, cluster,"/homerResults/", file), header = F, sep = "\t", 
                       nrows = 1, check.names = F)
    motif_name <- unlist(strsplit(motif$V2,"BestGuess:"))[[2]] 
    stats <- unlist(strsplit(motif$V6,","))
    homer_enrichment_summary <- rbind(homer_enrichment_summary, 
                                      data.frame(cluster = cluster,
                                                 description = motif_name,
                                                 motif_name = unlist(strsplit(motif_name, "[/]"))[[1]],
                                                 pval = as.numeric(unlist(strsplit(stats[3], ":"))[[2]])  ,
                                                 target = stats[1],
                                                 background = stats[2])
    )
  } 
} 

homer_enrichment_summary_signif <- homer_enrichment_summary[which(homer_enrichment_summary$pval <= 1e-12),] 
write.table(homer_enrichment_summary, file = paste0(homer_path, "homer_enrichment_results.txt"),
            quote = F, row.names = F, col.names = T)

# Retrieve significant results --------------------------------------------

homer_res <- homer_enrichment_summary_signif
homer_res$significance <- -log10(homer_res$pval)

motif_annotation <- read.table(motif_annotation_path, header = T, sep = "\t")
homer_res$motifID <- sapply(homer_res$motif_name, function(i) strsplit(i,"\\(")[[1]][1])
if(grepl("mouse", motif_annotation_path)){
  homer_res$gene_name <- sapply(homer_res$motifID, function(i) motif_annotation$mouse_gene_symbol[motif_annotation$name == i])
}else{
  homer_res$gene_name <- sapply(homer_res$motifID, function(i) motif_annotation$human_gene_symbol[motif_annotation$name == i])
}

signif_tfs <- unique(homer_res[homer_res$significance >=12, "gene_name"])
length(signif_tfs)
signif_tfs <- vector()
for(c in clusters_names){
  signif_tfs <- c(signif_tfs, homer_res[homer_res$significance >=12 & homer_res$cluster == c, "gene_name"])
}
signif_tfs <- unique(signif_tfs)

signif_matrix <- matrix(0, nrow = length(signif_tfs), ncol = length(clusters_names), dimnames = list(signif_tfs, clusters_names))

for(i in 1:length(signif_tfs)){
  for(j in clusters_names){
    signif_matrix[signif_tfs[i], j] <- ifelse(j %in% homer_res[homer_res$gene_name == signif_tfs[i],"cluster"],
                                               homer_res[homer_res$gene_name == signif_tfs[i] & homer_res$cluster == j,"significance"],
                                               0)
  }
}


# Get plot ----------------------------------------------------------------

pdf(paste0(output_path, "all_denovo_TFs.pdf"))

col_fun = colorRamp2(c(10, 12, 30, 50), c("white","#F4BA58","orange","#A03124"))

ht <-Heatmap(signif_matrix, name = "-log10(pval)", 
             cluster_rows = F, col = col_fun,
             rect_gp = gpar(col = "black", lwd = 1),
             column_names_rot = 0, 
             column_names_side = "top", row_names_side = "left",
             show_column_names = T, column_title = NULL, cluster_columns = F,
             show_row_names = T, row_dend_gp = gpar(lwd = 0.5))
draw(ht)
dev.off()
