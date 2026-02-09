
# Libraries ---------------------------------------------------------------

library(chipenrich)
library(ggplot2)

# Parameters --------------------------------------------------------------

peaks_file <- as.character(commandArgs(TRUE)[1])
output_path <- as.character(commandArgs(TRUE)[2])
geneset_param <- as.character(commandArgs(TRUE)[3])
nb_cpus <- as.numeric(commandArgs(TRUE)[4])
genome <- as.character(commandArgs(TRUE)[5])

dir.create(output_path, recursive = T)

# Load the peaks annotation -----------------------------------------------

peak_annotation_df <- read.table(peaks_file, 
                                 header = T, sep = "\t", check.names = F, quote = "")
peak_annotation_df$peak_id <- paste0(peak_annotation_df$seqnames, ":", peak_annotation_df$start, "-",peak_annotation_df$end)


# Run chipenrich ----------------------------------------------------------
locusDef = "nearest_tss"

anno_output = list()
anno_output_df = vector()
peaks_annotation_chipenrich = vector()
for(group in unique(peak_annotation_df$peak_cluster)){
  print(group)
  
  group_markers <- peak_annotation_df[which(peak_annotation_df$peak_cluster == group), "peak_id"]
  peaks_ids <- as.data.frame(peak_annotation_df[which(peak_annotation_df$peak_cluster == group),1:3])
  rownames(peaks_ids) <- group_markers
  colnames(peaks_ids) <- c("chr", "start", "end")
  
  results <- chipenrich::chipenrich(peaks = peaks_ids, genome = genome, genesets = geneset_param, locusdef = locusDef, 
                                    qc_plots = F, out_name = NULL, n_cores = nb_cpus)
  
  
  signif_pathways <- results$results#[which(results$results$FDR < 0.05),]
  peaks_assignment <- results$peaks
  peaks_assignment$peak_id <- paste0(peaks_assignment$chr, "-", peaks_assignment$peak_start, "-", peaks_assignment$peak_end)
  peaks_assignment$cell_type = group
  genes_names_matching <- unique(peaks_assignment[, c("gene_id", "gene_symbol")])
  rownames(genes_names_matching) <- genes_names_matching$gene_id 
  
  anno_output[[group]] = cbind(signif_pathways, data.frame(group = rep(group, nrow(signif_pathways))) )
  
  anno_output_df = rbind(anno_output_df, anno_output[[group]])
  
  peaks_annotation_chipenrich = rbind(peaks_annotation_chipenrich, peaks_assignment)
  
  gc()
}
write.table(anno_output_df,
            file = paste0(output_path, geneset_param, "_enrichment.txt"), 
            quote = T, row.names = F, col.names = T, sep = "\t", dec = ".",)


go_res_top = vector()
for(group in unique(anno_output_df$group)){
  subset = anno_output_df[which(anno_output_df$group == group),]
  subset = subset[order(subset$FDR, decreasing = F),]
  # subset = subset[]
  go_res_top = rbind(go_res_top, subset[1:min(20, nrow(subset)),] )

}
go_res_top$Significance <- -log10(go_res_top$FDR)
go_res_top$Description <- factor(go_res_top$Description, levels=unique((go_res_top$Description)[order(go_res_top$group)]))

pdf(paste0(output_path, geneset_param, "_pathways.pdf"), w = 15, h = 16)

p <- ggplot(go_res_top , aes(group, Description)) +
  geom_point(aes(size = Significance, color = Odds.Ratio), shape = 16, stroke = 1)
text_size=15
print(p + theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, size = text_size),
              axis.text.y = element_text( size = text_size),
              axis.title.y = element_text( size = text_size),
              text = element_text( size = text_size)))
dev.off()

go_res_top = vector()
for(group in unique(anno_output_df$group)){
  subset = anno_output_df[which(anno_output_df$group == group),]
  subset = subset[order(subset$FDR, decreasing = F),]
  go_res_top = rbind(go_res_top, subset[1:min(20, nrow(subset)),] )
  
}
go_res_top$Significance <- -log10(go_res_top$FDR)
go_res_top$Description <- factor(go_res_top$Description, levels=unique((go_res_top$Description)[order(go_res_top$group)]))

