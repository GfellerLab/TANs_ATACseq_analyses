
library(ComplexHeatmap)
library(circlize)
library(GenomicRanges)
library(ggplot2)
library(UpSetR)
library(dplyr)


called_peaks_path <- as.character(commandArgs(TRUE)[1])
output_path <- as.character(commandArgs(TRUE)[2])
peak_clustering <- as.character(commandArgs(TRUE)[3])
DAR_res <- as.character(commandArgs(TRUE)[4])
metadata_path <- as.character(commandArgs(TRUE)[5])
counts_path <- as.character(commandArgs(TRUE)[6])
fasta_file <- as.character(commandArgs(TRUE)[7])
homer_path <- as.character(commandArgs(TRUE)[8])
motif_annotation_path <- as.character(commandArgs(TRUE)[9])

# Functions ---------------------------------------------------------------

count_overlaps <- function(gr1, gr2) {
  overlaps <- findOverlaps(gr1, gr2, minoverlap = 0)
  return(length(unique(overlaps@from)))
}


find_unique_regions <- function(gr, other_gr_list) {
  # Combine all other GenomicRanges objects into one
  other_gr_combined <- unlist(as(other_gr_list, "GRangesList"))
  other_gr_combined <- reduce(other_gr_combined)
  
  overlaps <- findOverlaps(gr, other_gr_combined, minoverlap = 0)
  overlaping_pos <- unique(overlaps@from)
  
  # Return the number of unique regions
  return(length(gr) - length(unique(overlaps@from)))
}
find_overlaps <- function(gr1, gr2) {
  overlaps <- findOverlaps(gr1, gr2, minoverlap = 0)
  return(gr1[unique(overlaps@from)])
}

make_comb_mat_custom <- function (lt, mode, value_fun = length, top_n_sets = Inf, min_set_size = -Inf,
                                  universal_set = NULL, complement_size = NULL, set_on_rows = TRUE) 
{
  n = length(lt)
  nm = names(lt)
  lt = as.list(lt)
  
  union = getFromNamespace("union", ns = "BiocGenerics")
  intersect = getFromNamespace("intersect", ns = "BiocGenerics")
  setdiff = getFromNamespace("setdiff", ns = "BiocGenerics")
  
  set_size = sapply(lt, function(x) {
    value_fun(union(x, x[NULL]))
  })
  
  l = set_size >= min_set_size & rank(max(set_size) - set_size) <= top_n_sets
  set_size = set_size[l]
  lt = lt[l]
  n = length(lt)
  nm = names(lt)
  comb_mat = matrix(FALSE, nrow = n, ncol = sum(choose(n, 1:n)))
  rownames(comb_mat) = nm
  j = 1
  for (k in 1:n) {
    comb = combn(n, k)
    for (i in 1:ncol(comb)) {
      comb_mat[comb[, i], j] = TRUE
      j = j + 1
    }
  }
  get_comb_size = function(lt, mode, do = rep(TRUE, length(lt)), value_fun = length) {
    set1_index = which(do)
    set2_index = which(!do)
    if(length(set1_index) == 1){
      find_unique_regions(gr = lt[[set1_index]], other_gr_list = lt[set2_index])
    }else{
      s = lt[[set1_index[1]]]
      if (mode == "distinct") {
        for (i in set1_index[-1]) {
          s = intersect(s, lt[[i]])
        }
        for (i in set2_index) {
          s = setdiff(s, lt[[i]])
        }
      }else if (mode == "intersect") {
        for (i in set1_index[-1]) {
          s = find_overlaps(gr1 = s, gr2 = lt[[i]])
        }
      }
      length(s)
    }
    
  }
  comb_size = numeric(ncol(comb_mat))
  for (i in seq_len(ncol(comb_mat))) {
    comb_size[i] = get_comb_size(lt, mode = mode, do = comb_mat[, i], value_fun = value_fun)
  }
  comb_mat = comb_mat + 0
  l = comb_size > 0
  comb_mat = comb_mat[, l, drop = FALSE]
  comb_size = comb_size[l]
  if (!is.null(complement_size)) {
    comb_mat = cbind(rep(0, nrow(comb_mat)), comb_mat)
    comb_size = c(complement_size, comb_size)
  }
  attributes(comb_size) = NULL
  if (!set_on_rows) {
    comb_mat = t.default(comb_mat)
  }
  attr(comb_mat, "set_size") = unname(set_size)
  attr(comb_mat, "comb_size") = comb_size
  attr(comb_mat, "data") = lapply(lt, function(i) as.data.frame(i))
  param = list(mode = mode, value_fun = value_fun, universal_set = universal_set, 
               set_on_rows = set_on_rows)
  attr(comb_mat, "param") = param
  class(comb_mat) = c("comb_mat", "matrix")
  comb_mat = comb_mat[order.comb_mat(comb_mat)]
  return(comb_mat)
}


# Load called peaks in each condition -------------------------------------
all_grList <- readRDS(called_peaks_path)


# Get intersection sizes --------------------------------------------------

m_object <- make_comb_mat(all_grList)
m_subset <- m_object[comb_size(m_object) >= 500000]


intersection_size_df = as.data.frame(comb_size(m_object))
colnames(intersection_size_df) <- "Intersection_size"
intersection_size_df <- intersection_size_df[c("11000", "01100", "01010", "01001"),, drop = F]
intersection_size_df$Comparison <- c("Healthy_blood", "KP_blood", "Healthy_lung", "KP_lung_SiglecF_low")

pdf(paste0(output_path, "/peak_intersection_size.pdf"), w=8,h=7)
ggplot(data=intersection_size_df, aes(x=reorder(Comparison, -Intersection_size), y=Intersection_size)) +
  geom_bar(stat="identity") + theme_classic() + xlab("")+ ggtitle("Interstion between peaks called in KP_lung_SiglecF_high and in other conditions")
dev.off()


# Counts nb of overlaps ---------------------------------------------------

n <- length(all_grList)
overlap_matrix <- matrix(0, n, n, dimnames = list(names(all_grList), names(all_grList)))

for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    count <- count_overlaps(all_grList[[i]], all_grList[[j]])
    overlap_matrix[i, j] <- count
    overlap_matrix[j, i] <- count  # Symmetric matrix
  }
}

intersection_df <- data.frame(Nb_overlapping_peaks = overlap_matrix["KP_lung_SiglecF_high", c("Healthy_blood", "KP_blood", "Healthy_lung", "KP_lung_SiglecF_low")],#/length(all_grList$KP_lung_SiglecF_high), 
                              Comparison = c("Healthy_blood", "KP_blood", "Healthy_lung", "KP_lung_SiglecF_low"))

ggplot(data=intersection_df, aes(x=reorder(Comparison, -Nb_overlapping_peaks), y=Nb_overlapping_peaks)) +
  geom_bar(stat="identity") + theme_classic() + xlab("")+ ggtitle("Interstion between peaks called in KP_lung_SiglecF_high and in other conditions")

unique_regions <- sapply(1:length(all_grList), function(i) {
  find_unique_regions(gr = all_grList[[i]], other_gr_list = all_grList[-i])
})
names(unique_regions) <- names(all_grList)



df = data.frame(unique_regions = unique_regions, group = names(unique_regions))
df$prop <- df$unique_regions/sapply(1:length(all_grList), function(i) length(all_grList[[rownames(df)[i]]]))
pdf(paste0(output_path, "/unique_peaks.pdf"))
ggplot(data=df, aes(x=reorder(group, -prop), y=prop)) +
  geom_bar(stat="identity") + theme_classic() + xlab("")+ ylab("Proportions of unique regions") +
  geom_text(aes(label=unique_regions), position=position_dodge(width=0.9), vjust=-0.25)
ggplot(data=df, aes(x=reorder(group, -unique_regions), y=unique_regions)) +
  geom_bar(stat="identity") + theme_classic() + xlab("")+ ylab("Number of unique regions") +
  geom_text(aes(label=unique_regions), position=position_dodge(width=0.9), vjust=-0.25)
dev.off()


# Get nb overlapping peaks -------------------------------------------------


m_object <- make_comb_mat_custom(lt = all_grList,  mode = "intersect")

m_subset <- m_object[comb_name(m_object) %in% c("11111",
                                                "11000", "01100", "01010", "01001",
                                                "10001", "00101", "00011",
                                                "10010", "00110",
                                                "10100",
                                                "10000", "01000", "00100", "00010", "00001")]

pdf(paste0(output_path, "upset_plot.pdf"), w=10,h=5)
UpSet(m_subset, set_order = c("KP_lung_SiglecF_high", "KP_lung_SiglecF_low", "Healthy_lung", "Healthy_blood", "KP_blood" ),
      top_annotation = upset_top_annotation(m_subset, add_numbers = T),
      right_annotation = rowAnnotation(Called_peaks = anno_barplot(sapply(1:5, function(i) length(all_grList[[i]])),
                                                                   border = T, add_numbers = T, numbers_offset =  unit(-1, "cm")), 
                                       width = unit(3, "cm"))
)
dev.off()





# Peak annotations --------------------------------------------------------
anno_peak_clustering <- read.table(peak_clustering, header = T, sep = "\t",quote = "")

diff_res <- read.table(DAR_res, header = T, sep = "\t", quote = "")
colnames(diff_res)[1] = "IDs"

pval <- diff_res[, grepl("adj..p.value", colnames(diff_res))]
min_pval <- apply(pval, 1, min)
nonDARs <- which(min_pval > 0.05 )
  
peak_gr <- Signac::StringToGRanges(diff_res$IDs[nonDARs], sep = c(":", "-")) 

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

anno = ChIPseeker::annotatePeak(
  peak_gr,
  TxDb = txdb,
  addFlankGeneInfo = TRUE,
  overlap = "all", 
  annoDb = "org.Mm.eg.db"
)

anno_nonDARs = as.data.frame(anno)
head(anno_nonDARs)

anno_nonDARs$annotation_type = anno_nonDARs$annotation
anno_nonDARs$annotation_type[grep("Intron",anno_nonDARs$annotation_type)]="Intron"
anno_nonDARs$annotation_type[grep("Exon",anno_nonDARs$annotation_type)]="Exon"
anno_nonDARs$annotation_type[grep("Promoter",anno_nonDARs$annotation_type)]="Promoter"
anno_nonDARs$annotation_type[grep("Downstream",anno_nonDARs$annotation_type)]="Downstream"

table(anno_nonDARs$annotation_type)
anno_nonDARs$peak_cluster <- "non-DARs"
anno_all <- rbind(anno_peak_clustering, anno_nonDARs)

table(anno_all$peak_cluster)

proportions_table <- anno_all %>%
  group_by(peak_cluster, annotation_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(peak_cluster) %>%
  mutate(proportion = count / sum(count))

# View the table to ensure it looks correct
print(proportions_table)

pdf(paste0(output_path,"peak_clustering_annotations_withNonDARs.pdf"))

ggplot(proportions_table, aes(x = peak_cluster, y = proportion, fill = annotation_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "Proportion of Annotation Types by Peak Cluster",
       x = "Peak Cluster",
       y = "Proportion (%)",
       fill = "Annotation Type") + scale_fill_manual(values = c("Promoter"="indianred2","Downstream"="palegoldenrod", "3' UTR"="deepskyblue3", "Exon" = "green4", 
                                                                "Intron"="lightblue", "Distal Intergenic"="darkorange1", "5' UTR"="plum3", "thistle", "saddlebrown"))
dev.off()




# Peak cluster heatmap + summarized heatmaps ------------------------------

descriminative_colors = c("Healthy_blood"="#3288bdff",
                          "KP_lung_SiglecF_high"="#d90017ff",
                          "Healthy_lung"="#ffd92fff",
                          "KP_blood"="#70c4b4ff", 
                          "KP_lung_SiglecF_low" = "#ff7f00ff",
                          "batch1" = "#a0451fff", "batch2" = "black", "batch3" = "gray", "batch4" = "#984ea3ff", "batch5" = "yellow4")

# Read the raw counts 

counts_data <- read.table(counts_path, header = T, sep = "\t", row.names = 1)

# Load samples metadata 

metadata <- read.table(metadata_path, header = T, sep = "\t")
rownames(metadata) <- metadata$sample_name
metadata <- metadata[colnames(counts_data),] 

batch <- metadata$batch
condition <- metadata$groups

# Get diff_peaks 

diff_res <- read.table(DAR_res, header = T, sep = "\t", quote = "")
colnames(diff_res)[1] = "IDs"
summary(diff_res$IDs %in% rownames(counts_data))

pval <- diff_res[, grepl("adj..p.value", colnames(diff_res))]
min_pval <- apply(pval, 1, min)
logFC <- abs(diff_res[, grepl("Log2.Fold.Change", colnames(diff_res))])
max_logFC <- apply(logFC, 1, max)

diff_peaks <- which(min_pval < 0.05 )
length(unique(diff_res$IDs[diff_peaks]))

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
dataNorm <- gc_norm(input_matrix = as.matrix(counts_data), 
                    grouping_var = metadata$groups, 
                    fa_file = fasta_file, 
                    by_group = F)
norm_data <- log(dataNorm+1)
norm_data_subset <- norm_data[diff_res$IDs[diff_peaks] , ]

# Keep most variable peaks
rvdm = apply(norm_data_subset,1,var)
norm_data_subset = norm_data_subset[order(rvdm, decreasing = T)[1:5000] ,]

norm_data_subset <- norm_data_subset[paste0(anno_peak_clustering$seqnames, ":", anno_peak_clustering$start, "-", anno_peak_clustering$end),]
cluster_assignment <- anno_peak_clustering$peak_cluster
cl_colors <- c("#4DAF4A", "#984EA3", "#1B9E77","#6A3D9A",
               "#D95F02", "#7570B3", "#E7298A", "#66A61E",
               "#E6AB02", "#A6761D", "black", "gray")
names(cl_colors) <- paste0("C", 1:12)

groups_factor <- factor(metadata$groups, levels = c("Healthy_blood", "KP_blood", "Healthy_lung", 
                                                    "KP_lung_SiglecF_low", "KP_lung_SiglecF_high"))
sorted_indices <- order(groups_factor)
col_an = HeatmapAnnotation(Condition = metadata$groups[sorted_indices], 
                           annotation_name_rot = 45,
                           col = list(Condition = descriminative_colors
                           )
)
ha = rowAnnotation(Peak_cluster = cluster_assignment[order(cluster_assignment)] ,
                   col = list(Peak_cluster = cl_colors))


gene_names_to_add <- c("Cxcl2", "Thbs1", "Ptgs2", "Il1b", "Ccl3", "Vegfa", "Bcl2a1a", "Cdkn1a",
                       "Cxcr2", "Mmp9", "Trem1", "Gpr35", "Nod2", "Sell", "Dock2", "S100a9")
df <- norm_data_subset[order(cluster_assignment), sorted_indices]
rownames(anno_all) <- paste0(anno_all$seqnames, ":", anno_all$start, "-", anno_all$end)


panel_fun = function(index, nm) {
  grid.rect()
  grid.text(paste(gene_names_to_add[gene_names_to_add %in% anno_all[rownames(df)[index], "SYMBOL"]], 
                  collapse = "\n"), 0.5, 0.5,gp = gpar(fontface = "italic", fontsize = 10))
}
align_to = lapply(sort(unique(cluster_assignment)), function(x) which(cluster_assignment[order(cluster_assignment)] == x))
names(align_to) <- sort(unique(cluster_assignment))

pdf(paste0(output_path,"peak_clustering_heatmap_annotated.pdf"), w=8, h=10)

right_annotation <- rowAnnotation(foo = anno_block(
  align_to = align_to,
  panel_fun = panel_fun, 
  width = unit(3, "cm")
))
# col_fun = circlize::colorRamp2(c(3, 4, 5 ,6, 7), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124"))
Heatmap(df, show_row_names = F, show_column_names = F,
        left_annotation = ha, right_annotation = right_annotation, top_annotation = col_an, 
        show_row_dend = F, cluster_rows = F, cluster_columns = F)

dev.off()

col_an = HeatmapAnnotation(Condition = metadata$groups[sorted_indices], 
                           annotation_name_rot = 45,
                           col = list(Condition = descriminative_colors
                           )
)
data <- as.data.frame(norm_data_subset[order(cluster_assignment), sorted_indices])
df_avg <- data %>%
  mutate(Group = factor(cluster_assignment[order(cluster_assignment)])) %>%     
  group_by(Group) %>%              
  summarise(across(everything(), mean)) %>%  
  as.data.frame()

df_avg
rownames(df_avg) <- df_avg$Group


pdf(paste0(output_path,"peak_clustering_heatmap_averaged.pdf"), w=5, h=4)

col_fun = circlize::colorRamp2(quantile(unlist(df_avg)), c("#0000ffff", "#996efbff","white","#ff7a5aff","red"))
avg_heatmap <- Heatmap(as.matrix(df_avg[,2:ncol(df_avg)]), show_row_names = T, show_column_names = F, 
                       top_annotation = col_an,row_names_side = "left", 
                       row_names_gp = gpar(fontsize = 15),
                       heatmap_legend_param = list(direction = "horizontal"), col = col_fun,
                       show_row_dend = F, cluster_rows = F, cluster_columns = F, width = unit(3, "cm"))

print(avg_heatmap)

dev.off()


df_avg <- t(as.data.frame(t(df_avg[,2:ncol(df_avg)])) %>%
              mutate(Group = factor(metadata$groups[sorted_indices], levels = c("Healthy_blood", "KP_blood", "Healthy_lung", 
                                                                "KP_lung_SiglecF_low", "KP_lung_SiglecF_high"))) %>%     
              group_by(Group) %>%              
              summarise(across(everything(), mean)) %>%  
              as.data.frame())

colnames(df_avg) <- df_avg[1,]
df_avg <- as.data.frame(df_avg[-1,])
df_avg <- mutate_all(df_avg, function(x) as.numeric(as.character(x)))
col_an = HeatmapAnnotation(Condition = colnames(df_avg), 
                           annotation_name_rot = 45,
                           col = list(Condition = descriminative_colors)
)

pdf(paste0(output_path,"peak_clustering_heatmap_averaged.pdf"), w=5, h=4)

col_fun = circlize::colorRamp2(quantile(unlist(df_avg)), c("#1F5FA9", "#74C4EA","white","#F4BA58","#A03124"))
col_fun = circlize::colorRamp2(quantile(unlist(df_avg)), c("#0000ffff", "#996efbff","white","#ff7a5aff","red"))
avg_heatmap2 <- Heatmap(as.matrix(df_avg), show_row_names = T, show_column_names = F, #df_avg[,2:ncol(df_avg)]
                        top_annotation = col_an,row_names_side = "left", 
                        row_names_gp = gpar(fontsize = 15),
                        heatmap_legend_param = list(direction = "horizontal"), col = col_fun,
                        show_row_dend = F, cluster_rows = F, cluster_columns = F, width = unit(3, "cm"))

print(avg_heatmap2)

dev.off()

# Retrieve Homer results --------------------------------------------
homer_res <- read.table(homer_path, header = T)
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
clusters_names <- unique(homer_res$cluster)
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
signif_matrix <- t(signif_matrix)
signif_matrix <- signif_matrix[sort(rownames(signif_matrix)),]
# Get plot ----------------------------------------------------------------

pdf(paste0(output_path,"homer_res_heatmap.pdf"), 
    w=12, h=3)

col_fun = colorRamp2(c(min(signif_matrix), 12, 30, max(signif_matrix)), c("white","#F4BA58","orange","#A03124"))

ht <-Heatmap(signif_matrix, name = "-log10(pval)", 
             cluster_rows = F, col = col_fun,
             rect_gp = gpar(col = "black", lwd = 1),
             column_names_rot = 45, row_names_side = "left",
             heatmap_legend_param = list(direction = "horizontal"), 
             row_names_gp = gpar(fontface = "bold", fontsize = 15),
             column_names_gp = gpar(fontsize = 15),
             show_column_names = T, column_title = NULL, cluster_columns = F,
             show_row_names = T, row_dend_gp = gpar(lwd = 0.5))


draw(avg_heatmap + ht,  merge_legend = TRUE, heatmap_legend_side = "right")

dev.off()



pdf(paste0(output_path,"homer_res_heatmap2.pdf"), 
    w=12, h=3)

draw(avg_heatmap2 + ht,  merge_legend = TRUE, heatmap_legend_side = "right")

dev.off()

