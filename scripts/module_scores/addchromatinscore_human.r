library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GenomicFeatures)
library(EDASeq)
library(ChIPseeker)
library(msigdbr)
# Load required libraries
library(ggplot2)
library(readxl)
library(ggsignif)  
library(dplyr)
library(ggplot2)
library(tidyr)
library(limma)

library(Signac)
library(Seurat)
library(BiocParallel)
library(dplyr)
set.seed(2026)

# Functions ---------------------------------------------------------------
get_tss_info <- function(genes, mart_obj, window_size){
  tss_info <- getBM(attributes = c('external_gene_name', 'chromosome_name', 'transcription_start_site'),
                    filters = 'external_gene_name',
                    values = genes,
                    mart = mart_obj,
                    useCache   = FALSE)
  
  tss_info <- tss_info[!duplicated(tss_info$external_gene_name), ]
  tss_info <- tss_info[which(tss_info$chromosome_name != "MT"), ]
  tss_info$start <- tss_info$transcription_start_site - window_size
  tss_info$start <- ifelse(tss_info$start < 0, 0, tss_info$start)
  tss_info$end <- tss_info$transcription_start_site + window_size
  tss_info$chromosome_name <- paste0("chr", tss_info$chromosome_name)
  gr_ranges <- makeGRangesFromDataFrame(tss_info, keep.extra.columns = T)
  
  return(gr_ranges)
}
get_gene_body_info <- function(genes, mart_obj, window_size){
  gene_body_info <- getBM(attributes = c('external_gene_name',   # Gene symbol
                                         'chromosome_name',      # Chromosome name
                                         'start_position',       # Gene start position
                                         'end_position',
                                         'transcription_start_site'),        # tss
                          filters = 'external_gene_name',         # Filter based on gene symbols
                          values = genes,                         # Your list of genes
                          mart = mart_obj,useCache   = FALSE)
  gene_body_info <- gene_body_info[-which(duplicated(gene_body_info$external_gene_name)),]
  gene_body_info$diff1 <- gene_body_info$transcription_start_site-gene_body_info$start_position
  gene_body_info$diff2 <- gene_body_info$transcription_start_site-gene_body_info$end_position
  
  gene_body_info$start_position_updated <- gene_body_info$start_position
  gene_body_info$start_position_updated[gene_body_info$diff1 == 0] <- gene_body_info$start_position_updated[gene_body_info$diff1 == 0]  - window_size
  
  gene_body_info$end_position_updated <- gene_body_info$end_position
  gene_body_info$end_position_updated[gene_body_info$diff2 == 0] <- gene_body_info$end_position_updated[gene_body_info$diff2 == 0]  + window_size
  
  # test = gene_body_info[(gene_body_info$diff1 != 0 & gene_body_info$diff2 != 0),]
  # summary(test$transcription_start_site < test$end_position & test$transcription_start_site > test$start_position)
  
  gene_body_info <- gene_body_info[, c(1,2,8,9)]
  colnames(gene_body_info) <- c("gene", "chr", "start", "end")
  
  gene_body_info$chr <- paste0("chr", gene_body_info$chr)
  gr_ranges <- makeGRangesFromDataFrame(gene_body_info[,c(2:4,1)], keep.extra.columns = T)
  
  return(gr_ranges)
}


metadata_path <- as.character(commandArgs(TRUE)[1])
counts_path <- as.character(commandArgs(TRUE)[2])
output_path <- as.character(commandArgs(TRUE)[3])
senescence_gene_file <- as.character(commandArgs(TRUE)[4])
db_genesets <-  as.character(commandArgs(TRUE)[5])
mouse_human_map <-  as.character(commandArgs(TRUE)[6])


dir.create(output_path, recursive = T)

descriminative_colors = c("Blood"="#3288bdff",
                          "Tumor"="#d90017ff",
                          "Adjacent_lung"="#ffd92fff",
                          "batch_1" = "#a0451fff", "batch_2" = "black", "batch3" = "gray", "batch4" = "#984ea3ff", "batch5" = "yellow4")

# Read the raw counts -----------------------------------------------------

counts_data <- read.table(counts_path, header = T, sep = "\t", row.names = 1)
colnames(counts_data) <- gsub("[.]", "-",colnames(counts_data))

# Load samples metadata ---------------------------------------------------

metadata <- read.table(metadata_path, header = T, sep = "\t")
rownames(metadata) <- metadata$replicate

summary(colnames(counts_data) == metadata$replicate)



# Create ChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts  =  as(as.matrix(counts_data), "dgCMatrix"),
  sep     = c(":", "-"),
  genome  = "hg38",          
  fragments = NULL,          
  min.cells = 1,
  min.features = 1
)

# Create Seurat object with a ChromatinAssay
seurat.obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay  = "peaks",
  project = "ATAC"
)



# Load human genome data --------------------------------------------------
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# retrieve all genes info
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_genes <- genes(txdb)
gene_ids <- mcols(all_genes)$gene_id
gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
mcols(all_genes)$symbol <- gene_symbols


# Retrieve genesets -------------------------------------------------------
genesets_data = readRDS(db_genesets)
all_genesets = genesets_data[[1]]
gene_sets_origin = genesets_data[[2]]


# Map peaks to the genesets -----------------------------------------------
rowRanges <- Signac::StringToGRanges(rownames(counts_data), sep = c(":", "-"))
rowRanges$peak_id <- rownames(counts_data)

# # # match peaks to genesets 
window_size = 3000
tss_gr <- get_tss_info(all_genes$symbol, mart_obj = human, window_size = window_size )
length(unique(tss_gr$external_gene_name))

gene_body_gr <- get_gene_body_info(genes = all_genes$symbol, mart_obj = human, window_size = 5000)
length(unique(gene_body_gr$gene))

id_matching <- data.frame(entrez_ids = unique(unlist(all_genesets)),
                          symbols = mapIds(org.Hs.eg.db, keys = as.character(unique(unlist(all_genesets))), 
                                           column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
rownames(id_matching) <- id_matching$entrez_ids

total_genes <- id_matching[as.character(unique(unlist(all_genesets))), "symbols"]


overlapping_peaks <- findOverlaps(rowRanges,  tss_gr)
matching_peaks_tss <- data.frame(Genes = tss_gr[overlapping_peaks@to,]$external_gene_name,
                                 ATAC_regions = rowRanges[overlapping_peaks@from,]$peak_id)
matching_peaks_tss$both <- paste0(matching_peaks_tss$Genes, "|", matching_peaks_tss$ATAC_regions)
summary(unique(tss_gr$external_gene_name) %in% unique(matching_peaks_tss$Genes))


overlapping_peaks <- findOverlaps(rowRanges,  gene_body_gr)
matching_peaks_genebody <- data.frame(Genes = gene_body_gr[overlapping_peaks@to,]$gene, 
                                      ATAC_regions = rowRanges[overlapping_peaks@from,]$peak_id)
matching_peaks_genebody$both <- paste0(matching_peaks_genebody$Genes, "|", matching_peaks_genebody$ATAC_regions)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annot <- "org.Hs.eg.db"
anno <- ChIPseeker::annotatePeak(
  rowRanges,
  TxDb = txdb,
  addFlankGeneInfo = TRUE,
  overlap = "all", #any nearest gene is reported regardless of the overlap with the TSS
  annoDb = annot
)
anno_df = as.data.frame(anno)
anno_df = anno_df[abs(anno_df$distanceToTSS) < 100000,]
anno_df$peak <- paste0(anno_df$seqnames, ":",anno_df$start, "-", anno_df$end)

matching_peaks_nearestPeak <- anno_df[,c("peak", "SYMBOL", "distanceToTSS")]
colnames(matching_peaks_nearestPeak) <- c("ATAC_regions", "Genes", "distanceToTSS")
matching_peaks_nearestPeak <- matching_peaks_nearestPeak %>%
  group_by(Genes) %>%
  slice_min(order_by = abs(distanceToTSS), n = 1, with_ties = FALSE) %>%
  ungroup()

for(method in c("tss", "genebody", "nearestPeak")){
  all_genesets = genesets_data[[1]]
  gene_sets_origin = genesets_data[[2]]
  
  matching_peaks <- get(paste0("matching_peaks_", method))
  
  pathways_peaks <- list()
  df <- vector()
  for(i in 1:length(all_genesets)){
    print(i)
    genes <-  id_matching[as.character(unlist(all_genesets[i])), "symbols"]
    pathway_name = paste0(names(all_genesets[i]), " (", gene_sets_origin[i], ")")
    
    genes_with_match <- unique(matching_peaks[matching_peaks$Genes %in% genes, "Genes"])
    pct <- length(genes_with_match)/length(unique(genes))*100
    tmp_peaks <- unique(matching_peaks[matching_peaks$Genes %in% genes, "ATAC_regions"])
    pathways_peaks[[ pathway_name]] <- gsub(":", "-",tmp_peaks[!is.na(tmp_peaks)])
    df <- rbind(df, data.frame(Resource = gene_sets_origin[i], 
                               Pathway = pathway_name, 
                               GeneID = paste(genes, collapse = ","),
                               matching_pct = pct,
                               Matching_genes = paste(unique(matching_peaks[matching_peaks$Genes %in% genes, "Genes"]), collapse = ",")  )
    )
    
  }
  
  # convert mouse to human genes 
  genes_info <- read.table(mouse_human_map, header = T, sep = "\t")
  df_mouse <- genes_info %>%
    filter(Common.Organism.Name == "mouse, laboratory") %>%
    dplyr::select(DB.Class.Key, mouse_symbol = Symbol)
  df_human <- genes_info %>%
    filter(Common.Organism.Name == "human") %>%
    dplyr::select(DB.Class.Key, human_symbol = Symbol)
  df_map <- df_mouse %>%
    dplyr::inner_join(df_human, by = "DB.Class.Key")
  df_map
  
  Gungabeesoon <- list("Angiogenesis_GB" = c("Vegfa", "Snd1", "Mtdh", "Itga5", "Tnf", "Cxcl3", "Anxa3", "Hmgb1", "Hif1a", "Sema4d", "Lrg1", "Chil1"),
                       "Neutrophils_GB" = c("Abr", "Anxa3", "Cd177", "Itgam", "Itgb2", "Itgb2l", "Pikfyve", "Pram1", "Ptafr", "Spi1", "Stx11", "Syk"),
                       "Neutrophil_cytotoxicity_GB" = c("Ncf1", "Myd88", "Trem3", "Trem1", "Tusc2", "Cybb", "Cybc1", "Ncf2", "Ncf4", "Rac1", "Rac2"),
                       "Interferon_signaling_GB" = c("Adar", "Isg15", "Isg20", "Rsad2", "Ifit1", "Ifit3", "Ifitm1", "Ifitm3", "Irak1", "Oas3", "Stat1", "Stat2", "Irf7", "Cxcl10"),
                       "Myeloid_recruitment_GB" = c("Ccl3", "Mif", "Cxcl14", "Csf1", "Vegfa", "Ccl4", "Cxcl3"),
                       "ECM_remodeling_GB" = c("Adamdec1", "Ctsc", "Ctsb", "Rgcc", "Ctss", "Ctsz", "Adam17", "Adam10", "Adam8"),
                       "Tumor_proliferation_GB" = c("Tgfb1", "Tnf", "Il1a"),
                       "Immunosuppression_GB" = c("Havcr2", "Fcgr2b", "Il4ra", "Cd274", "Hif1a"))
  
  Gungabeesoon <- lapply(Gungabeesoon, function(x) df_map[df_map$mouse_symbol %in% x,"human_symbol"])
  all_genesets <- c(all_genesets, Gungabeesoon)
  gene_sets_origin <- c(gene_sets_origin,rep("Gungabeesoon", length(Gungabeesoon)))
  
  
  for(i in 1:length(Gungabeesoon)){
    print(i)
    genes <-  unlist(Gungabeesoon[i])
    pathway_name = paste0(names(Gungabeesoon)[i], " (Gungabeesoon)")
    
    genes_with_match <- unique(matching_peaks[matching_peaks$Genes %in% genes, "Genes"])
    pct <- length(genes_with_match)/length(unique(genes))*100
    tmp_peaks <- unique(matching_peaks[matching_peaks$Genes %in% genes, "ATAC_regions"])
    pathways_peaks[[ pathway_name]] <- gsub(":", "-",tmp_peaks[!is.na(tmp_peaks)])
    df <- rbind(df, data.frame(Resource = "Gungabeesoon", 
                               Pathway = pathway_name, 
                               GeneID = paste(genes, collapse = ","),
                               matching_pct = pct,
                               Matching_genes = paste(unique(matching_peaks[matching_peaks$Genes %in% genes, "Genes"]), collapse = ",")    )
    )
    
  }
  
  senescence_genes <- readxl::read_excel(senescence_gene_file,sheet = 1)
  table(senescence_genes$Classification)
  senescence_genes$Classification <- paste0(gsub(" ","-",senescence_genes$Classification), "_SEN" )
  
  senescence_genesets <- split(senescence_genes$`Gene(human)`, senescence_genes$Classification)
  senescence_genesets[["Senescence"]] <- senescence_genes$`Gene(human)`
  all_genesets <- c(all_genesets, senescence_genesets)
  gene_sets_origin <- c(gene_sets_origin, rep("senescence", length(senescence_genesets)))
  
  for(i in 1:length(senescence_genesets)){
    print(i)
    genes <-  unlist(senescence_genesets[i])
    pathway_name = paste0(names(senescence_genesets)[i], " (senescence)")
    
    genes_with_match <- unique(matching_peaks[matching_peaks$Genes %in% genes, "Genes"])
    pct <- length(genes_with_match)/length(unique(genes))*100
    tmp_peaks <- unique(matching_peaks[matching_peaks$Genes %in% genes, "ATAC_regions"])
    pathways_peaks[[ pathway_name]] <- gsub(":", "-",tmp_peaks[!is.na(tmp_peaks)])
    df <- rbind(df, data.frame(Resource = "senescence", 
                               Pathway = pathway_name, 
                               GeneID = paste(genes, collapse = ","),
                               matching_pct = pct,
                               Matching_genes = paste(unique(matching_peaks[matching_peaks$Genes %in% genes, "Genes"]), collapse = ",")  )
    )
    
  }
  
  register(MulticoreParam(workers = 3))    
  seurat.obj <- AddChromatinModule(seurat.obj, features = pathways_peaks, 
                                   genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                                   assay = "peaks", verbose = TRUE)
  
  
  saveRDS(pathways_peaks, file = paste0(output_path, method, "_peaks_sets.rds"))
  write.table(df, file = paste0(output_path, method, "_gene_peaks_association.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
  
  rownames(metadata) <- metadata$replicate
  
  metadata$Cell_type <- factor(
    metadata$groups,
    levels = c(
      "Blood" ,"Adjacent_lung","Tumor"  
    )
  )
  metadata$Sample_name <- rownames(metadata)
  # Merge module scores with metadata (transpose first to make samples rows)
  
  module_scores <- seurat.obj@meta.data[,4:ncol(seurat.obj@meta.data)]
  module_scores_save <- as.data.frame(t(module_scores))
  rownames(df) <- gsub("_", "-", df$Pathway)
  summary(rownames(module_scores_save) %in% rownames(df))
  names(pathways_peaks) <- gsub("_", "-", names(pathways_peaks))
  summary(rownames(module_scores_save) %in% names(pathways_peaks))
  
  module_scores_save$Pathway_genes <- sapply(rownames(module_scores_save), function(i) length(unlist(strsplit(df[i, "GeneID"], ","))))
  module_scores_save$Matching_genes <- sapply(rownames(module_scores_save), function(i) length(unlist(strsplit(df[i, "Matching_genes"], ","))))
  module_scores_save$matching_pct <- df[rownames(module_scores_save),"matching_pct"]
  module_scores_save$Pathway_regions <- sapply(rownames(module_scores_save), function(i) length(pathways_peaks[[i]]))
  module_scores_save$PathwayID <- rownames(module_scores_save)
  
  write.table(module_scores_save, file = paste0(output_path, method, "_samplesScores.txt"), 
              quote = F, row.names = F, col.names = T, sep = "\t")
  
  module_long <- as.data.frame(module_scores)
  summary(rownames(module_long) == rownames(metadata))
  module_long$Sample_name <- rownames(module_long)
  
  module_long <- module_long %>%
    left_join(metadata[, c("Sample_name", "Cell_type")], by = "Sample_name")
  
  # differential analysis on the signatures ---------------------------------
  conditions <- metadata$groups
  batch <- metadata$batch
  
  comp_list = list(c("Tumor", "Blood"), c("Tumor", "Adjacent_lung"), c("Adjacent_lung", "Blood"))
  groups <- unique(conditions)
  
  mm <- model.matrix( ~ 0 + conditions + batch)
  
  diff_summary <- vector()
  limma_fit <- lmFit(t(module_scores), mm)
  for(i in 1:length(comp_list)){
    
    cond1 = comp_list[[i]][1]
    cond2 = comp_list[[i]][2]
    print(cond1); print(cond2)
    
    contrast <- c(rep(0, length(groups)), rep(0, length(unique(batch))-1))
    names(contrast) <- c(sort(groups), paste0("batch", 2))
    
    contrast[cond1] <- 1
    contrast[cond2] <- -1
    contr <- makeContrasts(contrast, levels = colnames(coef(limma_fit)))
    tmp <- contrasts.fit(limma_fit, contr)
    tmp <- eBayes(tmp)
    
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    print(length(which(top.table$adj.P.Val < 0.05)))
    da_tests <- top.table
    da_tests$group1 <- rep(cond1, nrow(da_tests))
    da_tests$group2 <- rep(cond2, nrow(da_tests))
    da_tests$pathway <- rownames(da_tests)
    # da_tests$Mean_diff <- sapply(da_tests$pathway, function(j) mean(t(module_scores)[j, conditions == cond1]) - mean(t(module_scores)[j, conditions == cond2]))
    
    diff_summary <- rbind(diff_summary, da_tests)
  }
  head(diff_summary)
  
  diff_summary_subset <- diff_summary[diff_summary$adj.P.Val < 0.05, ]
  
  write.table(diff_summary, file = paste0(output_path, method, "_differential_analysis_res.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  
  
  
  pdf(paste0(output_path, method, "_boxplots.pdf" ))
  
  for(pathway_of_interest in c("inflammatory response (GOBP)",
                               "positive regulation of angiogenesis (GOBP)",
                               "canonical glycolysis (GOBP)", 
                               "negative regulation of apoptotic process (GOBP)",
                               "positive regulation of neutrophil migration (GOBP)",
                               "neutrophil degranulation (GOBP)",
                               "Angiogenesis-GB (Gungabeesoon)",
                               "Neutrophils-GB (Gungabeesoon)" ,
                               "Interferon-signaling-GB (Gungabeesoon)",
                               "Myeloid-recruitment-GB (Gungabeesoon)",
                               "ECM-remodeling-GB (Gungabeesoon)",
                               "Tumor-proliferation-GB (Gungabeesoon)",
                               "Immunosuppression-GB (Gungabeesoon)",
                               "Neutrophil-cytotoxicity-GB (Gungabeesoon)" ,
                               "positive regulation of inflammatory response (GOBP)", 
                               "regulation of angiogenesis (GOBP)",  
                               "angiogenesis (GOBP)",  
                               "positive regulation of response to oxidative stress (GOBP)",
                               "positive regulation of phagocytosis (GOBP)",
                               "regulation of type I interferon production (GOBP)",
                               "positive regulation of cytokine production involved in inflammatory response (GOBP)",
                               "positive regulation of acute inflammatory response (GOBP)",
                               "regulation of acute inflammatory response (GOBP)",   
                               "neutrophil chemotaxis (GOBP)",
                               "positive regulation of neutrophil chemotaxis (GOBP)",
                               "calcineurin-NFAT signaling cascade (GOBP)",
                               "positive regulation of calcineurin-NFAT signaling cascade (GOBP)"
                               # "ROS, RNS production in phagocytes (reactome)",
  )){
    module_plot <- module_long %>%
      select(Sample_name, Cell_type, all_of(pathway_of_interest))
    
    # Define cell type colors
    cell_type_colors <- c(
      "Blood" = "#3288BD",
      "Adjacent_lung" = "#FFD92F",
      "Tumor" = "#D90017"
    )
    
    # Define the pairwise comparisons you want to annotate
    comparisons <- list(
      c("Tumor", "Blood"),
      c("Tumor", "Adjacent_lung"),
      c("Adjacent_lung", "Blood")
    )
    
    # P-values from your external analysis (replace with your values)
    p_values <- unlist(lapply(comparisons, function(c) diff_summary$adj.P.Val[diff_summary$group1 == c[1] &
                                                                                diff_summary$group2 == c[2] & 
                                                                                diff_summary$pathway == pathway_of_interest]))
    
    # Convert p-values to labels (*, **, ***, ns)
    signif_labels <- ifelse(p_values < 0.0001, "****",
                            ifelse(p_values < 0.001, "***",
                                   ifelse(p_values < 0.01, "**",
                                          ifelse(p_values < 0.05, "*", "ns"))))
    
    # Choose automatic y-positions above current data range
    ymax <- max(module_plot[[pathway_of_interest]], na.rm = TRUE)
    
    # >>> Increased spacing between lines (0.12 instead of 0.06)
    step <- ifelse(is.finite(ymax) && ymax != 0, abs(ymax) * 0.10, 0.2) 
    
    y_positions <- ymax + seq(step, by = step, length.out = length(comparisons))
    # ------------------------------------------------------------------------------------
    
    # Box plot with border around axes and significance annotations
    print(ggplot(module_plot, aes(x = Cell_type, y = .data[[pathway_of_interest]], fill = Cell_type)) +
            geom_boxplot(alpha = 1) +
            scale_fill_manual(values = cell_type_colors) +
            labs(
              title = paste("module scores for", pathway_of_interest),
              x = "Cell type",
              y = "module score"
            ) +
            geom_signif(
              comparisons = comparisons,
              annotations = signif_labels,
              y_position = y_positions,
              tip_length = 0.02,
              textsize = 5,
              vjust = 0.5
            ) +
            # extra headroom so nothing is clipped or overlaps
            expand_limits(y = max(y_positions) + step * 1.2) +
            theme_minimal(base_size = 14) +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none",
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
            ))
    
    
  }
  
  dev.off()
}

