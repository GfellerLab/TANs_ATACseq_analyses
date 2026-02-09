############################################################
## Load packages
############################################################

library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

set.seed(123)

############################################################
## User settings
############################################################


# DEG_path <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/data/processed/proteomics/GSEA/Gene_list_for_GSEA.xlsx"          
# gmt_file <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/data/processed/proteomics/GSEA/m5.go.bp.v2025.1.Mm.symbols.gmt"  
# pathway_path <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/input_data/mouse_genesets_symbols.txt"
# output_path <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/results_proteomics/"

DEG_path <- as.character(commandArgs(TRUE)[1])
pathway_path <- as.character(commandArgs(TRUE)[2])
output_path <- as.character(commandArgs(TRUE)[3])

# Column names in the ranked file
gene_col <- "Gene"       
stat_col <- "t.statistic.NFATC1 KO-hi-WT-hi"     



############################################################
## Load and prepare the ranked gene list
############################################################

# Read Excel file
df <- read_excel(DEG_path)

# Keep only the columns of interest
df <- df %>%
  select(all_of(c(gene_col, stat_col)))

# Coerce statistic column to numeric
# (non-numeric values become NA)
df[[stat_col]] <- suppressWarnings(as.numeric(df[[stat_col]]))

# Drop rows with missing gene or stat
df <- df %>%
  filter(!is.na(.data[[gene_col]]),
         !is.na(.data[[stat_col]]))

# Optional: check how many rows remain
message("Rows after cleaning: ", nrow(df))

# Handle duplicate gene symbols:
# keep the row with the most extreme statistic (largest |t|)
df <- df %>%
  group_by(.data[[gene_col]]) %>%
  slice_max(order_by = abs(.data[[stat_col]]),
            n = 1,
            with_ties = FALSE) %>%
  ungroup()

# Build named numeric vector sorted by decreasing statistic
# names = gene, values = t-stat
gene_stats <- df %>%
  arrange(desc(.data[[stat_col]])) %>%
  transmute(gene = .data[[gene_col]],
            stat = .data[[stat_col]]) %>%
  deframe()

# Quick sanity check
message("Number of unique genes in ranked list: ", length(gene_stats))
message("First few entries of gene_stats:")
print(head(gene_stats))


############################################################
## GSEA using clusterProfiler 
############################################################


gene_sets_df <- read.table(pathway_path, header = T)
table(gene_sets_df$database)

head(gene_sets_df)

gs_rdata_code = sprintf("geneset.%s.%s", "GOMF", "mmu")
data(list = gs_rdata_code, package = "chipenrich.data",  envir = environment())
GOMF_genesets <- as.list(geneset.GOMF.mmu@set.gene)
names(GOMF_genesets) <- lapply(names(GOMF_genesets), function(X) as.character(mget(X, geneset.GOMF.mmu@set.name, ifnotfound = NA)))
GOMF_genesets <- Filter(function(x) length(x) >= 10 & length(x) < 1000, GOMF_genesets)

gs_rdata_code = sprintf("geneset.%s.%s", "GOCC", "mmu")
data(list = gs_rdata_code, package = "chipenrich.data",  envir = environment())
GOCC_genesets <- as.list(geneset.GOCC.mmu@set.gene)
names(GOCC_genesets) <- lapply(names(GOCC_genesets), function(X) as.character(mget(X, geneset.GOCC.mmu@set.name, ifnotfound = NA)))
GOCC_genesets <- Filter(function(x) length(x) >= 10 & length(x) < 1000, GOCC_genesets)

all_genesets <- c(GOMF_genesets, GOCC_genesets)
gene_sets_origin <- c(rep("GOMF", length(GOMF_genesets)), 
                      rep("GOCC", length(GOCC_genesets)))
all_genesets_symbols <- lapply(all_genesets, function(X) clusterProfiler::bitr(X, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")$SYMBOL)

names(all_genesets_symbols) <- paste0(names(all_genesets_symbols)," (",gene_sets_origin,')')
pathway_summary <- data.frame(
  database = gene_sets_origin,
  Pathway = names(all_genesets_symbols),
  Pathway_length = sapply(all_genesets_symbols, function(x) length(x)),
  Genes = sapply(all_genesets_symbols, function(x) paste(x, collapse = ","))
)
gene_sets_df <- rbind(gene_sets_df, pathway_summary)
table(gene_sets_df$database)



gsea_signif_summary <- vector()
for(db in c("GOBP", "GOMF", "GOCC", "kegg")){
  gene_sets_df_subset <- gene_sets_df[gene_sets_df$database %in% db,]
  genesets <- setNames(
    lapply(strsplit(gene_sets_df_subset$Genes, ","), trimws),  # split and trim spaces
    gene_sets_df_subset$Pathway                               # name each element by pathway
  )
  
  term2gene <- tibble::enframe(genesets, name = "TERM", value = "GENE") %>%
    unnest(GENE)
  
  gsea_cp <- GSEA(
    geneList     = gene_stats,
    TERM2GENE    = term2gene,
    pvalueCutoff = 1,
    verbose      = TRUE
  )
  
  # Inspect results
  print(gsea_cp@result)
  
  if(db == "GOBP"){
    pathway_name <- 'inflammatory response (GOBP)'
    # Enrichment plot for this pathway
    p_cp <- enrichplot::gseaplot2(
      gsea_cp,
      geneSetID = pathway_name,
      title     = paste0("GSEA – ", pathway_name)
    )
    
    pdf( paste0(output_path, "GOBP_INFLAMMATORY_RESPONSE.pdf"))
    print(p_cp)
    dev.off()
  }
  
  res_cp <- gsea_cp@result
  res_cp
  
  res_cp %>%
    dplyr::select(
      ID,
      setSize,
      enrichmentScore,  # ES
      NES,              # Normalized ES
      pvalue,           # nominal P
      p.adjust,         # FDR (BH)
      qvalue           # FDR q-value
    )
  
  res_cp_signif <- res_cp[res_cp$p.adjust < 0.05,]
  res_cp_signif$db <- db
  gsea_signif_summary <- rbind(gsea_signif_summary, res_cp_signif)
}


openxlsx::write.xlsx(gsea_signif_summary, file = paste0(output_path, "gsea_res.xlsx"))



