# library(ChIPseeker)
library(msigdbr)
library(biomaRt)
library(dplyr)
library(tidyr)

# Parameters --------------------------------------------------------------
senescence_gene_file <- as.character(commandArgs(TRUE)[1])
mouse_orthologs <- as.character(commandArgs(TRUE)[2])
output_path <- as.character(commandArgs(TRUE)[3])

############################################################################
############################################################################
###                                                                      ###
###                            FOR MOUSE DATA                            ###
###                                                                      ###
############################################################################
############################################################################

# Retrieve genesets 

gs_rdata_code = sprintf("geneset.%s.%s", "GOBP", "mmu")
data(list = gs_rdata_code, package = "chipenrich.data",  envir = environment())

GOBP_genesets <- as.list(geneset.GOBP.mmu@set.gene)
names(GOBP_genesets) <- lapply(names(GOBP_genesets), function(X) as.character(mget(X, geneset.GOBP.mmu@set.name, ifnotfound = NA)))
summary(sapply(GOBP_genesets, function(X) length(X)))
GOBP_genesets <- Filter(function(x) length(x) >= 10 & length(x) < 1000, GOBP_genesets)

gs_rdata_code = sprintf("geneset.%s.%s", "reactome", "mmu")
data(list = gs_rdata_code, package = "chipenrich.data",  envir = environment())

reactome_genesets <- as.list(geneset.reactome.mmu@set.gene)
names(reactome_genesets) <- lapply(names(reactome_genesets), function(X) as.character(mget(X, geneset.reactome.mmu@set.name, ifnotfound = NA)))
summary(sapply(reactome_genesets, function(X) length(X)))
reactome_genesets <- Filter(function(x) length(x) >= 10 & length(x) < 1000, reactome_genesets)

gs_rdata_code = sprintf("geneset.%s.%s", "kegg_pathway", "mmu")
data(list = gs_rdata_code, package = "chipenrich.data",  envir = environment())

kegg_genesets <- as.list(geneset.kegg_pathway.mmu@set.gene)
names(kegg_genesets) <- lapply(names(kegg_genesets), function(X) as.character(mget(X, geneset.kegg_pathway.mmu@set.name, ifnotfound = NA)))
summary(sapply(kegg_genesets, function(X) length(X)))
kegg_genesets <- Filter(function(x) length(x) >= 10 & length(x) < 1000, kegg_genesets)


hallmarks_genesets = msigdbr(species = "mouse", category = "H")
head(hallmarks_genesets)
hallmarks_genesets <- split(hallmarks_genesets$entrez_gene, hallmarks_genesets$gs_name)
hallmarks_genesets <- Filter(function(x) length(x) >= 10 & length(x) < 1000, hallmarks_genesets)

all_genesets <- c(GOBP_genesets, reactome_genesets, hallmarks_genesets, kegg_genesets)
gene_sets_origin <- c(rep("GOBP", length(GOBP_genesets)), 
                      rep("reactome", length(reactome_genesets)), 
                      rep("hallmarks", length(hallmarks_genesets)),
                      rep("kegg", length(kegg_genesets)))
gene_sets_origin <- gene_sets_origin[-which(duplicated(names(all_genesets)))]
all_genesets <- all_genesets[-which(duplicated(names(all_genesets)))]

all_genesets_symbols <- lapply(all_genesets, function(X) clusterProfiler::bitr(X, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")$SYMBOL)


Gungabeesoon <- list("Angiogenesis" = c("Vegfa", "Snd1", "Mtdh", "Itga5", "Tnf", "Cxcl3", "Anxa3", "Hmgb1", "Hif1a", "Sema4d", "Lrg1", "Chil1"),
                     "Neutrophils" = c("Abr", "Anxa3", "Cd177", "Itgam", "Itgb2", "Itgb2l", "Pikfyve", "Pram1", "Ptafr", "Spi1", "Stx11", "Syk"),
                     "Neutrophil_cytotoxicity" = c("Ncf1", "Myd88", "Trem3", "Trem1", "Tusc2", "Cybb", "Cybc1", "Ncf2", "Ncf4", "Rac1", "Rac2"),
                     "Interferon_signaling" = c("Adar", "Isg15", "Isg20", "Rsad2", "Ifit1", "Ifit3", "Ifitm1", "Ifitm3", "Irak1", "Oas3", "Stat1", "Stat2", "Irf7", "Cxcl10"),
                     "Myeloid_recruitment" = c("Ccl3", "Mif", "Cxcl14", "Csf1", "Vegfa", "Ccl4", "Cxcl3"),
                     "ECM_remodeling" = c("Adamdec1", "Ctsc", "Ctsb", "Rgcc", "Ctss", "Ctsz", "Adam17", "Adam10", "Adam8"),
                     "Tumor_proliferation" = c("Tgfb1", "Tnf", "Il1a"),
                     "Immunosuppression" = c("Havcr2", "Fcgr2b", "Il4ra", "Cd274", "Hif1a"))
all_genesets_symbols <- c(all_genesets_symbols, Gungabeesoon)
gene_sets_origin <- c(gene_sets_origin, rep("Gungabeesoon", length(Gungabeesoon)))


senescence_genes <- readxl::read_excel(senescence_gene_file,sheet = 2)
table(senescence_genes$Classification)
senescence_genes$Classification <- paste0(gsub(" ","-",senescence_genes$Classification), "_SEN" )

senescence_genesets <- split(senescence_genes$`Gene(murine)`, senescence_genes$Classification)
all_genesets_symbols <- c(all_genesets_symbols, senescence_genesets)
gene_sets_origin <- c(gene_sets_origin, rep("senescence", length(senescence_genesets)))

all_genesets_symbols[["Senescence"]] <- senescence_genes$`Gene(murine)`
gene_sets_origin <- c(gene_sets_origin, "senescence")

names(all_genesets_symbols) <- paste0(names(all_genesets_symbols)," (",gene_sets_origin,')')

pathway_summary <- data.frame(
  database = gene_sets_origin,
  Pathway = names(all_genesets_symbols),
  Pathway_length = sapply(all_genesets_symbols, function(x) length(x)),
  Genes = sapply(all_genesets_symbols, function(x) paste(x, collapse = ","))
)

openxlsx::write.xlsx(pathway_summary, file = paste0(output_path, "/mouse_genesets_symbols.xlsx"))
write.table(pathway_summary, file = paste0(output_path, "/mouse_genesets_symbols.txt"), 
            quote=T, row.names = F, col.names = T, sep = "\t")

pathway_summary_IDs <- data.frame(
  database = gene_sets_origin,
  Pathway = names(all_genesets_symbols),
  Pathway_length = sapply(all_genesets_symbols, function(x) length(x)),
  Genes = sapply(all_genesets_symbols, function(x) paste(x, collapse = ","))
)

openxlsx::write.xlsx(pathway_summary, file = paste0(output_path, "/mouse_genesets_IDs.xlsx"))




############################################################################
############################################################################
###                                                                      ###
###                            FOR HUMAN DATA                            ###
###                                                                      ###
############################################################################
############################################################################

gs_rdata_code = sprintf("geneset.%s.%s", "GOBP", "hsa")
data(list = gs_rdata_code, package = "chipenrich.data",  envir = environment())

GOBP_genesets <- as.list(geneset.GOBP.hsa@set.gene)
names(GOBP_genesets) <- lapply(names(GOBP_genesets), function(X) as.character(mget(X, geneset.GOBP.hsa@set.name, ifnotfound = NA)))
summary(sapply(GOBP_genesets, function(X) length(X)))
GOBP_genesets <- Filter(function(x) length(x) >= 10 & length(x) < 1000, GOBP_genesets)

gs_rdata_code = sprintf("geneset.%s.%s", "reactome", "hsa")
data(list = gs_rdata_code, package = "chipenrich.data",  envir = environment())

reactome_genesets <- as.list(geneset.reactome.hsa@set.gene)
names(reactome_genesets) <- lapply(names(reactome_genesets), function(X) as.character(mget(X, geneset.reactome.hsa@set.name, ifnotfound = NA)))
summary(sapply(reactome_genesets, function(X) length(X)))
reactome_genesets <- Filter(function(x) length(x) >= 10 & length(x) < 1000, reactome_genesets)

gs_rdata_code = sprintf("geneset.%s.%s", "kegg_pathway", "hsa")
data(list = gs_rdata_code, package = "chipenrich.data",  envir = environment())

kegg_genesets <- as.list(geneset.kegg_pathway.hsa@set.gene)
names(kegg_genesets) <- lapply(names(kegg_genesets), function(X) as.character(mget(X, geneset.kegg_pathway.hsa@set.name, ifnotfound = NA)))
summary(sapply(kegg_genesets, function(X) length(X)))
kegg_genesets <- Filter(function(x) length(x) >= 10 & length(x) < 1000, kegg_genesets)


hallmarks_genesets = msigdbr(species = "human", category = "H")
head(hallmarks_genesets)
hallmarks_genesets <- split(hallmarks_genesets$entrez_gene, hallmarks_genesets$gs_name)
hallmarks_genesets <- Filter(function(x) length(x) >= 10 & length(x) < 1000, hallmarks_genesets)

all_genesets <- c(GOBP_genesets, reactome_genesets, hallmarks_genesets, kegg_genesets)
gene_sets_origin <- c(rep("GOBP", length(GOBP_genesets)), 
                      rep("reactome", length(reactome_genesets)), 
                      rep("hallmarks", length(hallmarks_genesets)),
                      rep("kegg", length(kegg_genesets)))
gene_sets_origin <- gene_sets_origin[-which(duplicated(names(all_genesets)))]
all_genesets <- all_genesets[-which(duplicated(names(all_genesets)))]

all_genesets_symbols <- lapply(all_genesets, function(X) clusterProfiler::bitr(X, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")$SYMBOL)


# convert mouse to human genes 
genes_info <- read.table(mouse_orthologs, header = T, sep = "\t")
df_mouse <- genes_info %>%
  filter(Common.Organism.Name == "mouse, laboratory") %>%
  dplyr::select(DB.Class.Key, mouse_symbol = Symbol)
df_human <- genes_info %>%
  filter(Common.Organism.Name == "human") %>%
  dplyr::select(DB.Class.Key, human_symbol = Symbol)
df_map <- df_mouse %>%
  dplyr::inner_join(df_human, by = "DB.Class.Key")
df_map

Gungabeesoon <- list("Angiogenesis" = c("Vegfa", "Snd1", "Mtdh", "Itga5", "Tnf", "Cxcl3", "Anxa3", "Hmgb1", "Hif1a", "Sema4d", "Lrg1", "Chil1"),
                     "Neutrophils" = c("Abr", "Anxa3", "Cd177", "Itgam", "Itgb2", "Itgb2l", "Pikfyve", "Pram1", "Ptafr", "Spi1", "Stx11", "Syk"),
                     "Neutrophil_cytotoxicity" = c("Ncf1", "Myd88", "Trem3", "Trem1", "Tusc2", "Cybb", "Cybc1", "Ncf2", "Ncf4", "Rac1", "Rac2"),
                     "Interferon_signaling" = c("Adar", "Isg15", "Isg20", "Rsad2", "Ifit1", "Ifit3", "Ifitm1", "Ifitm3", "Irak1", "Oas3", "Stat1", "Stat2", "Irf7", "Cxcl10"),
                     "Myeloid_recruitment" = c("Ccl3", "Mif", "Cxcl14", "Csf1", "Vegfa", "Ccl4", "Cxcl3"),
                     "ECM_remodeling" = c("Adamdec1", "Ctsc", "Ctsb", "Rgcc", "Ctss", "Ctsz", "Adam17", "Adam10", "Adam8"),
                     "Tumor_proliferation" = c("Tgfb1", "Tnf", "Il1a"),
                     "Immunosuppression" = c("Havcr2", "Fcgr2b", "Il4ra", "Cd274", "Hif1a"))

Gungabeesoon <- lapply(Gungabeesoon, function(x) df_map[df_map$mouse_symbol %in% x,"human_symbol"])
all_genesets_symbols <- c(all_genesets_symbols, Gungabeesoon)
gene_sets_origin <- c(gene_sets_origin,rep("Gungabeesoon", length(Gungabeesoon)))


senescence_genes <- readxl::read_excel(senescence_gene_file,sheet = 1)
table(senescence_genes$Classification)
senescence_genes$Classification <- paste0(gsub(" ","-",senescence_genes$Classification), "_SEN" )

senescence_genesets <- split(senescence_genes$`Gene(human)`, senescence_genes$Classification)
all_genesets_symbols <- c(all_genesets_symbols, senescence_genesets)
gene_sets_origin <- c(gene_sets_origin, rep("senescence", length(senescence_genesets)))

all_genesets_symbols[["Senescence"]] <- senescence_genes$`Gene(human)`
gene_sets_origin <- c(gene_sets_origin, "senescence")

names(all_genesets_symbols) <- paste0(names(all_genesets_symbols)," (",gene_sets_origin,')')

pathway_summary <- data.frame(
  database = gene_sets_origin,
  Pathway = names(all_genesets_symbols),
  Pathway_length = sapply(all_genesets_symbols, function(x) length(x)),
  Genes = sapply(all_genesets_symbols, function(x) paste(x, collapse = ","))
)

openxlsx::write.xlsx(pathway_summary, file = paste0(output_path, "/human_genesets_symbols.xlsx"))
write.table(pathway_summary, file = paste0(output_path, "/human_genesets_symbols.txt"), 
            quote=T, row.names = F, col.names = T, sep = "\t")

pathway_summary_IDs <- data.frame(
  database = gene_sets_origin,
  Pathway = names(all_genesets_symbols),
  Pathway_length = sapply(all_genesets_symbols, function(x) length(x)),
  Genes = sapply(all_genesets_symbols, function(x) paste(x, collapse = ","))
)

openxlsx::write.xlsx(pathway_summary, file = paste0(output_path, "/human_genesets_IDs.xlsx"))


