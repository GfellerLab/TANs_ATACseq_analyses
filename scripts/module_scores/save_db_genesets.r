library(msigdbr)
library(biomaRt)
library(dplyr)
library(tidyr)

output_path = "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/input_data/"

# Mouse genesets ----------------------------------------------------------


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

genesets_data <- list(all_genesets, gene_sets_origin)
saveRDS(genesets_data, file = paste0(output_path, "peaks_db_genesets_mouse.rds"))


# Human genesets ----------------------------------------------------------
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

genesets_data <- list(all_genesets, gene_sets_origin)
saveRDS(genesets_data, file = paste0(output_path, "peaks_db_genesets_human.rds"))


