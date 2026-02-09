
# Libraries ---------------------------------------------------------------

library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2025)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(limma)

# Parameters --------------------------------------------------------------

metadata_path <- as.character(commandArgs(TRUE)[1])
counts_path <- as.character(commandArgs(TRUE)[2])
pfm_file <- as.character(commandArgs(TRUE)[3])
output_path <- as.character(commandArgs(TRUE)[4])
genome <- as.character(commandArgs(TRUE)[5])
# motif_annotation_path <- as.character(commandArgs(TRUE)[6])

if(genome == "hg38"){
  genome <- BSgenome.Hsapiens.UCSC.hg38
}else{
  genome <- BSgenome.Mmusculus.UCSC.mm10
}

dir.create(output_path, recursive = T)


# Load counts -------------------------------------------------------------
counts <- read.table(counts_path, header = T, sep = "\t", row.names = 1)
rowRanges <- Signac::StringToGRanges(rownames(counts), sep = c(":", "-"))
colnames(counts) <- gsub("[.]", "-",colnames(counts))

# Load metadata -----------------------------------------------------------
metadata <- read.table(metadata_path, header = T, sep = "\t")
rownames(metadata) <- metadata$sample_name
metadata <- metadata[colnames(counts),] 

# Get chromVAR deviation scores -------------------------------------------
counts <- SummarizedExperiment(assays=list(counts=as.matrix(counts)),
                               rowRanges=rowRanges, colData=metadata)

# add GC bias
counts <- addGCBias(counts, 
                    genome = genome)
head(rowData(counts))

# get PFM data and match motifs
pfm <- readRDS(pfm_file)
motif_ix <- matchMotifs(pfm, counts, 
                        genome = genome)

# register(MulticoreParam(3))
dev <- computeDeviations(object = counts, annotations = motif_ix, background_peaks = getBackgroundPeaks(counts,niterations=2000))
dev_matrix <- assays(dev)$deviations
rownames(dev_matrix) <- sapply(rownames(dev_matrix), function(i) pfm@listData[[i]]@name)


z_matrix <- assays(dev)$z
rownames(z_matrix) <- sapply(rownames(z_matrix), function(i) pfm@listData[[i]]@name)

# Save chromVAR deviation scores ------------------------------------------
saveRDS(dev_matrix, file = paste0(output_path, "dev_scores.rds"))
saveRDS(z_matrix, file = paste0(output_path, "z_scores.rds"))
z_matrix_toSave <- rbind(t(data.frame(metadata$groups)),z_matrix )
write.csv2(cbind(data.frame(row_names = rownames(z_matrix_toSave)),z_matrix_toSave), row.names = F, file = paste0(output_path, "z_scores_lung.csv"))
write.csv2(cbind(data.frame(row_names = rownames(dev_matrix)),dev_matrix), row.names = F, file = paste0(output_path, "dev_scores_lung.csv"))

