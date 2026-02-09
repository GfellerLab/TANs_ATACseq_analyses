library(Seurat)
library(Signac)
library(SuperCell)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)
library(getopt)
library(doParallel)
library(rhdf5)
spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputH5',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'fragmentFile', "f",1,  "character", "Fragment file path",
  'gamma', 'g', 1, "numeric", "gmaa for metacell identification",
  "pythonSeacellEnv", "p", 1, "character", "python path for seacell env to compute compactness/separation",
  "consensusPeakSet", "c", 1, "character", "path to consensus peak set file",
  "peakRDS", "a", 1, "character", "direct path to consensus macs2 peak set file (.rds)",
  "normalization", "n", 1, "character", "normalization method ('LogNormilsation' or 'SCT')",
  "macs2", "m", 1, "character", "path to macs2 to compue narrow peaks and use them instead of cellranger atac peak"
), byrow=TRUE, ncol=5)

opt = getopt(spec)


if (is.null(opt$normalization)) {
  opt$normalization <- "LogNormalization"
}



print(opt)

smp <- strsplit(opt$inputH5,"/")[[1]][length(strsplit(opt$inputH5,"/")[[1]])-1]


fragments_file <- opt$fragmentFile
x = opt$inputH5
dims <- as.integer(h5read(x, "matrix/shape"))
n_features <- dims[1]
n_cells    <- dims[2]

data_vals <- as.numeric(h5read(x, "matrix/data"))
row_indices <- as.integer(h5read(x, "matrix/indices"))
col_pointers <- as.integer(h5read(x, "matrix/indptr"))

barcodes <- h5read(x, "matrix/barcodes")

feature_ids   <- h5read(x, "matrix/features/id")       # length = n_features
feature_names <- h5read(x, "matrix/features/name")     # length = n_features
feature_intervals <- h5read(x, "matrix/features/interval")  # length = n_features
feature_type  <- h5read(x, "matrix/features/feature_type")

full_counts <- new(
  "dgCMatrix",
  Dim      = c(n_features, n_cells),
  Dimnames = list(feature_names, barcodes),
  x        = data_vals,
  i        = row_indices,
  p        = col_pointers
)

rna_rows  <- which(feature_type == "Gene Expression")
atac_rows <- which(feature_type == "Peaks")

rna_counts  <- full_counts[rna_rows, , drop = FALSE]
rownames(rna_counts) <- feature_names[rna_rows]
atac_counts <- full_counts[atac_rows, , drop = FALSE]
rownames(atac_counts) <- feature_names[atac_rows]



rna_counts <- rna_counts[-which(duplicated(rownames(rna_counts))),]
obj <- CreateSeuratObject(CreateAssayObject(rna_counts))
obj$Sample <- smp

if (!is.null(opt$peakRDS)) {
  consensus_peakset <- readRDS(opt$peakRDS)
  
  peaks <- FeatureMatrix(
    fragments = CreateFragmentObject(opt$fragmentFile),
    features = consensus_peakset,
    cells = colnames(obj)
  )
  
  print("macs2 peak loaded")
  

  annotations.peakset.inter <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations.peakset.inter) <- 'UCSC'
  genome(annotations.peakset.inter) <- "mm10"
  
  chrom_assay<- CreateChromatinAssay(
    counts = peaks,
    #sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1, # peaks have been identified after a first analysis
    annotation=annotations.peakset.inter,
    fragments = opt$fragmentFile
  )
  print("chrom assay created")
} else {
  print("using cell ranger peaks")

  annotations<- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "mm10"
  
  # same code as in the original study
  chrom_assay<- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = opt$fragmentFile,
    min.cells = 10,
    annotation = annotations
  )
  
}


if (!is.null(opt$consensusPeakSet)) {
  print("Loading consensus peak set table...")
  consensus_peakset <- read.csv(opt$consensusPeakSet,row.names = 1)
  consensus_peakset <- paste(consensus_peakset$seqnames,consensus_peakset$start,consensus_peakset$end,sep = "-")
  
  
  peaks <- FeatureMatrix(
    fragments = Fragments(chrom_assay),
    features = consensus_peakset,
    cells = colnames(chrom_assay)
  )
  
  grange.counts.peakset.inter <- StringToGRanges(rownames(peaks), sep = c(":", "-"))
  grange.use.peakset.inter <- seqnames(grange.counts.peakset.inter) %in% standardChromosomes(grange.counts.peakset.inter)
  atac_counts.peakset.inter <- peaks[as.vector(grange.use.peakset.inter), ]
  annotations.peakset.inter <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations.peakset.inter) <- 'UCSC'
  genome(annotations.peakset.inter) <- "mm10"
  
  chrom_assay<- CreateChromatinAssay(
    counts = peaks,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = Fragments(chrom_assay)
  )
}

if (!is.null(opt$macs2)) {
  print("Call peaks with macs2")
  peaks <- CallPeaks(chrom_assay,macs2.path = opt$macs2)
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
  # quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(chrom_assay),
    features = peaks,
    cells = colnames(chrom_assay)
  )
  #atac

  
  #overwrite ATAC assay
  new_chrom_assay<- CreateChromatinAssay(
    counts = macs2_counts,
    genome = 'mm10',
    min.cells = 1,
    fragments = Fragments(chrom_assay),
    annotation = annotations
  )
  obj[["ATAC"]] <- new_chrom_assay
  rm(chrom_assay)
  rm(new_chrom_assay)
  
} else {
  obj[["ATAC"]] <- chrom_assay
  rm(chrom_assay)
}



DefaultAssay(obj) <- "RNA"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

rm(inputdata.10x)
gc()

# RNA processing


if (opt$normalization == "LogNormalization") {
obj <- NormalizeData(obj, verbose = FALSE) %>% FindVariableFeatures(nfeatures = 2000,verbose = FALSE) %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
rnaAssay <- "RNA"
} else {
  #SCT
  obj <- SCTransform(obj, verbose = FALSE,vst.flavor = "v2")
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = c(1:50),reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  rnaAssay <- "SCT"
}

# ATAC processing

DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)

# We exclude the first dimension as this is typically correlated with sequencing depth
png(paste0(opt$outdir,"/depthCor.png"))
DepthCor(obj)
dev.off()


obj <- RunUMAP(obj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
obj <- FindClusters(obj, graph.name = "wknn", algorithm = 1, verbose = FALSE)


p1 <- DimPlot(obj, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(obj, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(obj, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

pdf(paste0(opt$outdir,"/umaps.pdf"))
p1 + NoLegend()
p2 + NoLegend()
p3 + NoLegend()
dev.off()

# Metacell identification with SuperCell
outputDirMcFragment <- paste0(opt$outdir,"/aggregated_fragment_file/")
fragmentFiles <- list()
fragmentFiles[["ATAC"]] <- GetFragmentData(object = Fragments(obj)[[1]], slot = "path")

if(!is.null(opt$peakRDS)) {
  saveRDS(obj,paste0(opt$outdir,"/sc_sobj.rds"))
}

obj.mc <- SuperCell::SCimplify_for_Seurat_v5(seurat = obj,
                               k.knn = 30,
                               kernel = T,
                               gamma = opt$gamma,
                               assay = c(rnaAssay,"ATAC"),
                               dims = list(1:50,2:50),
                               reduction = list("pca","lsi"),
                               fragmentFiles = fragmentFiles,
                               prefixMC = paste0("Metacell_"),
                               tmpPath =paste0(opt$outdir,"/tmp/"),
                               outputDirMcFragment = outputDirMcFragment,
                               nb_cl = 10, 
                               graph.name = "knn")



DefaultAssay(obj.mc) <- "ATAC"
obj.mc <- FindTopFeatures(obj.mc, min.cutoff = 'q0') 
obj.mc <- RunTFIDF(obj.mc)

DefaultAssay(obj.mc) <- "RNA"

if (opt$normalization == "LogNormalization") {
  obj.mc <- NormalizeData(obj.mc, verbose = FALSE) %>% FindVariableFeatures(nfeatures = 2000,verbose = FALSE) %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
} else {
  obj.mc <- SCTransform(obj.mc, verbose = FALSE,vst.flavor = "v2")
  obj.mc <- RunPCA(obj.mc, verbose = FALSE)
}


pdf(paste0(opt$outdir,"/wnn.umap.metacells.pdf"))
DimPlotSC(seurat.mc = obj.mc,seurat = obj,metacell.col = "seurat_clusters",sc.col = "seurat_clusters",reduction = "wnn.umap") + theme_classic()
dev.off()

# Save seurat metacell object
saveRDS(obj.mc,paste0(opt$outdir,"/metacells.rds"))


