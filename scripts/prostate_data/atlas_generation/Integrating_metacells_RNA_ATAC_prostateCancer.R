{library(getopt)
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(EnsDb.Mmusculus.v79)
  library(Signac)
  library(S4Vectors)
  library(patchwork)
  set.seed(1234)}

options(future.globals.maxSize= 15000*1024^2)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "RNAnormalization", "n", 1, "character", "normalization method ('LogNormalization' or 'SCT')",
  "RNAcomp", "p", 1, "character", "range of components to consider for wnn analysis (eg 1:40 for RNA pca)",
  "ATACcomp", "q", 1, "character", "range of components to consider for wnn analysis (eg 1:50 for ATAC lsi)",
  'ATACintegration', "i",1,"character","ATAC integration method, ('harmony' default or 'rlsi' Signac)",
  'topFeaturesATAC', "t",1,"character","ATAC top feature param",
  "filterResolution", "f", 1, "numeric", "filtering resolution (default 0.5)."
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if (is.null(opt$RNAcomp)) {
  opt$RNAcomp <- c(1:50)
} else {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
}

if (is.null(opt$ATACcomp)) {
  opt$ATACcomp <- c(2:50)
} else {
  ci <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][2])
  opt$ATACcomp <- c(ci:cf)
}

if (is.null(opt$RNAnormalization)) {opt$RNAnormalization = "LogNormalize"}

if (is.null(opt$ATACintegration)) {opt$ATACintegration = "harmony"}


if (is.null(opt$filterResolution)) {
  opt$filterResolution <- 0.5
} 

if(is.null(opt$topFeaturesATAC)){
  opt$topFeaturesATAC = "q0"
}


print(opt)


# Get the files
metacells.files <- list.files(gsub("_harmony_q0|_harmony_q75|_signac_q0|_signac_q75","",opt$outdir),recursive = T,pattern = "metacells.rds",full.names = T,include.dirs = T)
metacells.files <- metacells.files[grep(pattern = "combined.metacells.rds|4_5M_1", metacells.files, invert = T)]
print(metacells.files)

metacell.objs <- lapply(metacells.files, FUN = function(x){ 
  x <- readRDS(x)
  x <- RenameCells(x,add.cell.id = unique(x$Sample))
  x@misc$new.membership <- paste0(unique(x$Sample),"_Metacell_", x@misc$membership)
  names(x@misc$new.membership) <- names(x@misc$membership)
  DefaultAssay(x) <- "RNA"
  x[["ATAC"]] <- NULL
  return(x)
}
)

#Same strategy as in the original study
#1 - integrate RNA,
#2 - discard low qual cluster
#3 - integrate ATAC/RNA on remaining cells
#4 - multimodal wnn

#Integrating RNA assay
metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
  x <- FindVariableFeatures(x)
})


## Integration of RNA data

if (opt$RNAnormalization != "SCT") {
  
  rnaAssay = "RNA"	
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) { 
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x,normalization.method = "LogNormalize",assay = "RNA")
    x <- FindVariableFeatures(x)
    return(x)})
  features <- SelectIntegrationFeatures(object.list = metacell.objs,verbose = F)
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
    x <- ScaleData(x,features = features);
    x <- RunPCA(x, features = features, verbose = FALSE);
    return(x)
  })
  
  new.memberships <- unlist(lapply(metacell.objs,FUN = function(x) {x@misc$new.membership}))
  
  print("Finding anchors")
  
  anchors <- FindIntegrationAnchors(object.list = metacell.objs, 
                                    reduction = "rpca",
                                    anchor.features = features,
                                    verbose = F,
                                    dims = opt$RNAcomp)
  
  print("Integrating data")
  rna.All.combined <- IntegrateData(anchorset = anchors, dims = opt$RNAcomp,verbose = F)
  rna.All.combined@misc$membership <- new.memberships
  
  rm(metacell.objs,anchors)
  gc()
  
  DefaultAssay(rna.All.combined) <- "integrated"
  rna.All.combined <- ScaleData(rna.All.combined, verbose = FALSE)
  
} else {
  rnaAssay = "SCT"
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) { DefaultAssay(x) <- "RNA";x <- SCTransform(x, verbose = FALSE,vst.flavor = "v2"); return(x)})
  features <- SelectIntegrationFeatures(object.list = metacell.objs,nfeatures = 3000,verbose = F)
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  new.memberships <- unlist(lapply(metacell.objs,FUN = function(x) {x@misc$new.membership}))
  
  metacell.objs <- PrepSCTIntegration(object.list = metacell.objs, anchor.features = features,verbose = F)
  
  print("Finding anchors")
  
  anchors <- FindIntegrationAnchors(object.list = metacell.objs, 
                                    reduction = "rpca", 
                                    normalization.method = "SCT",
                                    anchor.features = features,
                                    verbose = F,
                                    dims = opt$RNAcomp)
  
  print("Integrating data")
  rna.All.combined <- IntegrateData(anchorset = anchors, 
                                    dims = opt$RNAcomp,
                                    verbose = F,
                                    normalization.method = "SCT")
  rna.All.combined@misc$membership <- unlist(lapply(metacell.objs,FUN = function(x) {x@misc$new.membership}))
  
  rm(metacell.objs,anchors)
  gc()
  
  DefaultAssay(rna.All.combined) <- "integrated"
}



#saveRDS(rna.All.combined,"rna.all.combined.rds")

rna.All.combined <- RunPCA(rna.All.combined, verbose = FALSE)
rna.All.combined <- RunUMAP(rna.All.combined, dims = opt$RNAcomp,reduction.name = "rna.umap")

rna.All.combined <- FindNeighbors(rna.All.combined, dims = opt$RNAcomp)
rna.All.combined <- FindClusters(rna.All.combined, resolution = c(1:10)*0.1)

genes <- c("Epcam","Krt8","Chga","Krt5","Svs5","Col1a2","Pecam1","Ptprc","S100a9","C1qa","Cd3e","Cd19","Plp1")
DefaultAssay(rna.All.combined) <- rnaAssay
for (r in c(1:10)*0.1) {
  pdf(paste0(opt$outdir,'res',r,"_all_rna_mcs.pdf"))
  Idents(rna.All.combined) <- paste0('integrated_snn_res.',r)
  print(DimPlot(rna.All.combined,label = T,reduction = "rna.umap"))
  print(VlnPlot(rna.All.combined,features = c("nFeature_RNA"),pt.size = 0.001,log = T))
  print(VlnPlot(rna.All.combined,features = c("nCount_RNA"),pt.size = 0.001,log = T))
  print(FeaturePlot(rna.All.combined,features = genes[1:7]))
  print(FeaturePlot(rna.All.combined,features = genes[8:13]))
  
  
  dev.off()
}
DefaultAssay(rna.All.combined) <- "integrated"

filter.res <- paste0("integrated_snn_res.",opt$filterResolution)
# Remove low quality cluster
df <- rna.All.combined@meta.data %>%
  group_by_at(filter.res) %>%
  summarise(medNF = median(nFeature_RNA))

print(df)

lowQualClust <- df[which.min(df$medNF),filter.res]
print(lowQualClust)
Idents(rna.All.combined) <- filter.res
rna.filtered.combined <- subset(rna.All.combined,idents = lowQualClust[[filter.res]],invert = T)
rna.filtered.combined <- rna.filtered.combined[, rna.filtered.combined$nFeature_RNA > min(df$medNF)] 
rna.filtered.combined@misc$membership <- new.memberships[new.memberships %in% colnames(rna.filtered.combined)]

remove(rna.All.combined)
gc()

# integrate RNA clean data
print("Integrate RNA clean data")
DefaultAssay(rna.filtered.combined) <- rnaAssay
metacell.objs <- SplitObject(rna.filtered.combined,split.by = "Sample")

if (opt$RNAnormalization != "SCT") {
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) { DefaultAssay(x) <- "RNA";NormalizeData(x,normalization.method = "LogNormalize",assay = "RNA") %>% FindVariableFeatures(); return(x)})
  features <- SelectIntegrationFeatures(object.list = metacell.objs,verbose = F)
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
    x <- ScaleData(x,features = features);
    x <- RunPCA(x, features = features, verbose = FALSE);
    return(x)
  })
  
  print("Finding anchors")
  
  anchors <- FindIntegrationAnchors(object.list = metacell.objs, 
                                    reduction = "rpca",
                                    anchor.features = features,verbose = F,
                                    dims = opt$RNAcomp)
  
  print("Integrating data")
  rna.filtered.combined <- IntegrateData(anchorset = anchors, dims = opt$RNAcomp,verbose = F)
  rm(metacell.objs,anchors)
  gc()
  
  DefaultAssay(rna.filtered.combined) <- "integrated"
  rna.filtered.combined <- ScaleData(rna.filtered.combined, verbose = FALSE)
  
} else {
  rnaAssay = "SCT"
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) { DefaultAssay(x) <- "RNA";x <- SCTransform(x, verbose = FALSE,vst.flavor = "v2"); return(x)})
  features <- SelectIntegrationFeatures(object.list = metacell.objs,verbose = F)
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  metacell.objs <- PrepSCTIntegration(object.list = metacell.objs, anchor.features = features,verbose = F)
  
  print("Finding anchors")
  
  anchors <- FindIntegrationAnchors(object.list = metacell.objs, 
                                    reduction = "rpca", 
                                    normalization.method = "SCT",
                                    anchor.features = features, 
                                    verbose = F,
                                    dims = opt$RNAcomp)
  
  print("Integrating data")
  rna.filtered.combined <- IntegrateData(anchorset = anchors, 
                                         dims = opt$RNAcomp, 
                                         verbose = F,  
                                         normalization.method = "SCT")
  rm(metacell.objs,anchors)
  gc()
  
  DefaultAssay(rna.filtered.combined) <- "integrated"
}


rna.filtered.combined <- RunPCA(rna.filtered.combined, verbose = FALSE,reduction.name = "integrated_pca")
rna.filtered.combined <- RunUMAP(rna.filtered.combined, dims = opt$RNAcomp,reduction = "integrated_pca",reduction.name = "rna.umap")
rna.filtered.combined <- FindNeighbors(rna.filtered.combined,reduction = "integrated_pca", dims = opt$RNAcomp)
rna.filtered.combined <- FindClusters(rna.filtered.combined, resolution = 0.2)

# saveRDS(rna.filtered.combined,"tmp.rds")

kept.cells <- as.vector(colnames(rna.filtered.combined))
print("Loading ATAC")

metacell.objs.atac <- lapply(metacells.files,function(x) {
  x <- readRDS(x)
  DefaultAssay(x) <- "ATAC"
  smp <- unique(x$Sample)
  print(smp)
  x <- RenameCells(x,add.cell.id = unique(x$Sample))
  x <- x[,colnames(x) %in% kept.cells]
  x[["RNA"]] <- NULL
  gc()
  return(x)
})

if (sum(rownames(metacell.objs.atac[[1]]) %in% rownames(metacell.objs.atac[[2]]) )  < 50000) {
  
  print("Get a consensus peak set...")
  peak.list <- lapply(X =metacell.objs.atac,FUN = function(x) {
    GenomicRanges::GRanges(sub(x=rownames(x),pattern="-",replacement = ":"))
  })
  
  combined.peaks <- GenomicRanges::reduce(x=c(peak.list[[1]],peak.list[[2]]))
  for (p in peak.list[2:length(peak.list)]) {
    combined.peaks <- GenomicRanges::reduce(x=c(combined.peaks,p))
  }
  
  peakwidths <- GenomicRanges::width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
  combined.peaks
  
  print("Recompute peak count matrix...")
  
  library(future.apply)
  
  metacell.objs.atac <- future_lapply(X = metacell.objs.atac,FUN = function(x) {
    counts <- FeatureMatrix(
      fragments = Fragments(x),
      features = combined.peaks,
      cells = colnames(x)
    )
    
    x[["ATAC"]] <- CreateChromatinAssay(counts,genome = 'mm10', fragments = Fragments(x))
  }
  )
  
  
}



atac.filtered.combined <- merge(metacell.objs.atac[[1]],metacell.objs.atac[-1])

print("Creating combined object")
#################add ATAC assay to RNA assay
rna.filtered.combined[["ATAC"]]<- atac.filtered.combined[["ATAC"]]
rm(atac.filtered.combined)
gc()
DefaultAssay(rna.filtered.combined) <- "ATAC"
rna.filtered.combined <- FindTopFeatures(rna.filtered.combined, min.cutoff = opt$topFeaturesATAC) 
rna.filtered.combined <- RunTFIDF(rna.filtered.combined)
rna.filtered.combined <- RunSVD(rna.filtered.combined)
rna.filtered.combined <-  RunUMAP(rna.filtered.combined, reduction = "lsi",dims = opt$ATACcomp)



if (opt$ATACintegration == "harmony") {

  
  DefaultAssay(rna.filtered.combined) <- "ATAC"
  rna.filtered.combined <- harmony::RunHarmony(
    object = rna.filtered.combined,
    group.by.vars = 'Sample',
    reduction = 'lsi',reduction.save = 'integrated_lsi',
    assay.use = 'ATAC',
    project.dim = FALSE
  )
  
} else {
  
  metacell.objs.atac <- lapply(metacell.objs.atac,function(x) {
    DefaultAssay(x) <- "ATAC"
    x <- FindTopFeatures(x, min.cutoff = opt$topFeaturesATAC) 
    x <- RunTFIDF(x)
    x <- RunSVD(x)
    return(x)
  })

  ################# batch correction ATAC data using Signac
  
  depth.cor.list <- lapply(metacell.objs.atac, function(x){
    depth.cor <- DepthCor(x, n = max(opt$ATACcomp),reduction ="lsi")
    depth.cor <- abs(depth.cor$data$nCount_ATAC)
    return(depth.cor)
  })
  depth.cor.df <- do.call(rbind, depth.cor.list)
  dims.lsi.cusanovich.smp <- which(colMeans(depth.cor.df) < 0.5)
  print("Nb of variable features for ATAC")
  print(length(VariableFeatures(rna.filtered.combined)))
  

  integration.anchors <- FindIntegrationAnchors(
    object.list = metacell.objs.atac,
    anchor.features = VariableFeatures(rna.filtered.combined), #rownames(rna.filtered.combined), #VariableFeatures(rna.filtered.combined)
    reduction = "rlsi",
    dims = dims.lsi.cusanovich.smp #opt$ATACcomp #dims.lsi.cusanovich.smp
  )
  

  success <- FALSE
  k.weight <- 100
  # Loop until success or k.weight becomes too small
  while (!success && k.weight > 0) {
    # Try running the IntegrateEmbeddings function
    result <- tryCatch({
      atac.filtered.combined <- IntegrateEmbeddings(
        anchorset = integration.anchors,
        reductions = rna.filtered.combined[["lsi"]],
        new.reduction.name = "integrated_lsi",
        dims.to.integrate = 1:max(opt$ATACcomp),
        k.weight = k.weight
      )
      list(TRUE,NA)  # If successful, return TRUE
    }, error = function(e) {
      cat("Error encountered with k.weight =", k.weight, ": ", e$message, "\n")
      new.k.weight <- strsplit(e$message,split = "k.weight to less than ")[[1]][2]
      new.k.weight <- as.numeric(strsplit(new.k.weight, split = " or increase k.anchor.")[[1]][1]) - 10
      
      cat("New k.weight =", new.k.weight, "\n")
      
      return(list(FALSE, new.k.weight))  # If there's an error, return FALSE
    })
    
    # Check if the function was successful
    if (result[[1]]) {
      success <- TRUE
      cat("integrated_lsi", " Integration successful with k.weight =", k.weight, "\n")
    } else {
      # Decrement k.weight
      k.weight <- result[[2]]
    }
  }
  
  rna.filtered.combined@reductions[['lsi']] <- atac.filtered.combined@reductions$lsi
  rna.filtered.combined@reductions[['integrated_lsi']] <- atac.filtered.combined@reductions$integrated_lsi
  remove(atac.filtered.combined)
  gc()
}

#########re-compute the UMAP using corrected LSI embeddings
rna.filtered.combined <- RunUMAP(rna.filtered.combined, dims = opt$ATACcomp, reduction = 'integrated_lsi',reduction.name = "atac.umap")
DefaultAssay(rna.filtered.combined) <- "ATAC"
rna.filtered.combined <- FindMultiModalNeighbors(rna.filtered.combined, reduction.list = list("integrated_pca", "integrated_lsi"), dims.list = list(opt$RNAcomp, opt$ATACcomp))
rna.filtered.combined <- FindClusters(rna.filtered.combined, resolution = 0.1,graph.name = "wsnn") 
rna.filtered.combined <- RunUMAP(rna.filtered.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

saveRDS(new.memberships,"membership.rds")
rna.filtered.combined@misc$membership <- new.memberships[new.memberships %in% colnames(rna.filtered.combined)]

saveRDS(rna.filtered.combined,paste0(opt$outdir,"/combined.metacells.rds"))


p1 <- DimPlot(rna.filtered.combined, reduction = "rna.umap",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(rna.filtered.combined, reduction = "atac.umap",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(rna.filtered.combined,reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

pdf(file=paste0(opt$outdir,"/umaps.pdf"))
p1
p2
p3

print(DimPlot(rna.filtered.combined,label = T,reduction = "wnn.umap"))
print(VlnPlot(rna.filtered.combined,features = c("nFeature_RNA"),pt.size = 0.001,log = T))
print(VlnPlot(rna.filtered.combined,features = c("nCount_RNA"),pt.size = 0.001,log = T))
DefaultAssay(rna.filtered.combined) <- "RNA"
print(FeaturePlot(rna.filtered.combined,features = genes[1:7],reduction = "wnn.umap"))
print(FeaturePlot(rna.filtered.combined,features = genes[8:13],reduction = "wnn.umap"))

dev.off()
## 


















