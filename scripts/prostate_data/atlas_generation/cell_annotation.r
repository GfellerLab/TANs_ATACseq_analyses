
# Libraries ---------------------------------------------------------------

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(EnsDb.Mmusculus.v79)
library(Signac)
library(S4Vectors)
library(patchwork)
library(SuperCell)
library(reshape2)
library(colorspace)
library(getopt)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'input_rds',  'i', 1, "character", "REQUIRED : path to the h5 files",
  'output',     'o',1, "character", 'Outdir path (default ./)'
), byrow=TRUE, ncol=5)

opt = getopt(spec)


combined.metacells <- readRDS(opt$input_rds)
output_path <- opt$output

# Get clustering ----------------------------------------------------------

combined.metacells <- FindClusters(combined.metacells, resolution = c(c(1:7)*0.01),graph.name = "wsnn",algorithm = 3) 


Idents(combined.metacells) <- "wsnn_res.0.07" 

# Annotate cells ----------------------------------------------------------
df <- combined.metacells@meta.data %>%
  group_by_at("wsnn_res.0.07") %>%
  summarise(medNF = median(nFeature_RNA))

print(df)


pdf(paste0(output_path, 'markers_major_types2.pdf'))
print(VlnPlot(combined.metacells,features = c("nFeature_RNA"),pt.size = 0.001,log = T))
DimPlot(combined.metacells,reduction = "wnn.umap", label = T)
print(FeaturePlot(combined.metacells,features = "nFeature_RNA"))


#Luminal cluster 1 #1 and maybe 12

DefaultAssay(combined.metacells) <- "RNA"
VlnPlot(combined.metacells,features = c("Krt19","Nupr1","Lmo7","Ceacam1"),ncol = 2,pt.size = 0,assay = "RNA")
VlnPlot(combined.metacells,features = c("Fcgbp","Slc12a2","Lmo7","Ceacam1"),ncol = 2,pt.size = 0,assay = "RNA")

# Basal 6
VlnPlot(combined.metacells,features = c("Krt14","Krt15","Krt5","Krt17"),ncol = 2,pt.size = 0,assay = "RNA")

# Neuroendocrine 2 10? #11 2 

VlnPlot(combined.metacells,features = c("Chga","Nrxn1","Hcn1","Fgf14"),ncol = 2,pt.size = 0,assay = "RNA")
VlnPlot(combined.metacells,features = c("Cadm2","Kcnb2","Kcnip4","Cntn4"),ncol = 2,pt.size = 0,assay = "RNA")

# Seminal vesicle 3 #8/7
VlnPlot(combined.metacells,features = c("Svs5","Svs2","Cdo1","Fyb2"),ncol = 2,pt.size = 0,assay = "RNA")

# Mesenchymal 1 0
VlnPlot(combined.metacells,features = c("Col1a2","Apod","Dcn","Serping1"),ncol = 2,pt.size = 0,assay = "RNA")

# Mesenchymal 2 8 #5
VlnPlot(combined.metacells,features = c("Col1a2","Apod","Flt1","Serping1"),ncol = 2,pt.size = 0,assay = "RNA")

# Endothelial 8 #9
VlnPlot(combined.metacells,features = c("Pecam1","Emcn","Selp","Flt1"),ncol = 2,pt.size = 0,assay = "RNA")

# Neutrophil 4 #4
VlnPlot(combined.metacells,features = c("S100a9","Clec4d","Cxcl2","Acod1"),ncol = 2,pt.size = 0,assay = "RNA")

# Macrophages 7 #3
VlnPlot(combined.metacells,features = c("Cd74","Cd86","Cybb"),ncol = 2,pt.size = 0,assay = "RNA")

#T cell 9 #10
VlnPlot(combined.metacells,features = c("Skap1","Itk","Ptpn22","Cd3e"),ncol = 2,pt.size = 0,assay = "RNA")

# B cell 11
VlnPlot(combined.metacells,features = c("Igkc","Ebf1","Mef2c","Bank1"),ncol = 2,pt.size = 0,assay = "RNA")

# Neuron 13
VlnPlot(combined.metacells,features = c("Plp1"),ncol = 2,pt.size = 0,assay = "RNA") + NoLegend()
FeaturePlot(combined.metacells,c("Plp1"),reduction = "wnn.umap")

combined.metacells <- FindSubCluster(combined.metacells,cluster = 9,graph.name = "wsnn",resolution = 0.05)
DimPlot(combined.metacells[,combined.metacells$wsnn_res.0.07 == 9],group.by = "sub.cluster",reduction = "wnn.umap")
DefaultAssay(combined.metacells) <- "RNA"

Idents(combined.metacells) <- "sub.cluster"
# VlnPlot(combined.metacells,features = c("Igkc","Ebf1","Mef2c","Bank1"),ncol = 2,pt.size = 0,assay = "RNA")
VlnPlot(combined.metacells,features = c("Plp1"),ncol = 2,pt.size = 0,assay = "RNA")

dev.off()

combined.metacells$celltype <- combined.metacells$sub.cluster
Idents(combined.metacells) <- "celltype"
combined.metacells <-  RenameIdents(combined.metacells, '1' = 'Luminal_1', '10' = 'Luminal_2')
combined.metacells <- RenameIdents(combined.metacells, '6' = 'Basal')
combined.metacells <- RenameIdents(combined.metacells, '2' = 'Neuroendocrine')
combined.metacells <- RenameIdents(combined.metacells, '9_2' = 'Neuron')
combined.metacells <- RenameIdents(combined.metacells, '3' = 'Seminal_vesicle')
combined.metacells <- RenameIdents(combined.metacells, '9_0' = 'T_cell')
combined.metacells <- RenameIdents(combined.metacells, '0' = "Mesenchymal_1",  '7' = "Mesenchymal_2")
combined.metacells <- RenameIdents(combined.metacells, '8' = "Endothelial")
combined.metacells <- RenameIdents(combined.metacells, '4' = "Neutrophil")
combined.metacells <- RenameIdents(combined.metacells, '9_1' = "T_cell")
combined.metacells <- RenameIdents(combined.metacells, '5' = "Macrophage")
combined.metacells <- RenameIdents(combined.metacells, '11' = "B_cell")

combined.metacells$celltype <- Idents(combined.metacells)


# Visualize annotation ----------------------------------------------------


my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

colors <- c('Luminal_1'='#E5D2DD',
            "Luminal_2" = '#CCC9E6',
            'Neuroendocrine'='#53A85F',
            'Basal' = '#F1BB72',
            "Seminal_vesicle"= '#F3B1A0',
            'Mesenchymal_1'='#D6E7A3',
            'Mesenchymal_2'='#57C3F3',
            "Endothelial"='#476D87',
            "Neutrophil"='#E95C59',
            "Macrophage"='#AB3282',
            "T_cell" = '#23452F',
            "B_cell" = '#BD956A',
            "Neuron" = '#8C549C')


Idents(combined.metacells) <- "celltype"

plot1 <- DimPlot(combined.metacells,reduction = "wnn.umap",cols = colors,label = T,repel = T)
plot1 + NoLegend()


combined.metacells@meta.data$sample <- factor(combined.metacells$Sample,
                                              levels = c("WT_1","WT_2",
                                                         "2W_1","2W_2",
                                                         "1M_1","1M_2",
                                                         "2_5M_1","2_5M_2","2_5M_3","2_5M_4",
                                                         "3_5M_1","3_5M_2",
                                                         "4_5M_1","4_5M_2",
                                                         "6M"))


combined.metacells$time <- combined.metacells$Sample
combined.metacells$time[startsWith(combined.metacells$Sample,"1M")] <- "1M"
combined.metacells$time[startsWith(combined.metacells$Sample,"2_5M")] <- "2_5M"
combined.metacells$time[startsWith(combined.metacells$Sample,"2W")] <- "2W"
combined.metacells$time[startsWith(combined.metacells$Sample,"3_5M")] <- "3_5M"
combined.metacells$time[startsWith(combined.metacells$Sample,"4_5M")] <- "4_5M"
combined.metacells$time[startsWith(combined.metacells$Sample,"6M")] <- "6M"
combined.metacells$time[startsWith(combined.metacells$Sample,"WT")] <- "WT"

combined.metacells@meta.data$time <- factor(combined.metacells$time,levels = c("WT","2W","1M",
                                                                               "2_5M","3_5M","4_5M",
                                                                               "6M"))
pdf(paste0(output_path, "celltype_in_time.pdf"))

ggplot(combined.metacells@meta.data,aes(x=sample,fill = time)) + geom_bar(stat = "count")  +
  scale_x_discrete(guide = guide_axis(angle = 45)) + ylab("metacell counts")
dev.off()


# Save major types --------------------------------------------------------

combined.metacells$major_type <- combined.metacells$celltype
combined.metacells$major_type <- as.character(gsub(x=combined.metacells$major_type, pattern = "_[0-9]",replacement = ""))


colors <- c('Luminal'='#E5D2DD',
            'Neuroendocrine'='#53A85F',
            'Basal' = '#F1BB72',
            "Seminal_vesicle"= '#F3B1A0',
            'Mesenchymal'='#D6E7A3',
            "Endothelial"='#476D87',
            "Neutrophil"='#E95C59',
            "Macrophage"='#AB3282',
            "T_cell" = '#23452F',
            "B_cell" = '#BD956A',
            "Neuron" = '#8C549C')

combined.metacells@meta.data$major_type <- factor(combined.metacells$major_type,                                                
                                                  levels = (c("Luminal","Neuroendocrine","Basal",
                                                              "Seminal_vesicle","Mesenchymal","Endothelial",
                                                              "Neutrophil","Macrophage",
                                                              "T_cell","B_cell","Neuron"))) 

pdf(paste0(output_path, "major_celltypes.pdf"))
DimPlot(combined.metacells,reduction = "wnn.umap",group.by = "major_type",cols = colors)
dev.off()



# Check cell types frequencies --------------------------------------------

smpCounts <- aggregate(combined.metacells$size, by=list(sample = combined.metacells$sample,
                                                        major_type = combined.metacells$major_type,
                                                        celltype = combined.metacells$celltype,
                                                        time = combined.metacells$time), FUN=sum)


darker_colors <- darken(colors, 0.2)

contingencyTable <- xtabs(x ~ major_type+time,data = smpCounts)
colSums(contingencyTable)

freqMatrix <- apply(contingencyTable,1,FUN = function(x){x/colSums(contingencyTable)})

freqMatrix_df <- melt(freqMatrix)
freqMatrix_df$major_type <- factor(freqMatrix_df$major_type ,
                                   levels = c("Basal","Endothelial","Luminal","Neuroendocrine","B_cell","T_cell","Macrophage","Neutrophil","Seminal_vesicle","Mesenchymal","Neuron"))
freqMatrix_df$group <- freqMatrix_df$major_type
freqMatrix_df$group[freqMatrix_df$time == "WT"] <- NA

pdf(paste0(output_path, "frequency_plots.pdf"))
ggplot(freqMatrix_df,aes(x = time,y=value*100,color=major_type,group = group)) + geom_point() + geom_line() + facet_wrap("~major_type") +
  scale_x_discrete(guide = guide_axis(angle = 45)) + theme_bw() + scale_color_manual(values = darker_colors) + ylab("% total cells (mean of samples)")

ggplot(freqMatrix_df[freqMatrix_df$major_type %in% c("Luminal","Neuroendocrine","Neutrophil"),],aes(x = time,y=value*100,color=major_type,group = group)) + geom_point() + geom_line() +
  scale_x_discrete(guide = guide_axis(angle = 45)) + theme_bw() + scale_color_manual(values = darker_colors) + ylab("% total cells")
dev.off()

# Save immune cells ---------------------------------------


dir.create(paste0(output_path, "/tmp/"), recursive = T)
DefaultAssay(combined.metacells) <- "ATAC"

# add fragments files

new.paths <- list.files(output_path,
                        pattern='MC_atac_fragments.tsv.gz',
                        recursive = T,
                        full.names = T)
new.paths <- new.paths[grepl(x = new.paths,pattern = "MC_atac_fragments.tsv.gz")]
new.paths <- new.paths[!grepl(x = new.paths,pattern = "tbi|matching")]

immune.cells <- combined.metacells[,combined.metacells$major_type %in% c("T_cell",'B_cell',"Macrophage","Neutrophil")]

frags <- list() 
for (i in seq_along(new.paths)) {
  sample <- basename(dirname(new.paths[i]))
  # frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[i]) # update path
  cell_names <- gsub(paste0(sample,"_"), "", colnames(immune.cells)[immune.cells$Sample == sample])
  names(cell_names) <- colnames(immune.cells)[immune.cells$Sample == sample]
  frags[[i]] <- CreateFragmentObject(new.paths[i], cells = cell_names)
}
Fragments(immune.cells) <- frags
saveRDS(immune.cells,
        paste0(output_path, "immune.combined.metacells.rds"))


peaks <- CallPeaks(
  object = immune.cells,
  outdir = paste0(output_path, "/tmp/"),
  group.by = "major_type"
)


peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

saveRDS(peaks, paste0(output_path,"/peaks.rds"))


# Extract counts for the new peaks ----------------------------------------

macs2_counts <- FeatureMatrix(
  fragments = Fragments(immune.cells),
  features = peaks,
  cells = colnames(immune.cells)
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

immune.cells[["ATAC2"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  genome = 'mm10',
  min.cells = 1,
  fragments = Fragments(immune.cells),
  annotation = annotations
)

saveRDS(immune.cells,
        paste0(output_path, "immune.combined.metacells.peakCalling2.rds"))


# Annotate immune cells ---------------------------------------------------


immune.cells <- FindMultiModalNeighbors(immune.cells, 
                                        reduction.list = list("integrated_pca", "integrated_lsi"),
                                        dims.list = list(1:50, 2:50))

immune.cells <- FindClusters(immune.cells, 
                                   resolution = c(c(1:9)*0.1),
                                   graph.name = "wsnn") 

immune.cells <- RunUMAP(immune.cells,
                              nn.name = "weighted.nn", 
                              reduction.name = "wnn.umap", 
                              reduction.key = "wnnUMAP_")



pdf(paste0(output_path, "immune_subtypes_markers.pdf"))

DimPlot(immune.cells,reduction = "wnn.umap",group.by = "major_type",cols = colors)

DefaultAssay(immune.cells) <- "RNA"
Idents(immune.cells) <- "wsnn_res.0.4"

DimPlot(immune.cells,reduction = "wnn.umap",group.by = "wsnn_res.0.4")

VlnPlot(immune.cells, features = c("Cd8a", "Cd4", "Cd19",
                                   "Clec10a", "Ciita", "Flt3",
                                   "Trem2", "Gpnmb", "Lrp12",
                                   "Fn1", "Cd86",
                                   "Cxcr2", "Mreg"), pt.size = 0.001, ncol = 3)
FeaturePlot(immune.cells,
            features = c("Nfkb1","Hdc","Il1b","Siglecf","Trem1","Cxcr2","P2rx7","Mreg"),reduction = "wnn.umap")


# Neutrophils
VlnPlot(immune.cells[,immune.cells$wsnn_res.0.4 %in% c(0,2,4)], features = c("Nfkb1","Hdc","Il1b","Siglecf","Trem1","Cxcr2"),group.by = "wsnn_res.0.4")



Idents(immune.cells) <- 'wsnn_res.0.4'
immune.cells <- RenameIdents(immune.cells,
                             "0" = "Neutro_inf_low",
                             "2" = "Neutro_inf_high",
                             "1" = "Macro",
                             "3" = "Macro",
                             "5" = "Macro",
                             "6" = "CD8/NK",
                             "4" = "Neutro_Siglecf_high",
                             "7" = "CD4",
                             "8" = "B cells"
)

mycolors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
              '#E95C59', '#E59CC4', '#AB3282', '#BD956A', '#8C549C',
              '#9FA3A8', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
              '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
              '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
              '#968175')


immune.cells$fine_celltype = Idents(immune.cells)
Idents(immune.cells) <- "fine_celltype"
DimPlot(immune.cells,reduction = "wnn.umap",label = T,repel = T, cols = mycolors)
dev.off()

neutrophils <- immune.cells[,immune.cells$fine_celltype %in% c("Neutro_inf_low", "Neutro_inf_high", "Neutro_Siglecf_high")]

saveRDS(neutrophils,
        paste0(output_path, "neutrophils.rds"))

