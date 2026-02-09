
######################################################
# Get consensus peaks
######################################################
# Libraries ---------------------------------------------------------------
require(readxl)


# Source functions  -------------------------------------------------------
source("/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/workflow/scripts/preprocessing/consensus_peaks_functions.r")

# Parameters --------------------------------------------------------------

min_score <- as.numeric(commandArgs(TRUE)[1])
organism <- as.character(commandArgs(TRUE)[2])
pep <- as.character(commandArgs(TRUE)[3])
chrom_sizes_path <- as.character(commandArgs(TRUE)[4])
blacklist_file <- as.character(commandArgs(TRUE)[5])
data_path <- as.character(commandArgs(TRUE)[6])
sampleList <- as.character(commandArgs(TRUE)[7])
output_path <- as.character(commandArgs(TRUE)[8])
bampe_mode <- as.character(commandArgs(TRUE)[9])

group_id <- "groups"

dir.create(output_path, recursive = T)

if(organism != "mouse"){
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
}else{
  genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
}


# Read samples metadata ---------------------------------------------------
samples <- read.table(sampleList, header = F)
colnames(samples) <- "sample_name"
samples$patient <- gsub("-B|-L|-T|-Bm|-Tm|-1|-2", "", samples$sample_name)
samples$patient <- gsub("GC", "G-LC", samples$patient)
samples$groups <- gsub(".*-(B|T|L|Bm|Tm).*", "\\1", samples$sample_name)
samples$groups <- plyr::revalue(samples$groups, c("B"="Blood","Bm"="Blood", "L"="Adjacent_lung", "T"="Tumor","Tm"="Tumor"))
samples <- samples[1:19,]

stats_files = list.files(
  data_path,
  pattern = "stats.tsv",
  recursive = T,
  full.names = T
)

summary_df = data.frame(sample=vector(), NRF=vector(), PBC1=vector(), PBC2=vector(), TSS_score=vector(), FRiP=vector(), Peak_count=vector(), Alignment_rate=vector(), Alignment_rate_chrM2x= vector())
for(i in stats_files){
  sample = unlist(strsplit(i, "/"))
  sample = sample[length(sample)-1]
  stats_df = read.table(i, header = F, sep = "\t")
  summary_df = rbind(summary_df,data.frame(sample=sample,
                                           NRF=as.numeric(stats_df$V2[which(stats_df$V1 == "NRF")]),
                                           PBC1=as.numeric(stats_df$V2[which(stats_df$V1 == "PBC1")]),
                                           PBC2=as.numeric(stats_df$V2[which(stats_df$V1 == "PBC2")]),
                                           TSS_score=as.numeric(stats_df$V2[which(stats_df$V1 == "TSS_score")]),
                                           FRiP=as.numeric(stats_df$V2[which(stats_df$V1 == "FRiP")]),
                                           Peak_count=as.numeric(stats_df$V2[which(stats_df$V1 == "Peak_count")]),
                                           Alignment_rate=as.numeric(stats_df$V2[which(stats_df$V1 == "Alignment_rate")])))
}

samples_metadata <- merge(samples, summary_df, by.x = "sample_name", by.y = "sample")
if(bampe_mode != "bampe"){
  samples_metadata$summits <- paste0(data_path, "/",samples_metadata$sample_name, "/peak_calling_hg38/", samples_metadata$sample_name,"_peaks_normalized.narrowPeak")
}else{
  samples_metadata$summits <- paste0(data_path, "/",samples_metadata$sample_name, "/peak_calling_bampe_hg38/", samples_metadata$sample_name,"_peaks_fixedWidth.narrowPeak")
}
samples_metadata$organism <- organism
samples_metadata$biological_rep <- paste0(samples_metadata$patient, ".", samples_metadata$groups)
samples_metadata$technical_rep <- c(rep(1,7), 2, 1,1,2,1,2,1,1,2,1,1,1)
write.table(samples_metadata, file = paste0(output_path, "metadata.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


# Extract reproducible set of peaks for each cell protocol --------------------

samples_metadata <- samples_metadata[!samples_metadata$sample_name %in% c("G-LC-09-B","G-LC-09-L","G-LC-09-T"),]

all_grList <- list()
group_id <- "biological_rep"
group_names <- unique(samples_metadata[,group_id])


for(group in group_names){
  print(paste0("group: ",group))

  samples <- samples_metadata$sample_name[which(samples_metadata[,group_id] == group)]
  write.table(samples_metadata[which(samples_metadata[, group_id] == group), ],
              file = paste0(output_path,"samples_annotations.csv"), row.names = F, col.names = T, quote = F, sep = ",")

  # Change project name in the pepatac config
  system(command = paste0('sed \'2s/.*/name: ',group,'/\' ',pep,' > ', output_path, group,'_config.yaml') )
  pep_group <- paste0(output_path, group,'_config.yaml')

  # Load the project metadata
  prj <- pepr::Project(pep_group)
  project_name    <- pepr::config(prj)$name
  project_samples <- pepr::sampleTable(prj)$sample_name
  sample_table    <- data.table::data.table(sample_name = pepr::sampleTable(prj)$sample_name,
                                            groups = sapply(pepr::sampleTable(prj)$sample_name, function(i)
                                              samples_metadata[which(samples_metadata$sample == i), group_id]),
                                            genome = pepr::sampleTable(prj)$genome,
                                            summits = sapply(pepr::sampleTable(prj)$sample_name, function(i)
                                              samples_metadata[which(samples_metadata$sample == i), "summits"])
  )

  print(project_samples)
  print(sample_table)

  sample_table$samples    <- sample_table$sample_name

  gr_merged_across_study  <- extendedPeakSet(df = as.data.frame(sample_table),
                                             BSgenome = genome,
                                             blacklist = blacklist_file,
                                             scorePerMillion = min_score,
                                             selectionRules = length(sample_table$samples), #ifelse(length(sample_table$samples) != 2,"round(n/2)", 2),
                                             extend = 250,
                                             final_converge = F
  )

  all_grList[[group]] <- gr_merged_across_study$groupGRlist[[group]]
  group_consensus_peaks <- data.frame(GeneID = paste0( seqnames(x = gr_merged_across_study$groupGRlist[[group]] ),":",
                                                       start(gr_merged_across_study$groupGRlist[[group]]),"-",
                                                       end(gr_merged_across_study$groupGRlist[[group]])),
                                      Chr = chrom(gr_merged_across_study$groupGRlist[[group]]),
                                      Start = start(gr_merged_across_study$groupGRlist[[group]]),
                                      End = end(gr_merged_across_study$groupGRlist[[group]]),
                                      Strand = rep(".",length(gr_merged_across_study$groupGRlist[[group]])))
  colnames(group_consensus_peaks) <- c("GeneID","Chr","Start","End","Strand")


}

group_id <- "groups"
group_names <- unique(samples_metadata[,group_id])

peaks_by_group <- list()
for(group in group_names){
	grList_subset <- all_grList[unique(samples_metadata[samples_metadata[,group_id] == group, "biological_rep"])]
	peaks_by_group[[group]] <- collapse_peaks_from_GRList(grList = grList_subset,
                                              scorePerMillion = min_score,
                                              group_name = group, selectionRules = 2)

}
# Final list of peaks (from all groups groups) ---------------------------
final_peaks     <- collapse_peaks_from_GRList(grList = peaks_by_group,
                                              scorePerMillion = 0,
                                              group_name = "grouped", selectionRules = 1)
names(final_peaks) <- paste0(seqnames(final_peaks), ":", start(final_peaks), "-", end(final_peaks))
final_peaks_df = as.data.frame(final_peaks)
final_peaks_df$GeneID = rownames(final_peaks_df)
final_peaks_df$Strand = "."
colnames(final_peaks_df)[1:3] = c("Chr","Start","End")

summary(final_peaks_df$score)
write.table(final_peaks_df[,c("Chr","Start","End","Strand")] , file = paste0(output_path, "peaks.bed"), quote = F, row.names = F, col.names = F, sep = "\t")
write.table(final_peaks_df[,c("GeneID","Chr","Start","End","Strand")] , file = paste0(output_path, "peaks.saf"), quote = F, row.names = F, col.names = T, sep = "\t")

saveRDS(all_grList, file = paste0(output_path,"all_peaks.rds"))
