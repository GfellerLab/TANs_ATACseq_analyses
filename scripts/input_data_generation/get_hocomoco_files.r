
library(TFBSTools)
library(jsonlite)
library(dplyr)

#release data retrieved from here: https://zenodo.org/records/10012938 
hocomoco_path <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/SuperCellMultiomicsAnalyses/HOCOMOCO_v12_release/"
data_path <- "/mnt/curnagl/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/input_data/"

# Read annotation file ----------------------------------------------------

# Read the file line by line
lines <- readLines(paste0(hocomoco_path, "H12CORE_annotation.jsonl"))
if (length(lines) > 0 && nchar(trimws(lines[length(lines)])) == 0) {
  lines <- lines[-length(lines)]  # Remove the last empty line if present
}
json_data_list <- lapply(lines, fromJSON)

print(json_data_list[[1]])

jsob_subset <- lapply(json_data_list, function(x) {
  data.frame(name = x$name,
             tf = x$tf,
             collection = x$collection,
             subtype_order = x$subtype_order,
             datatype = x$datatype,
             quality = x$quality,
             motif_datatype = x$original_motif$motif_datatype,
             origin = x$original_motif$origin,
             original.motif.species = x$original_motif$species,
             tfclass_id = x$masterlist_info$tfclass_id,
             tfclass_superclass = x$masterlist_info$tfclass_superclass,
             tfclass_class = x$masterlist_info$tfclass_class,
             tfclass_family = x$masterlist_info$tfclass_family,
             tfclass_subfamily = ifelse(is.null(x$masterlist_info$tfclass_subfamily), "",x$masterlist_info$tfclass_subfamily ),
             human_gene_symbol = x$masterlist_info$species$HUMAN$gene_symbol,
             mouse_gene_symbol = ifelse(is.null(x$masterlist_info$species$MOUSE$gene_symbol), "",x$masterlist_info$species$MOUSE$gene_symbol)
  )
})
json_subset <- do.call(rbind, jsob_subset)

# mouse pfm list ----------------------------------------------------------

mouse_filtered <- as.data.frame(json_subset[json_subset$mouse_gene_symbol != "",] %>%
                                  group_by(tf) %>%                  
                                  filter(subtype_order == 0 | n() == 1) %>%  # Keep only rows where subtype_order is 0 or if there is only one row per group
                                  ungroup())   
dim(mouse_filtered)
openxlsx::write.xlsx(mouse_filtered, file = paste0(data_path, "HOCOMOCO_v12_annotation_mouse.xlsx"))
write.table(mouse_filtered, file = paste0(data_path, "HOCOMOCO_v12_annotation_mouse.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

pwf.objs <- list()
for(i in 1:nrow(mouse_filtered)){
  motif <- mouse_filtered$name[i]
  lines <- readLines(paste0(hocomoco_path, "jaspar/",motif,"_jaspar_format.txt"))
  
  # Split the file into motifs using the '>' symbol
  motif_blocks <- split(lines, cumsum(grepl("^>", lines)))
  
  # Iterate over motifs and create PFMatrix objects
  for (block in motif_blocks) {
    if (length(block) > 1) {
      motif_info <- strsplit(sub(">", "", block[1]), "\\s+")[[1]]
      motif_name <- motif_info[1]
      motif_matrix <- do.call(rbind, lapply(block[-1], function(row) as.numeric(strsplit(row, "\\s+")[[1]])))
      rownames(motif_matrix) <- c("A", "C", "G", "T")
      pwf.objs[[motif_name]] <- TFBSTools::PFMatrix(profileMatrix = motif_matrix, 
                                                    ID = motif_name, 
                                                    name = mouse_filtered$mouse_gene_symbol[i],
                                                    tags = list(tf =  mouse_filtered$mouse_gene_symbol[i],
                                                                family = mouse_filtered$tfclass_family[i],
                                                                sub.family = mouse_filtered$tfclass_subfamily[i])
      )
    }
  }
  
  pwf.objs <- Filter(Negate(is.null), pwf.objs)
}

pfmatrix_list <- do.call(TFBSTools::PFMatrixList, pwf.objs)

saveRDS(pfmatrix_list, file = paste0(data_path,"/H12CORE_mouse_pfm.rds"))


# Create mouse jaspar file for TOBIAS ------------------------------

output_file <- paste0(data_path, 'tobias_motifs_file.txt')
file.create(output_file)

for (i in 1:nrow(mouse_filtered)) {
  motif <- mouse_filtered$name[i]
  f <- paste0(hocomoco_path, "jaspar/",motif,"_jaspar_format.txt")
  lines <- readLines(f)
  con <- file(output_file, open = "a")  # append mode
  writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
  # writeLines("\n",con = con)  
  close(con)
}

# Create mouse homer file ------------------------------

output_file <- paste0(data_path, 'homer_motifs_file.txt')
file.create(output_file)

for (i in 1:nrow(mouse_filtered)) {
  motif <- mouse_filtered$name[i]
  f <- paste0(hocomoco_path, "homer/pvalue_0.001/",motif,"_homer_format_0.001.motif")
  lines <- readLines(f)
  con <- file(output_file, open = "a")  # append mode
  writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
  # writeLines("\n",con = con)  
  close(con)
}

# human pfm list ----------------------------------------------------------
human_filtered <- as.data.frame(json_subset %>%
                                  group_by(tf) %>%                  
                                  filter(subtype_order == 0 | n() == 1) %>%  # Keep only rows where subtype_order is 0 or if there is only one row per group
                                  ungroup())   
dim(human_filtered)
openxlsx::write.xlsx(human_filtered, file = paste0(data_path, "HOCOMOCO_v12_annotation_human.xlsx"))
write.table(human_filtered, file = paste0(data_path, "HOCOMOCO_v12_annotation_human.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


pwf.objs <- list()
for(i in 1:nrow(human_filtered)){
  motif <- human_filtered$name[i]
  lines <- readLines(paste0(hocomoco_path, "jaspar/",motif,"_jaspar_format.txt"))
  
  # Split the file into motifs using the '>' symbol
  motif_blocks <- split(lines, cumsum(grepl("^>", lines)))
  
  # Iterate over motifs and create PFMatrix objects
  for (block in motif_blocks) {
    if (length(block) > 1) {
      motif_info <- strsplit(sub(">", "", block[1]), "\\s+")[[1]]
      motif_name <- motif_info[1]
      motif_matrix <- do.call(rbind, lapply(block[-1], function(row) as.numeric(strsplit(row, "\\s+")[[1]])))
      rownames(motif_matrix) <- c("A", "C", "G", "T")
      pwf.objs[[motif_name]] <- TFBSTools::PFMatrix(profileMatrix = motif_matrix, 
                                                    ID = motif_name, 
                                                    name = human_filtered$human_gene_symbol[i],
                                                    tags = list(tf = human_filtered$human_gene_symbol[i],
                                                                family = human_filtered$tfclass_family[i],
                                                                sub.family = human_filtered$tfclass_subfamily[i])
      )
    }
  }
  
  pwf.objs <- Filter(Negate(is.null), pwf.objs)
}

pfmatrix_list <- do.call(TFBSTools::PFMatrixList, pwf.objs)

saveRDS(pfmatrix_list, file = paste0(data_path,"/H12CORE_human_pfm.rds"))

# Create mouse homer file ------------------------------

output_file <- paste0(data_path, 'homer_motifs_file_human.txt')
file.create(output_file)

for (i in 1:nrow(human_filtered)) {
  motif <- human_filtered$name[i]
  f <- paste0(hocomoco_path, "homer/pvalue_0.001/",motif,"_homer_format_0.001.motif")
  lines <- readLines(f)
  con <- file(output_file, open = "a")  # append mode
  writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
  # writeLines("\n",con = con)  
  close(con)
}
