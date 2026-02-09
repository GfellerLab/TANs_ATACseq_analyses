
# Parameters --------------------------------------------------------------
dar_res_path <- as.character(commandArgs(TRUE)[1])
comparison <- as.character(commandArgs(TRUE)[2])
direction <- as.character(commandArgs(TRUE)[3])
log_thr <- as.numeric(commandArgs(TRUE)[4])
output_path <- as.character(commandArgs(TRUE)[5])

dir.create(paste0(output_path, "/"))

# Read diff peaks ---------------------------------------------------------
dar_res <- read.table(dar_res_path, header = T, sep = "\t", quote = "")
colnames(dar_res)[1] = "peak_id" 

groups <- unlist(strsplit(comparison, "_VS_"))
cond1 <- groups[1]
cond2 <- groups[2]
print(cond1); print(cond2)

colnames <- colnames(dar_res)[which(grepl(cond1, colnames(dar_res)) & grepl(cond2, colnames(dar_res)))] 
adj_pval_col <- grep("adj..p.value", colnames, value = T)
logFC_col <- grep("Log2.Fold.Change", colnames, value = T)
conditions <- unlist(strsplit(logFC_col, ".vs.."))
reverse <- ifelse(cond1 == conditions[1], T, F)

if(reverse){
  if(direction == "up"){
    dar_subset <- dar_res[which(dar_res[, adj_pval_col] < 0.05 &
                                  dar_res[, logFC_col] <= -log_thr),]
    dim(dar_subset)
  }else {
    dar_subset <- dar_res[which(dar_res[, adj_pval_col] < 0.05 &
                                  dar_res[, logFC_col] >= log_thr), ]
    dim(dar_subset)
  }
  

  output_data <- data.frame(log2FC = dar_subset[, logFC_col],
                            p_val_adj= dar_subset[, adj_pval_col],
                            group1 = rep(cond2, nrow(dar_subset)),
                            group2 = rep(cond1, nrow(dar_subset)),
                            peak = dar_subset$peak_id)
}else{
  if(direction == "up"){
    dar_subset <- dar_res[which(dar_res[, adj_pval_col] < 0.05  &
                                  dar_res[, logFC_col] >= log_thr), ]
    dim(dar_subset)
  }else {
    
    dar_subset <- dar_res[which(dar_res[, adj_pval_col] < 0.05  &
                                  dar_res[, logFC_col] <= -log_thr),]
    dim(dar_subset)
    
  }
  
  output_data <- data.frame(log2FC = dar_subset[, logFC_col],
                            p_val_adj= dar_subset[, adj_pval_col],
                            group1 = rep(cond1, nrow(dar_subset)),
                            group2 = rep(cond2, nrow(dar_subset)),
                            peak = dar_subset$peak_id)
} 

write.table(output_data, file = paste0(output_path, comparison, "_", direction, "_logFC", log_thr, ".txt"), row.names = F, col.names = T, sep = "\t", quote = F)


