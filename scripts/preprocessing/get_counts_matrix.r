counts_path <- as.character(commandArgs(TRUE)[1])
output_path <- as.character(commandArgs(TRUE)[2])

counts <- read.table(counts_path, skip = 1, header = T,check.names = F)
rownames_save <- counts$Geneid
counts <- counts[, c(7:ncol(counts))]
rownames(counts) <- rownames_save
colnames(counts) <- sapply(colnames(counts), function(i) unlist(strsplit(basename(i), '_sort_dedup.bam'))[[1]])
dim(counts)

counts = counts[grep("chrY|chrUn|random",rownames(counts),invert = T),]
dim(counts)

df_output <- cbind(data.frame(IDs = rownames(counts)), counts)
write.table(df_output, file=paste0(output_path, "raw_counts.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

