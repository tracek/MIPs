if (!require("plyr")){ 
  install.packages("plyr") 
} 

txt_files <- list.files(pattern = "\\.txt$")

input_file_name <- txt_files[1]

# Read the data
mips <- read.table(input_file_name, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Remove extension from input file name
input_file_name <- sub("^([^.]*).*", "\\1", input_file_name) 

# Sort the data
mips_sorted <- mips[order(mips$chr, mips$feature_start_position),]

to_be_copied <- c("Coverage.ext.probe", "Coverage.lig.probe", "Average.probe.coverage")
rest <- setdiff(colnames(mips), to_be_copied)

col_base <- mips_sorted[rest]

for (txt_file in txt_files) {
  mips <- read.table(txt_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  mips_sorted <- mips[order(mips$chr, mips$feature_start_position),]
  input_file_name <- sub("^([^.]*).*", "\\1", txt_file)
  fileid <- unlist(strsplit(txt_file, "_"))[1]
  cat("Processing", fileid, "\n")
  col_ext <- mips_sorted[to_be_copied]
  
  ext_probe_name_id <- paste(fileid, "Coverage.ext.probe", sep="_")
  lig_probe_name_id <- paste(fileid, "Coverage.lig.probe", sep="_")
  avg_probe_name_id <- paste(fileid, "Average.probe.coverage", sep="_")
  colnames(col_ext) <- c(ext_probe_name_id, lig_probe_name_id, avg_probe_name_id)
  
  col_base <- cbind(col_base, col_ext)
}

write.table(col_base, file="mips_merge.txt", sep="\t", row.names=FALSE)
