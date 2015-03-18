library("plyr")

go_to_bed <- function(exom, filename)
{
  filename <- paste(filename, "bed", sep=".")
  cat("track name=MIP_plus itemRgb=on\n", file=filename, append=FALSE)
  exom_plus_condition <- exom$probe_strand == "+"
  exom_plus <- exom[exom_plus_condition, ]
  mips_plus_total <- nrow(exom_plus)
  
  chr_plus <- paste("chr", exom_plus$chr, sep="")
  ext_probe_start_minus_1 <- exom_plus$ext_probe_start - 1
  lig_probe_stop <- exom_plus$lig_probe_stop
  mip_pick_count_rank_score <- paste(exom_plus$X.mip_pick_count, exom_plus$rank_score, sep="_")
  score <- rep(800, mips_plus_total)
  probe_strand <- rep("+", mips_plus_total)
  mip_target_start_position_minus_1 <- exom_plus$mip_target_start_position - 1
  mip_target_stop_position <- exom_plus$mip_target_stop_position
  RGB_plus <- rep("85,107,47", mips_plus_total)
  
  df_plus <- data.frame(chr_plus,
                        ext_probe_start_minus_1, 
                        lig_probe_stop,
                        mip_pick_count_rank_score,
                        score,
                        probe_strand,
                        mip_target_start_position_minus_1,
                        mip_target_stop_position,
                        RGB_plus)
  write.table(df_plus, file=filename, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  
  cat("\ntrack name=MIP_minus itemRgb=on\n", file=filename, append=TRUE)
  
  exom_minus <- exom[!exom_plus_condition, ]
  mips_minus_total <- nrow(exom_minus)
  
  chr_minus <- paste("chr", exom_minus$chr, sep="")
  lig_probe_start_plus_1 <- exom_minus$lig_probe_start + 1
  ext_probe <- exom_minus$ext_probe_stop
  mip_pick_count_rank_score <- paste(exom_minus$X.mip_pick_count, exom_minus$rank_score, sep="_")
  score <- rep(600, mips_minus_total)
  probe_strand <- rep("-", mips_minus_total)
  mip_target_start_position_minus_1 <- exom_minus$mip_target_start_position - 1
  mip_target_stop_position <- exom_minus$mip_target_stop_position
  RGB_minus <- rep("65,105,225", mips_minus_total)
  
  df_minus <- data.frame(chr_minus,
                         lig_probe_start_plus_1, 
                         ext_probe,
                         mip_pick_count_rank_score,
                         score,
                         probe_strand,
                         mip_target_start_position_minus_1,
                         mip_target_stop_position,
                         RGB_minus)
  write.table(df_minus, file=filename, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")    
}

find_missing_exon_cover <- function(exon) 
{
  exon <- exon[order(exon$mip_target_start_position),]
  
  start <- numeric()
  end <- numeric()
  
  mips_total <- nrow(exon)
  exon_start <- exon[1,]$feature_start_position
  exon_stop <- exon[mips_total,]$feature_stop_position
  
  mip_first <- exon[1,]
  mip_last  <- exon[mips_total,]
  chr <- mip_first$chr
  
  # Check start
  
  starts_of_all_mips_in_exon <- exon[1:mips_total, "mip_target_start_position"]
  ends_of_all_mips_in_exon <- exon[1:mips_total, "mip_target_stop_position"]
  
  if (all(exon_start < starts_of_all_mips_in_exon))
  {
    start <- c(start, exon_start)
    end_temp <- exon_start + min(starts_of_all_mips_in_exon - exon_start)
    end <- c(end, end_temp)
  }
  if (all(exon_stop > ends_of_all_mips_in_exon))
  {
    end <- c(end, exon_stop)
    start_temp <- exon_stop - min(exon_stop - ends_of_all_mips_in_exon)
    start <- c(start, start_temp)
  }
  
  if (mips_total > 1)
  {   
    for(mip_no in 1:(mips_total - 1))
    {
      mip <- exon[mip_no,]
      
      starts_of_all_subsequent_mips <- exon[(mip_no + 1):mips_total, "mip_target_start_position"]
      
      if (all(mip$mip_target_stop_position < starts_of_all_subsequent_mips))
      {
        if (mip$mip_target_stop_position == 17256698)
        {
          print(mip$mip_target_stop_position)
        }
        
        start <- c(start, mip$mip_target_stop_position)
        end_temp <- mip$mip_target_stop_position + min(starts_of_all_subsequent_mips - mip$mip_target_stop_position)
        end <- c(end, end_temp)
      }
    }
  }
  
  missing_regions_total <- length(start)
  if (missing_regions_total > 0) # We found some without coverege
  {
    chr_column <- rep(paste("chr", chr, sep=""), missing_regions_total) # Make column with chr
    df <- data.frame(chr_column, start, end)
    filename <- paste(input_file_name, "exons_not_covered.txt", sep="_")
    write.table(df, file=filename, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  }
}

input_file_name <- file.choose()

mips <- read.table(input_file_name, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Remove extension from input file name
input_file_name <- sub("^([^.]*).*", "\\1", input_file_name) 

# Sort the data
mips_sorted <- mips[order(mips$chr, mips$feature_start_position),]

# Report exons that are not covered
ddply(mips_sorted,
      c("feature_start_position", "feature_stop_position"),
      find_missing_exon_cover,
      .inform = TRUE)

go_to_bed(mips_sorted, input_file_name)
