if (!require(plyr)){ 
    install.packages(plyr) 
} 

input_file_name <- "test.txt"
too_intronic_param <- 0

# Functions' definition

check_probe_strands <- function(mip, mip_next, mip_next2)
{
    if (mip$probe_strand == '-')
    {
        return (mip_next$probe_strand == '+' & mip_next2$probe_strand == '+')
    }
    else if (mip$probe_strand == '+')
    {
        return (mip_next$probe_strand == '-' & mip_next2$probe_strand == '-')
    }
    else if (mip$probe_strand == '+')
    {
        return (mip_next$probe_strand == '+' & mip_next2$probe_strand == '-')
    }
    else if (mip$probe_strand == '-')
    {
        return (mip_next$probe_strand == '-' & mip_next2$probe_strand == '+')
    }
    else
    {
        stop("Unknown probe strand")
    }
}

find_missing_exon_cover <- function(exon) 
{
    # result <- numeric()
    start <- numeric()
    end <- numeric()
    
    mips_total <- nrow(exon)
    exon_start <- exon[1,]$feature_start_position
    exon_stop <- exon[mips_total,]$feature_stop_position
    
    mip_first <- exon[1,]
    mip_last  <- exon[mips_total,]
    chr <- mip_first$chr
    
    # Chech start
    if (exon_start < mip_first$mip_target_start_position)
    {
        start <- c(start, exon_start)
        end <- c(end, mip_first$mip_target_start_position)
        # result <- c(result, exon_start, mip_first$mip_target_start_position)
        # write_missing_mips(chr, exon_start, mip_first$mip_target_start_position)
    }
    if (exon_stop > mip_last$mip_target_stop_position)
    {
        start <- c(start, exon_stop)
        end <- c(end, mip_last$mip_target_stop_position)
        # result <- c(result, exon_stop, mip_last$mip_target_stop_position)
        # write_missing_mips(chr, exon_stop, mip_last$mip_target_stop_position)
    }
    
    if (mips_total > 1)
    {   
        for(mip_no in 1:(mips_total - 1))
        {
            mip <- exon[mip_no,]
            mip_next <- exon[mip_no + 1,]
            
            if (mip$mip_target_stop_position < mip_next$mip_target_start_position)
            {
                start <- c(start, mip$mip_target_stop_position)
                end <- c(end, mip_next$mip_target_start_position)
                # result <- c(result, mip$mip_target_stop_position, mip_next$mip_target_start_position)
                # write_missing_mips(chr, mip$mip_target_stop_position, mip_next$mip_target_start_position)
            }
        }
    }
    
    missing_regions_total <- length(start)
    if (missing_regions_total > 0) # We found some without coverege
    {
        chr_column <- rep(paste("chr", chr, sep=""), missing_regions_total) # Make column with chr
        df <- data.frame(chr_column, start, end)
        write.table(df, file="exons_not_covered.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    }
}

check_exon_covered <- function(exon) 
{
    mips_total <- nrow(exon)
    if (mips_total > 1)
    {
        mip_first <- exon[1,]
        mip_last  <- exon[mips_total,]
        
        if (mip_first$mip_target_start_position > mip_first$feature_start_position |
                mip_last$mip_target_stop_position < mip_first$feature_stop_position)
        {
            return(FALSE)
        }
        
        for(mip_no in 1:(mips_total - 1)) {
            mip <- exon[mip_no,]
            mip_next <- exon[mip_no + 1,]
            
            if (mip$mip_target_stop_position < mip_next$mip_target_start_position)
            {
                return(FALSE)
            }
        }
    }
    
    return(TRUE)
}

check_excessive_mips <- function(mip_1, mip_2, mip_3, end) 
{
    if (mip_1$mip_target_stop_position >= mip_2$mip_target_start_position &
            mip_1$mip_target_stop_position >= mip_3$mip_target_start_position)
    {
        if (mip_2$mip_target_stop_position >= end &
                mip_3$mip_target_stop_position >= end)
        {
            if (mip_2$rank_score > mip_3$rank_score) 
            {
                excessive_mip <- mip_3$X.mip_pick_count
            }
            else if (mip_2$rank_score > mip_3$rank_score)
            {
                excessive_mip <- mip_2$X.mip_pick_count
            }
            else # Situation is ambivalent - take out the second one
            {
                excessive_mip <- mip_3$X.mip_pick_count
            }
            # cat("Excessive mip: ", excessive_mip, "\n")
            excessive_mips <<- c(excessive_mips, excessive_mip)
            return(TRUE)
        }
    }
    
    return(FALSE)
}

find_excessive_mips_in_exon <- function(exon)
{
    if (check_exon_covered(exon) == TRUE) 
    {
        mips_total <- nrow(exon)
        if (mips_total >= 3)
        {
            for(row_no in 1:(mips_total - 2))
            {
                mip <- exon[row_no,]
                mip_2 <- exon[row_no + 1,]
                mip_3 <- exon[row_no + 2,]
                if (check_probe_strands(mip, mip_2, mip_3) == TRUE) 
                {
                    found <- FALSE
                    if ((row_no + 2) == mips_total)
                    {
                        found <- check_excessive_mips(mip, mip_2, mip_3, mip$feature_stop_position)
                    }
                    else
                    {
                        mip_4 <- exon[row_no + 3,]
                        found <- check_excessive_mips(mip, mip_2, mip_3, mip_4$mip_target_start_position)
                    }
                    
                    if (found)
                    {
                        if (check_exon_covered(exon) == FALSE)
                        {
                            stop("Something went terribly wrong - we removed a much needed MIP! We are DOOMED!")
                        }
                    }
                }
            }
        }
    }
}

save_mips <- function(mips, name)
{
    write.table(mips, file=name, sep="\t", row.names=FALSE)
}

#############################################################################


cat("track name=MIP_candidates description=\"MIP_candidates\" visibility=1 color=0,255,0\n", file="exons_not_covered.txt")
excessive_mips <- numeric()

# Read the data
mips <- read.table(input_file_name, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Step 1: Remove duplicates 
condition_duplicated <- duplicated(mips[,c('ext_probe_start','ext_probe_stop','lig_probe_start', 'lig_probe_stop')])
mips_duplicated <- mips[condition_duplicated,]
save_mips(mips_duplicated, "01_mips_duplicated.txt")
mips_no_dup <- mips[!condition_duplicated,]

# Step 2: Exclude high copy count
condition_too_high_copy_count <-  (mips$ext_copy_count > 100 | mips$lig_copy_count > 100) |
    (mips$ext_copy_count > 5 & mips$lig_copy_count > 5)
mips_too_high_copy_count <- mips[condition_too_high_copy_count,]
save_mips(mips_too_high_copy_count, "02_mips_too_high_copy_count.txt")
mips_no_high_copy <- mips[!condition_too_high_copy_count,]

# Step 3: Exclude too far intronic
condition_too_intronic <- (mips_no_high_copy$feature_start_position) > mips_no_high_copy$mip_target_stop_position |
                          (mips_no_high_copy$feature_stop_position) < mips_no_high_copy$mip_target_start_position
mips_too_intronic <- mips_no_high_copy[condition_too_intronic,]
save_mips(mips_too_intronic, "03_mips_too_intronic.txt")
mips_exonic <- mips_no_high_copy[!condition_too_intronic,] 

ddply(mips_exonic)

# Remove excessive MIPs
ddply(mips_exonic,                                           # Apply a function to mips_exonic
      c("feature_start_position", "feature_stop_position"),  # such that 'start' and 'end' are the same
      find_excessive_mips_in_exon,                           # Name of the function
      .inform = TRUE)                                        # Detailed error reporting - off for better performance

# 'find_excessive_mips_in_exon' saved results in 'excessive_mips' vector
# Now remove all excessive MIPs
condition_mips_excessive <- mips_exonic$X.mip_pick_count %in% excessive_mips
mips_excessive <- mips_exonic[condition_mips_excessive,]
save_mips(mips_excessive, "04_mips_excessive.txt")
mips_excessive_removed <- mips_exonic[!condition_mips_excessive,]

# Save final result
save_mips(mips_excessive_removed, "mips_final.txt")

# Report exons that are not covered
ddply(mips_excessive_removed,
      c("feature_start_position", "feature_stop_position"),
      find_missing_exon_cover,
      .inform = TRUE)

print("DONE!")

