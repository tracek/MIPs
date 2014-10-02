if (!require("plyr")){ 
    install.packages("plyr") 
} 

input_file_name <- file.choose()
too_intronic_param <- 0
overlap <- 200

# Functions' definition

check_probe_strands <- function(mip, mip_2, mip_3)
{
    if (mip$probe_strand == '-')
    {
        return (mip_3$probe_strand == '+')
    }
    else if (mip$probe_strand == '+')
    {
        return (mip_3$probe_strand == '-')
    }
    else
    {
        stop("Unknown probe strand")
    }
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

check_exon_covered <- function(exon) 
{   
    mips_total <- nrow(exon)
    
    if (mips_total > 1)
    {
        exon <- exon[order(exon$mip_target_start_position),]
        exon_start <- exon[1,]$feature_start_position
        exon_stop <- exon[mips_total,]$feature_stop_position        
        
        starts_of_all_mips_in_exon <- exon[1:mips_total, "mip_target_start_position"]
        ends_of_all_mips_in_exon <- exon[1:mips_total, "mip_target_stop_position"]
        
        if (all(exon_start < starts_of_all_mips_in_exon))
        {
            return(FALSE)
        }
        
        if (all(exon_stop > ends_of_all_mips_in_exon))
        {
            return(FALSE)
        }        
        
        
        for(mip_no in 1:(mips_total - 1)) 
        {
            mip <- exon[mip_no,]
            
            starts_of_all_subsequent_mips <- exon[(mip_no + 1):mips_total,"mip_target_start_position"]
            
            if (all(mip$mip_target_stop_position < starts_of_all_subsequent_mips))
            {
                return(FALSE)
            }
        }
    }
    
    return(TRUE)
}


eliminate_mip_based_on_rank <- function(mip_1, mip_2)
{
    if (mip_1$rank_score > mip_2$rank_score) 
    {
        return(mip_2$X.mip_pick_count)
    }
    else if (mip_1$rank_score > mip_2$rank_score)
    {
        return(mip_1$X.mip_pick_count)
    }
    else # Situation is ambivalent - take out the second one
    {
        return(mip_2$X.mip_pick_count)
    }       
}

# Assign numbers to signs to find which form majority
transform_sign_into_number <- function(sign)
{
    if (sign == "+")
    {
        return(1)
    }
    else if (sign == "-")
    {
        return(-1)
    }
    else if (sign == 0)
    {
        return(0)
    }
    else
    {
        stop("Unknown sign!")
    }
}

# Sum the signs
sum_signs <- function(...) # Takes variable list of arguments
{
    signs <- list(...) # make a list of the input arguments
    list_of_numbers <- lapply(signs, transform_sign_into_number) # Transform signs into numbers
    total <- do.call(sum, list_of_numbers) # Sum the numbers
    return(total)
}

 
eliminate_mip <- function(mip_1, mip_2, mip_0, mip_3)
{
    # Find leading and lagging probe strand
    lagging_probe_strand <- mip_0$probe_strand
    if (missing(mip_3)) 
    {
        # It's the END! of an exon.
        leading_probe_strand <- 0
    }
    else
    {
        leading_probe_strand <- mip_3$probe_strand
    }
    
    # If they have identical signs, then only rank matters
    if (mip_1$probe_strand == mip_2$probe_strand)
    {
        mip_to_eliminate <- eliminate_mip_based_on_rank(mip_1, mip_2)
        return(mip_to_eliminate)
    }
    else # Signs are different - try to get opposite signs overlap 
    {
        # The sign that forms majority is the one to eliminate
        probe_strands_sum <- sum_signs(lagging_probe_strand, mip_2$probe_strand, mip_0$probe_strand, leading_probe_strand)
        if (probe_strands_sum > 0) # majority is "+"
        {
            if (mip_1$probe_strand == "+")
            {
                return(mip_1$X.mip_pick_count)
            }
            else
            {
                return(mip_2$X.mip_pick_count)
            }
        }
        else if (probe_strands_sum < 0) # majority is "-"
        {
            if (mip_1$probe_strand == "-")
            {
                return(mip_1$X.mip_pick_count)
            }
            else
            {
                return(mip_2$X.mip_pick_count)
            }            
        }
        else # same amount of "+" and "-"
        {
            mip_to_eliminate <- eliminate_mip_based_on_rank(mip_1, mip_2)
            return(mip_to_eliminate)        
        }
    }
}


check_excessive_mips <- function(mip_1, mip_2, mip_3, mip_4) 
{
    if (mip_1$mip_target_stop_position >= mip_2$mip_target_start_position &
        mip_1$mip_target_stop_position >= mip_3$mip_target_start_position)
    {
        if (mip_2$mip_target_stop_position >= mip_4$mip_target_start_position)
        {
            if (mip_3$mip_target_stop_position >= mip_4$mip_target_start_position)
            {
                excessive_mip <<- eliminate_mip(mip_2, mip_3, mip_1, mip_4)
                excessive_mips <<- c(excessive_mips, excessive_mip)
                return(excessive_mip)
            }
            else
            {
                excessive_mips <<- c(excessive_mips, mip_3$X.mip_pick_count)
                return(excessive_mip)              
            }
        }
    }
    
    return(NULL)
}

check_excessive_mips_end <- function(mip_1, mip_2, mip_3) 
{
    if (mip_1$mip_target_stop_position >= mip_2$mip_target_start_position &
        mip_1$mip_target_stop_position >= mip_3$mip_target_start_position)
    {
        exon_stop <- mip_3$feature_stop_position
        if (mip_2$mip_target_stop_position > exon_stop &
            mip_3$mip_target_stop_position > exon_stop)
        {
            excessive_mip <<- eliminate_mip(mip_2, mip_3, mip_1)
            excessive_mips <<- c(excessive_mips, excessive_mip)
            return(excessive_mip)            
        }
        else if (mip_2$mip_target_stop_position > exon_stop)
        {
            excessive_mip <- mip_3$X.mip_pick_count
            excessive_mips <<- c(excessive_mips, excessive_mip)
            return(excessive_mip)                    
        }
        else if (mip_3$mip_target_stop_position > exon_stop)
        {
            excessive_mip <- mip_2$X.mip_pick_count
            excessive_mips <<- c(excessive_mips, excessive_mip)
            return(excessive_mip)                    
        }
        else
        {
            return(NULL)
        }
    }
    
    return(NULL)
}

find_excessive_mips_in_exon <- function(exon)
{
    # Sort first by mip_target_start_position
    exon <- exon[order(exon$mip_target_start_position),]

    # Check if exon is covered - otherwise skip removal of excessive MIPs
    if (check_exon_covered(exon) == TRUE) 
    {
        mips_total <- nrow(exon) # Number of MIPs in exon
        if (mips_total >= 3)     # 
        {
            # Eliminate lagging excessive mips
            if (exon[1,"mip_target_start_position"] < exon[1,"feature_start_position"] &
                exon[2,"mip_target_start_position"] < exon[2,"feature_start_position"])
            {
                excessive_mips <<- c(excessive_mips, exon[1,"X.mip_pick_count"])
                exon <- exon[!(exon$X.mip_pick_count == exon[1,"X.mip_pick_count"]),]
                mips_total <- nrow(exon)
            }
            
            # Eliminate leading excessive mips
            if (exon[mips_total - 1, "mip_target_stop_position"] > exon[mips_total,"feature_stop_position"] &
                exon[mips_total    , "mip_target_stop_position"] > exon[mips_total,"feature_stop_position"])
            {
                excessive_mips <<- c(excessive_mips, exon[mips_total,"X.mip_pick_count"])
                exon <- exon[!(exon$X.mip_pick_count == exon[mips_total,"X.mip_pick_count"]),]
                mips_total <- nrow(exon)
            }            
            
            for(row_no in 1:(mips_total - 2))
            {
                if (row_no > mips_total - 2)
                {
                    break
                }
                mip <- exon[row_no,]
                mip_2 <- exon[row_no + 1,]
                mip_3 <- exon[row_no + 2,]
                
                if (check_probe_strands(mip, mip_2, mip_3) == TRUE) 
                {
                    if ((row_no + 2) == mips_total)
                    {
                        found <- check_excessive_mips_end(mip, mip_2, mip_3)
                    }
                    else
                    {
                        mip_4 <- exon[row_no + 3,]
                        found <- check_excessive_mips(mip, mip_2, mip_3, mip_4)
                    }
                    
                    if (length(found) > 0)
                    {
                        exon <- exon[!(exon$X.mip_pick_count == found),]
                        mips_total <- nrow(exon)
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

merge_exons <- function(exons)
{   
    mips_total <- nrow(exons)
    for (row_no in 1:(mips_total - 1))
    {
        mip <- exons[row_no,]
        mip_2 <- exons[row_no + 1,]
        
        
        if (mip$feature_start_position != mip_2$feature_start_position) # different exon 
        {
            if ((mip$feature_stop_position + overlap) > mip_2$feature_start_position)
            {
                # We have sufficient overlap - merge!
                new_exon_start <- mip$feature_start_position
                new_exon_stop <- mip_2$feature_stop_position
                
                exons[exons$feature_start_position == new_exon_start, "feature_stop_position"] <- new_exon_stop
                exons[exons$feature_stop_position == new_exon_stop, "feature_start_position"] <- new_exon_start
#                 cat("Merge ", mip$feature_stop_position, "id:", mip$X.mip_pick_count,
#                     " and", mip_2$feature_start_position, "id:", mip_2$X.mip_pick_count, "\n")
            }
        }
    }
    
    return(exons)
}
      
#############################################################################

excessive_mips <- numeric()

# Read the data
mips <- read.table(input_file_name, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Remove extension from input file name
input_file_name <- sub("^([^.]*).*", "\\1", input_file_name) 

filename <- paste(input_file_name, "exons_not_covered.txt", sep="_")
cat("track name=exons_not_covered description=\"not_covered\" visibility=1 color=0,255,0\n", file=filename)

# Sort the data
mips_sorted <- mips[order(mips$chr, mips$feature_start_position),]
filename <- paste(input_file_name, "mips_sorted.txt", sep="_")
save_mips(mips_sorted, filename)

# Step 1: Remove duplicates 
condition_duplicated <- duplicated(mips_sorted[,c('ext_probe_start','ext_probe_stop','lig_probe_start', 'lig_probe_stop')])
mips_duplicated <- mips_sorted[condition_duplicated,]
filename <- paste(input_file_name, "mips_duplicated.txt", sep="_")
save_mips(mips_duplicated, filename)
mips_no_dup <- mips_sorted[!condition_duplicated,]

# Step 2: Exclude high copy count 
condition_too_high_copy_count <-  (mips_no_dup$ext_copy_count > 10 | mips_no_dup$lig_copy_count > 10) |
  (mips_no_dup$ext_copy_count > 1 & mips_no_dup$lig_copy_count > 3)
mips_too_high_copy_count <- mips_no_dup[condition_too_high_copy_count,]
filename <- paste(input_file_name, "mips_too_high_copy_count.txt", sep="_")
save_mips(mips_too_high_copy_count, filename)
mips_no_high_copy <- mips_no_dup[!condition_too_high_copy_count,]

# Step 3: Exclude too far intronic
condition_too_intronic <- (mips_no_high_copy$feature_start_position) > mips_no_high_copy$mip_target_stop_position |
                          (mips_no_high_copy$feature_stop_position) < mips_no_high_copy$mip_target_start_position
mips_too_intronic <- mips_no_high_copy[condition_too_intronic,]
filename <- paste(input_file_name, "mips_too_intronic.txt", sep="_")
save_mips(mips_too_intronic, filename)
mips_exonic <- mips_no_high_copy[!condition_too_intronic,] 

# Merge exons that are close to each other within a chromosome
mips_exonic_merged <- ddply(mips_exonic, "chr", merge_exons, .inform = TRUE)

# Remove excessive MIPs
ddply(mips_exonic_merged,                                    # Apply a function to mips_exonic_merged
      c("feature_start_position", "feature_stop_position"),  # such that 'start' and 'end' are the same
      find_excessive_mips_in_exon,                           # Name of the function
      .inform = TRUE)                                        # Detailed error reporting - off for better performance

# 'find_excessive_mips_in_exon' saved results in 'excessive_mips' vector
# Now remove all excessive MIPs
condition_mips_excessive <- mips_exonic_merged$X.mip_pick_count %in% excessive_mips
mips_excessive <- mips_exonic_merged[condition_mips_excessive,]
filename <- paste(input_file_name, "mips_excessive.txt", sep="_")
save_mips(mips_excessive, filename)
mips_excessive_removed <- mips_exonic_merged[!condition_mips_excessive,]

condition_mips_duplicated <- duplicated(mips_excessive_removed[,c('mip_target_start_position','mip_target_stop_position')])
filename <- paste(input_file_name, "mips_duplicated.txt", sep="_")
save_mips(mips_excessive_removed[condition_mips_duplicated,], filename)
mips_duplicated_removed <- mips_excessive_removed[!condition_mips_duplicated,]

# Save final result
filename <- paste(input_file_name, "mips_final.txt", sep="_")
save_mips(mips_duplicated_removed, filename)

# Report exons that are not covered
ddply(mips_duplicated_removed,
      c("feature_start_position", "feature_stop_position"),
      find_missing_exon_cover,
      .inform = TRUE)

removed_total <- nrow(mips) - nrow(mips_duplicated_removed)
cat("Removed mips: ", removed_total, "\n")

go_to_bed(mips_duplicated_removed, input_file_name)

print("DONE!")

