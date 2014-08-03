library(plyr)


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
    else
    {
        stop("Unknown probe strand")
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
            cat("Excessive mip: ", excessive_mip, "\n")
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
    else
    {
        # print(exon$X.mip_pick_count)
    }

    return(exon)
}

too_intronic_param <- 0

excessive_mips <- numeric()

# Read the data
mips <-  read.table('test.txt', sep="\t", header=TRUE, stringsAsFactors=FALSE)

# mips  <- head(mips,20)

# Step 1: Remove duplicates 
mips_no_dup <- mips[!duplicated(mips[,c('ext_probe_start','ext_probe_stop','lig_probe_start', 'lig_probe_stop')]),]

# Step 2: Exclude high copy count
condition_too_high_copy_count <-  (mips_no_dup$ext_copy_count > 100 | mips_no_dup$lig_copy_count > 100) |
    (mips_no_dup$ext_copy_count > 5 & mips_no_dup$lig_copy_count > 5)
excluded_too_high_copy_count <- mips_no_dup[condition_too_high_copy_count,]
mips_no_high_copy <- mips_no_dup[!condition_too_high_copy_count,]

# Step 3: Exclude too far intronic
condition_too_intronic <- (mips_no_high_copy$feature_start_position + too_intronic_param) > mips_no_high_copy$mip_target_stop_position |
                          (mips_no_high_copy$feature_stop_position - too_intronic_param) < mips_no_high_copy$mip_target_start_position
excluded_too_intronic <- mips_no_high_copy[condition_too_intronic,]
mips_exonic <- mips_no_high_copy[!condition_too_intronic,] 

# Remove excessive MIPs
ddply(mips_exonic,                                           # Apply a function to mips_exonic
      c("feature_start_position", "feature_stop_position"),  # such that 'start' and 'end' are the same
      find_excessive_mips_in_exon,                           # Name of the function
           .inform = TRUE)                                   # Detailed error reporting - off for better performance

# 'find_excessive_mips_in_exon' saved results in 'excessive_mips' vector
# Now remove all excessive MIPs
mips_excessive_removed <- mips_exonic[!mips_exonic$X.mip_pick_count %in% excessive_mips,]



