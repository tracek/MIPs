library(plyr)


check_same_exom <- function(mip1, mip2) {
    start_equal <-  mip1$feature_start_position == mip2$feature_start_position
    end_equal <- mip1$feature_stop_position == mip2$feature_stop_position
    return (start_equal & end_equal)
}

check_probe_strands <- function(mip, mip_next, mip_next2) {
    if (mip$probe_strand == '-')
    {
        return (mip_next == '+' & mip_next2 == '+')
    }
    else if (mip$probe_strand == '+')
    {
        return (mip_next == '-' & mip_next2 == '-')
    }
    else
    {
        stop("Unknown probe strand")
    }
}

check_preconditions <- function(mip, mip_next, mip_next2) {
    same_exom <- check_same_exom(mip, mip_next) & check_same_exom(mip, mip_next2)
    if (same_exom == TRUE)
    {
        return(check_probe_strands(mip, mip_next, mip_next2))
    }
    else
    {
        return(FALSE)
    }
}

check_overlap <- function(mip, mip_next, mip_next2) {
    if (check_preconditions(mip, mip_next, mip_next2))
    {
        
    }
}

too_intronic_param <- 0

# Read the data
mips <-  read.table('test.txt', sep="\t", header=TRUE)

mips  <- head(mips,17)

# Step 1: Remove duplicates 
mipsNoDup <- mips[!duplicated(mips[,c('ext_probe_start','ext_probe_stop','lig_probe_start', 'lig_probe_stop')]),]

# Step 2: Exclude high copy count
condition_too_high_copy_count <-  (mipsNoDup$ext_copy_count > 100 | mipsNoDup$lig_copy_count > 100) |
    (mipsNoDup$ext_copy_count > 5 & mipsNoDup$lig_copy_count > 5)
excluded_too_high_copy_count <- mipsNoDup[condition_too_high_copy_count,]
mipsNoHighCopy <- mipsNoDup[!condition_too_high_copy_count,]

# Step 3: Exclude too far intronic
condition_too_intronic <- (mipsNoHighCopy$feature_start_position + too_intronic_param) > mipsNoHighCopy$mip_target_stop_position |
    (mipsNoHighCopy$feature_stop_position - too_intronic_param) < mipsNoHighCopy$mip_target_start_position
excluded_too_intronic <- mipsNoHighCopy[condition_too_intronic,]
mipsTooIntronic <- mipsNoHighCopy[!condition_too_intronic,] 

a <- ddply(mipsTooIntronic,c("feature_start_position", "feature_stop_position"),
    function(x) {
        print(nrow(x))
    })

