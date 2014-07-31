
check_exom_covered <- function(exom) {
    mips_total <- nrow(exom)
    if (mips_total > 1)
    {
        mip_first <- exom[1,]
        mip_last  <- exom[mips_total,]
        
        if (mip_first$mip_target_start_position > mip_first$feature_start_position |
                mip_last$mip_target_stop_position < mip_first$feature_end_position)
        {
            stop("Exom not covered - cannot recover")
        }
        
        for(mip_no in 1:(mips_total - 1)) {
            mip <- exom[mip_no,]
            mip_next <- exom[mip_no + 1,]
            
            if (mip$mip_target_stop_position < mip_next$mip_target_stop_position)
            {
                return(FALSE)
            }
        }
    }
    
    return(TRUE)
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

a <- ddply(mips,c("feature_start_position", "feature_stop_position"),
           function(exom) {
               print(nrow(exom))
               if (check_exom_covered(exom) == FALSE)
               {
                   stop("Exom not covered - cannot recover")
               }
               return(exom)
           })