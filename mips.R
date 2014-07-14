too_intronic_param <- 0

# Read the data
mips <-  read.table('MIP_animal_candidates_new.70mers.txt', sep="\t", header=TRUE)

# Step 1: Remove duplicates 
mipsNoDup <- mips[!duplicated(mips[,c('ext_probe_start','ext_probe_stop','lig_probe_start', 'lig_probe_stop')]),]

# Step 2: Exclude high copy count
condition_too_high_copy_count = (mipsNoDup$ext_copy_count > 100 | mipsNoDup$lig_copy_count > 100) |
    (mipsNoDup$ext_copy_count > 5 & mipsNoDup$lig_copy_count > 5)
excluded_too_high_copy_count <- mipsNoDup[condition_too_high_copy_count,]
mipsNoHighCopy <- mipsNoDup[!condition_too_high_copy_count,]

# Step 3: Exclude too far intronic
condition_too_intronic <- (mipsNoHighCopy$feature_start_position + too_intronic_param) > mipsNoHighCopy$mip_target_stop_position |
    (mipsNoHighCopy$feature_stop_position - too_intronic_param) < mipsNoHighCopy$mip_target_start_position
excluded_too_intronic <- mipsNoHighCopy[condition_too_intronic,]
mipsTooIntronic <- mipsNoHighCopy[!condition_too_intronic,]