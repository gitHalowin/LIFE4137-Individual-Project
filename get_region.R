leaf_As75_peaks <- read.csv("filtered1/Anchored/leaf_As75/top_10_sig_SNP_positions.csv", header = F)
# make 500 bp windows around the peak locus
leaf_As75_peaks$start <- leaf_As75_peaks$V2-500
leaf_As75_peaks$end <- leaf_As75_peaks$V2+500
leaf_As75_peaks <- leaf_As75_peaks[,-2]
# standardize the variable "Chr" for intersecting with gff
leaf_As75_peaks$V1 <- paste0("Chr",leaf_As75_peaks$V1)
write.table(leaf_As75_peaks,"filtered1/Anchored/bed/leaf_As75_peaks.txt", row.names = F, col.names = F, quote = F)

leaf_Zn66_peaks <- read.csv("filtered1/Anchored/leaf_Zn66/top_10_sig_SNP_positions.csv", header = F)
leaf_Zn66_peaks$start <- leaf_Zn66_peaks$V2-500
leaf_Zn66_peaks$end <- leaf_Zn66_peaks$V2+500
leaf_Zn66_peaks <- leaf_Zn66_peaks[,-2]
leaf_Zn66_peaks$V1 <- paste0("Chr",leaf_Zn66_peaks$V1)
write.table(leaf_Zn66_peaks,"filtered1/Anchored/bed/leaf_Zn66_peaks.txt", row.names = F, col.names = F, quote = F)

leaf_K39_peaks <- read.csv("filtered2/Anchored/leaf_K39/top_10_sig_SNP_positions.csv", header = F)
leaf_K39_peaks$start <- leaf_K39_peaks$V2-500
leaf_K39_peaks$end <- leaf_K39_peaks$V2+500
leaf_K39_peaks <- leaf_K39_peaks[,-2]
leaf_K39_peaks$V1 <- paste0("Chr",leaf_K39_peaks$V1)
write.table(leaf_K39_peaks,"filtered2/Anchored/leaf_K39/leaf_K39_peaks.txt", row.names = F, col.names = F, quote = F)

leaf_Mn55_peaks <- read.csv("filtered2/Anchored/leaf_Mn55/top_10_sig_SNP_positions.csv", header = F)
leaf_Mn55_peaks$start <- leaf_Mn55_peaks$V2-500
leaf_Mn55_peaks$end <- leaf_Mn55_peaks$V2+500
leaf_Mn55_peaks <- leaf_Mn55_peaks[,-2]
leaf_Mn55_peaks$V1 <- paste0("Chr",leaf_Mn55_peaks$V1)
write.table(leaf_Mn55_peaks,"filtered2/Anchored/leaf_Mn55/leaf_Mn55_peaks.txt", row.names = F, col.names = F, quote = F)

seed_Mn55_peaks <- read.csv("filtered2/Anchored/seed_Mn55/top_10_sig_SNP_positions.csv", header = F)
seed_Mn55_peaks$start <- seed_Mn55_peaks$V2-500
seed_Mn55_peaks$end <- seed_Mn55_peaks$V2+500
seed_Mn55_peaks <- seed_Mn55_peaks[,-2]
seed_Mn55_peaks$V1 <- paste0("Chr",seed_Mn55_peaks$V1)
write.table(seed_Mn55_peaks,"filtered2/Anchored/seed_Mn55/seed_Mn55_peaks.txt", row.names = F, col.names = F, quote = F)


seed_Ni60_peaks <- read.csv("filtered2/Anchored/seed_Ni60/top_10_sig_SNP_positions.csv", header = F)
seed_Ni60_peaks$start <- seed_Ni60_peaks$V2-500
seed_Ni60_peaks$end <- seed_Ni60_peaks$V2+500
seed_Ni60_peaks <- seed_Ni60_peaks[,-2]
seed_Ni60_peaks$V1 <- paste0("Chr",seed_Ni60_peaks$V1)
write.table(seed_Ni60_peaks,"filtered2/Anchored/seed_Ni60/seed_Ni60_peaks.txt", row.names = F, col.names = F, quote = F)
