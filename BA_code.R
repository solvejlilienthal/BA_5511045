################################################################################

# Bachelor thesis - R script
# Script by Solvej Lilienthal
# last edited: 19.11.2023

################################################################################
################################################################################

# - unique refers to the clonotype approach
# - duplicate refers to the clones approach

# For protection of privacy purposes are the file names, donor sample ids and 
# directory pathways anonymized.

################################################################################
################################################################################

# used libraries
library(Seurat) 
library(BiocParallel)
library(tidyverse)
library(dplyr)
library(data.table)
library(plyr)

library(circlize)
library(eulerr)
library(car)
library(ggbeeswarm)
library(RColorBrewer)

# setting the work environment
setwd("/R")
loc_name = "/R/" 
load(file = paste0(loc_name,"MH_filtered_noXCR.simple.SO_all.data.RData"))

################################################################################

# Data preparation

################################################################################

# load the contig annotation file
isotype_all <- read.csv("/scratch/testuser34/MH_filtered_noXCR/BCR/outs_dummy_donor/vdj_b/filtered_contig_annotations.csv")

# order the data first by barcode, then umis (descending) and then reads (descending)
isotype_all <- isotype_all[order(isotype_all$barcode, -isotype_all$umis, -isotype_all$reads), ]  

# kepp only the cells with a complete BCR annonation (IgH and IgL chain)
igh <- filter(isotype_all, isotype_all$chain == "IGH")
igl <- filter(isotype_all, isotype_all$chain != "IGH")
fullBCR <- igh$barcode[igh$barcode %in% igl$barcode]

isotype <- colnames(SO_all)
clonotype <- colnames(SO_all)

# determine the clonotype of cells with a complete BCR 
# (clonotype of IgL == clonotype of IgH)
clonotype <- lapply(clonotype, function(clonotype) 
  ifelse(clonotype %in% fullBCR,
         isotype_all$raw_clonotype_id[match(clonotype, isotype_all$barcode)], 0))

# determine the isotype of the cell (C gene of the IgH contig)
isotype_filtered <- filter(isotype_all, (barcode %in% fullBCR) & (chain == "IGH"))
isotype <- lapply(isotype, function(isotype) 
  ifelse(isotype %in% fullBCR,
         isotype_filtered$c_gene[match(isotype, isotype_filtered$barcode)], 0))

clonotype <- unlist(clonotype)
isotype <- unlist(isotype)

# save both isotype and clonotype in the SO
SO <- SO_all
SO$clonotype <- clonotype
SO$isotype <- isotype

# save the isotype and barcode of each cell to load it into the 
# 10x Genomics Loupe Browser
# mark cells without an isotype that are found in the SO as contamination (NA)
iso_df <- data.frame(as.vector(colnames(SO)), SO$isotype)
colnames(iso_df) = c("barcode", "isotype")
iso_df = iso_df[sample(1:nrow(iso_df)),]
iso_df$isotype[iso_df$isotype ==  "" | iso_df$isotype == 0] = NA

fwrite(iso_df, file = "MH_isotype.csv")

################################################################################

# Overlap

################################################################################

# load the cell type categorization made in the 10x Genomics Loupe Browser
# add the experiment and clonotype of each cell to the data frame
forOverlap <- read_csv("ForOverlaps.csv")
forOverlap$Donor <- SO$experiment
forOverlap$Clonotype <- SO$clonotype

# filter out cells without complete BCR annotation and without GEX
forOverlap <- filter(forOverlap, forOverlap$Barcode %in% fullBCR)
forOverlap <- filter(forOverlap, forOverlap$ForOverlaps != "Contamination")

# split the samples by donors (according to the experiment/ sample id)
donor2 <- filter(forOverlap, forOverlap$Donor == "BM_PCs_2_83" | forOverlap$Donor == "BM_Bsm_2" | forOverlap$Donor == "Blood_Bsm_2")
donor3 <- filter(forOverlap, forOverlap$Donor == "BM_PCs_3_83" | forOverlap$Donor == "BM_Bsm_3" | forOverlap$Donor == "Blood_Bsm_3")
donor5 <- filter(forOverlap, forOverlap$Donor == "BM_PCs-Bmem_3_555" | forOverlap$Donor == "BM_PCs-Bmem_3_556-blood555")
donor6 <- filter(forOverlap, forOverlap$Donor == "Blood6_PCs_Bmem" | forOverlap$Donor == "BM_6_PC_Bmem")
donor7 <- filter(forOverlap, forOverlap$Donor == "BM_PCs_Bm_7" | forOverlap$Donor == "Blood_PCs_Bm_7")
donor8 <- filter(forOverlap, forOverlap$Donor == "BM_PCs_8" | forOverlap$Donor == "BM_Bmem_8" | forOverlap$Donor == "Blood_Bmem_8")
donor9 <- filter(forOverlap, forOverlap$Donor == "BM_PCs_9" | forOverlap$Donor == "BM_Bmem_9" | forOverlap$Donor == "Blood_Bmem_9")
donor10 <- filter(forOverlap, forOverlap$Donor == "BM_PCs_10" | forOverlap$Donor == "BM_Bmem_10" | forOverlap$Donor == "Blood_Bmem_10")

donorList <- list(donor2, donor3, donor5, donor6, donor7, donor8, donor9, donor10)
donorNumber <- c("2", "3", "5", "6", "7", "8", "9", "10")

# calculate the true overlap for both approaches
donorIndex = 1
for (donor in donorList) {
  BM_Bmem <- filter(donor, ForOverlaps == "BM_Bmem")
  BL_Bmem <- filter(donor, ForOverlaps == "Blood_Bmem")
  BM_PC <- filter(donor, ForOverlaps == "BM_PC")
  BL_PB <- filter(donor, ForOverlaps == "Blood_PB")
  
  # exclude cells with clonotype in Bmem overlap
  Bmem_Overlap <- intersect(BM_Bmem$Clonotype, BL_Bmem$Clonotype)
  
  BM_Bmem <- filter(BM_Bmem, !(Clonotype %in% Bmem_Overlap))
  BL_Bmem <- filter(BL_Bmem, !(Clonotype %in% Bmem_Overlap))
  BM_PC <- filter(BM_PC, !(Clonotype %in% Bmem_Overlap))
  BL_PB <- filter(BL_PB, !(Clonotype %in% Bmem_Overlap))
  
  
  # clonotype approach
  Bmem_Overlap <- intersect(BM_Bmem$Clonotype, BL_Bmem$Clonotype)
  PCB_Overlap <- intersect(BM_PC$Clonotype, BL_PB$Clonotype)
  BM_Overlap <- intersect(BM_Bmem$Clonotype, BM_PC$Clonotype)
  BL_Overlap <- intersect(BL_Bmem$Clonotype, BL_PB$Clonotype)
  BMPC_BLBmem_Overlap <- intersect(BM_PC$Clonotype, BL_Bmem$Clonotype)
  BMBmem_BLPB_Overlap <- intersect(BM_Bmem$Clonotype, BL_PB$Clonotype)
  
  originalOverlap_unique <- rbind(c(length(unique(BM_PC$Clonotype)), length(BM_Overlap), length(PCB_Overlap), length(BMPC_BLBmem_Overlap)),
                                  c(length(BM_Overlap), length(unique(BM_Bmem$Clonotype)), length(BMBmem_BLPB_Overlap), length(Bmem_Overlap)),
                                  c(length(PCB_Overlap), length(BMBmem_BLPB_Overlap), length(unique(BL_PB$Clonotype)), length(BL_Overlap)),
                                  c(length(BMPC_BLBmem_Overlap), length(Bmem_Overlap), length(BL_Overlap), length(unique(BL_Bmem$Clonotype))))
  
  colnames(originalOverlap_unique) <- c("BMPC", "BMBmem", "BLPB", "BLBmem")
  rownames(originalOverlap_unique) <- c("BMPC", "BMBmem", "BLPB", "BLBmem")
  write.table(originalOverlap_unique,  file = paste0("originalOverlaps_unique_D", donorNumber[donorIndex], ".txt", sep = ""), quote = FALSE, sep ="\t" ,row.names = TRUE, col.names = TRUE)
  
  
  # clones approach
  Bmem_Overlap_r <- BM_Bmem$Clonotype[BM_Bmem$Clonotype %in% BL_Bmem$Clonotype]
  PCB_Overlap_r <- BM_PC$Clonotype[BM_PC$Clonotype %in% BL_PB$Clonotype]
  BM_Overlap_r <- BM_Bmem$Clonotype[BM_Bmem$Clonotype %in% BM_PC$Clonotype]
  BL_Overlap_r <- BL_Bmem$Clonotype[BL_Bmem$Clonotype %in% BL_PB$Clonotype]
  BMPC_BLBmem_Overlap_r <- BM_PC$Clonotype[BM_PC$Clonotype %in% BL_Bmem$Clonotype]
  BMBmem_BLPB_Overlap_r <- BM_Bmem$Clonotype[BM_Bmem$Clonotype %in% BL_PB$Clonotype]
  
  Bmem_Overlap_l <- BL_Bmem$Clonotype[BL_Bmem$Clonotype %in% BM_Bmem$Clonotype]
  PCB_Overlap_l <- BL_PB$Clonotype[BL_PB$Clonotype %in% BM_PC$Clonotype]
  BM_Overlap_l <- BM_PC$Clonotype[BM_PC$Clonotype %in% BM_Bmem$Clonotype]
  BL_Overlap_l <- BL_PB$Clonotype[BL_PB$Clonotype %in% BL_Bmem$Clonotype]
  BMPC_BLBmem_Overlap_l <- BL_Bmem$Clonotype[BL_Bmem$Clonotype %in% BM_PC$Clonotype]
  BMBmem_BLPB_Overlap_l <- BL_PB$Clonotype[BL_PB$Clonotype %in% BM_Bmem$Clonotype]
  
  originalOverlap_duplicate <- rbind(c(length(BM_PC$Barcode), length(BM_Overlap_l), length(PCB_Overlap_r), length(BMPC_BLBmem_Overlap_r)),
                                     c(length(BM_Overlap_r), length(BM_Bmem$Barcode), length(BMBmem_BLPB_Overlap_r), length(Bmem_Overlap_r)),
                                     c(length(PCB_Overlap_l), length(BMBmem_BLPB_Overlap_l), length(BL_PB$Barcode), length(BL_Overlap_l)),
                                     c(length(BMPC_BLBmem_Overlap_l), length(Bmem_Overlap_l), length(BL_Overlap_r), length(BL_Bmem$Barcode)))
  
  colnames(originalOverlap_duplicate) <- c("BMPC", "BMBmem", "BLPB", "BLBmem")
  rownames(originalOverlap_duplicate) <- c("BMPC", "BMBmem", "BLPB", "BLBmem")
  write.table(originalOverlap_duplicate, file = paste0("originalOverlaps_duplicate_D", donorNumber[donorIndex], ".txt", sep = ""), quote = FALSE, sep ="\t", col.names = TRUE)
  
  donorIndex = donorIndex + 1
}

# calculate the randomized overlap for both approaches
donorIndex = 1
for (donor in donorList) {
  BM_Bmem_r <- filter(donor, ForOverlaps == "BM_Bmem")
  BL_Bmem_r <- filter(donor, ForOverlaps == "Blood_Bmem")
  BM_PC_r <- filter(donor, ForOverlaps == "BM_PC")
  BL_PB_r <- filter(donor, ForOverlaps == "Blood_PB")
  
  Bmem_Overlap <- intersect(BM_Bmem_r$Clonotype, BL_Bmem_r$Clonotype)
  
  BM_Bmem_r <- filter(BM_Bmem_r, !(Clonotype %in% Bmem_Overlap))
  BL_Bmem_r <- filter(BL_Bmem_r, !(Clonotype %in% Bmem_Overlap))
  BM_PC_r <- filter(BM_PC_r, !(Clonotype %in% Bmem_Overlap))
  BL_PB_r <- filter(BL_PB_r, !(Clonotype %in% Bmem_Overlap))
  
  originalData <- c(BM_Bmem_r$Clonotype, BL_Bmem_r$Clonotype, 
                    BM_PC_r$Clonotype, BL_PB_r$Clonotype)
  originalData <- originalData[!is.na(originalData)]
  
  randomisedOverlap_unique <- data.frame(matrix(ncol = 10, nrow = 0))
  randomisedOverlap_r <- data.frame(matrix(ncol = 10, nrow = 0))
  randomisedOverlap_l <- data.frame(matrix(ncol = 10, nrow = 0))
  
  set.seed(42)
  for (i in 1:1000) {
    randomOverlap <- sample(originalData) # randomize by clonotypes shuffling
    
    BM_Bmem_temp <- randomOverlap[1:length(BM_Bmem_r$Clonotype)]
    BL_Bmem_temp <- randomOverlap[(length(BM_Bmem_r$Clonotype) + 1):
                                    (length(BM_Bmem_r$Clonotype) + length(BL_Bmem_r$Clonotype))]
    BM_PC_temp <- randomOverlap[(length(BM_Bmem_r$Clonotype) + length(BL_Bmem_r$Clonotype) + 1):
                                  (length(originalData) - length(BL_PB_r$Clonotype))]
    BL_PB_temp <- randomOverlap[(length(originalData) - length(BL_PB_r$Clonotype) + 1):
                                  length(originalData)]
    
    
    # clonotype approach
    Bmem_Overlap <- intersect(BM_Bmem_temp, BL_Bmem_temp)
    PCB_Overlap <- intersect(BM_PC_temp, BL_PB_temp)
    BM_Overlap <- intersect(BM_PC_temp, BM_Bmem_temp)
    BL_Overlap <- intersect(BL_PB_temp, BL_Bmem_temp)
    BMPC_BLBmem_Overlap <- intersect(BM_PC_temp, BL_Bmem_temp)
    BMBmem_BLPB_Overlap <- intersect(BM_Bmem_temp, BL_PB_temp)
    
    randomisedOverlap_unique <- rbind(randomisedOverlap_unique, 
                                      c(length(Bmem_Overlap), length(PCB_Overlap), length(BM_Overlap), length(BL_Overlap),
                                        length(BMPC_BLBmem_Overlap), length(BMBmem_BLPB_Overlap)))
    
    # clones approach - calculate both sided overlaps (first AB; then BA)
    Bmem_Overlap_r <- BM_Bmem_temp[BM_Bmem_temp %in% BL_Bmem_temp]
    PCB_Overlap_r <- BM_PC_temp[BM_PC_temp %in% BL_PB_temp]
    BM_Overlap_r <- BM_Bmem_temp[BM_Bmem_temp %in% BM_PC_temp]
    BL_Overlap_r <- BL_Bmem_temp[BL_Bmem_temp %in% BL_PB_temp]
    BMPC_BLBmem_Overlap_r <- BM_PC_temp[BM_PC_temp %in% BL_Bmem_temp]
    BMBmem_BLPB_Overlap_r <- BM_Bmem_temp[BM_Bmem_temp %in% BL_PB_temp]
    
    randomisedOverlap_r <- rbind(randomisedOverlap_r, 
                                 c(length(Bmem_Overlap_r), length(PCB_Overlap_r), length(BM_Overlap_r), length(BL_Overlap_r),
                                   length(BMPC_BLBmem_Overlap_r), length(BMBmem_BLPB_Overlap_r)))
    
    Bmem_Overlap_l <- BL_Bmem_temp[BL_Bmem_temp %in% BM_Bmem_temp]
    PCB_Overlap_l <- BL_PB_temp[BL_PB_temp %in% BM_PC_temp]
    BM_Overlap_l <- BM_PC_temp[BM_PC_temp %in% BM_Bmem_temp]
    BL_Overlap_l <- BL_PB_temp[BL_PB_temp %in% BL_Bmem_temp]
    BMPC_BLBmem_Overlap_l <- BL_Bmem_temp[BL_Bmem_temp %in% BM_PC_temp]
    BMBmem_BLPB_Overlap_l <- BL_PB_temp[BL_PB_temp %in% BM_Bmem_temp]
    
    randomisedOverlap_l <- rbind(randomisedOverlap_l, 
                                 c(length(Bmem_Overlap_l), length(PCB_Overlap_l), length(BM_Overlap_l), length(BL_Overlap_l),
                                   length(BMPC_BLBmem_Overlap_l), length(BMBmem_BLPB_Overlap_l)))
  }
  
  colnames(randomisedOverlap_r) <- c("Bmem_Overlap", "PCB_Overlap", "BM_Overlap", "BL_Overlap", "BMPC_BLBmem_Overlap", "BMBmem_BLPB_Overlap")
  colnames(randomisedOverlap_l) <- c("Bmem_Overlap", "PCB_Overlap", "BM_Overlap", "BL_Overlap", "BMPC_BLBmem_Overlap", "BMBmem_BLPB_Overlap") 
  colnames(randomisedOverlap_unique) <- c("Bmem_Overlap", "PCB_Overlap", "BM_Overlap", "BL_Overlap", "BMPC_BLBmem_Overlap", "BMBmem_BLPB_Overlap") 
  
  file_name <- paste("randomisedOverlaps_r_D", donorNumber[donorIndex], ".txt", sep = "")
  write.table(randomisedOverlap_r,  file = file_name, quote = FALSE, sep ="\t", col.names = TRUE)
  
  file_name <- paste("randomisedOverlaps_l_D", donorNumber[donorIndex], ".txt", sep = "")
  write.table(randomisedOverlap_l,  file = file_name, quote = FALSE, sep ="\t", col.names = TRUE)
  
  file_name <- paste("randomisedOverlaps_unique_D", donorNumber[donorIndex], ".txt", sep = "")
  write.table(randomisedOverlap_unique,  file = file_name, quote = FALSE, sep ="\t", col.names = TRUE)
  
  donorIndex = donorIndex + 1
}

################################################################################

# Quality control - comparison

################################################################################

# Significance of difference between original and randomized overlaps
p_value_overlaps = as.data.frame(matrix(ncol = 4, nrow = 0))
for (donorIndex in 1:length(donorNumber)) {
  originalOverlap <- read.delim(paste0("originalOverlaps_unique_D", donorNumber[donorIndex], ".txt", sep = ""))
  randomisedOverlaps <- read.delim(paste0("randomisedOverlaps_unique_D", donorNumber[donorIndex], ".txt", sep = ""))
  
  smaller = c(length(randomisedOverlaps$BM_Overlap[randomisedOverlaps$BM_Overlap < originalOverlap[1,2]]),
              length(randomisedOverlaps$BMPC_BLBmem_Overlap[randomisedOverlaps$BMPC_BLBmem_Overlap < originalOverlap[1,4]]))
  
  greater = c(length(randomisedOverlaps$BM_Overlap[randomisedOverlaps$BM_Overlap > originalOverlap[1,2]]),
              length(randomisedOverlaps$BMPC_BLBmem_Overlap[randomisedOverlaps$BMPC_BLBmem_Overlap > originalOverlap[1,4]]))
  
  p_value_overlaps = rbind(p_value_overlaps, c(1 - abs(smaller - greater)/1000))
}

# Significance of difference between original and randomized overlaps (FC)
for (donorIndex in 1:length(donorNumber)) {
  originalOverlap <- read.delim(paste0("originalOverlaps_unique_D", donorNumber[donorIndex], ".txt", sep = ""))
  randomisedOverlaps <- read.delim(paste0("randomisedOverlaps_unique_D", donorNumber[donorIndex], ".txt", sep = ""))
  
  Bmem_Overlap <- ceiling(mean(randomisedOverlaps$Bmem_Overlap))
  PCB_Overlap <- ceiling(mean(randomisedOverlaps$PCB_Overlap))
  BM_Overlap <- ceiling(mean(randomisedOverlaps$BM_Overlap))
  BL_Overlap <- ceiling(mean(randomisedOverlaps$BL_Overlap))
  BMPC_BLBmem_Overlap <- ceiling(mean(randomisedOverlaps$BMPC_BLBmem_Overlap))
  BMBmem_BLPB_Overlap <- ceiling(mean(randomisedOverlaps$BMBmem_BLPB_Overlap))
  
  BMPC <- c(originalOverlap[1,1], BM_Overlap, PCB_Overlap, BMPC_BLBmem_Overlap)
  BMBmem <- c(BM_Overlap, originalOverlap[2,2], BMBmem_BLPB_Overlap, Bmem_Overlap)
  BLPB <- c(PCB_Overlap, BMBmem_BLPB_Overlap, originalOverlap[3,3], BL_Overlap)
  BLBmem <- c(BMPC_BLBmem_Overlap, Bmem_Overlap, BL_Overlap, originalOverlap[4,4])
  
  matrix = rbind(BMPC, BMBmem, BLPB, BLBmem)
  colnames(matrix) <- c("BMPC", "BMBmem", "BLPB", "BLBmem")
  
  write.table(matrix, file = paste0("randomisedOverlaps_mean_unique_D", donorNumber[donorIndex], ".txt", sep = ""), quote = FALSE, sep ="\t", col.names = TRUE)
  write.table(matrix/originalOverlap, file = paste0("FC_Overlap_unique_ro_D", donorNumber[donorIndex], ".txt", sep = ""), quote = FALSE, sep ="\t", col.names = TRUE)
}

################################################################################

# Visualization

################################################################################

# proportional Venn/ Euler diagram - unique (original & randomized)
for (donorIndex in 1:length(donorNumber)) {
  randomisedOverlaps_ <- read.delim(paste0("originalOverlaps_unique_D", donorNumber[donorIndex], ".txt", sep = ""))
  #randomisedOverlaps_ <- read.delim(paste0("randomisedOverlaps_mean_unique_D", donorNumber[donorIndex], ".txt", sep = ""))
  
  # percentage of clonotypes in organ A that are also in organ B (unique)
  fit <- euler(c(BM_Bmem = randomisedOverlaps_[2,2],
                 BM_PC = randomisedOverlaps_[1,1],
                 "BM_Bmem&BM_PC" = randomisedOverlaps_[1,2]),
               input = "union", shape = "ellipse")
  
  png(file = paste0("BM_Overlap_unique_D", donorNumber[donorIndex], ".png", sep = ""), width = 600, height = 350)
  print(plot(fit, 
             key = TRUE, counts = TRUE, 
             quantities = TRUE, list(type = c("counts", "percent"), font=3, round=1, cex=0.8), 
             fills = c("darkolivegreen2", "dodgerblue3"), 
             alpha = 0.6))
  dev.off()
  
  
  fit <- euler(c(BM_PC =  randomisedOverlaps_[1,1],
                 BL_Bmem =  randomisedOverlaps_[4,4],
                 "BM_PC&BL_Bmem" =  randomisedOverlaps_[1,4]),
               input = "union", shape = "ellipse")
  
  png(file = paste0("BMPC_BLBmem_Overlap_unique_D", donorNumber[donorIndex], ".png", sep = ""), width = 600, height = 350)
  print(plot(fit, 
             key = TRUE, counts = TRUE, 
             quantities = TRUE, list(type = c("counts", "percent"), font=3, round=1, cex=0.8), 
             fills = c("dodgerblue3", "orange"), 
             alpha = 0.6))
  dev.off()
  
  donorIndex = donorIndex + 1
}

# source: https://www.researchgate.net/post/Plotting_the_Euler_diagram_with_percentages_with_fractions_decimals_using_eulerr


# chord diagram - duplicated (proportional & normalized) & ordered - 100 normed
col = adjustcolor(c("lightblue2", "#8DB6CD", "#63B8FF", "#FFA07A", "orange2", "darkolivegreen3", "violetred1"), alpha.f = 0.3)
col_border = c("lightblue2", "#8DB6CD", "#63B8FF", "#FFA07A", "orange2", "darkolivegreen3", "violetred1")
for (donorIndex in 1:length(donorNumber)) {
  #overlap_matrix <- read.delim(paste0("originalOverlaps_duplicate_D", donorNumber[donorIndex], ".txt", sep = ""))
  overlap_matrix <- read.delim(paste0("randomisedOverlaps_mean_duplicate_D", donorNumber[donorIndex], ".txt", sep = ""))
  norm_matrix = rbind(c(rep(overlap_matrix[1,1], 4)), c(rep(overlap_matrix[2,2], 4)), c(rep(overlap_matrix[3,3], 4)), c(rep(overlap_matrix[4,4], 4)))
  overlap_matrix = overlap_matrix / norm_matrix
  
  png(file = paste0("chordDiagram_100_randomised_D", donorNumber[donorIndex], ".png", sep = ""), width = 2300, height = 2050, res = 300)
  
  # chord diagram
  circos.clear()
  sectors = c("BMPC", "BMBmem", "BLBmem")
  circos.par(cell.padding = c(0, 0, 0, 0), gap.after = c(1, 1, 1), gap.degree = 5, "track.height"=0.8)
  circos.initialize(factors = sectors, xlim = c(0, 1))
  circos.track(ylim = c(0,1),
               track.height = 0.06, bg.border = c("lightblue", "darkolivegreen3", "orange"),
               bg.col = c(BMBmem = "lightblue", BLBmem = "darkolivegreen3", BMPC = "orange"),
               panel.fun=function(x, y) {
                 chr=CELL_META$sector.index
                 xlim=CELL_META$xlim
                 ylim=CELL_META$ylim
                 circos.text(mean(xlim), mean(ylim), chr, cex=1, facing="bending.inside", niceFacing = TRUE)
                 circos.axis(h="top")})
  
  # BMPC - BMBmem
  circos.link(sectors[1], c(overlap_matrix[1,4], overlap_matrix[1,4] + overlap_matrix[1,2]), 
              sectors[2], c(0, overlap_matrix[2,1]), 
              col = col[6], border = col_border[6])
  
  # BMPC - BLBmem
  circos.link(sectors[1], c(0, overlap_matrix[1,4]), sectors[3], c(0, overlap_matrix[4,1]), 
              col = col[5], border = col_border[5])
  dev.off()
}


# histogram - unique
for (donorIndex in 1:length(donorNumber)) {
  originalOverlap_ <- read.delim(paste0("originalOverlaps_unique_D", donorNumber[donorIndex], ".txt", sep = ""))
  randomisedOverlap_ <- read.delim(paste0("randomisedOverlaps_unique_D", donorNumber[donorIndex], ".txt", sep = ""))
  
  
  png(file = paste0("BM_Overlap_histogram_unique_D", donorNumber[donorIndex], ".png", sep = ""), width = 600, height = 350)
  print(ggplot(randomisedOverlap_, aes(x = BM_Overlap)) + 
          geom_histogram(color="skyblue3", fill="lightblue") +
          geom_vline(xintercept = originalOverlap_[1,2], color = "red", linetype = "dashed", size = 1) +
          labs(title = paste0("BM Overlap: unique randomised vs original - Donor ", donorIndex),
               x = "BM Overlap", y = "Count"))
  dev.off()

  
  png(file = paste0("BMPC_BLBmem_Overlap_histogram_unique_D", donorNumber[donorIndex], ".png", sep = ""), width = 600, height = 350)
  print(ggplot(randomisedOverlap_, aes(x = BMPC_BLBmem_Overlap)) + 
          geom_histogram(color="skyblue3", fill="lightblue") +
          geom_vline(xintercept = originalOverlap_[1,2], color = "red", linetype = "dashed", size = 1) +
          labs(title = paste0("BMPC - BLBmem Overlap: unique randomised vs original - Donor ", donorIndex), 
               x = "BMPC_BLBmem Overlap", y = "Count"))
  dev.off()
  
  donorIndex = donorIndex + 1 
}

## FC beeplot
FC = as.data.frame(matrix(ncol = 2, nrow = 0))
for (donorIndex in 1:length(donorNumber)) {
    FC_temp <- read.delim(paste0("FC_Overlap_unique_or_D", donorNumber[donorIndex], ".txt", sep = ""))
    
    FC = rbind(FC, c(FC_temp[1,2], "BMPC-BMBmem"))
    FC = rbind(FC, c(FC_temp[1,4], "BMPC-BLBmem"))
}
colnames(FC) = c("fraction", "overlap")
FC$donor = as.factor(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8))
FC$fraction = as.numeric(FC$fraction)

FC_log = FC
FC_log$fraction = log(FC_log$fraction, 10)

FC = filter(FC, donor != 3 & donor != 4 & donor != 5)

png(file = "FC_overlap_beeplot_woBmem.png", width = 2300, height = 2050, res = 300)
ggplot(FC, aes(overlap, fraction, color = donor)) + 
  geom_beeswarm(cex = 3, size = 3) + 
  scale_color_manual(values=c("olivedrab3", "black", "#912CEE", "#EE9A00", "#7AC5CD"))
dev.off()


# calculate a Whitney Wilcoxon test
# p-value = 0.625 (original/random)
# p-value = 0.4375 (random/original)
wilcox.test(FC$fraction[FC$overlap == "BMPC-BMBmem"], 
            FC$fraction[FC$overlap == "BMPC-BLBmem"], 
            paired = TRUE, exact = TRUE, correct = TRUE, conf.int = TRUE) 

# without donor 1
# p-value = 0.25 (original/random)
# p-value = 0.4375 (random/original)
wilcox.test(FC$fraction[FC$overlap == "BMPC-BMBmem" & FC$donor != 1], 
            FC$fraction[FC$overlap == "BMPC-BLBmem" & FC$donor != 1], 
            paired = TRUE, exact = TRUE, correct = TRUE, conf.int = TRUE) 

################################################################################

# Precursor analysis

################################################################################

# merging the different donors in one data frame
mutationsDonor <- rbind(read_delim("8_BM_PCs_2.txt"), read_delim("8_BM_Bsm_2.txt"), read_delim("8_Blood_Bsm_2.txt"),
                        read_delim("8_BM_PCs_3.txt"), read_delim("8_BM_Bsm_3.txt"), read_delim("8_Blood_Bsm_3.txt"),
                        read_delim("8_BM_PCs_8.txt"), read_delim("8_BM_Bmem_8.txt"), read_delim("8_Blood_Bmem_8.txt"), 
                        read_delim("8_BM_PCs_9.txt"), read_delim("8_BM_Bmem_9.txt"), read_delim("8_Blood_Bmem_9.txt"), 
                        read_delim("8_BM_PCs_10.txt"), read_delim("8_BM_Bmem_10.txt"), read_delim("8_Blood_Bmem_10.txt"))

# merge the correct contigs for one barcode together based on prior contig selection (UMI and read counts)
mutationsDonor$barcode <- substr(mutationsDonor$`Sequence ID`, 1, nchar(mutationsDonor$`Sequence ID`) - 10)
mutationsDonor$sampleID <- c(rep("1", dim(read_delim("8_BM_PCs_2.txt"))[1]), rep("2", dim(read_delim("8_BM_Bsm_2.txt"))[1]), rep("3", dim(read_delim("8_Blood_Bsm_2.txt"))[1]),
                             rep("4", dim(read_delim("8_BM_PCs_3.txt"))[1]), rep("5", dim(read_delim("8_BM_Bsm_3.txt"))[1]), rep("6", dim(read_delim("8_Blood_Bsm_3.txt"))[1]),
                             rep("13", dim(read_delim("8_BM_PCs_8.txt"))[1]), rep("14", dim(read_delim("8_BM_Bmem_8.txt"))[1]), rep("15", dim(read_delim("8_Blood_Bmem_8.txt"))[1]),
                             rep("16", dim(read_delim("8_BM_PCs_9.txt"))[1]), rep("17", dim(read_delim("8_BM_Bmem_9.txt"))[1]), rep("18", dim(read_delim("8_Blood_Bmem_9.txt"))[1]),
                             rep("19", dim(read_delim("8_BM_PCs_10.txt"))[1]), rep("20", dim(read_delim("8_BM_Bmem_10.txt"))[1]), rep("21", dim(read_delim("8_Blood_Bmem_10.txt"))[1]))
mutationsDonor$barcode <- paste0(mutationsDonor$barcode, mutationsDonor$sampleID, sep = "")
mutationsDonor$contig <- paste(mutationsDonor$barcode, substr(mutationsDonor$`Sequence ID`, nchar(mutationsDonor$`Sequence ID`) - 8, nchar(mutationsDonor$`Sequence ID`)), sep = "")

mutationsDonor <- filter(mutationsDonor, `V-DOMAIN Functionality` == "productive")

# chose the mutation numbers with possible mutations in D or J segment for FR3 and CDR3
# numbers < 100 --> deleting the () symbols at the end and the beginning of non (x) numbers -> "" -> as.numeric("") = NA therefore it works
mutationsDonor$`CDR3-IMGT Nb of mutations` <- ifelse(is.na(as.numeric(str_sub(sub('.* ()', '', mutationsDonor$`CDR3-IMGT Nb of mutations`), 2, -2))), 
                                                     as.numeric(mutationsDonor$`CDR3-IMGT Nb of mutations`),
                                                     as.numeric(str_sub(sub('.* ()', '', mutationsDonor$`CDR3-IMGT Nb of mutations`), 2, -2)))
mutationsDonor$`CDR3-IMGT Nb of nucleotides` <- ifelse(is.na(as.numeric(str_sub(sub('.* ()', '', mutationsDonor$`CDR3-IMGT Nb of nucleotides`), 2, -2))), 
                                                       as.numeric(mutationsDonor$`CDR3-IMGT Nb of nucleotides`),
                                                       as.numeric(str_sub(sub('.* ()', '', mutationsDonor$`CDR3-IMGT Nb of nucleotides`), 2, -2)))

mutationsDonor$`FR3-IMGT Nb of mutations` <- ifelse(is.na(as.numeric(str_sub(sub('.* ()', '', mutationsDonor$`FR3-IMGT Nb of mutations`), 2, -2))),
                                                    as.numeric(mutationsDonor$`FR3-IMGT Nb of mutations`),
                                                    as.numeric(str_sub(sub('.* ()', '', mutationsDonor$`FR3-IMGT Nb of mutations`), 2, -2)))

# numbers > 100 --> deleting the () symbols at the end and the beginning of non (x) numbers -> "1" -> as.numeric("1") != NA therefore it does not work
mutationsDonor$`FR3-IMGT Nb of nucleotides` <- ifelse(as.numeric(str_sub(sub('.* ()', '', mutationsDonor$`FR3-IMGT Nb of nucleotides`), 2, -2)) < 10, 
                                                      as.numeric(mutationsDonor$`FR3-IMGT Nb of nucleotides`),
                                                      as.numeric(str_sub(sub('.* ()', '', mutationsDonor$`FR3-IMGT Nb of nucleotides`), 2, -2)))


# select the same contigs as before
firstIgL_iso = filter(isotype_all, chain != "IGH" & barcode %in% fullBCR)
firstIgH_iso = filter(isotype_all, chain == "IGH" & barcode %in% fullBCR)
firstIgL = firstIgL_iso$contig_id[match(mutationsDonor$barcode, firstIgL_iso$barcode)]
firstIgH = firstIgH_iso$contig_id[match(mutationsDonor$barcode, firstIgH_iso$barcode)]

mutationsDonor <- filter(mutationsDonor, contig %in% c(firstIgL, firstIgH))

# merge contigs with the same barcode together by adding up the number of mutations
# and the length of the reagion for all FR1-3 and CDR1-3
mutationsDonor <- ddply(mutationsDonor,.(barcode), summarize,
                        CDR1_m = sum(`CDR1-IMGT Nb of mutations`), CDR1_n = sum(`CDR1-IMGT Nb of nucleotides`),
                        CDR2_m = sum(`CDR2-IMGT Nb of mutations`), CDR2_n = sum(`CDR2-IMGT Nb of nucleotides`),
                        CDR3_m = sum(`CDR3-IMGT Nb of mutations`), CDR3_n = sum(`CDR3-IMGT Nb of nucleotides`),
                        FR1_m = sum(`FR1-IMGT Nb of mutations`), FR1_n = sum(`FR1-IMGT Nb of nucleotides`),
                        FR2_m = sum(`FR2-IMGT Nb of mutations`), FR2_n = sum(`FR2-IMGT Nb of nucleotides`),
                        FR3_m = sum(`FR3-IMGT Nb of mutations`), FR3_n = sum(`FR3-IMGT Nb of nucleotides`))

# determine the donor number for each cell
mutationsDonor$donor <- ifelse(substr(mutationsDonor$barcode, 18, 20) %in% c("1", "2", "3"), 2, 
                               ifelse(substr(mutationsDonor$barcode, 18, 20) %in% c("4", "5", "6"), 3,
                                      ifelse(substr(mutationsDonor$barcode, 18, 20) %in% c("13", "14", "15"), 8,
                                             ifelse(substr(mutationsDonor$barcode, 18, 20) %in% c("16", "17", "18"), 9, 10))))

# calculate the mutation rate in each region
mutationsDonor$CDR1 <- mutationsDonor$CDR1_m / mutationsDonor$CDR1_n
mutationsDonor$CDR2 <- mutationsDonor$CDR2_m / mutationsDonor$CDR2_n
mutationsDonor$CDR3 <- mutationsDonor$CDR3_m / mutationsDonor$CDR3_n
mutationsDonor$FR1 <- mutationsDonor$FR1_m / mutationsDonor$FR1_n
mutationsDonor$FR2 <- mutationsDonor$FR2_m / mutationsDonor$FR2_n
mutationsDonor$FR3 <- mutationsDonor$FR3_m / mutationsDonor$FR3_n

# align the isotypes to it's corresponding position in the genome
isotype_ordered = c("IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2") 
mutationsDonor$isotype <- isotype_filtered$c_gene[match(mutationsDonor$barcode, isotype_filtered$barcode)]
mutationsDonor$isotype_ <- match(mutationsDonor$isotype, isotype_ordered)
mutationsDonor$isotype_[is.na(mutationsDonor$isotype_) | mutationsDonor$isotype == ""] = 0
mutationsDonor = filter(mutationsDonor, isotype_ != 0)

# add the cell categorization
mutationsDonor$celltype <- forOverlap$ForOverlaps[match(mutationsDonor$barcode, forOverlap$Barcode)]

# add the clonotype
mutationsDonor$clonotype <- forOverlap$Clonotype[match(mutationsDonor$barcode, forOverlap$Barcode)] 
mutationsDonor <- filter(mutationsDonor, !is.na(clonotype))

# selecting clonotypes that are in either one of the considered overlaps (1 or 2)
clonotype_BMPC <- unique(mutationsDonor$clonotype[mutationsDonor$celltype == "BM_PC"])
clonotype_BMBmem <- unique(mutationsDonor$clonotype[mutationsDonor$celltype == "BM_Bmem"])
clonotype_BLBmem <- unique(mutationsDonor$clonotype[mutationsDonor$celltype == "Blood_Bmem"])

Bmem <- intersect(clonotype_BMBmem, clonotype_BLBmem)
BM <- intersect(clonotype_BMPC, clonotype_BMBmem)[!(intersect(clonotype_BMPC, clonotype_BMBmem) %in% Bmem)]
BMPCBLBmem <- intersect(clonotype_BMPC, clonotype_BLBmem)[!(intersect(clonotype_BMPC, clonotype_BLBmem) %in% Bmem)]

mutationsDonor_filtered <- filter(mutationsDonor, clonotype %in% c(BM, BMPCBLBmem) & !(clonotype %in% Bmem))

sameMR <- function(mD_fil_Bmem, BMPC, i) {
  return(mD_fil_Bmem$FR1 == BMPC$FR1[i] & mD_fil_Bmem$FR2 == BMPC$FR2[i] & mD_fil_Bmem$FR3 == BMPC$FR3[i] & 
         mD_fil_Bmem$CDR1 == BMPC$CDR1[i] & mD_fil_Bmem$CDR2 == BMPC$CDR2[i] & mD_fil_Bmem$CDR3 == BMPC$CDR3[i])
}
smallerMR <- function(mD_fil_Bmem, BMPC, i) {
  return(mD_fil_Bmem$FR1 <= BMPC$FR1[i] & mD_fil_Bmem$FR2 <= BMPC$FR2[i] & mD_fil_Bmem$FR3 <= BMPC$FR3[i] & 
         mD_fil_Bmem$CDR1 <= BMPC$CDR1[i] & mD_fil_Bmem$CDR2 <= BMPC$CDR2[i] & mD_fil_Bmem$CDR3 <= BMPC$CDR3[i])
}

# count the number of possible extrafollicular and interfollicular precursors
# in the clones approach
for (donorIndex in unique(mutationsDonor_filtered$donor)) {
  mD_fil_BMPC = filter(mutationsDonor_filtered, celltype == "BM_PC" & donor == donorIndex & isotype_ != 0)
  mD_fil_BMBmem = filter(mutationsDonor_filtered, celltype == "BM_Bmem" & donor == donorIndex & clonotype %in% BM & isotype_ != 0)
  mD_fil_BLBmem = filter(mutationsDonor_filtered, celltype == "Blood_Bmem" & donor == donorIndex & clonotype %in% BMPCBLBmem & isotype_ != 0)
  
  straight_BM = 0
  straight_BMBL = 0
  down_BM = 0
  down_BMBL = 0
  
  for (i in 1:length(mD_fil_BMPC$barcode)) {
    # BM PC - Bm Bmem (overlap 1)
    straight_BM = straight_BM + length(mD_fil_BMBmem$barcode[mD_fil_BMBmem$clonotype == mD_fil_BMPC$clonotype[i] & 
                                                             sameMR(mD_fil_BMBmem, mD_fil_BMPC, i) == 1 &
                                                             mD_fil_BMBmem$isotype_ == mD_fil_BMPC$isotype_[i]])
    
    down_BM = down_BM + length(mD_fil_BMBmem$barcode[mD_fil_BMBmem$clonotype == mD_fil_BMPC$clonotype[i] & 
                                                     smallerMR(mD_fil_BMBmem, mD_fil_BMPC, i) &
                                                     mD_fil_BMBmem$isotype_ <= mD_fil_BMPC$isotype_[i]])
    
    # BM PC - BL Bmem (overlap 2)
    straight_BMBL = straight_BMBL + length(mD_fil_BLBmem$barcode[mD_fil_BLBmem$clonotype == mD_fil_BMPC$clonotype[i] & 
                                                                 sameMR(mD_fil_BLBmem, mD_fil_BMPC, i) == 1 &
                                                                 mD_fil_BLBmem$isotype_ == mD_fil_BMPC$isotype_[i]])
    
    down_BMBL = down_BMBL + length(mD_fil_BLBmem$barcode[mD_fil_BLBmem$clonotype == mD_fil_BMPC$clonotype[i] & 
                                                         smallerMR(mD_fil_BLBmem, mD_fil_BMPC, i) &
                                                         mD_fil_BLBmem$isotype_ <= mD_fil_BMPC$isotype_[i]])
  }
  
  df_ = as.data.frame(matrix(ncol = 0, nrow = 4))
  df_$fraction <- c(c(straight_BM, down_BM - straight_BM)/ length(mD_fil_BMPC$barcode[mD_fil_BMPC$clonotype %in% BM]),
                    c(straight_BMBL, down_BMBL - straight_BMBL)/ length(mD_fil_BMPC$barcode[mD_fil_BMPC$clonotype %in% BMPCBLBmem]))
  df_$overlap <- c("BMPC-BMBmem", "BMPC-BMBmem", "BMPC-BLBmem", "BMPC-BLBmem")
  df_$mR <- c("Bmem = PC", "Bmem < PC", "Bmem = PC", "Bmem < PC")
  
  png(file = paste0("mR_fraction_woBmem_D", donorIndex, ".png", sep = ""), width = 600, height = 350)
  print(ggplot(df_, aes(y = fraction, x = overlap, fill = mR)) + 
    geom_col(position = "dodge2") + 
    scale_fill_brewer(palette = "Pastel1") +
    ylim(c(0, 0.34)))
  dev.off()
  
  write.table(df_, file = paste0("mutationRates_original_woBmem_D", donorIndex, ".txt", sep = ""), quote = FALSE, sep ="\t", col.names = TRUE)
}



# count the number of the randomized extrafollicular and interfollicular precursors
# in the clones approach
for (donorIndex in unique(mutationsDonor_filtered$donor)) {
  mutationsDonor_randomised = as.data.frame(matrix(nrow = 0, ncol = 4))
  set.seed(42)
  for (counter in 1:1000) {
    mutationsDonor_ftemp = mutationsDonor_filtered
    mutationsDonor_ftemp$celltype = sample(mutationsDonor_filtered$celltype)
    
    mD_fil_BMPC = filter(mutationsDonor_ftemp, celltype == "BM_PC" & donor == donorIndex & isotype_ != 0)
    mD_fil_BMBmem = filter(mutationsDonor_ftemp, celltype == "BM_Bmem" & donor == donorIndex & clonotype %in% BM & isotype_ != 0)
    mD_fil_BLBmem = filter(mutationsDonor_ftemp, celltype == "Blood_Bmem" & donor == donorIndex & clonotype %in% BMPCBLBmem & isotype_ != 0)
    
    straight_BM_r = 0
    straight_BMBL_r = 0
    down_BM_r = 0
    down_BMBL_r = 0
    
    for (i in 1:length(mD_fil_BMPC$barcode)) {
      # BM PC - BM Bmem (overlap 1)
      straight_BM_r = straight_BM_r + length(mD_fil_BMBmem$barcode[mD_fil_BMBmem$clonotype == mD_fil_BMPC$clonotype[i] &
                                                                   sameMR(mD_fil_BMBmem, mD_fil_BMPC, i) == 1 & 
                                                                   mD_fil_BMBmem$isotype_ == mD_fil_BMPC$isotype_[i]])
      
      down_BM_r = down_BM_r + length(mD_fil_BMBmem$barcode[mD_fil_BMBmem$clonotype == mD_fil_BMPC$clonotype[i] &
                                                           smallerMR(mD_fil_BMBmem, mD_fil_BMPC, i) == 1 & 
                                                           mD_fil_BMBmem$isotype_ <= mD_fil_BMPC$isotype_[i]])

      # BM PC - BL Bmem (overlap 2)
      straight_BMBL_r = straight_BMBL_r + length(mD_fil_BLBmem$barcode[mD_fil_BLBmem$clonotype == mD_fil_BMPC$clonotype[i] &
                                                                       sameMR(mD_fil_BLBmem, mD_fil_BMPC, i) == 1 &
                                                                       mD_fil_BLBmem$isotype_ == mD_fil_BMPC$isotype_[i]])
      
      down_BMBL_r = down_BMBL_r + length(mD_fil_BLBmem$barcode[mD_fil_BLBmem$clonotype == mD_fil_BMPC$clonotype[i] &
                                                               smallerMR(mD_fil_BLBmem, mD_fil_BMPC, i) == 1 &
                                                               mD_fil_BLBmem$isotype_ <= mD_fil_BMPC$isotype_[i]])
    }
    
    mutationsDonor_randomised <- rbind(mutationsDonor_randomised, 
                                       c(straight_BM_r/ length(mD_fil_BMPC$barcode[mD_fil_BMPC$clonotype %in% BM]),
                                         (down_BM_r - straight_BM_r)/ length(mD_fil_BMPC$barcode[mD_fil_BMPC$clonotype %in% BM]),
                                         straight_BMBL_r/ length(mD_fil_BMPC$barcode[mD_fil_BMPC$clonotype %in% BMPCBLBmem]),
                                         (down_BMBL_r - straight_BMBL_r)/ length(mD_fil_BMPC$barcode[mD_fil_BMPC$clonotype %in% BMPCBLBmem])))
  }
  colnames(mutationsDonor_randomised) = c("straight_BM", "down_BM", "straight_BMBL", "down_BMBL")
  write.table(mutationsDonor_randomised, file = paste0("mutationRates_randomised_woBmem_D", donorIndex, ".txt", sep = ""), quote = FALSE, sep ="\t", col.names = TRUE)
  
  # bar plot with percentage y-axis scale
  mR_df = as.data.frame(matrix(nrow = 4, ncol = 0))
  mR_df$fraction = c(mean(mutationsDonor_randomised$straight_BM), mean(mutationsDonor_randomised$down_BM), 
                     mean(mutationsDonor_randomised$straight_BMBL), mean(mutationsDonor_randomised$down_BMBL))
  mR_df$overlap = c("BMPC-BMBmem", "BMPC-BMBmem", "BMPC-BLBmem", "BMPC-BLBmem")
  mR_df$mR = c("Bmem = PC", "Bmem < PC", "Bmem = PC", "Bmem < PC")
  
  png(file = paste0("mR_fraction_random_woBmem_D", donorIndex, ".png", sep = ""), width = 600, height = 350)
  print(ggplot(mR_df, aes(y = fraction, x = overlap, fill = mR)) + 
          geom_col(position = "dodge2") + 
          scale_fill_brewer(palette = "Pastel1") +
          ylim(c(0, 0.34)))
  dev.off()
}

# generate a data frame containing all informations necessary for the bee plot
i = 0
precursor_df = as.data.frame(matrix(nrow = 0, ncol = 5))
for (donorIndex in unique(mutationsDonor_filtered$donor)) {
  original <- read.delim(paste0("mutationRates_original_woBmem_D", donorIndex, ".txt", sep=""))
  randomised <- read.delim(paste0("mutationRates_randomised_woBmem_D", donorIndex, ".txt", sep=""))
  
  df_temp = as.data.frame(matrix(nrow = 8, ncol = 0))
  df_temp$fraction = c(original$fraction, mean(randomised$straight_BM), 
                       mean(randomised$down_BM), mean(randomised$straight_BMBL),
                       mean(randomised$down_BMBL))
  df_temp$donor = c(rep(donorIndex, 8))
  df_temp$randomized = c(rep("original", 4), rep("randomized", 4))
  df_temp$order = c(i, i+1, i+2, i+3, i, i+1, i+2, i+3)
  df_temp$split = c("BMPC == BMBmem", "BMPC < BMBmem", "BMPC == BLBmem", "BMPC < BLBmem",
                    "r(BMPC == BMBmem)", "r(BMPC < BMBmem)", "r(BMPC == BLBmem)", "r(BMPC < BLBmem)")
  df_temp = df_temp[order(df_temp$order), ]
  precursor_df = rbind(precursor_df, df_temp)
  i = i + 8
}

precursor_df$facet = ifelse(precursor_df$split == "BMPC == BMBmem" | precursor_df$split == "r(BMPC == BMBmem)", "BMPC == BMBmem", 
                            ifelse(precursor_df$split == "BMPC == BLBmem" | precursor_df$split == "r(BMPC == BLBmem)", "BMPC == BLBmem", 
                                   ifelse(precursor_df$split == "BMPC < BMBmem" | precursor_df$split == "r(BMPC < BMBmem)", "BMPC < BMBmem", "BMPC < BLBmem")))


# Bee plot
png(file = paste0("precursor_beeplot.png", sep = ""), width = 2300, height = 2050, res = 300)
precursor_df %>%
  mutate(split = fct_relevel(split, 
                             "BMPC == BMBmem", "r(BMPC == BMBmem)",
                             "BMPC < BMBmem", "r(BMPC < BMBmem)",
                             "BMPC == BLBmem", "r(BMPC == BLBmem)",
                             "BMPC < BLBmem", "r(BMPC < BLBmem)")) %>%
  ggplot(aes(as.factor(randomized), fraction, color = as.factor(donor))) + 
  geom_beeswarm(cex = 4, size = 3) +
  scale_color_manual(values=c("olivedrab3", "black", "#912CEE", "#EE9A00", "#7AC5CD")) +
  facet_wrap("facet",scales="free")
dev.off()


# calculate all 4 Whitney Wilcoxon tests - randomized vs true observed ratios

# extrafollicular precursor in the BM
# p value = 0.4375 --> not significant
wilcox.test(precursor_df$fraction[precursor_df$split == "BMPC == BMBmem"], 
            precursor_df$fraction[precursor_df$split == "r(BMPC == BMBmem)"], 
            paired = TRUE, exact = TRUE, correct = TRUE, conf.int = TRUE)

# interfollicular precursor in the BM
# p value = 0.0625 --> not significant
wilcox.test(precursor_df$fraction[precursor_df$split == "BMPC < BMBmem"], 
            precursor_df$fraction[precursor_df$split == "r(BMPC < BMBmem)"], 
            paired = TRUE, exact = TRUE, correct = TRUE, conf.int = TRUE)

# different GC emigration time points or BL Bmem reentry in the GC
# p value = 0.0625 --> not significant
wilcox.test(precursor_df$fraction[precursor_df$split == "BMPC == BLBmem"], 
            precursor_df$fraction[precursor_df$split == "r(BMPC == BLBmem)"], 
            paired = TRUE, exact = TRUE, correct = TRUE, conf.int = TRUE)

# interfollicular precursor in the blood
# p value = 0.0625 --> not significant
wilcox.test(precursor_df$fraction[precursor_df$split == "BMPC < BLBmem"], 
            precursor_df$fraction[precursor_df$split == "r(BMPC < BLBmem)"], 
            paired = TRUE, exact = TRUE, correct = TRUE, conf.int = TRUE)




# calculate all 2 Whitney Wilcoxon tests - overlap 1 vs overlap 2

# extrafollicular differentiation in overlap 1 versus different GC emigration 
# time points in overlap 2
# p value = 0.825 --> not significant
wilcox.test(precursor_df$fraction[precursor_df$split == "BMPC == BMBmem"], 
            precursor_df$fraction[precursor_df$split == "BMPC == BLBmem"], 
            paired = TRUE, exact = TRUE, correct = TRUE, conf.int = TRUE)

# interfollicular differentiation between overlap 1 and 2
# p value = 1 --> not significant
wilcox.test(precursor_df$fraction[precursor_df$split == "BMPC < BMBmem"], 
            precursor_df$fraction[precursor_df$split == "BMPC < BLBmem"], 
            paired = TRUE, exact = TRUE, correct = TRUE, conf.int = TRUE)












