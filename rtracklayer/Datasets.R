# datasets.R
# by Jie Na
# 1-15-12
# to get the annotation datasets
##
# cd /proj/genomics/Projects/shared/GMI/Ting_Methylation/.common/DPT_Projects/Prostate_DPT/rtracklayer
source("getTable1.R")

#################################
## 1. CpGisland
CpG <- getTableFromGB(sessionNM = "hg18", TrackNM = "cpgIslandExt", vnm = "CpGi")
save(CpG, file = "Datasets/hg18_CpGislands.RData")
#28226x10

#################################
### 2. CTCF
CTCF <- getTableFromGB(sessionNM = "hg18", TrackNM = "wgEncodeBroadChipSeq", vnm = "CTCF")
save(CTCF, file = "Datasets/hg18_CTCF.RData")
#102340

#################################
### 3. RNA
RNA <- getTableFromGB(sessionNM = "hg18", TrackNM = "wgRna", vnm = "miRNA")
save(RNA, file = "Datasets/hg18_miRNA.RData")
#1120

#################################
### 4. Repeated famlily
REP <- getTableFromGB(sessionNM = "hg18", TrackNM = "rmsk", vnm = "rmsk")
save(REP, file = "Datasets/hg18_Repeatfamily.RData")
#5,057,005

#################################
### 5. expression/ Riken Cage loc
EXP <- getTableFromGB(sessionNM = "hg18", TrackNM = "wgEncodeRikenCage", vnm = "RickenCage")
save(EXP, file = "Datasets/Consolidated_Ricken_CAGE.RData")
#str(EXP)
#'data.frame':   21,482,388  obs. of  5 variables:

#################################
### 6. refFlat
FLAT.tx <- getTableFromGB(sessionNM = "hg18", TrackNM = "refGene", reftype = "tx", vnm = "Gene")
save(FLAT.tx, file = "Datasets/refFlat_tx.RData")

FLAT.cds <- getTableFromGB(sessionNM = "hg18", TrackNM = "refGene", reftype = "cd", vnm = "CodingRegion")
save(FLAT.cds, file = "Datasets/refFlat_cds.RData")

FLAT.exon <- getTableFromGB(sessionNM = "hg18", TrackNM = "refGene", reftype = "exon", vnm = "Exon")
save(FLAT.exon, file = "Datasets/refFlat_exon.RData")

#str(FLAT)
#'data.frame':   40,019  obs. of 11 variables:

#################################
### 7. Mapability
MAP <- getTableFromGB(sessionNM = "hg18", TrackNM = "wgEncodeMapability", vnm = "Map")
save(MAP, file = "Datasets/Mapability.RData")


#################################
### 8. DNaseI
DNaseI <- getTableFromGB(sessionNM = "hg18", TrackNM = "wgEncodeRegDnaseClustered", vnm = "DNaseI")
save(DNaseI, file = "Datasets/DNaseI.RData")
#str(DNaseI)
#'data.frame':   956991  obs. of  5 variables:


