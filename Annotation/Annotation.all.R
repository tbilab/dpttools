# IRanges annotation in genereal

library(IRanges)
source("Annotation.fun.R")

# load the results
test.report.sigs <- read.csv(paste("../", list.files("../")[grep("sites.csv",list.files("../"))], sep = ""), as.is = T, stringsAsFactors = FALSE, colClasses=c("character","numeric","numeric","numeric","character","character","numeric","numeric","numeric","character","character"))

test.report.sigs$id <- paste(test.report.sigs$chr, test.report.sigs$start, test.report.sigs$end, sep = "_")

# make it Rangeddata
test.report.sigs$space <- test.report.sigs$chr
test <- as(test.report.sigs,"RangedData") 

test0 <- ranges(test)

#########
# 1. CpG 
load("../Datasets/hg18_CpGislands.RData")
CpG <- remove_(CpG)

map2result(test0, CpG, "CpGi", cal.rate = FALSE) # 9182

# 12979 on 7/18/12
# 11402 on 8/10/12

# test
load(paste("result", "CpGi", "map.RData", sep = "_"))
dim(temp)
temp[1,]

#########
# 2 promoter

load("../Datasets/refFlat_tx.RData")
prom0 <- remove_(FLAT.tx)

# create promoter regions by defn of TSS - 1000 to TSS + 500 for +, TES - 500 to TES +1000 for -
rd <- as.data.frame(prom0)
rd$start[rd$strand == "+"] <- rd$txStart[rd$strand == "+"] - 1000
rd$end[rd$strand == "+"] <- rd$txStart[rd$strand == "+"] + 500
rd$start[rd$strand == "-"] <- rd$txEnd[rd$strand == "-"] - 500
rd$end[rd$strand == "-"] <- rd$txEnd[rd$strand == "-"] + 1000

# check uniqueness of geneName and Name
dim(rd)
length(unique(rd$geneName)) # not unique
length(unique(rd$name)) # not unique

start(prom0) <- rd$start
end(prom0) <- rd$end

## 
map2result(test0, prom0, "prom", cal.rate = FALSE) #6408
# 11036 on 7/18/12
#  9634 on 8/10/12

load(paste("result", "prom", "map.RData", sep = "_"))
dim(temp)
temp[1,]

#########
# 2.1 transcript ends

load("../Datasets/refFlat_tx.RData")
trans <- remove_(FLAT.tx)

# create promoter regions by defn of TSS - 1000 to TSS + 500 for +, TES - 500 to TES +1000 for -
rd <- as.data.frame(trans)
rd$start[rd$strand == "+"] <- rd$txEnd[rd$strand == "+"] - 500
rd$end[rd$strand == "+"] <- rd$txEnd[rd$strand == "+"] + 1000
rd$start[rd$strand == "-"] <- rd$txStart[rd$strand == "-"] - 1000
rd$end[rd$strand == "-"] <- rd$txStart[rd$strand == "-"] + 500

# check uniqueness of geneName and Name
dim(rd)
length(unique(rd$geneName)) # not unique
length(unique(rd$name)) # not unique

start(trans) <- rd$start
end(trans) <- rd$end

## 
map2result(test0, trans, "txend", cal.rate = FALSE)
# 7721 on 8/9/12
#  6194 on 8/10/12

load(paste("result", "txend", "map.RData", sep = "_"))
dim(temp)
temp[1,]

#########
# 2.2 Exon

load("../Datasets/refFlat_exon.RData")
exon <- remove_(FLAT.exon)

map2result(test0, exon, "exon", cal.rate = FALSE) #180994
#  142226 on 8/10/12

load(paste("result", "exon", "map.RData", sep = "_"))
dim(temp)
temp[1,]


#########
# 3 CAGE

load("../Datasets/Consolidated_Ricken_CAGE.RData")
CAGE.all <- CAGE[["rd"]]
CAGE.all <- remove_(CAGE.all)
CAGE.all$value <- "CAGE"
map2result(test0, CAGE.all, "CAGE", cal.rate = FALSE)
load(paste("result", "CAGE", "map.RData", sep = "_"))
dim(temp) #  142226 on 8/10/12

temp[1,]

#########
# 4 miRNA
load("../Datasets/hg18_miRNA.RData")
RNA <- remove_(RNA)

map2result(test0, RNA, "miRNA", cal.rate = FALSE) #56
load(paste("result", "miRNA", "map.RData", sep = "_"))
dim(temp)
temp[1,]

#########
# 5 CTCF
load("../Datasets/hg18_CTCF.RData")
CTCF.all <- CTCF[["rd"]]
CTCF.all <- remove_(CTCF.all)


map2result(test0, CTCF.all, "CTCF", cal.rate = FALSE) #5831
load(paste("result", "CTCF", "map.RData", sep = "_"))
dim(temp)
temp[1,]

#########
# 6 DNaseI
load("../Datasets/DNaseI.RData")
DNaseI <- remove_(DNaseI)


map2result(test0, DNaseI, "DNaseI", cal.rate = FALSE)# 22159
load(paste("result", "DNaseI", "map.RData", sep = "_"))
dim(temp)
temp[1,]


#########
# 7 RepeatFamily
load("../Datasets/hg18_Repeatfamily.RData")
REP <- remove_(REP)

map2result(test0, REP, "REP", cal.rate = FALSE)# 30329
load(paste("result", "REP", "map.RData", sep = "_"))
dim(temp)
temp[1,]

#
########
########
# summarize the annotation results
########

# load in the anno results first
# CpGi

resultpath <- "." 
load(paste(resultpath, paste("result", "CpGi", "map.RData", sep = "_"), sep = "/"))
length(unique(temp$id)) 
result_allChr_CpGi <- temp[, c("id", "chr", "start", "end", "CpGi.value")]
result_allChr_CpGi <- result_allChr_CpGi[!duplicated(result_allChr_CpGi),]

#promoter genes
load(paste(resultpath, paste("result", "prom", "map.RData", sep = "_"), sep = "/")) 
length(unique(temp$id)) 
result_allChr_prom <- temp[, c("id", "chr", "start", "end", "prom.value")]
result_allChr_prom <- result_allChr_prom[!duplicated(result_allChr_prom),] 
# collapse the gene names for duplicate record
out <- result_allChr_prom[!duplicated(result_allChr_prom$id),]
out <- out[order(out$id),]
Test <- tapply(result_allChr_prom$prom.value, result_allChr_prom$id, paste, collapse = ";")
out$prom.value<- Test

result_allChr_prom <- out 

#txend genes
load(paste(resultpath, paste("result", "txend", "map.RData", sep = "_"), sep = "/")) 
length(unique(temp$id)) 
result_allChr_txend <- temp[, c("id", "chr", "start", "end", "txend.value")]
result_allChr_txend <- result_allChr_txend[!duplicated(result_allChr_txend),] 
# collapse the gene names for duplicate record
out <- result_allChr_txend[!duplicated(result_allChr_txend$id),]
out <- out[order(out$id),]
Test <- tapply(result_allChr_txend$txend.value, result_allChr_txend$id, paste, collapse = ";")
out$txend.value<- Test

result_allChr_txend <- out

#exon
load(paste(resultpath, paste("result", "exon", "map.RData", sep = "_"), sep = "/")) 
length(unique(temp$id)) 
result_allChr_exon <- temp[, c("id", "chr", "start", "end", "exon.value")]
result_allChr_exon <- result_allChr_exon[!duplicated(result_allChr_exon),] 
# collapse the gene names for duplicate record
out <- result_allChr_exon[!duplicated(result_allChr_exon$id),]
out <- out[order(out$id),]
Test <- tapply(result_allChr_exon$exon.value, result_allChr_exon$id, paste, collapse = ";")
out$exon.value<- Test

result_allChr_exon <- out

# CAGE
load(paste(resultpath, paste("result", "CAGE", "map.RData", sep = "_"), sep = "/")) 
temp$CAGE.value <- "CAGE"
length(unique(temp$id)) 
result_allChr_CAGE <- temp[, c("id", "chr", "start", "end", "CAGE.value")]
result_allChr_CAGE <- result_allChr_CAGE[!duplicated(result_allChr_CAGE),]

#miRNA
load(paste(resultpath, paste("result", "miRNA", "map.RData", sep = "_"), sep = "/"))
length(unique(temp$id)) 
result_allChr_miRNA <- temp[, c("id", "chr", "start", "end", "miRNA.value")]
result_allChr_miRNA <- result_allChr_miRNA[!duplicated(result_allChr_miRNA),] 
# collapse the gene names for duplicate record
out <- result_allChr_miRNA[!duplicated(result_allChr_miRNA$id),]
out <- out[order(out$id),]
Test <- tapply(result_allChr_miRNA$miRNA.value, result_allChr_miRNA$id, paste, collapse = ";")
out$miRNA.value<- Test

result_allChr_miRNA <- out 

# CTCF
load(paste(resultpath, paste("result", "CTCF", "map.RData", sep = "_"), sep = "/")) 
length(unique(temp$id)) 
result_allChr_CTCF <- temp[, c("id", "chr", "start", "end", "CTCF.value")]
result_allChr_CTCF <- result_allChr_CTCF[!duplicated(result_allChr_CTCF),]

#DNaseI
load(paste(resultpath, paste("result", "DNaseI", "map.RData", sep = "_"), sep = "/")) 
length(unique(temp$id)) 
result_allChr_DNaseI <- temp[, c("id", "chr", "start", "end", "DNaseI.value")]
result_allChr_DNaseI <- result_allChr_DNaseI[!duplicated(result_allChr_DNaseI),]

# REP
load(paste(resultpath, paste("result", "REP", "map.RData", sep = "_"), sep = "/")) 
length(unique(temp$id)) 
result_allChr_REP <- temp[, c("id", "chr", "start", "end", "REP.value")]
result_allChr_REP <- result_allChr_REP[!duplicated(result_allChr_REP),]


#summ$CpG <- ifelse(is.na(summ$CpGi.start), 0, 1)
#summ$Prom <- ifelse(is.na(summ$prom.start), 0, 1)
#summ$CAGE<- ifelse(is.na(summ$CAGE.start), 0, 1)
#summ$miRNA<- ifelse(is.na(summ$miRNA.start), 0, 1)
#summ$CTCF<- ifelse(is.na(summ$CTCF.start), 0, 1)
#summ$DNaseI<- ifelse(is.na(summ$DNaseI.start), 0, 1)
#summ$REP<- ifelse(is.na(summ$REP.start), 0, 1)

summ <- merge(test.report.sigs, result_allChr_CpGi[, c("id", "CpGi.value")], by.x = "id", by.y = "id", all.x = TRUE)
summ <- merge(summ, result_allChr_prom[, c("id", "prom.value")], by.x = "id", by.y = "id", all.x = TRUE)
summ <- merge(summ, result_allChr_txend[, c("id", "txend.value")], by.x = "id", by.y = "id", all.x = TRUE)
summ <- merge(summ, result_allChr_exon[, c("id", "exon.value")], by.x = "id", by.y = "id", all.x = TRUE)
summ <- merge(summ, result_allChr_CAGE[, c("id", "CAGE.value")], by.x = "id", by.y = "id", all.x = TRUE)
summ <- merge(summ, result_allChr_miRNA[, c("id", "miRNA.value")], by.x = "id", by.y = "id", all.x = TRUE)
summ <- merge(summ, result_allChr_CTCF[, c("id", "CTCF.value")], by.x = "id", by.y = "id", all.x = TRUE)
summ <- merge(summ, result_allChr_DNaseI[, c("id", "DNaseI.value")], by.x = "id", by.y = "id", all.x = TRUE)
summ <- merge(summ, result_allChr_REP[, c("id", "REP.value")], by.x = "id", by.y = "id", all.x = TRUE)


#remove duplicates
summ1 <- subset(summ, !duplicated(summ$id)) 
summ1 <- summ1[, -which(colnames(summ1) == "space")]
#output the csv file
#
out <- summ1[,c("site.id", "chr", "start", "end", "width", "pattern", "contrast", "ratio", "p.value", "qvalue", "region.id", "CpGi.value", "prom.value", "txend.value", "exon.value", "CAGE.value", "miRNA.value", "CTCF.value", "DNaseI.value", "REP.value")]
write.csv(out, paste("MyTable_S1_", substr(Sys.time(), 1, 10), ".csv", sep = ""), na = "", row.names = FALSE) 
