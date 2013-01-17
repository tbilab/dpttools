## getTableFromGB.R

# functions to get the table from UCSC genome brower
# by : Jeanie Na
# 12-16-11
# modified on 1-18-12 for newer version of rtracklayer


# start of the function
getTableFromGB <- function(sessionNM = "hg18", 
                           TrackNM = "cpgIslandExt",
                           vnm = "CpGi",
                           reftype = NULL, starti=1, endi=1){
  
  library(IRanges)
  library(rtracklayer)
  library(gdata)
  
  session <- browserSession()
  genome(session) <- sessionNM 
  
  
  if(TrackNM == "rmsk"){
    seqnames <- c("chr1","chr1_random","chr2","chr2_random","chr3","chr3_random",
                  "chr4","chr4_random","chr5","chr5_h2_hap1","chr5_random","chr6","chr6_cox_hap1","chr6_qbl_hap2", "chr6_random",
                  "chr7","chr7_random","chr8","chr8_random","chr9","chr9_random","chr10","chr10_random","chr11","chr11_random",
                  "chr12","chr13","chr13_random","chr14","chr15","chr15_random","chr16","chr16_random","chr17","chr17_random",
                  "chr18","chr18_random","chr19","chr19_random","chr20","chr21","chr21_random","chr22","chr22_h2_hap1","chr22_random", 
                  "chrX","chrX_random","chrY","chrM" )
    
    query <- ucscTableQuery(session, TrackNM, range = GRangesForUCSCGenome(genome =sessionNM, chrom = seqnames[1]))
    tn <- tableNames(query)
    tableName(query) <- tn 	
    table.df <- getTable(query)
    for ( i in 2:length(seqnames)){
      query <- ucscTableQuery(session, TrackNM, range = GRangesForUCSCGenome(genome =sessionNM, chrom = seqnames[i]))
      tn <- tableNames(query)
      tableName(query) <- tn 	
      table.df <- rbind(table.df, getTable(query))
    }
    # get the IRangeData
    table.df$value <- vnm
    rn <- table.df
    rn$space <- table.df$genoName
    rn$start <- table.df$genoStart
    rn$end <- table.df$genoEnd
    rd <- as(rn, "RangedData")
    
  } else {
    query <- ucscTableQuery(session, TrackNM)
    
    ## list the table names
    tn <- tableNames(query)
    
    if(length(tn) == 1){
      ## get the table
      tableName(query) <- tn 
      ## get a data.frame
      table.df <- getTable(query)
      table.df$value <- vnm
      if(vnm == "miRNA") table.df$value <- table.df$name

      # get the IRangeData			
      rn <- table.df
      rn$space <- table.df$chrom
      rn$start <- table.df$chromStart
      rn$end <- table.df$chromEnd
        
      rd <- as(rn, "RangedData")
    }
    else{
      if ( TrackNM == "wgEncodeBroadChipSeq") {
        ctcf <- grep("Ctcf", tn)
        tn1 <- tn[ctcf][grep("Peaks", tn[ctcf])]
      }
      if ( TrackNM == "wgEncodeRikenCage") {
        tn1 <- tn[-grep("Alig", tn)]
      }
      if ( TrackNM == "refGene") {
        tn1 <- "refFlat"
      }
      if ( TrackNM == "wgEncodeMapability") {
        #mapi <- grep("36", tn)
        #mapi <- c(mapi, grep("50", tn))
	  mapi <- grep("Broad", tn)
        
        tn1 <- tn[mapi]
      }
      if (TrackNM == "refGene") {
      	## get the table
        tableName(query) <- tn1 
        ## get a data.frame
        table.df <- getTable(query)
	
	# get the IRangeData        
        rn <- table.df
        
        rn$space <- rn$chrom
        if(reftype == "tx"){
        rn$start <- rn$txStart
        rn$end <- rn$txEnd
        rn$value <- paste(rn$geneName, rn$name, sep = "/+/")
        }
        if(reftype == "exon"){
          rn <- NULL
          #for(i in 1:dim(table.df)[1]){
          for(i in starti:endi){
              
            print(i)
                                        # exon part
            table.df.temp <- table.df[rep(i, table.df$exonCount[i]),]
            table.df.temp$exonStartsi <- as.integer(unlist(strsplit(as.character(table.df$exonStarts[i]), ","))[1])
            table.df.temp$exonEndsi <- as.integer(unlist(strsplit(as.character(table.df$exonEnds[i]), ","))[1])
            
            for(j in 1: table.df$exonCount[i]){
              table.df.temp$exonStartsi[j] <- as.integer(unlist(strsplit(as.character(table.df$exonStarts[i]), ","))[j])
              table.df.temp$exonEndsi[j] <- as.integer(unlist(strsplit(as.character(table.df$exonEnds[i]), ","))[j])
              table.df.temp$exonCounti[j] <- j
            }
            table.df.temp.exon <- table.df.temp
            table.df.temp.exon$value <- paste( table.df.temp.exon$geneName, table.df.temp.exon$name, paste("exon",  table.df.temp.exon$exonCounti, sep = ""), sep = "/")
            
            if(table.df$exonCount[i] > 1) {
                                        # intron part
              table.df.temp <- table.df[rep(i, table.df$exonCount[i]),]
              x0 <- IRanges(start = table.df.temp.exon$exonStartsi, end = table.df.temp.exon$exonEndsi)
              
              intron <- gaps(x0)
              if(length(intron) >0){
                for(j in 1: length(start(intron))){
                  table.df.temp$exonStartsi[j] <- start(intron)[j]
                  table.df.temp$exonEndsi[j] <- end(intron)[j]
                  table.df.temp$exonCounti[j] <- j
                }
                table.df.temp.intron <- table.df.temp
                table.df.temp.intron$value <-  paste( table.df.temp.exon$geneName, table.df.temp.exon$name, paste("intron",  table.df.temp.exon$exonCounti, sep = ""), sep = "/")
                temp <- interleave(table.df.temp.exon, table.df.temp.intron, drop=FALSE)
                temp <- temp[-dim(temp)[1],]
                
                rn <- rbind(rn,temp)
              }else{
                rn <- rbind(rn, table.df.temp.exon)
              }
            }else{
              rn <- rbind(rn, table.df.temp.exon)
            }
            #print(paste("   ", dim(rn), sep = ""))
          }
          rn$start <- rn$exonStartsi
          rn$end <- rn$exonEndsi
          rn$space <- rn$chrom
          
        }
        if(reftype == "cd"){
        rn$start <- rn$cdsStart
        rn$end <- rn$cdsEnd
        rn$value <- paste(rn$geneName, rn$name, sep = "/+/")
        }
        
        rd <- as(rn, "RangedData")
        
      } else {
        rds <- RangedDataList()
      	for (i in 1:length(tn1)){
          ## get the table
          tableName(query) <- tn1[i] 
          ## get a data.frame
          table.df <- getTable(query)
          table.df$value <- vnm
        
          # get the IRangeData
          rn <- table.df
          rn$space <- table.df$chrom
          rn$start <- table.df$chromStart
          rn$end <- table.df$chromEnd
          rd <- as(rn, "RangedData")
          rds[[tn1[i]]]<- rd
        }
        if( TrackNM != "wgEncodeMapability") {
      # consolidate the data
          cons <- as.data.frame(rds[[1]])
          for( i in 2:length(tn1)){
            cons <- rbind(cons, as.data.frame(rds[[i]]))
          }
          value <- rep(vnm, dim(cons)[1])                      
          rd <- reduce(RangedData(IRanges(start = cons$start, end = cons$end), value , space = cons$space))
        }
      }
    }
    
  }

  if( TrackNM == "wgEncodeBroadChipSeq") {
    return(list(raw.full = rds, rd = rd))
  }
  else if( TrackNM == "wgEncodeMapability") {
    return(raw.full = rds)
  }
  else{
    return(rd)
  }
}


