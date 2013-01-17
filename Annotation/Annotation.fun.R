# IRanges annotation in genereal

library(IRanges)
#/home/naj/MyProjects_SSGB/shared/DPTAnalyses/Ting_ProstateClinical/COMMON/iDPT/analyses/BW/Ann

# functions

# remove chr with _ X and Y
# note: change to _ and M later based on the result file

remove_ <- function(rd){
  if(length(grep("[_M]", space(rd))) != 0){
    rd[-grep("[_M]", space(rd)),]
  }
  else {
    rd
  }
#	rd[-grep("[_XY]", space(rd)),]
}

# overlap with the IRanges annotation
getOverlap <- function(chrnm = "chr10", query = test0, subject = rd, mymaxgap = 50, cal.rate = TRUE, rev = FALSE){
  if(chrnm %in% space(subject)){
  anno <- subject[chrnm]
  
  # get the space chrnm of result and annotation file
  chri.anno <- ranges(anno)@unlistData#chri.anno <- subject[[chrnm]]
  chri.result <- query[[chrnm]]
  anno <- anno@values@unlistData@listData
    
  if(!is.null(chri.anno) & !is.null(chri.result)){ # If both subsets are not empty
    # find overlap of the query with the subjects
    if (rev){
		chri.anno <- gaps(chri.anno)
	}
    ol <- findOverlaps(chri.result, chri.anno, maxgap = mymaxgap)
    ind <- as.matrix(ol)

	#par(mfrow = c(2, 1))
	#plotRanges(chri.anno[ind[,2]])
	#plotRanges(chri.result[ind[,1]])

    # get info of the overlap if there is any overlap
    if(dim(ind)[1]>0){ # If there is any overlap
      # inititate result data frame
      results <- data.frame(query.start =  rep(0, dim(ind)[1]), query.end = 0, subject.start=0, subject.end = 0, query.rate = 0, subject.rate = 0, value = as.character(rep("", dim(ind)[1])))
      results$chr <- chrnm
      results$value <- ""
      for( i in 1: dim(ind)[1]){ # for each overlapping record
	  # Save the start and end of the overlapping segments 
        results$query.start[i] <- start(chri.result[ind[i,1]])
        results$query.end[i] <- end(chri.result[ind[i,1]])
        results$subject.start[i] <- start(chri.anno[ind[i,2]])
        results$subject.end[i] <- end(chri.anno[ind[i,2]])
        results$value[i] <- as.character(anno$value[ind[i,2]])

        # if there is any gap and cal.rate = TRUE
        if (cal.rate){
	  #Get the gap of the overlapping segments and initiate the signs
          gap <- intersect(chri.anno[ind[i,2]], chri.result[ind[i,1]])
          anno.sign <- 1 # case 3.1 4.1
          result.sign <- 1
          
          if(length(gap) != 0) { # if the overlap is positive, calculate deltas, case1-4
            delta1 <- abs(start(gap) - start(chri.anno[ind[i,2]]))
            delta2 <- abs(end(gap) - end(chri.anno[ind[i,2]]))
            delta3 <- abs(start(gap) - start(chri.result[ind[i,1]]))
            delta4 <- abs(end(gap) - end(chri.result[ind[i,1]]))
            
	    if (delta1 > delta2 & !(delta3 ==0 & delta4 ==0)) anno.sign <- -1 # case 1
            if (delta1 < delta2 & !(delta3 ==0 & delta4 ==0)) result.sign <- -1 # case 2
            if(width(gap) == width(chri.anno[ind[i,2]]) | width(gap) == width(chri.result[ind[i,1]])){ # case 3 or 4
              if(delta3 == delta4 & delta3 == 0 & delta1 > delta2) anno.sign <- -1 #case 3.2
              if(delta1 == delta2 & delta1 == 0 & delta3 > delta4) result.sign <- -1 #case 4.2
            }
            
            results$query.rate[i] <- result.sign * width(gap)/width(chri.result[ind[i,1]])* 100
            results$subject.rate[i] <- anno.sign * width(gap)/width(chri.anno[ind[i,2]]) * 100
          }else{# if the overlap is really a gap, case 5 6
            gap <- pgap(chri.anno[ind[i,2]], chri.result[ind[i,1]])
	    delta1 <- abs(start(gap) - start(chri.anno[ind[i,2]]))
            delta2 <- abs(end(gap) - end(chri.anno[ind[i,2]]))
            delta3 <- abs(start(gap) - start(chri.result[ind[i,1]]))
            delta4 <- abs(end(gap) - end(chri.result[ind[i,1]]))
            
            if(delta2 == delta3) result.sign <- -1 # caes 5
            if(delta1 == delta4) anno.sign <- -1 # case 6
            results$query.rate[i] <- result.sign * (width(gap) + width(chri.result[ind[i,1]]))/width(chri.result[ind[i,1]])
            results$subject.rate[i] <- anno.sign * (width(gap) + width(chri.anno[ind[i,2]]))/width(chri.anno[ind[i,2]])
          }
        }
#ol1 <- c(chri.anno[ind[i,2]], chri.result[ind[i,1]], gap)
#plotRanges(ol1, col = c(1, 2, 3))
#browser()
      }
      results
    }
    else{
      return(NULL)
    }
  }
  else{
    return(NULL)
  }
  }
  else{
	return(NULL)
  }
}

map2result <- function(result.rd, anno, nm, mymaxgap = 50, cal.rate = TRUE, rev = FALSE){
  
  temp <- NULL
  for (chr in  unique(space(result.rd)) ) {
    temp <- rbind(temp, getOverlap(chrnm = chr, query = result.rd, subject = anno, mymaxgap, cal.rate, rev))
  }

  temp <- temp[, c(dim(temp)[2], 1:(dim(temp)[2]-1))]
  colnames(temp) <- c("chr", "start", "end", paste(nm, "start", sep = "."), paste(nm, "end", sep = "."),
                          paste(nm, "O_result.rate", sep = "_"), paste(nm, "O", paste(nm, "rate", sep = "."), sep ="_"),  paste(nm, "value", sep ="."))
  temp$id <- paste(temp$chr, temp$start, temp$end, sep = "_")
  save(temp, file = paste("result", nm, "map.RData", sep = "_"))
}
