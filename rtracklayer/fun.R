GenomicRanges <- function (start = integer(), end = integer(), chrom = NULL, genome = NULL) 
{
    ir <- IRanges(start, end)
    if (!is.null(chrom)) {
        if (!is.factor(chrom)) 
            chrom <- factor(chrom, unique(chrom))
        if (!length(ir)) 
            chrom <- chrom[FALSE]
        if (length(chrom) != length(ir)) 
            stop("length of 'chrom' must equal that of 'start' and 'end'", 
                " unless 'start' and 'end' are zero length")
        rl <- split(ir, chrom)
    }
    else rl <- RangesList(ir)
    universe(rl) <- genome
    rl
}
