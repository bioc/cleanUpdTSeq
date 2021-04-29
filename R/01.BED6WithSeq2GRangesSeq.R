#' Covert (extended) BED6 file to a GRanges object
#' 
#' Convert to a GRanges object from a (extended) BED6 file with at least six 
#' columns: chrom, chromStart, strEnd, name, score and strand, and optional 
#' upstream sequences (including pA sites) and downstream sequences of pA sites
#' 
#' @param file A character(1) vector, representing a path to a extended BED file
#' containing at least six columns in the order of chrom, chromStart, strEnd, 
#' name, score and strand. The strand information must be designated as "+",
#' or "-". Optional fields--upstream sequences (including pA sites) and 
#' downstream sequences of pA sites--are allowed. For more details about the BED
#' format, see https://genome.ucsc.edu/FAQ/FAQformat.html#format1.
#' @param skip A integer(1) vector, indicating how many rows (header lines) to 
#' skip when the BED file is read into R.
#' @param withSeq A logical(1) vector, indicating that upstream and downstream 
#' sequences flanking pA sites are included in the file
#' @param upstream.seq.ind An integer(1),vector delineating the column 
#' location of upstream sequences of the putative pA site
#' @param downstream.seq.ind An integer(1),vector delineating the column 
#' location of downstream sequences of the putative pA site
#' @return An object of GRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom utils read.delim
#' @author  Haibo Liu, Lihua J. Zhu
#' @examples
#' 
#' testFile <- system.file("extdata", "test.bed",
#'                         package = "cleanUpdTSeq")
#' peaks <- BED6WithSeq2GRangesSeq(file = testFile, 
#'                                skip = 1L, withSeq = TRUE)
#' 
#' @export
BED6WithSeq2GRangesSeq <- function(file, 
                                   skip = 1L,
                                   withSeq = TRUE,
                                   upstream.seq.ind = 7L,
                                   downstream.seq.ind = 8L) 
{
    if (missing(file) || !file.exists(file)){
        stop("file is required and must be readable!")
    }
    if (!is.numeric(skip)){stop("skip must be an integer")}
    if (withSeq && (!is.numeric(upstream.seq.ind) || 
        !is.numeric(downstream.seq.ind))) {
        stop("upstream.seq.ind and downstream.seq.ind must be single ",
             "integers when withSeq = TRUE")
    }
    test_line <- readLines(con = file, n = skip + 2, ok = TRUE, 
                           warn = TRUE, encoding = "unknown",
                           skipNul = FALSE)[skip + 2]
    if (!grepl("\t", test_line)){
        stop("file is not tab-delimited")
    }
    myPeak <- read.delim(file, skip = skip, header = FALSE,
                         sep = "\t", stringsAsFactors = FALSE)
    if (ncol(myPeak) < 6 || !is.numeric(myPeak[, 2]) ||
        !is.numeric(myPeak[, 3]) || !is.character(myPeak[, 4]) ||
        !is.character(myPeak[, 6])) {
        stop("file ", file, "is not a valid file. file must contain at ", 
             "least 6 fields in the order of chrom, chromStart, strEnd, ",
             "name, score and strand.")
    }
    if (length(unique(myPeak[, 4])) != nrow(myPeak)) {
        stop("unique names must be specified in column 4 of the BED6 file")
    }
    if (!any(myPeak[, 6] %in% c("+", "-"))) {
        stop("strand info must be '+' or '-'")
    }
    if (withSeq && (upstream.seq.ind < 6 || upstream.seq.ind > ncol(myPeak) ||
        downstream.seq.ind < 6 || downstream.seq.ind > ncol(myPeak) ||
        upstream.seq.ind == downstream.seq.ind)) {
        stop("upstream.seq.ind and downstream.seq.ind are not correctly specified")    
    }
    
    pA_grs <- GRanges(seqnames = myPeak[, 1], 
                      ranges = IRanges(start = myPeak[, 2] + 1, 
                                       end = myPeak[, 3],
                                       names = myPeak[, 4]),
                      strand = myPeak[, 6])
    if (withSeq){
        pA_grs$upstream.seq <- myPeak[, upstream.seq.ind]
        pA_grs$downstream.seq <-  myPeak[, downstream.seq.ind]
    }
    # here GRanges is 1-based
    pA_grs
}