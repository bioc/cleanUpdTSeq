#' Retrieve upstream and downstream sequences
#' 
#' Retrieve upstream and downstream sequences of pA sites from a BSgenome object
#' based on a GRanges object
#'
#' @param peaks An object of GRanges representing pA sites
#' @param upstream An integer(1) vector, length of upstream sequence of pA 
#'   sites, including pA site.
#' @param downstream An integer(1) vector, length of downstream sequences of 
#'   pA sites
#' @param genome An object of BSgenome.
#' @return A data.frame containing sequences upstream and downstream pA sites:
#' \describe{
#'     \item{upstream.seq}{sequence upstream pA site, including pA site}
#'     \item{downstream.seq}{sequence downstream pA site}
#' }
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Biostrings getSeq
#' @author Haibo Liu
#' @export
#'
#' @examples
#' library(BSgenome.Drerio.UCSC.danRer7)
#' testFile <- system.file("extdata", "test.bed",
#'                         package = "cleanUpdTSeq")
#' peaks <- BED6WithSeq2GRangesSeq(file = testFile, 
#'                                skip = 1L, withSeq = FALSE)
#' peaks_seq <- getContextSequences(peaks, 
#'                                  upstream = 40L, 
#'                                  downstream = 30L,
#'                                  genome = Drerio) 
#'                                                         
getContextSequences <- function(peaks,
                                upstream = 40L, 
                                downstream = 30L,
                                genome)
{
    if (missing(genome) || !is(genome, "BSgenome")) {
        stop("genome must be a BSgenome object!")
    }
    if (missing(peaks) || !is(peaks, "GRanges")) {
        stop("peaks must be a GRanges object!")
    }
    if (!is.numeric(upstream) || !is.numeric(downstream) ||
        length(upstream) != 1 || length(downstream) != 1 ||
        any(grepl("\\.", c(upstream, downstream))))
    {
        stop("upstream and downstream must be single integers")
    }
    if(!all(unique(seqnames(peaks)) %in% seqnames(genome))){
        stop("Not all the sequence names are in the genome.", 
             "\nPlease doulbe check the seqnames.")
    }
    peaks$pA.site <- start(peaks)
    ## watch for out of bound cases
    start(peaks) <- ifelse(strand(peaks) == "+", 
                           start(peaks) - upstream + 1, 
                           start(peaks) - downstream)
    start(peaks) <- ifelse(start(peaks) < 1, 1, start(peaks))
    
    end(peaks) <- ifelse(strand(peaks) == "+", 
                         end(peaks) + downstream, 
                         end(peaks) + upstream -1)
    chr_length <- seqlengths(genome)[as.character(seqnames(peaks))]
    end(peaks) <- ifelse(end(peaks) > chr_length, chr_length, end(peaks))
    peaks$pA.site <- ifelse(strand(peaks) == "+", 
                            peaks$pA.site - start(peaks) + 1,
                            end(peaks) - peaks$pA.site +1)
    seq <- getSeq(genome, peaks, as.character=TRUE)
    data.frame(upstream.seq = substr(seq, 1, peaks$pA.site),
               downstream.seq = substr(seq, peaks$pA.site + 1, 
                                       upstream + downstream), 
               stringsAsFactors = FALSE)
}
