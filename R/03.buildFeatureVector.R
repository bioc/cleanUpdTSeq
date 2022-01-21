#' build Feature Vector_2
#' 
#' This function creates a data frame. Fields include peak name, upstream
#' sequence, downstream sequence, and features to be used in classifying the
#' putative polyadenylation site.
#' 
#' 
#' @param peaks An object of GRanges that may contain the upstream and
#' downstream sequence information. This item is created by the function
#' \link{BED6WithSeq2GRangesSeq}.
#' @param genome Name of the genome to get sequences from. To find out a
#' list of available genomes, please type BSgenome::available.genomes() in R.
#' @param upstream An integer(1) vector, length of upstream sequence to retrieve.
#' @param downstream An integer(1) vector, length of downstream sequence to 
#' retrieve.
#' @param wordSize An integer(1) vector,  size of the kmer feature for the
#' upstream sequence. wordSize = 6 should always be used.
#' @param alphabet A character(1) vector, a string containing DNA bases.
#'  By default, "ACTG".
#' @param sampleType A character(1) vector, indicating type of sequences for 
#' building feature vectors. Options are TP (true positive) and TN 
#' (true negative) for training data, or unknown for test data.
#' @param replaceNAdistance An integer(1) vector, specifying an number for
#' avg.distanceA2PeakEnd, the average distance of As to the putative pA site,
#' when there is no A in the downstream sequence.
#' @param method A character(1) vector, specifying a machine learning method 
#' to  to use. Currently, only "NaiveBayes" is implemented.
#' @param fetchSeq A logical (1), indicating whether upstream and downstream 
#' sequences should be retrieved from the BSgenome object at this step or not.
#' @param return_sequences A logical(1) vector, indicating whether upstream and
#'   downstream sequences should be included in the output
#' @return An object of "\code{\link{featureVector}}"
#' @author Sarah Sheppard, Haibo Liu, Jianhong Ou, Nathan Lawson, Lihua J. Zhu
#' @importFrom Biostrings DNAStringSet alphabetFrequency
#' @importFrom stringr str_locate_all
#' @importFrom GenomicRanges strand start end mcols seqnames
#' @importFrom seqinr s2c count
#' @importFrom S4Vectors metadata
#' @import BSgenome.Drerio.UCSC.danRer7
#' @importMethodsFrom GenomicRanges end<- start<-
#' 
#' @examples
#' library(BSgenome.Drerio.UCSC.danRer7)
#' testFile <- system.file("extdata", "test.bed",
#'                         package = "cleanUpdTSeq")
#' peaks <- BED6WithSeq2GRangesSeq(file = testFile, 
#'                                skip = 1L, withSeq = TRUE)
#' ## build the feature vector for the test set with sequence information 
#' testSet.NaiveBayes = buildFeatureVector(peaks,
#'                                         genome = Drerio, 
#'                                         upstream = 40L,
#'                                         downstream = 30L, 
#'                                         wordSize = 6L, 
#'                                         alphabet = "ACGT",
#'                                         sampleType = "unknown",
#'                                         replaceNAdistance = 30, 
#'                                         method = "NaiveBayes", 
#'                                         fetchSeq = FALSE,
#'                                         return_sequences = TRUE)
#' 
#' ## convert the test set to GRanges without upstream and downstream 
#' ## sequence information
#' peaks <- BED6WithSeq2GRangesSeq(file = testFile, 
#'                                skip = 1L, withSeq = FALSE)
#' #build the feature vector for the test set without sequence information
#' testSet.NaiveBayes = buildFeatureVector(peaks,
#'                                         genome = Drerio, 
#'                                         upstream = 40L,
#'                                         downstream = 30L, 
#'                                         wordSize = 6L,
#'                                         alphabet = "ACGT",
#'                                         sampleType = "unknown",
#'                                         replaceNAdistance = 30,
#'                                         method = "NaiveBayes", 
#'                                         fetchSeq = TRUE,
#'                                         return_sequences = TRUE)
#' 
#' @export
#' 

buildFeatureVector <- function(peaks, 
                                 genome = Drerio, 
                                 upstream = 40L, 
                                 downstream = 30L, 
                                 wordSize = 6L,
                                 alphabet = "ACGT", 
                                 sampleType = c("TP", "TN", "unknown"),
                                 replaceNAdistance = 30L, 
                                 method = c("NaiveBayes", "SVM"), 
                                 fetchSeq = FALSE,
                                 return_sequences = FALSE)
{
    if((!class(upstream) %in% c("integer", "numeric")) ||
           (!class(downstream) %in% c("integer", "numeric")) ||
           (!class(wordSize) %in% c("integer", "numeric")))
        stop("upstream, downstream and wordSize must be integers")
    if(!is(genome, "BSgenome"))
        stop("genome must be an object BSgenome")
    ## use peak end as the separator to get 50bp upstream and 40 bp 
    ## downstream for peak at negative strand, it is the reverse 
    ## complement of the plus strand
    if (fetchSeq == TRUE)
    {
        sequences <- getContextSequences(peaks,
                                         upstream = upstream, 
                                         downstream = downstream,
                                         genome = genome)
        upstream.seq <- sequences$upstream.seq
        downstream.seq <- sequences$downstream.seq
    } else {
        upstream.seq <- peaks$upstream.seq
        end.up <- vapply(upstream.seq, nchar, numeric(1))
        start.up <- end.up - upstream + 1
        start.up <- vapply(start.up, function(i) {max(i, 1)}, 
                           numeric(1))
        upstream.seq <- substr(upstream.seq, start.up, end.up) 
        downstream.seq <- peaks$downstream.seq
        end.down <- vapply(downstream.seq, function(i) {
                           min(nchar(i), downstream)}, numeric(1))
        downstream.seq <- substr(downstream.seq, 1, end.down)
    }
    
    ## upstream hexamer presence/absence features
    if (method == "NaiveBayes") {
        upstream_hexamer <- do.call(rbind, 
             lapply(1:length(upstream.seq), 
                    function(i) {
                    t(count(s2c(upstream.seq[i]), 
                      wordsize = wordSize, 
                      alphabet = s2c(alphabet)))}))
        upstream_hexamer[upstream_hexamer > 0] <- "1"
        upstream_hexamer <- as.data.frame(upstream_hexamer,
                                          stringsAsFactors = TRUE)
    } else {
        upstream_hexamer <- cbind(upstream.seq, downstream.seq)
    }
    
    ## Downstream features: mono-, dinucleotide frequency, mean distances 
    ## between A and pA sites.
    mononuc_freq <- 
        as.data.frame(alphabetFrequency(DNAStringSet(downstream.seq), 
                                        baseOnly = TRUE)[, -5, drop=FALSE])
    A_pos <- str_locate_all(pattern = "A", downstream.seq)
    avg.distanceA2PeakEnd <- 
        vapply(A_pos, function(.x) {mean(.x[,1])}, numeric(1))
    avg.distanceA2PeakEnd[is.na(avg.distanceA2PeakEnd)] <- replaceNAdistance 
    
    dinuc_freq <- 
        do.call(rbind, 
                lapply(1:length(downstream.seq),
                               function(i) {
                               t(count(s2c(downstream.seq[i]), 
                                       wordsize = 2, 
                                       alphabet = s2c(alphabet)))}))
    
    feature_df <- cbind(mononuc_freq,
                        avg.distanceA2PeakEnd, 
                        dinuc_freq,
                        upstream_hexamer)
                                   
    feature_names <- c("y","n.A.Downstream", 
                        "n.C.Downstream",
                        "n.G.Downstream",
                        "n.T.Downstream",
                        "avg.distanceA2PeakEnd")
    if (sampleType ==  "TP")
    {
        feature_df <- cbind(rep(1, nrow(upstream_hexamer)), feature_df)
        colnames(feature_df)[1:6] <- feature_names
    } else if (sampleType == "TN") {
        feature_df <- cbind(rep(0, nrow(upstream_hexamer)), feature_df)
        colnames(feature_df)[1:6] <- feature_names
    } else {
        colnames(feature_df)[1:5] <- feature_names[-1]	
    }
    if (return_sequences)
    {
        if (method == "NaiveBayes")
        {
            feature_df <- cbind(feature_df, upstream.seq, downstream.seq)
        } else {
            nParam <- ncol(feature_df)
            colnames(feature_df)[(nParam - 1):nParam] <- 
                c("upstream.seq", "downstream.seq")
        }
    } 
    rownames(feature_df) <- names(peaks)
    
    md <- metadata(genome)
    feature_vectors <- new("featureVector", 
                           data = feature_df,
                           info = new("modelInfo", 
                                    upstream = as.integer(upstream), 
                                    downstream = as.integer(downstream), 
                                    wordSize = as.integer(wordSize), 
                                    alphabet = alphabet,
                                    genome = metadata(genome)$genome,
                                    metadata = md))
}
