#' Build a Naive Bayes Classifier
#' 
#' Computes the conditional a-posterior probabilities of a categorical class
#' variable given independent predictor variables using the Bayes rule.
#' 
#' @param Ndata.NaiveBayes A data.frame, containing features for the negative 
#' training data, described further in \code{\link[cleanUpdTSeq]{data.NaiveBayes}}.
#' @param Pdata.NaiveBayes A data.frame, containing features for the positive 
#' training data, described further in \code{\link[cleanUpdTSeq]{data.NaiveBayes}}.
#' @param genome Name of the genome to get sequences from. To find out a
#' list of available genomes, please type BSgenome::available.genomes() in R.
#' @param upstream An integer(1) vector, length of upstream sequence to retrieve.
#' @param downstream An integer(1) vector, length of downstream sequence to 
#' retrieve.
#' @param wordSize An integer(1) vector,  size of the kmer feature for the
#' upstream sequence. wordSize = 6 should always be used.
#' @param alphabet A character(1) vector, a string containing DNA bases.
#'  By default, "ACTG".
#' @return An object of class "naiveBayes".
#' @author Jianhong Ou
#' @seealso \code{\link[e1071]{naiveBayes}}
#' @importFrom e1071 naiveBayes 
#' @importFrom stats as.formula predict
#' @examples
#' 
#' if (interactive()){
#'     data(data.NaiveBayes)
#'     classifier <- buildClassifier(data.NaiveBayes$Negative, 
#'                                   data.NaiveBayes$Positive)
#' }
#' 
#' @export 

buildClassifier <- function(Ndata.NaiveBayes, 
                            Pdata.NaiveBayes, 
                            upstream = 40L,
                            downstream = 30L, 
                            wordSize = 6L,
                            genome = Drerio, 
                            alphabet = c("ACGT")){
    if ((!class(upstream) %in% c("integer", "numeric")) ||
           (!class(downstream) %in% c("integer", "numeric")) ||
           (!class(wordSize) %in% c("integer", "numeric"))) {
        stop("upstream, downstream, and wordSize must be objects of integer")
    }
    i <- ncol(Ndata.NaiveBayes) - 2
    xnam <- colnames(Ndata.NaiveBayes)[2:i]
    fmla <- as.formula(paste("y ~ ", paste(xnam, collapse = "+")))
    trainingData <- rbind(Pdata.NaiveBayes, Ndata.NaiveBayes)
    classifier <- naiveBayes(fmla, data = trainingData, laplace = 1)
    new("PASclassifier", 
        classifier = classifier, 
        info = new("modelInfo", 
                   upstream = as.integer(upstream), 
                   downstream = as.integer(downstream), 
                   wordSize = as.integer(wordSize),
                   alphabet = alphabet,
                   genome = genome))
}
