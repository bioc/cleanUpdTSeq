#' predict authenticity of putative pA sites
#' 
#' classify putative pA sites into true and false bins.
#' 
#' @param Ndata.NaiveBayes A data.frame, containing features for the negative 
#' training data, which is built using the function 
#' \code{\link{buildFeatureVector}}. It is described further
#' in\code{\link[cleanUpdTSeq]{data.NaiveBayes}}.
#' @param Pdata.NaiveBayes A data.frame, containing features for the positive 
#' training data, which is built using the function 
#' \code{\link{buildFeatureVector}}. It is described further
#' in\code{\link[cleanUpdTSeq]{data.NaiveBayes}}.
#' @param classifier An object of class \linkS4class{PASclassifier}.
#' @param testSet.NaiveBayes An object of \code{\linkS4class{featureVector}} for
#' test data built for Naive Bayes analysis using the function 
#' \code{\link{buildFeatureVector}}.
#' @param outputFile A character(1) vector, file name for outputting prediction
#' results. The prediction output is written to the file, tab separated.
#' @param assignmentCutoff A numeric(1) vector, specifying the cutoff for
#' classifying a putative pA site into a true or false pA class. It should be 
#' any number between 0 and 1. For example, assignmentCutoff = 0.5 will assign
#' an putative pA site with prob_true_pA > 0.5 to the True class (1),
#' and any putative pA site with prob_true_pA < = 0.5 as False (0).
#' @return  A data.frame including all info as described below. The upstream 
#' and downstream sequence used in assessing the putative pA site might be 
#' included when return_sequences = TRUE. 
#' \item{peak_name}{the name of the putative pA site (originally from the 4th field 
#' in the bed file).} 
#' \item{prob_fake_pA}{the probability that the putative pA site is false}
#' \item{prob_true_pA}{the probability that the putative pA site is true} 
#' \item{pred_class}{the predicted class of the putative pA site, 
#'   based on the assignment cutoff. 0 = Falsee/oligo(dT) internally primed,
#'   1 = True} 
#' \item{upstream_seq}{the upstream sequence of the putative pA site  used in 
#' the analysis}
#' \item{downstream_seq}{the downstream sequence of the putative pA site
#'   used in the analysis.}
#' @param return_sequences A logical(1) vector, indicating whether upstream and
#'   downstream sequences should be included in the output
#' @author Sarah Sheppard, Haibo Liu, Jianhong Ou, Nathan Lawson, Lihua J. Zhu
#' @importFrom utils write.table
#' @references Sheppard S, Lawson ND, Zhu LJ. Accurate identification of
#'   polyadenylation sites from 3' end deep sequencing using a naive Bayes
#'   classifier. Bioinformatics. 2013;29(20):2564-2571.
#'   
#' @examples
#' library(BSgenome.Drerio.UCSC.danRer7)
#' testFile <- system.file("extdata", "test.bed",
#'                         package = "cleanUpdTSeq")
#' ## convert the test set to GRanges without upstream and downstream sequence
#' ## information
#' peaks <- BED6WithSeq2GRangesSeq(file = testFile, 
#'                                skip = 1L, withSeq = TRUE)
#' ## build the feature vector for the test set without sequence information
#' testSet.NaiveBayes = buildFeatureVector(peaks,
#'                                         genome = Drerio, 
#'                                         upstream = 40L,
#'                                         downstream = 30L, 
#'                                         wordSize = 6L, 
#'                                         alphabet = c("ACGT"),
#'                                         sampleType = "unknown",
#'                                         replaceNAdistance = 30,
#'                                         method = "NaiveBayes", 
#'                                         fetchSeq = TRUE,
#'                                         return_sequences = TRUE)
#' data(data.NaiveBayes)
#' ## sample the test data for code testing, DO NOT do this for real data
#' samp <- c(1:22, sample(23:4118, 50), 4119, 4120)
#' Ndata.NaiveBayes <- data.NaiveBayes$Negative[, samp]
#' Pdata.NaiveBayes <- data.NaiveBayes$Positive[, samp]
#' testSet.NaiveBayes@data <- testSet.NaiveBayes@data[, samp[-1]-1]
#'     
#' test_out <- predictTestSet(Ndata.NaiveBayes, 
#'                            Pdata.NaiveBayes,
#'                            testSet.NaiveBayes,
#' 	                          outputFile = tempfile(), 
#'                            assignmentCutoff = 0.5)
#' 
#' @export
predictTestSet <- function(Ndata.NaiveBayes = NULL,
                           Pdata.NaiveBayes = NULL, 
                           testSet.NaiveBayes,
                           classifier = NULL,
                           outputFile = "test-predNaiveBayes.tsv", 
                           assignmentCutoff = 0.5,
                           return_sequences = FALSE){
    if (missing(testSet.NaiveBayes)) {stop("testSet.NaiveBayes is required")}
    if (!is(testSet.NaiveBayes, "featureVector")){
        stop("testSet.NaiveBayes must be an object of class ",
             "\"featureVector\"")
    }
    i <- ncol(testSet.NaiveBayes@data) - 1
    if (!is.null(classifier)) {
        if (!is(classifier, "PASclassifier")) {
            stop("classifier must be an object of class \"PASclassifier\"")
        }
        if (identical(classifier@info, testSet.NaiveBayes@info)) {
            classifier <- classifier@classifier
        } else {
            stop("upstream, downstream, wordSize, and alphabet of ", 
                 "classifier and testSet.NaiveBayes must be the same")
        }
    } else {
        if (is.null(Ndata.NaiveBayes) || is.null(Pdata.NaiveBayes)) {
            stop("Positive and negative data for training a Naive Bayes ",
            "classifier is needed")
        }
        xnam <- colnames(Ndata.NaiveBayes)[2:i]
        fmla <- as.formula(paste("y ~ ", paste(xnam, collapse = "+")))
        trainingData <- rbind(Pdata.NaiveBayes, Ndata.NaiveBayes)
        classifier <- naiveBayes(fmla, data = trainingData, laplace = 1)
    }
    testSet.NaiveBayes <- testSet.NaiveBayes@data
    pred.prob.test <- as.data.frame(predict(classifier, type = "raw", 
                                        newdata = testSet.NaiveBayes))
    pred.class.test <- as.numeric(pred.prob.test[,2] > assignmentCutoff)
    
    col_names <- c("peak_name", "prob_fake_pA", 
                   "prob_true_pA", "pred_class",
                   "upstream_seq", "downstream_seq")
    if (return_sequences){
        pred.names.test <- cbind(rownames(testSet.NaiveBayes), 
                                 pred.prob.test, 
                                 pred.class.test, 
                                 testSet.NaiveBayes[,i], 
                                 testSet.NaiveBayes[,i+1], 
                                 stringsAsFactors = FALSE)
        colnames(pred.names.test) <- col_names
    } else {
        pred.names.test <- cbind(rownames(testSet.NaiveBayes), 
                                 pred.prob.test, 
                                 pred.class.test, 
                                 stringsAsFactors = FALSE)
        colnames(pred.names.test) <- col_names[1:4]
    }
    if (!is.null(outputFile))
    {
        if (!file.create(outputFile)){
            outputFile <-  paste0(basename(tempfile()),".prediction.results.txt")
            message("outputFile is not a valid file name. output will be ",
                    "written to ", outputFile)
        }
        write.table(pred.names.test, file = outputFile, quote = FALSE, 
                    sep = "\t", row.names = FALSE)
    }
    pred.names.test
}
