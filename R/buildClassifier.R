buildClassifier <- function(Ndata.NaiveBayes, Pdata.NaiveBayes, 
                            upstream=40L, downstream=30L, wordSize=6L,
                            genome=Drerio, alphabet=c("ACGT")){
    if((!class(upstream) %in% c("integer", "numeric")) ||
           (!class(downstream) %in% c("integer", "numeric")) ||
           (!class(wordSize) %in% c("integer", "numeric")))
        stop("upstream, downstream and wordSize must be objects of integer")
    upstream <- as.integer(upstream)
    downstream <- as.integer(downstream)
    wordSize <- as.integer(wordSize)
    i <- length(colnames(Ndata.NaiveBayes)) - 2
    xnam <- colnames(Ndata.NaiveBayes)[2:i]
    fmla <- as.formula(paste("y ~ ", paste(xnam, collapse = "+")))
    trainingData <- rbind(Pdata.NaiveBayes, Ndata.NaiveBayes)
    classifier <- naiveBayes(fmla, data = trainingData, laplace = 1)
    new("PASclassifier", classifier=classifier, 
        info=new("modelInfo", upstream=upstream, downstream=downstream, 
                 wordSize=wordSize, genome=genome, alphabet=alphabet))
}