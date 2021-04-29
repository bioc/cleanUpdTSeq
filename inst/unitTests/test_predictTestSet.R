test_predictTestSet <- function(){
    library(BSgenome.Drerio.UCSC.danRer7)
    tmpfile <- tempfile()
    load(system.file("extdata", "test-predNaiveBayes.rda", 
                     package = "cleanUpdTSeq"))
    out <- predictTestSet(Ndata.NaiveBayes,
                          Pdata.NaiveBayes,
                          testSet.NaiveBayes,
                          outputFile = tmpfile, 
                          assignmentCutoff = 0.5,
                          return_sequences = TRUE)

    x <- read.delim(system.file("extdata", 
                                "test-predNaiveBayes.txt",
                                package = "cleanUpdTSeq"), 
                    stringsAsFactors =  FALSE)
    y <- read.delim(tmpfile)
    checkIdentical(format(x, digits = 7), format(y, digits = 7))

    classifier <- buildClassifier(Ndata.NaiveBayes, 
                                  Pdata.NaiveBayes)
    out <- predictTestSet(testSet.NaiveBayes = testSet.NaiveBayes,
                   classifier = classifier,
                   outputFile = tmpfile, 
                   assignmentCutoff  =  0.5,
                   return_sequences = TRUE)
    z <- read.delim(tmpfile)
    checkIdentical(format(x, digits = 7), format(z, digits = 7))
}