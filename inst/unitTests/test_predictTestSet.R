test_predictTestSet <- function(){
    tmpfile <- tempfile()
    testFile = system.file("extdata", "test.bed", package="cleanUpdTSeq")
    load(system.file("extdata", "test-predNaiveBayes.rds", package="cleanUpdTSeq"))
    predictTestSet(Ndata.NaiveBayes, 
                   Pdata.NaiveBayes,
                   testSet.NaiveBayes,
                   outputFile=tmpfile, 
                   assignmentCutoff = 0.5)
    x <- read.delim(system.file("extdata", "test-predNaiveBayes.xls", package="cleanUpdTSeq"))
    y <- read.delim(tmpfile)
    checkIdentical(x, y)
    unlink(tmpfile)
    classifier <- buildClassifier(Ndata.NaiveBayes, 
                                  Pdata.NaiveBayes)
    predictTestSet(testSet.NaiveBayes=testSet.NaiveBayes,
                   classifier=classifier,
                   outputFile=tmpfile, 
                   assignmentCutoff = 0.5)
    z <- read.delim(tmpfile)
    checkIdentical(x, z)
    unlink(tmpfile)
}