predictTestSet <-
    function (Ndata.NaiveBayes, Pdata.NaiveBayes, 
              testSet.NaiveBayes=testSet.NaiveBayes, 
              outputFile="test-predNaiveBayes.tsv", 
              assignmentCutoff=0.5)
        
    {
        i <- length(colnames(Ndata.NaiveBayes)) - 2
        xnam <- colnames(Ndata.NaiveBayes)[2:i]
        fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
        trainingData <- rbind(Pdata.NaiveBayes, Ndata.NaiveBayes)
        n.t <- dim(trainingData)[1]
        
        testSet.NaiveBayes <- cbind(y=0, testSet.NaiveBayes)
        trainingData <- rbind(trainingData, testSet.NaiveBayes)
        classifier <- naiveBayes(fmla, data=trainingData[1:n.t, ], laplace=1)
        
        test.start <- n.t + 1
        test.end <- dim(trainingData)[1]
        pred.prob.test <- predict(classifier, type="raw", 
                                 newdata=trainingData[test.start:test.end,2:i])
##        pred.prob.test[1:2,]
        
        pred.class.test <- unlist(lapply(pred.prob.test[,2], function(p) {
            if (p >assignmentCutoff) {1} else {0} }))
        comp.class.test <- table(pred.class.test, trainingData[test.start:test.end,1])
##        comp.class.test
        pred.names.test <- cbind(as.character(rownames(testSet.NaiveBayes)), 
                                pred.prob.test, 
                                pred.class.test, 
                                as.character(testSet.NaiveBayes[,i+1]), 
                                as.character(testSet.NaiveBayes[,i+2]))
##        dim(pred.names.test)
##        pred.names.test[1,]
        colnames(pred.names.test) <- c( "PeakName", 
                                       "prob False/oligodT internally primed", 
                                       "prob True", 
                                       "pred.class",
                                       "UpstreamSeq", 
                                       "DownstreamSeq")
        write.table(pred.names.test, file=outputFile, sep="\t", row.names=FALSE)  
    }
