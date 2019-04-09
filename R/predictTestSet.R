predictTestSet <-
    function (Ndata.NaiveBayes, Pdata.NaiveBayes, 
              testSet.NaiveBayes,
              classifier=NULL,
              outputFile="test-predNaiveBayes.tsv", 
              assignmentCutoff=0.5)
        
    {
        if(missing(testSet.NaiveBayes))
            stop("testSet.NaiveBayes is required")
        if(!is(testSet.NaiveBayes, "featureVector"))
            stop("testSet.NaiveBayes must be an object of class \"featureVector\"")
        i <- length(colnames(testSet.NaiveBayes@data)) - 1
        if(!is.null(classifier)){
            if(!is(classifier, "PASclassifier"))
                stop("classifier must be an object of class \"PASclassifier\"")
            if(classifier@info@upstream==testSet.NaiveBayes@info@upstream &&
                   classifier@info@downstream==testSet.NaiveBayes@info@downstream &&
                   classifier@info@wordSize==testSet.NaiveBayes@info@wordSize &&
                   classifier@info@alphabet==testSet.NaiveBayes@info@alphabet){
                if(organism(classifier@info@genome)!=organism(testSet.NaiveBayes@info@genome))
                    message("genome of classifier is different from testSet.NaiveBayes.")
                classifier <- classifier@classifier
            }else{
                stop("upstream, downstream wordSize and alphabet of classifier and testSet.NaiveBayes must be same")
            }
        }else{
            xnam <- colnames(Ndata.NaiveBayes)[2:i]
            fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
            trainingData <- rbind(Pdata.NaiveBayes, Ndata.NaiveBayes)
            classifier <- naiveBayes(fmla, data=trainingData, laplace=1)
        }
        
        testSet.NaiveBayes <- testSet.NaiveBayes@data
        
        pred.prob.test <- predict(classifier, type="raw", 
                                 newdata=testSet.NaiveBayes)
        pred.class.test <- as.numeric(pred.prob.test[,2] > assignmentCutoff)
        pred.names.test <- cbind(as.character(rownames(testSet.NaiveBayes)), 
                                pred.prob.test, 
                                pred.class.test, 
                                as.character(testSet.NaiveBayes[,i]), 
                                as.character(testSet.NaiveBayes[,i+1]))

        colnames(pred.names.test) <- c( "PeakName", 
                                       "prob False/oligodT internally primed", 
                                       "prob True", 
                                       "pred.class",
                                       "UpstreamSeq", 
                                       "DownstreamSeq")
        if((!is.null(outputFile)) && (!is.na(outputFile)) && nchar(outputFile)>0) 
            write.table(pred.names.test, file=outputFile, sep="\t", row.names=FALSE)
        
        pred.names.test <- as.data.frame(pred.names.test, stringsAsFactors=FALSE)
        mode(pred.names.test[,2]) <- "numeric"
        mode(pred.names.test[,3]) <- "numeric"
        return(invisible(pred.names.test))
    }
