buildFeatureVector <- function(peaks, BSgenomeName=Drerio, 
                               upstream=50, downstream=40, 
                               wordSize=6, alphabet=c("ACGT"), 
                               sampleType=c("TP", "TN",  "unknown"),
                               replaceNAdistance=40, 
                               method=c("NaiveBayes", "SVM"), 
                               ZeroBasedIndex=1, fetchSeq=FALSE)
{
    #### use peak end as the separator to get 50bp upstream and 40 bp downstream
    #### for peak at negative strand, it is the reverse complement of the plus strand
    if (fetchSeq == TRUE)
    {
        if (ZeroBasedIndex)
        {
            start(peaks) <- start(peaks) +1
            end(peaks) <- end(peaks) + 1
        }
        seq <- getUpstreamSequence(peaks, upstream=upstream, genome=BSgenomeName)
        upstream.seq <- seq$sequence
        seq <- getDownstreamSequence(peaks , downstream=downstream, 
                                    genome=BSgenomeName)
        downstream.seq <- seq$sequence
    }
    else
    {
        upstream.seq <- peaks$upstream.seq
        end.up <- unlist(lapply(upstream.seq, nchar))
        start.up <- end.up - upstream + 1
        start.up <- unlist(lapply(start.up, function(i) {max(i,1)}))
        upstream.seq <- substr(upstream.seq,start.up, end.up) 
        downstream.seq <- peaks$downstream.seq
        end.down <- unlist(lapply(downstream.seq, function(i) {
            min(nchar(i),downstream)
        }))	
        downstream.seq <- substr(downstream.seq, 1, end.down)
    }
    if (method !="SVM")
    {
        newData <- do.call(rbind, lapply(1:length(upstream.seq), 
                                        function(i) {
                                            t(count(s2c(as.character(upstream.seq[i])), 
                                                    wordsize=wordSize, 
                                                    alphabet=s2c(alphabet)))
                                        }
        ))
        newData[newData>0] <- 1
        newData <- as.data.frame(newData)
        newData[,1:ncol(newData)] <- lapply(newData[,1:ncol(newData)], factor)
    }
    else
    {
        newData <- cbind(as.character(upstream.seq),as.character(downstream.seq))
    }
    n.A <- unlist(lapply(1:length(downstream.seq), function(i) {
        table(factor(s2c(as.character(downstream.seq[i])), levels=c("A")))
    }))
    n.C <- unlist(lapply(1:length(downstream.seq), function(i) {
        table(factor(s2c(as.character(downstream.seq[i])), levels=c("C")))
    }))
    n.T <- unlist(lapply(1:length(downstream.seq), function(i) {
        table(factor(s2c(as.character(downstream.seq[i])), levels=c("T")))
    }))
    n.G <- unlist(lapply(1:length(downstream.seq), function(i) {
        table(factor(s2c(as.character(downstream.seq[i])), levels=c("G")))
    }))	
    
    avg.distanceA2PeakEnd <- unlist(lapply(1:dim(newData)[1], function(i) {
        temp <- mean(grep("A", s2c(as.character(downstream.seq[i]))))
        if (is.na(temp)) {replaceNAdistance} else {temp}
    }))
    
    newData2 <- do.call(rbind, lapply(1:length(downstream.seq),
                                     function(i) {
                                         t(count(s2c(as.character(downstream.seq[i])), 
                                                 wordsize=2, 
                                                 alphabet=s2c(alphabet)))
                                     }
    ))
    newData2 <- as.data.frame(newData2)
    
    if (sampleType == "TP")
    {
        
        newData1 <- data.frame(cbind(rep(1,dim(newData)[1]), 
                                    n.A, 
                                    n.C, 
                                    n.T, 
                                    n.G,
                                    avg.distanceA2PeakEnd,
                                    newData2, 
                                    newData))
        colnames(newData1)[1] <- "y"
        colnames(newData1)[2] <- "n.A.Downstream"
        colnames(newData1)[3] <- "n.C.Downstream"
        colnames(newData1)[4] <- "n.T.Downstream"
        colnames(newData1)[5] <- "n.G.Downstream"
        colnames(newData1)[6] <- "avg.distanceA2PeakEnd"
        id <- 23:ncol(newData1)
    }
    else if (sampleType == "TN")
    {
        
        newData1 <- data.frame(cbind(rep(0,dim(newData)[1]), 
                                    n.A, 
                                    n.C, 
                                    n.T, 
                                    n.G, 
                                    avg.distanceA2PeakEnd, 
                                    newData2,
                                    newData))
        colnames(newData1)[1] <- "y"
        colnames(newData1)[2] <- "n.A.Downstream"
        colnames(newData1)[3] <- "n.C.Downstream"
        colnames(newData1)[4] <- "n.T.Downstream"
        colnames(newData1)[5] <- "n.G.Downstream"
        colnames(newData1)[6] <- "avg.distanceA2PeakEnd"
        id <- 23:ncol(newData1)
    }
    else
    {
        newData1 <- data.frame(cbind(n.A, n.C, n.T, n.G,
                                    avg.distanceA2PeakEnd, 
                                    newData2,
                                    newData))
        colnames(newData1)[1] <- "n.A.Downstream"
        colnames(newData1)[2] <- "n.C.Downstream"
        colnames(newData1)[3] <- "n.T.Downstream"
        colnames(newData1)[4] <- "n.G.Downstream"
        colnames(newData1)[5] <- "avg.distanceA2PeakEnd"	
        id <- 22:ncol(newData1)
    }
    
    newData1$n.A.Downstream <- as.numeric(as.character(newData1$n.A.Downstream))
    newData1$n.C.Downstream <- as.numeric(as.character(newData1$n.C.Downstream))
    newData1$n.T.Downstream <- as.numeric(as.character(newData1$n.T.Downstream))
    newData1$n.G.Downstream <- as.numeric(as.character(newData1$n.G.Downstream))
    newData1$avg.distanceA2PeakEnd <- as.numeric(as.character(newData1$avg.distanceA2PeakEnd))
    if (method != "SVM")
    {
        newData1 <- cbind(newData1, upstream.seq, downstream.seq)
    }
    else
    {
        nParam <- dim(newData1)[2]
        u.index <- nParam-1
        colnames(newData1)[u.index:nParam] <- c("upstream.seq", "downstream.seq")
    }
    if (fetchSeq)
    {
        rownames(newData1) <- as.character(names(seq))
    }
    else
    {
        rownames(newData1) <- as.character(names(peaks))
    }
    newData1
}
