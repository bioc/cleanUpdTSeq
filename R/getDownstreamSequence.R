getDownstreamSequence <-
    function(peaks, downstream=20, genome)
    {
        if (missing(genome) || class(genome) != "BSgenome") {
            stop("genome is required BSgenome object!")
        }
        if (missing(peaks) || class(peaks) != "GRanges")
        {
            stop("peaks is required as GRanges object!")
        }        
        
        peaksStrand <- strand(peaks)
        plus.peaks <- peaks[peaksStrand %in% c("*", "+", "1")]
        minus.peaks <- peaks[peaksStrand %in% c("-", "-1")]
        
        Start <- c(end(plus.peaks) + 1, 
                  start(minus.peaks)-1- as.numeric(downstream))
        End <- Start + as.numeric(downstream)
        
        strand <- c(rep("+", length(start(plus.peaks))), 
                   rep("-", length(start(minus.peaks))))
        
        chr <- as.character(c(seqnames(plus.peaks), seqnames(minus.peaks)))
        
        ##clean up the chromosome names.
        if(!all(chr %in% seqnames(genome))){
            message("Not all the sequence names are in the genome.\nTry to convert.")
            if(any(grep("^chr", seqnames(genome)))){
                pattern <- "^(\\d+|V?I{0,3}|IV|MT|M|X|Y)$"
                chr[grepl(pattern, chr)] <- paste("chr",
                                                  chr[grepl(pattern, chr)],
                                                  sep="")
            } else {
                chr <- gsub("^chr", chr, ignore.case=TRUE)
            }
        }
        if(!all(chr %in% seqnames(genome))){
            stop("Not all the sequence names are in the genome.\nPlease doulbe check the seqnames.")
        }
        
        ##extract sequence from genome
        seqgr <- GRanges(seqnames=chr, 
                         ranges=IRanges(start=Start, 
                                        end=ifelse(End < seqlengths(genome)[chr],
                                                   End,
                                                   seqlengths(genome)[chr])),
                         strand=strand)
        seq <- getSeq(genome, seqgr, as.character=TRUE)
        
        ##export GRanges
        seqinfo <- seqinfo(genome)
        seqinfo <- seqinfo[unique(chr)[order(unique(chr))]]
        GRanges(seqnames=chr,
                ranges=IRanges(start=c(start(plus.peaks), start(minus.peaks)),
                               end=c(end(plus.peaks), end(minus.peaks)),
                               names=c(names(plus.peaks), names(minus.peaks))),
                strand=strand,
                sequence=seq,
                seqinfo=seqinfo
        )
    }
