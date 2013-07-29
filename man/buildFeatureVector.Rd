\name{buildFeatureVector}
\alias{buildFeatureVector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
buildFeatureVector
}
\description{
This function creates a data frame. Fields include peak name, upstream sequence, downstream sequence, and features to be used in classifying the putative polyadenylation site. 
}
\usage{
buildFeatureVector(peaks, BSgenomeName = Drerio, upstream = 50, 
	downstream = 40, wordSize = 6, alphabet = c("ACGT"),
 sampleType = c("TP", "TN", "unknown"), replaceNAdistance = 40,
 method = c("NaiveBayes", "SVM"), ZeroBasedIndex = 1, fetchSeq = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{
An object of GRanges that may contain the upstream and downstream sequence information. This item is created by the function \link{BED2GRangesSeq}.
}
  \item{BSgenomeName}{
Name of the genome to use to get sequence. To find out a list of available genomes, please type available.genomes() in R.
}
  \item{upstream}{
This is the length of upstream sequence to use in the analysis.
}
  \item{downstream}{
This is the length of downstream sequence to use in the analysis.
}
  \item{wordSize}{
This is the size of the word to use as a feature for the upstream sequence. wordSize = 6 should always be used. }
  \item{alphabet}{
These are the bases to use, for example DNA bases ACTG.
}
  \item{sampleType}{
This is the type of sample. For example TP (true positive) or TN (true negative) for training data or unknown for test data.
}
  \item{replaceNAdistance}{
If there is no A in the downstream sequence, then use this for  the average distance of As to the putative polyadenylation site.
}
  \item{method}{
This is which machine learning method to specify. For this release, method should always be set to NaiveBayes.
}
  \item{ZeroBasedIndex}{
If the coordinates are set using Zero Based indexing, set this = 1.
}
  \item{fetchSeq}{
Boolean, for getting upstream and downstream sequence at this step or not. 
}
}
\value{
  An object of data frame. Fields include peak name, upstream sequence, downstream sequence, and features to be used in classifying the putative polyadenylation site
}
\author{
Sarah Sheppard, Jianhong Ou, Nathan Lawson, Lihua Julie Zhu
}
\examples{
testFile = system.file("extdata", "test.bed", package="cleanUpdTSeq")
testSet = read.table(testFile, sep = "\t", header = TRUE)
	
#convert the test set to GRanges with upstream and downstream sequence information
peaks = BED2GRangesSeq(testSet[1:10, ], upstream.seq.ind = 7, downstream.seq.ind = 8, withSeq=TRUE)
#build the feature vector for the test set with sequence information 
testSet.NaiveBayes = buildFeatureVector(peaks,BSgenomeName = Drerio, upstream = 40,
	   downstream = 30, wordSize = 6, alphabet=c("ACGT"),
	   sampleType = "unknown",replaceNAdistance = 30, 
	   method = "NaiveBayes", ZeroBasedIndex = 1, fetchSeq = FALSE)
    
#convert the test set to GRanges without upstream and downstream sequence information
peaks = BED2GRangesSeq(testSet[1:10, ],withSeq=FALSE)
        
#build the feature vector for the test set without sequence information
testSet.NaiveBayes = buildFeatureVector(peaks,BSgenomeName = Drerio, upstream = 40,
        downstream = 30, wordSize = 6, alphabet=c("ACGT"),
        sampleType = "unknown",replaceNAdistance = 30,
        method = "NaiveBayes", ZeroBasedIndex = 1, fetchSeq = TRUE)
}

\keyword{ misc }