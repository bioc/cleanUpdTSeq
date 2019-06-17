\name{cleanUpdTSeq-package}
\alias{cleanUpdTSeq-package}
\alias{cleanUpdTSeq}
\docType{package}
\title{
This package classifies putative polyadenylation sites.
}
\description{
3'ends of transcripts have generally been poorly annotated. With the advent of deep sequencing, many methods have been developed to identify 3'ends. The majority of these methods use an oligodT primer which can bind to internal adenine-rich sequences, and lead to artifactual identification of polyadenylation sites. Heuristic filtering methods rely on a certain number of As downstream of a putative polyadenylation site to classify the site as true or oligodT primed. This package provides a robust method to classify putative polyadenylation sites using a Naive Bayes classifier.
}
\details{
\tabular{ll}{
Package: \tab cleanUpdTSeq\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-07-22\cr
License: \tab GPL-2\cr
}
}
\author{
Sarah Sheppard, Jianhong Ou, Nathan Lawson, Lihua Julie Zhu
Maintainer: Sarah Sheppard <Sarah.Sheppard@umassmed.edu>, Jianhong Ou <Jianhong.Ou@umassmed.edu>, Lihua Julie Zhu <Julie.Zhu@umassmed.edu>
}
\references{
1. Meyer, D., et al., e1071: Misc Functions of the Department of Statistics (e1071), TU Wien. 2012.

2. Pages, H., BSgenome: Infrastructure for Biostrings-based genome data packages.

3. Sheppard, S., Lawson, N.D. and Zhu, L.J., 2013. Accurate identification of polyadenylation sites from 3' end deep sequencing using a naive Bayes classifier. Bioinformatics, 29(20), pp.2564-2571.
}

\examples{
#read in a test set
#### first install the package using the following command
#### BiocManager::install("cleanUpdTSeq")
if (interactive())
{
	library(cleanUpdTSeq)
	testFile = system.file("extdata", "test.bed", package="cleanUpdTSeq")
	testSet = read.table(testFile, sep = "\t", header = TRUE)
	
	#convert the test set to GRanges with upstream and downstream sequence information
	peaks = BED2GRangesSeq(testSet,upstream.seq.ind = 7, downstream.seq.ind = 8, withSeq=TRUE)
	#build the feature vector for the test set with sequence information 
	library(BSgenome.Drerio.UCSC.danRer7)
	testSet.NaiveBayes = buildFeatureVector(peaks,BSgenomeName = Drerio, upstream = 40,
	 downstream = 30, wordSize = 6, alphabet=c("ACGT"),
	 sampleType = "unknown",replaceNAdistance = 30, 
	method = "NaiveBayes", ZeroBasedIndex = 1, fetchSeq = FALSE)
	
	#convert the test set to GRanges without upstream and downstream sequence information
        peaks = BED2GRangesSeq(testSet,withSeq=FALSE)
        
	#build the feature vector for the test set without sequence information
	testSet.NaiveBayes = buildFeatureVector(peaks,BSgenomeName = Drerio, upstream = 40,
         downstream = 30, wordSize = 6, alphabet=c("ACGT"),
         sampleType = "unknown",replaceNAdistance = 30,
        method = "NaiveBayes", ZeroBasedIndex = 1, fetchSeq = TRUE)

	#predict the test set
	data(data.NaiveBayes)
	predictTestSet(data.NaiveBayes$Negative, data.NaiveBayes$Positive, testSet.NaiveBayes,
	outputFile = "test-predNaiveBayes.tsv", assignmentCutoff = 0.5)
}
}

\keyword{ package }
