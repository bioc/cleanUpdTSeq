\name{predictTestSet}
\alias{predictTestSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
predictTestSet
}
\description{
This function can be used to predict the probabilities for a set of putative pA sites.
}
\usage{
predictTestSet(Ndata.NaiveBayes, Pdata.NaiveBayes, testSet.NaiveBayes, classifier=NULL,
 outputFile = "test-predNaiveBayes.tsv", assignmentCutoff = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ndata.NaiveBayes}{
This is the negative training data, described further in \link{data.NaiveBayes}.
}
  \item{Pdata.NaiveBayes}{
This is the positive training data, described further in \link{data.NaiveBayes}.
}
  \item{classifier}{An object of class "naiveBayes".}
  \item{testSet.NaiveBayes}{
This is the test data, a feature vector that has been built for Naive Bayes analysis using the function "buildFeatureVector".
}
  \item{outputFile}{
This is the name of the file the output will be written to.
}
  \item{assignmentCutoff}{
This is the cutoff used to assign whether a putative pA is true or false. This can be any floating point number between 0 and 1. For example, assignmentCutoff = 0.5 will assign an putative pA site with prob.1 > 0.5 to the True class (1), and any putative pA site with prob.1 <= 0.5 as False (0).
}
}
\value{
The output is written to a tab separated file containing fields for peak name, the probability of the putative pA site being false (prob.0), the probability of the putative pA site being true (prob.1), the predicted class (0/False or 1/True) depending on the assignment cutoff, and the upstream and downstream sequence used in assessing the putative pA site.
	\item{PeakName }{This is the name of the putative pA site (originally from the 4th field in the bed file).}
	\item{prob False/oligodT internally primed}{This is the probability that the putative pA site is false. Values range from 0-1, with 1 meaning the site is False/oligodT internally primed.}
	\item{prob True}{This is the probability that the putative pA site is true. Values range from 0-1, with 1 meaning the site is True.}
	\item{pred.class }{This is the predicted class of the putative pA site, based on the assignment cutoff. 0= Falsee/oligodT internally primed, 1 = True}
	\item{UpstreamSeq}{This is the upstream sequence of the putative pA site used in the analysis.}
	\item{DownstreamSeq}{This is the downstream sequence of the putative pA site used in the analysis.}
    
The function also return an invisible matrix including all info as decribed above.
}
\references{
Sarah Sheppard, Nathan D. Lawson, and Lihua Julie Zhu.  2013. Accurate identification of polyadenylation sites from 3' end deep sequencing using a na\"{i}ve Bayes classifier. Bioinformatics. Under revision 
}
\author{
Sarah Sheppard, Jianhong Ou, Nathan Lawson, Lihua Julie Zhu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
    testFile = system.file("extdata", "test.bed", package="cleanUpdTSeq")
    testSet = read.table(testFile, sep = "\t", header = TRUE)
		
	#convert the test set to GRanges without upstream and downstream sequence information
        peaks = BED2GRangesSeq(testSet,withSeq=FALSE)
        
	#build the feature vector for the test set without sequence information
	testSet.NaiveBayes = buildFeatureVector(peaks,BSgenomeName = Drerio, upstream = 40,
         downstream = 30, wordSize = 6, alphabet=c("ACGT"),
         sampleType = "unknown",replaceNAdistance = 30,
        method = "NaiveBayes", ZeroBasedIndex = 1, fetchSeq = TRUE)
        
    data(data.NaiveBayes)
    
    ## sample the test data for code testing, DO NOT do this for real data
    ## START SAMPLING
    samp <- c(1:22, sample(23:4119, 50), 4119, 4120)
    Ndata.NaiveBayes <- data.NaiveBayes$Negative[,samp]
    Pdata.NaiveBayes <- data.NaiveBayes$Positive[,samp]
    testSet.NaiveBayes@data <- testSet.NaiveBayes@data[, samp-1]
    ## END SAMPLING
    
	predictTestSet(Ndata.NaiveBayes, 
                   Pdata.NaiveBayes,
                   testSet.NaiveBayes,
	               outputFile="test-predNaiveBayes.xls", 
                   assignmentCutoff = 0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }