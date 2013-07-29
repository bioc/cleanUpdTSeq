\name{getUpstreamSequence}
\alias{getUpstreamSequence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get upstream sequences of the putative pA site 
}
\description{
This function gets the sequence upstream of a putative pA site (including the site)
}
\usage{
getUpstreamSequence(peaks, upstream = 40, genome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{
GRanges containing putative pA sites 
}
  \item{upstream}{
upstream = int. This is the length of the sequence to get.
}
  \item{genome}{
BSgenomeName
}
}
\value{
Returns an object of GRanges with uptream sequences.
}
\author{
Sarah Sheppard, Jianhong Ou, Nathan Lawson, Lihua Julie Zhu
}

\examples{
        testFile <- system.file("extdata", "test.bed", package="cleanUpdTSeq")
        testSet <- read.table(testFile, sep="\t", header=TRUE)
        peaks <- BED2GRangesSeq(testSet[1:10, ], withSeq=FALSE)
        seq = getUpstreamSequence(peaks, upstream=40, genome=Drerio)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }