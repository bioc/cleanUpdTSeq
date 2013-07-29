\name{getDownstreamSequence}
\alias{getDownstreamSequence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getDownstreamSequence}
\description{
This function gets the sequence upstream of a putative pA site (including the site)
}
\usage{
getDownstreamSequence(peaks, downstream = 20, genome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{
GRanges containing putative pA sites 
}
  \item{downstream}{
downstream = int. This is the length of the sequence to get.
}
  \item{genome}{
BSgenomeName
}
}
\value{
Returns an object of GRanges with downtream sequences.
}
\author{
Sarah Sheppard, Jianhong Ou, Nathan Lawson, Lihua Julie Zhu
}
\examples{
        testFile <- system.file("extdata", "test.bed", package="cleanUpdTSeq")
        testSet <- read.table(testFile, sep="\t", header=TRUE)
        peaks <- BED2GRangesSeq(testSet[1:10, ], withSeq=FALSE)
        seq = getDownstreamSequence(peaks, downstream=30, genome=Drerio)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
