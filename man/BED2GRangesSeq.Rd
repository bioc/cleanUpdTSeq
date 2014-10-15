\name{BED2GRangesSeq}
\alias{BED2GRangesSeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BED2GRangesSeq
}
\description{
This function converts an object of data.frame from a bed file with sequence information to GRanges with sequence information. 
}
\usage{
BED2GRangesSeq(data.BED, header = FALSE, 
               upstream.seq.ind = 7, downstream.seq.ind = 8, 
               withSeq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.BED}{
An object of data.frame from a bed file. The data.frame should at least contains 3 required fields: chrome, chromStart, chromend. The fourth field for "name" is suggested for keeping track of the putative polyadenylation site from the input to the output. The sixth field for "strand" is suggested, as this will affect the classification. For this function, the bed data.frame may also contain two additional fields for the sequence upstream and downstream of the putative pA site. If these are not supplied, the sequence may be obtained when the feature vector is built. Please see http://genome.ucsc.edu/FAQ/FAQformat.html#format1 for more information about the bed file format.
}
  \item{header}{
header = Boolean
TRUE if the first row is the header
FALSE if the first row is data
}
  \item{upstream.seq.ind}{
upstream.seq.ind = int, to delineate the column location containing the sequence upstream of the putative pA site
}
  \item{downstream.seq.ind}{
downstream.seq.ind = int, to delineate the column location containing the sequence downstream of the putative pA site
}
 \item{withSeq}{TRUE indicates that there are upstream and downstream sequences in the file, FALSE indicates that there is no upstream or downstream sequences in the file}
}

\value{
Returns object of GRanges
}

\author{
Sarah Sheppard, Jianhong Ou, Nathan Lawson, Lihua Julie Zhu
}

\examples{
testFile <- system.file("extdata", "test.bed", package="cleanUpdTSeq")
testSet <- read.table(testFile, sep = "\t", header = TRUE)
peaks <- BED2GRangesSeq(testSet,withSeq=TRUE)
}

\keyword{ misc }