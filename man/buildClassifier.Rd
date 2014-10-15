\name{buildClassifier}
\alias{buildClassifier}
\title{
get Naive Bayes Classifier
}
\description{
Computes the conditional a-posterior probabilities of a categorical class variable given independent predictor variables using the Bayes rule.
}
\usage{
buildClassifier(Ndata.NaiveBayes, Pdata.NaiveBayes, 
                upstream=40L, downstream=30L, wordSize=6L,
                genome=Drerio, alphabet=c("ACGT"))
}
\arguments{
  \item{Ndata.NaiveBayes}{
This is the negative training data, described further in \code{\link[cleanUpdTSeq]{data.NaiveBayes}}.
}
  \item{Pdata.NaiveBayes}{
This is the positive training data, described further in \code{\link[cleanUpdTSeq]{data.NaiveBayes}}.
}
  \item{upstream}{
This is the length of upstream sequence to use in the analysis.
}
  \item{downstream}{
This is the length of downstream sequence to use in the analysis.
}
  \item{wordSize}{
This is the size of the word to use as a feature for the upstream sequence. wordSize = 6 should always be used. 
}
  \item{genome}{
Name of the genome to use to get sequence. To find out a list of available genomes, please type available.genomes() in R.
}
  \item{alphabet}{
These are the bases to use, for example DNA bases ACTG.
}
}
\value{
An object of class "naiveBayes".
}
\author{
Jianhong Ou
}
\seealso{
\code{\link[e1071]{naiveBayes}}
}
\examples{
if (interactive()){
    data(data.NaiveBayes)
    classifier <- buildClassifier(data.NaiveBayes$Negative, data.NaiveBayes$Positive)
}
}
\keyword{ misc }
