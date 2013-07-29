\name{data.NaiveBayes}
\alias{data.NaiveBayes}
\docType{data}
\title{
Training Data
}
\description{
This is the negative and positive training data.
}
\usage{data(data.NaiveBayes)}
\format{
  A list with 2 data frame, "Negative" and "Positive". Negative has 9219 observations on the following 4120 variables. And Positive is a data frame with 22770 observations on the following 4120 variables. 
  The format is:
  \describe{
    \item{\code{Negative}}{'data.frame':    9219 obs. of  4120 variables:}
    \item{\code{Positive}}{'data.frame':   22770 obs. of  4120 variables:}
  }
  
  Both of them have same structure. 
  \describe{
        \item{\code{y}}{a numeric vector}
        \item{\code{n.A.Downstream}}{a numeric vector}
        \item{\code{n.C.Downstream}}{a numeric vector}
        \item{\code{n.T.Downstream}}{a numeric vector}
        \item{\code{n.G.Downstream}}{a numeric vector}
        \item{\code{avg.distanceA2PeakEnd}}{a numeric vector}
        \item{\code{dimer}: such as AA, AC, AG, AT, CA, ... etc.}{a numeric vector}
        \item{\code{heximer}: such as AAAAAA, ACGTAC, ... etc.}{a factor with levels \code{0} \code{1}}
        \item{\code{ upstream.seq }}{a vector of sequence string}
        \item{\code{ downstream.seq }}{a vector of sequence string}
    }
}
\examples{
data(data.NaiveBayes)
head(str(data.NaiveBayes$Negative))
head(str(data.NaiveBayes$Positive))
}
\keyword{datasets}