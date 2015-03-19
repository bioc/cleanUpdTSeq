\name{featureVector-class}
\docType{class}
\alias{featureVector-class}
\alias{featureVector}
\alias{$,featureVector-method}
\alias{$<-,featureVector-method}

\title{Class \code{"featureVector"}}
\description{
  An object of class \code{"featureVector"} represents the output of \code{\link{buildFeatureVector}}
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("featureVector", data, info)}.
}
\section{Slots}{
  \describe{
    \item{\code{data}}{
      An object of data frame. Fields include peak name, upstream sequence, downstream sequence, and features to be used in classifying the putative polyadenylation site
    }
    \item{info}{ Object of class \code{\link{modelInfo}}}
  }
}
\section{Methods}{
    \describe{
        \item{$, $<-}{Get or set the slot of \code{\link{featureVector}}}
    }
}
\keyword{classes}

