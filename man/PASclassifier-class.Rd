\name{PASclassifier-class}
\docType{class}
\alias{PASclassifier-class}
\alias{PASclassifier}
\alias{$,PASclassifier-method}
\alias{$<-,PASclassifier-method}

\title{Class \code{"PASclassifier"}}
\description{
  An object of class \code{"PASclassifier"} represents the output of \code{\link{buildClassifier}}
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("PASclassifier", classifier, info)}.
}
\section{Slots}{
  \describe{
    \item{\code{classifier}}{Object of class \code{"naiveBayes"} The output of \code{\link[e1071]{naiveBayes}}.
      An object of class "naiveBayes" including components:
        \describe{
            \item{\code{apriori}}{Class distribution for the dependent variable.}
            \item{\code{tables}}{A list of tables, one for each predictor variable. For each categorical variable a table giving, for each attribute level, the conditional probabilities given the target class. For each numeric variable, a table giving, for each target class, mean and standard deviation of the (sub-)variable.}
        }
    }
    \item{info}{ Object of class \code{\link{modelInfo}}}
  }
}
\section{Methods}{
    \describe{
        \item{$, $<-}{Get or set the slot of \code{\link{PASclassifier}}}
    }
}
\examples{
data(classifier)
classifier$info$upstream
classifier$info$wordSize
classifier$info$alphabet
}
\keyword{classes}
