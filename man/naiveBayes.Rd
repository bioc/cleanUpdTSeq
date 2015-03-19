\name{naiveBayes-class}
\docType{class}
\alias{naiveBayes-class}
\alias{naiveBayes}
\alias{$,naiveBayes-method}
\alias{$<-,naiveBayes-method}

\title{Class \code{"naiveBayes"}}
\description{
  An object of class \code{"naiveBayes"} represents the conditional a-posterior 
  probabilities of a categorical class variable given independent predictor 
  variables using the Bayes rule.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("naiveBayes", apriori,
                          tables,
                          levels,
                          call)}.
}
\section{Slots}{
  \describe{
    \item{\code{apriori}}{Class distribution for the dependent variable.}
    \item{\code{tables}}{A list of tables, one for each predictor variable. 
    For each categorical variable a table giving, for each attribute level, 
    the conditional probabilities given the target class. For each numeric 
    variable, a table giving, for each target class, mean and standard 
    deviation of the (sub-)variable.}
  }
}
\section{Methods}{
    \describe{
        \item{$, $<-}{Get or set the slot of \code{\link{naiveBayes}}}
    }
}
\keyword{classes}
