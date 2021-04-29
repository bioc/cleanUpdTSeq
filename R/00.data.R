#' Training Data
#' 
#' A RData containing negative and positive training data
#'
#' @name data.NaiveBayes
#' @docType data
#' @format A list with 2 data frame, "Negative" and "Positive". Negative has
#' 9219 observations on the following 4120 variables. And Positive is a data
#' frame with 22770 observations on the following 4120 variables.  The format
#' is: 
#' \describe{\item{list("Negative")}{'data.frame': 9219 obs. of 4120 variables:}
#'       \item{list("Positive")}{'data.frame': 22770 obs. of 4120 variables:}}
#' 
#' Both of them have same structure.
#' \describe{\item{list("y")}{a numeric vector} 
#'           \item{list("n.A.Downstream")}{a numeric vector}
#'           \item{list("n.C.Downstream")}{a numeric vector}
#'           \item{list("n.T.Downstream")}{a numeric vector}
#'           \item{list("n.G.Downstream")}{a numeric vector}
#'           \item{list("avg.distanceA2PeakEnd")}{a numeric vector}
#'           \item{list("dimer")}{a numeric vector}
#'           \item{: such as AA, AC, AG, AT, CA, ... etc.}{a numeric vector}
#'           \item{list("heximer")}{a factor with levels \code{0} \code{1}}
#'           \item{: such as AAAAAA, ACGTAC, ... etc.}{a factor with levels 
#'                \code{0} \code{1}} 
#'           \item{list("upstream.seq")}{a vector of sequence string} 
#'           \item{list("downstream.seq")}{a vector of sequence string}}
#' @keywords datasets
#' @examples
#' 
#' data(data.NaiveBayes)
#' head(str(data.NaiveBayes$Negative))
#' head(str(data.NaiveBayes$Positive))
"data.NaiveBayes"


#' NaiveBayes classifier
#' 
#' An object of class "naiveBayes" generated from data.NaiveBayes
#' 
#' @name classifier
#' @docType data
#' @format An object of class "\code{\link{PASclassifier}}" including
#' components:
#' @keywords datasets
#' @examples
#' 
#' data(classifier)
#' names(classifier)
#' 
"classifier"