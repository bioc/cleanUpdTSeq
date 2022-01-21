#' Class \code{"modelInfo"}
#' 
#' An object of class \code{"modelInfo"} represents the information of sequence
#' to use in the analysis
#' 
#' @name modelInfo-class
#' @aliases modelInfo-class modelInfo $,modelInfo-method $<-,modelInfo-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("modelInfo", upstream, downstream, wordSize, alphabe, genome, metadata)}.
#' @importClassesFrom BSgenome BSgenome
#' @exportClass modelInfo PASclassifier featureVector
#' @keywords classes
#' 
setClass("modelInfo", slots = c(upstream="integer",
                                downstream="integer", 
                                wordSize="integer",
                                alphabet="character",
                                genome="character",
                                metadata="list"))


setMethod("$", "modelInfo", function(x, name) slot(x, name))
setReplaceMethod("$", "modelInfo", 
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x})

#' Class \code{"naiveBayes"}
#' 
#' An object of class \code{"naiveBayes"} represents the conditional
#' a-posterior probabilities of a categorical class variable given independent
#' predictor variables using the Bayes rule.
#' 
#' 
#' @name naiveBayes-class
#' @aliases naiveBayes-class naiveBayes $,naiveBayes-method
#' $<-,naiveBayes-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("naiveBayes", apriori, tables, levels, call)}.
#' @keywords classes
#' @import methods
#' @exportMethod $ $<-
setClass("naiveBayes", slots = c(apriori="table",
                                 tables="list",
                                 levels="NULL",
                                 call="call"))


setMethod("$", "naiveBayes", function(x, name) slot(x, name))
setReplaceMethod("$", "naiveBayes", 
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x})

#' Class \code{"PASclassifier"}
#' 
#' An object of class \code{"PASclassifier"} represents the output of
#' \code{\link{buildClassifier}}
#' 
#' 
#' @name PASclassifier-class
#' @aliases PASclassifier-class PASclassifier $,PASclassifier-method
#' $<-,PASclassifier-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("PASclassifier", classifier, info)}.
#' @keywords classes
#' @examples
#' 
#' data(classifier)
#' classifier$info$upstream
#' classifier$info$wordSize
#' classifier$info$alphabet
#' 
setClass("PASclassifier", slots = c(classifier="naiveBayes", 
                                    info="modelInfo"))

setMethod("$", "PASclassifier", function(x, name) slot(x, name))
setReplaceMethod("$", "PASclassifier", 
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x})

#' Class \code{"featureVector"}
#' 
#' An object of class \code{"featureVector"} represents the output of
#' \code{\link{buildFeatureVector}}
#' 
#' 
#' @name featureVector-class
#' @aliases featureVector-class featureVector $,featureVector-method
#' $<-,featureVector-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("featureVector", data, info)}.
#' @keywords classes
setClass("featureVector", slots = c(data="data.frame", 
                                    info="modelInfo"))

setMethod("$", "featureVector", function(x, name) slot(x, name))
setReplaceMethod("$", "featureVector", 
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x})