setClass("modelInfo", representation(upstream="integer", 
                                     downstream="integer", 
                                     wordSize="integer",
                                     alphabet="character",
                                     genome="BSgenome"))

setMethod("$", "modelInfo", function(x, name) slot(x, name))
setReplaceMethod("$", "modelInfo", 
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })

setClass("naiveBayes", representation(apriori="table",
                                      tables="list",
                                      levels="NULL",
                                      call="call"))

setMethod("$", "naiveBayes", function(x, name) slot(x, name))
setReplaceMethod("$", "naiveBayes", 
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })

setClass("PASclassifier", representation(classifier="naiveBayes", 
                                         info="modelInfo"))

setMethod("$", "PASclassifier", function(x, name) slot(x, name))
setReplaceMethod("$", "PASclassifier", 
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })

setClass("featureVector", representation(data="data.frame", 
                                         info="modelInfo"))

setMethod("$", "featureVector", function(x, name) slot(x, name))
setReplaceMethod("$", "featureVector", 
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })