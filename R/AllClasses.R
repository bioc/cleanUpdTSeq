setClass("modelInfo", representation(upstream="integer", 
                                     downstream="integer", 
                                     wordSize="integer",
                                     alphabet="character",
                                     genome="BSgenome"))

setClass("naiveBayes", representation(apriori="table",
                                      tables="list",
                                      levels="NULL",
                                      call="call"))

setClass("PASclassifier", representation(classifier="naiveBayes", 
                                         info="modelInfo"))

setClass("featureVector", representation(data="data.frame", 
                                         info="modelInfo"))