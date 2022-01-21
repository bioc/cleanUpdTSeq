#' This package classifies putative polyadenylation sites.
#'
#' 3'ends of transcripts have generally been poorly annotated. With the advent 
#' of deep sequencing, many methods have been developed to identify 3'ends. 
#' The majority of these methods use an oligodT primer which can bind to 
#' internal adenine-rich sequences, and lead to artifactual identification 
#' of polyadenylation sites. Heuristic filtering methods rely on a certain 
#' number of As downstream of a putative polyadenylation site to classify 
#' the site as true or oligodT primed. This package provides a robust method 
#' to classify putative polyadenylation sites using a Naive Bayes classifier.
#'
#' @docType package
#' @author Sarah Sheppard, Haibo Liu, Jianhong Ou, Nathan Lawson, Lihua Julie Zhu
#' @name cleanUpdTSeq
#' @aliases cleanUpdTSeq-package
globalVariables(c("Drerio"))
