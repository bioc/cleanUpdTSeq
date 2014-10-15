\name{modelInfo-class}
\docType{class}
\alias{modelInfo-class}
\alias{modelInfo}

\title{Class \code{"modelInfo"}}
\description{
  An object of class \code{"modelInfo"} represents the information of sequence to use in the analysis
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("modelInfo", upstream, downstream, wordSize, alphabe, genome)}.
}
\section{Slots}{
  \describe{
    \item{genome}{
Name of the genome to use to get sequence. To find out a list of available genomes, please type available.genomes() in R.
}
  \item{upstream}{
This is the length of upstream sequence to use in the analysis.
}
  \item{downstream}{
This is the length of downstream sequence to use in the analysis.
}
  \item{wordSize}{
This is the size of the word to use as a feature for the upstream sequence. wordSize = 6 should always be used. }
  \item{alphabet}{
These are the bases to use, for example DNA bases ACTG.
}
  }
}
\keyword{classes}

