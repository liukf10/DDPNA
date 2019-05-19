\name{ME_inf}
\alias{ME_inf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
module eigengenes information
}
\description{
put sample names as rownames in WGCNA module eigenvalue data.frame.
}
\usage{
ME_inf(MEs, data, intensity.type = "LFQ", rowname = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MEs}{
  module eigenvalue which is calculated in WGCNA package.
}
  \item{data}{
  protein quantification data. column is sample. row is protein ID.
}
  \item{intensity.type}{
  quantification data type, which can help extract sample name.
  This must be (an abbreviation of) one of the strings "\code{LFQ}","\code{intensity}","\code{iBAQ}","\code{none}".
}
  \item{rowname}{
  sample names when "\code{intensity.type}" is "\code{none}", rowname will be used.
}
}


\author{
Kefu Liu
}

\examples{
rm(list = ls())
library(DDPNA)
data(net)
data(imputedData)
data <- imputedData
logD <- data$log2_value
MEs <- ME_inf(net$MEs, logD)
}

