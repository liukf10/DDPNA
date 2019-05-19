\name{groupmean}
\alias{groupmean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
groupmean
}
\description{
mean of sample group
}
\usage{
groupmean(data, group, method = c("mean", "median"), name = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  protein quantification data. column is sample. row is protein ID.
}
  \item{group}{
  sample group information
}
  \item{method}{
  Arithmetic mean of sample group or median of sample group.
  This must be (an abbreviation of) one of the strings "\code{mean}","\code{median}".
}
  \item{name}{
  a logical value indicated whether add "mean" or "median" in sample group name.
}
}


\author{
Kefu Liu
}



\examples{
rm(list = ls())
library(DDPNA)
data(imputedData)
data <- imputedData
logD <- data$log2_value
group <- gsub("[0-9]+","", colnames(logD))
datamean <- groupmean(logD, group, name = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the


