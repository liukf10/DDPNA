\name{fc.pos}
\alias{fc.pos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fc.pos
}
\description{
Pick up proteins based on foldchange and return proteins position in data.
}
\usage{
fc.pos(fc, vs.set2, vs.set1 = "WT",
       cutoff = 1, datatype = c("none", "log2"),
       fctype = "all", order = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fc}{
proteomic data of mean value in groups.
}
  \item{vs.set2}{
compared group 2 name
}
  \item{vs.set1}{
compared group 1 name
}
  \item{cutoff}{
a numberic value indicated foldchange threshold.
}
  \item{datatype}{
The quantification data is normal data or log2 data.
This must be (an abbreviation of) one of the strings "\code{none}","\code{log2}".
}
  \item{fctype}{
foldchange is ordered by up-regulated or down-regulated or changed
}
  \item{order}{
a logical value indicated that whether ordered by foldchange.
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
relative <- data$relative_value
rownames(relative) <- data$inf$ori.ID
group <- gsub("[0-9]+", "", colnames(relative))
datamean <- groupmean(relative, group, name = FALSE)
fc_1vs2 <- fc.pos(datamean, vs.set2 = "ad", vs.set1 = "ctl",
                  cutoff = 1, datatype = "none",
                  fctype = "up", order = TRUE)
fc_ID <- rownames(relative)[fc_1vs2]


}

