\name{moduleID}
\alias{moduleID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
extract intersection ID between dataset and module
}
\description{
extract intersection ID between dataset and one of module
}
\usage{
moduleID(inf, module, num, coln = "new.ID")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inf}{
dataset protein ID information. a vector of protein IDs.
}
  \item{module}{
module information which is getted in \code{Module_inf} function.
}
  \item{num}{
module number which will extract to compared with dataset ID information.
}
  \item{coln}{
column names of module protein IDs.
}
}
\details{
 column \code{coln} information in \code{module} when module number is \code{num} intersect with inf.
}


\author{
Kefu Liu
}



\examples{
data(net)
data(imputedData)
data <- imputedData
logD <- data$log2_value
rownames(logD) <- data$inf$ori.ID
group <- gsub("[0-9]+","", colnames(logD))
Module <- Module_inf(net, data$inf)
up <- changedID(logD, group, vs.set2 = "ad",vs.set1 = "ctl",
              rank = "foldchange",anova = FALSE, Padj = "none",cutoff = 1,
              datatype = "log2",fctype = "up")
intersection <- moduleID(up, Module, 5, coln = "ori.ID")
}

