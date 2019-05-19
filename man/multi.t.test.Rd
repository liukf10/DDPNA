\name{multi.t.test}
\alias{multi.t.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
multi.t.test
}
\description{
multiple comparisons t test and choose significant proteins in proteomic data.
}
\usage{
multi.t.test(data, group,
             sig = 0.05, Adj.sig = TRUE,
             grpAdj = "bonferroni",
             geneAdj = "fdr", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
protein quantification data. column is sample. row is protein ID.
}
  \item{group}{
sample group information
}
  \item{sig}{
siginificant P value threshold. The default is 0.05.
}
  \item{Adj.sig}{
a logical value indicated that whether adjust P-values for multiple proteins comparisons in each two groups.
}
  \item{grpAdj}{
adjust multiple groups comparisions P-value in each two groups. The default is \code{"bonferroni"}.
it can seen in \code{p.adjust.methods}.
}
  \item{geneAdj}{
adjust multiple proteins comparisions P-value in each group. The default is \code{"fdr"}.
it can seen in \code{p.adjust.methods}.
}
  \item{\dots}{
Other arguments.
}
}


\author{
Kefu Liu
}

\examples{
rm(list = ls())
library(DDPNA)
data(Dforimpute)
data <- Data_impute(Dforimpute, miss.value = 0, distmethod = "manhattan")
logD <- data$log2_value
rownames(logD) <- data$inf$ori.ID
group <- gsub("[0-9]+", "", colnames(logD))
Tsig_P <- multi.t.test(logD[1:100,], group, Adj.sig = FALSE, geneAdj = "fdr")

}

