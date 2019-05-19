\name{single_mod_enrichplot}
\alias{single_mod_enrichplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
single_mod_enrichplot
}
\description{
FCS enrichment analysis of a sets of proteins in one module.
}
\usage{
single_mod_enrichplot(module, Mod_Nam, classifiedID,
                      coln = "new.ID", datainf = NULL,
                      plot = TRUE, filename = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{module}{
module information which is getted in \code{Module_inf} function.
}
  \item{Mod_Nam}{
the module name which module will be calculate.
}
  \item{classifiedID}{
a sets of protein IDs which is ordered by change value/ p value and so on.
}
  \item{coln}{
column name of module contains protein IDs. it could be matched with "\code{classifiedID}"
}
  \item{datainf}{
proteomic data protein ID information.
The default value is "\code{NULL}". which is means that the "\code{classifiedID}" come from proteomic information is the same as the module construction proteomic information.
If they are different, proteomic data information should be given.
}
  \item{plot}{
a logical value indicating whether draw enrichment variation trend plot.
}
  \item{filename}{the filename of plot. The default value is \code{NULL} which means no file saving.
  The plot will be saved to "plot" folder and saved in pdf format. }
  \item{\dots}{
Other arguments.
}
}

\author{
Kefu Liu
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
rm(list = ls())
library(DDPNA)
data(net)
data(imputedData)
data <- imputedData
logD <- data$log2_value
rownames(logD) <- data$inf$ori.ID
group <- gsub("[0-9]+","", colnames(logD))
Module <- Module_inf(net, data$inf)
up <- changedID(logD, group, vs.set2 = "ad",vs.set1 = "ctl",
                rank = "foldchange",anova = FALSE, Padj = "none", cutoff = 1,
                datatype = "log2", fctype = "up")
m5enrich <- single_mod_enrichplot(Module, 5, up, coln="ori.ID")

}

