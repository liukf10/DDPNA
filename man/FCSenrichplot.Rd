\name{FCSenrichplot}
\alias{FCSenrichplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FCSenrichplot
}
\description{
plot of FCS enrichment analysis
}
\usage{
FCSenrichplot(FCSenrich, count = 1, p = 0.05, filter = "p",
              plot = TRUE, filename = NULL,filetype = "pdf", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FCSenrich}{
FCS enrichment information which is getted in \code{module_enrich} function.
}
  \item{count}{
a numeric value. Module will choosed when countnumber is larger than count value .
}
  \item{p}{
a numeric value. Module will choosed when any Fisher's extract test p value is less than count value .
}
  \item{filter}{
filter methods.
This must be (an abbreviation of) one of the strings "\code{p}","\code{p.adj}", "\code{none}".
}
  \item{plot}{
a logical value indicating whether draw enrichment variation trend plot.
}
 \item{filename}{
the filename of plot. The default value is \code{NULL} which means no file saving.
The plot will be saved to "plot" folder.
}
 \item{filetype}{
the file type of plot.  the type should be one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).
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
data(net)
data <- Data_impute(Dforimpute, miss.value=0, distmethod="manhattan")
logD <- data$log2_value
rownames(logD) <- data$inf$ori.ID
group <- gsub("[0-9]+","", colnames(logD))
Module <- Module_inf(net, data$inf)
pos<-which(Module$moduleNum \%in\% c(11:22))
up <- changedID(logD[-pos,], group, vs.set2 = "ad",vs.set1 = "ctl",
              rank = "foldchange",anova = FALSE, Padj = "none",cutoff = 1,
              datatype = "log2",fctype = "up")
FCSenrich <- Module_Enrich(Module[-pos,], up, coln="ori.ID")
FCSenrich <- FCSenrichplot(FCSenrich)

}
