\name{modpcomp}
\alias{modpcomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
modpcomp
}
\description{
extract module pca component
}
\usage{
modpcomp(data, colors, nPC = 2,
         plot = FALSE, filename = NULL, group = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  protein quantification data. column is sample. row is protein ID.
}
  \item{colors}{
  protein and module information. which is calculated in WGCNA package.
}
  \item{nPC}{
  how many PCA component will saved.
}
  \item{plot}{
  a logical value indicating whether draw PCA plot. This function need load \pkg{ggfortify} first.
}
  \item{filename}{
  The filename of plot. The default value is \code{NULL} which means no file saving.
  The plot will be saved to "plot" folder and saved in pdf format.
}
  \item{group}{
  sample group information.
}
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
Module_PCA <- modpcomp(logD, net$colors)

# if plot PCA and plot module 6 PCA
group <- gsub("[0-9]+", "", colnames(logD))
pos <- which(net$colors == 6)
if (requireNamespace("ggfortify", quietly = TRUE)){
require("ggfortify")
Module_PCA <- modpcomp(logD[pos,], net$colors[pos], plot = TRUE, group = group)
}

}

