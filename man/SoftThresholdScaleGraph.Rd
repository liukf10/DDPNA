\name{SoftThresholdScaleGraph}
\alias{SoftThresholdScaleGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SoftThresholdScaleGraph
}
\description{
pick soft thresholding powers for WGCNA analysis and plot
}
\usage{
SoftThresholdScaleGraph(data,
                        xlab = "Soft Threshold (power)",
                        ylab = "Scale Free Topology Model Fit, signed R^2",
                        main = "Scale independence",
                        filename = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  protein quantification data. row is sample. column is protein ID.
}
  \item{xlab}{
  x axis label
}
  \item{ylab}{
  y axis label
}
  \item{main}{
  plot title
}
  \item{filename}{
  the filename of plot. The default value is \code{NULL} which means no file saving.
  The plot will be saved to "plot" folder and saved in pdf format.
  }
}
\details{
pick soft thresholding powers for WGCNA analysis and plot. The function is also can replaced by "\code{pickSoftThreshold}" function in WGCNA package.
}
\value{
A list with the following components:
  \item{powerEstimate }{the lowest power fit for scale free topology.}
  \item{fitIndices }{a data frame containing the fit indices for scale free topology. }
%% ...
}

\author{
Kefu Liu
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{pickSoftThreshold} in WGCNA package.
}
\examples{
rm(list = ls())
library(DDPNA)
data(Dforimpute)
data <- Data_impute(Dforimpute, miss.value = 0, distmethod = "manhattan")
logD <- data$log2_value
rownames(logD) <- data$inf$ori.ID
sft <- SoftThresholdScaleGraph(t(logD))
}


