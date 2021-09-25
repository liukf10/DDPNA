\name{Data_impute}
\alias{Data_impute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data_impute
}
\description{
data clean process: detect and remove outlier sample and impute missing value.
The process is following:
1. Remove some genes which the number of missing value larger than maxNAratio.
2. Outlier sample detect and remove these sample.
3. Repeat Steps 1-2 untile meet the iteration times or no outlier sample can be detected.
4. impute the missing value.
The function also can only do gene filter or remove outlier or impute missing value.
}
\usage{
Data_impute(data, inf = "inf", intensity = "LFQ", miss.value = NA,
            splNExt = TRUE, maxNAratio = 0.5,
            removeOutlier = TRUE,
            outlierdata = "intensity", iteration = NA, sdout = 2,
            distmethod = "manhattan", A.IAC = FALSE,
            dohclust = FALSE, treelabels = NA,
            plot = TRUE, filename = NULL,
            text.cex = 0.7, text.col = "red", text.pos = 1,
            text.labels = NA, abline.col = "red", abline.lwd = 2,
            impute = TRUE, verbose = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  MaxQconvert data or a list Vector which contain two data.frame:ID information and quantification data
  }
  \item{inf}{
  the data.frame name contain protein ID information
  }
  \item{intensity}{
  the data.frame name only contain quantification data
  }
  \item{miss.value}{
  the type of miss.value showed in quantificaiton data.
  The default value is \code{NA}. The miss.value usually can be \code{NA} or  \code{0}.}
  \item{splNExt}{
  a logical value whether extract sample name.(suited for MaxQuant quantification data)
  }
  \item{maxNAratio}{
  The maximum percent missing data allowed in any row (default 50\%).For any rows with more than maxNAratio\% missing will deleted.
  }
  \item{removeOutlier}{
  a logical value indicated whether remove outlier sample.
  }
  \item{outlierdata}{
  The value is deprecated.
  which data will be used to analysis outlier sample detect.This must be (an abbreviation of) one of the strings "\code{intensity}","\code{relative_value}","\code{log2_value}".
  }
  \item{iteration}{
  a numberic value indicating how many times it go through the outlier sample detect and remove loop.\code{NA} means do loops until no outlier sample.
  }
  \item{sdout}{
  a numberic value indicating the threshold to judge the outlier sample. The default \code{2} means 0.95 confidence intervals
  }
  \item{distmethod}{
  The distance measure to be used. This must be (an abbreviation of) one of the strings "\code{manhattan}","\code{euclidean}", "\code{canberra}","\code{correlation}"
  }
  \item{A.IAC}{
  a logical value indicated whether decreasing \code{correlation} variance.
  }
  \item{dohclust}{
  a logical value indicated whether doing hierarchical clustering and plot dendrograms.
  }
  \item{treelabels}{
  labels of dendrograms
  }
  \item{plot}{
  a logical value indicated whether plot numbersd scatter diagrams.
  }
  \item{filename}{
  the filename of plot. The number and plot type information will added automatically. The default value is \code{NULL} which means no file saving.
  all the plot will be saved to "plot" folder and saved in pdf format.
  }
  \item{text.cex}{
  outlier sample annotation text size(scatter diagrams parameters)
  }
  \item{text.col}{
  outlier sample annotation color(scatter diagrams parameters)
  }
  \item{text.pos}{
  outlier sample annotation position(scatter diagrams parameters)
  }
  \item{text.labels}{
  outlier sample annotation (scatter diagrams parameters)
  }
  \item{abline.col}{
  the threshold line color (scatter diagrams parameters)}
  \item{abline.lwd}{the threshold line width (scatter diagrams parameters)
  }
  \item{impute}{
  a logical value indicated whether do knn imputation.
  }
  \item{verbose}{
  integer level of verbosity. Zero means silent, 1 means have some Diagnostic Messages.
}
  \item{\dots}{
  Other arguments.
  }
}
\details{
detect and remove outlier sample and impute missing value.

}
\value{
a list of proteomic data.
  \item{inf}{Portein information included protein IDs and other information.}
  \item{intensity}{Quantification informaton.}
  \item{relative_value}{intensity divided by geometric mean}
  \item{log2_value}{log2 of relative_value}
}

\author{
Kefu Liu
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(Dforimpute)
data <- Data_impute(Dforimpute,distmethod="manhattan")
}
