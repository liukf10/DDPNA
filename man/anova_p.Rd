\name{anova_p}
\alias{anova_p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
anova_p
}
\description{
anova analysis in proteomic data.
}
\usage{
anova_p(data, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
protein quantification data. column is sample. row is protein ID.
}
  \item{group}{
sample group information
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
anova_P <- anova_p(logD[1:100,], group)

}

