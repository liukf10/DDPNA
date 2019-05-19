\name{Module_inf}
\alias{Module_inf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Module and protein information.
}
\description{
module and protein information match
}
\usage{
Module_inf(net, inf, inftype = "Convert", IDname = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
  module network which is calculated in WGCNA package.
}
  \item{inf}{
  proteome quantification data information which contains protein IDs.
}
  \item{inftype}{
  data information type.
  This must be (an abbreviation of) one of the strings "\code{Convert}","\code{MaxQ}","\code{none}".
  "\code{Convert}" means protein ID is converted by \code{MaxquantDataconvert} function.
  "\code{MaxQ}" means original Maxquant software quantification data information.

}
  \item{IDname}{
  IDname is "\code{inf}" column names which will extract.
}
  \item{\dots}{
  other argument.
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
data <- Data_impute(Dforimpute, miss.value = 0, distmethod = "manhattan")
Module <- Module_inf(net, data$inf)
}

