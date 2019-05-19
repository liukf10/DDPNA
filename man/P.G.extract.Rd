\name{P.G.extract}
\alias{P.G.extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Protein Groups information extract.
}
\description{
uniprot ID, ENTRYNAME and status information extract.(only fit for MaxQuant data.)
}
\usage{
P.G.extract(inf, ncol = 4, justID = FALSE,
            status1 = FALSE, ENTRY1 = FALSE, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inf}{
  protein groups IDs information.
}
  \item{ncol}{
  column numbers of output result.
}
  \item{justID}{
  a logical value indicated whether only extract uniprot ID.
}
  \item{status1}{
  a logical value indicated whether extract the first ID status.
}
  \item{ENTRY1}{
  a logical value indicated whether extract the first ID ENTRY NAME.
}
  \item{verbose}{
  integer level of verbosity. Zero means silent, 1 means have Diagnostic Messages.
}
}


\author{
Kefu Liu
}


\examples{
rm(list = ls())
library(DDPNA)
data(ProteomicData)
MaxQdata <- MaxQprotein(ProteomicData$MaxQ)
inf <- P.G.extract(MaxQdata$protein_IDs, justID = TRUE, status = TRUE, ENTRY = TRUE)

}

