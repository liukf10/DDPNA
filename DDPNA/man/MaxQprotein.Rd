\name{MaxQprotein}
\alias{MaxQprotein}

\title{
read proteomic quantification data and seperate the protein information and quantification information.
}
\description{
The function will seperate data into 4 parts: protein information, intensity, iBAQ and LFQ (iBAQ and LFQ only fit for 'MaxQuant' software result). For MaxQ data, it can remove the contaminant and reverse protein.
}
\usage{
MaxQprotein(proteinGroups, IDname = "Majority.protein.IDs",
            IDtype = "MaxQ", remove = TRUE, QuanCol = NULL,
            verbose = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{proteinGroups}{
  the proteomic quantification data
}
  \item{IDname}{
  The column name of uniprot ID. The default value is "\code{Majority.protein.IDs}"" which is the column name in MaxQuant data.
}
  \item{IDtype}{
  "\code{MaxQ}" means \code{proteinGroups} is Maxquant quantification data, "\code{none}" means other type data. This must be (an abbreviation of) one of the strings: "\code{MaxQ}","\code{none}".
}
  \item{remove}{
  a logical value indicated whether remove contaminant and reverse ID.
}
  \item{QuanCol}{
  The quantification data columns. It's only needed when \code{IDtype} is "\code{none}".
  When \code{IDtype} is "\code{none}" and \code{QuanCol} is not given, the intensity will auto extract all columns except \code{IDname} as quantification data. It may have error in next analysis.
}
  \item{verbose}{
  integer level of verbosity. Zero means silent, 1 means have Diagnostic Messages.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
a list of proteomic information.
  \item{protein_IDs}{Portein IDs which is \code{IDname} column information.}
  \item{intensity}{Quantification intensity informaton. When  \code{IDtype} is "\code{none}", it is the  \code{QuanCol} columns information.}
  \item{iBAQ}{Quantification iBAQ intensity informaton.(only for \code{IDtype} is "\code{MaxQ}")}
  \item{LFQ}{Quantification LFQ intensity informaton.(only for \code{IDtype} is "\code{MaxQ}")}
%% ...
}

\author{
Kefu Liu
}

\examples{
data(ProteomicData)
# example for MaxQ Data
MaxQdata <- MaxQprotein(ProteomicData$MaxQ)
# example for other type Data
otherdata <- MaxQprotein(ProteomicData$none, IDname = "Protein",
                         IDtype = "none", QuanCol = 2:9)

}

