\name{rename_dupnewID}
\alias{rename_dupnewID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rename_dupnewID
}
\description{
rename the duplicated newID in moduleinf and renew the ID in DEPstat
}
\usage{
rename_dupnewID(DEPstat, moduleinf, DEPfromMod = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DEPstat}{
a dataframe contains columns:"new.ID" and "ori.ID". it can get from \code{dataStatInf}.
}
  \item{moduleinf}{
a dataframe contains columns:"new.ID" and "ori.ID". it can get from \code{Module_inf}.
}
  \item{DEPfromMod}{
a logical value indicated that whether DEPstat and moduleinf is getted from the same datasets.
The default value is FALSE.
}
}

\value{
a data.frame contains DEPstat information and renewed the new.ID column.
}

\author{
Kefu Liu
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(net)
data(imputedData)
Module <- Module_inf(net, imputedData$inf)
group <- gsub("[0-9]+","", colnames(imputedData$intensity))
data <- imputedData
data$inf <- data$inf[1:100,]
data$intensity <- data$intensity[1:100,]
stat <- dataStatInf(data, group, meanmethod = "median",
                    T.test = "pairwise", Aadj = "fdr",
                    Tadj = "fdr", cutoff = FALSE)
stat <- rename_dupnewID(stat, Module, DEPfromMod = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the

