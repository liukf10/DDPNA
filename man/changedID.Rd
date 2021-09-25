\name{changedID}
\alias{changedID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
changedID
}
\description{
extract significant differential protein
}
\usage{
changedID(relative_value, group, vs.set2, vs.set1 = "WT",
          rank = "none", anova = TRUE, anova.cutoff = 0.05,
          T.cutoff = 0.05, Padj = "fdr",
          cutoff = 1.5, datatype = c("none","log2"), fctype = "all",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{relative_value}{protein quantification data}
  \item{group}{sample group information}
  \item{vs.set2}{compared group 2 name}
  \item{vs.set1}{compared group 1 name}
  \item{rank}{order by which type. This must be (an abbreviation of) one of the strings "\code{none}","\code{foldchange}", "\code{anova}","\code{t}"}
  \item{anova}{a logical value indicating whether do anova analysis.}
  \item{anova.cutoff}{a numberic value indicated that anova test p value upper limit.}
  \item{T.cutoff}{a numberic value indicated that t.test p value upper limit.}
  \item{Padj}{
  p adjust methods of multiple comparisons.
  it can seen in \code{p.adjust.methods}.
}
  \item{cutoff}{a numberic value indicated that foldchange lower limit.}
  \item{datatype}{The quantification data is normal data or log2 data.}
  \item{fctype}{foldchange is ordered by up-regulated or down-regulated or changed}
  \item{\dots}{
Other arguments.
}
  }
\details{
extract significant differential protein ID based on foldchange, t.test p value, anova p value.
}
\value{
a vector of protein ID information.
}
\author{
Kefu Liu
}



\examples{
data(imputedData)
data <- imputedData
logD <- data$log2_value
rownames(logD) <- data$inf$ori.ID
group <- gsub("[0-9]+","", colnames(logD))
up <- changedID(logD[201:260,], group, vs.set2 = "ad", vs.set1 = "ctl",
              rank = "foldchange",anova = FALSE, Padj = "none", cutoff = 1,
              datatype = "log2", fctype = "up")

}
