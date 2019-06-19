\name{dataStatInf}
\alias{dataStatInf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
dataStatInf
}
\description{
summrize the statistics information of data
}
\usage{
dataStatInf(prodata, group, intensity = "intensity",
            Egrp = NULL, Cgrp = "ctl",
            meanmethod = "mean", datatype = c("none", "log2"),
            anova = TRUE, T.test = c("pairwise", "student", "none"),
            Aadj = "none", Tadj = "none", cutoff = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prodata}{
  proteome data. a list Vector which contain two data.frame: ID information and quantification data
}
  \item{intensity}{
  the data.frame name only contain quantification data
%%     ~~Describe \code{intensity} here~~
}
  \item{group}{
  sample group information
}
  \item{Egrp}{
  experiment group name. It must be assigned when use Student T.test.
%%     ~~Describe \code{Egrp} here~~
}
  \item{Cgrp}{
  control group name. It must be assigned. The default value is "ctl".
}
  \item{meanmethod}{
  Arithmetic mean of sample group or median of sample group.
  This must be (an abbreviation of) one of the strings "\code{mean}","\code{median}".
}
  \item{datatype}{
  The quantification data is normal data or log2 data.
%%     ~~Describe \code{datatype} here~~
}
  \item{anova}{
  a logical value indicating whether do anova analysis.
%%     ~~Describe \code{anova} here~~
}
  \item{T.test}{
  T.test method. "none" means not running t.test.
  "pairwise" means calculate pairwise comparisons between group levels with corrections for multiple testing
  "student" means student t test.
  This must be (an abbreviation of) one of the strings "\code{pairwise}","\code{student} and "\code{none}"".
}
  \item{Aadj}{
  anova P value adjust methods. it can seen in \code{p.adjust.methods}.
%%     ~~Describe \code{Aadj} here~~
}
  \item{Tadj}{
  T test P value adjust methods. it can seen in \code{p.adjust.methods}.
%%     ~~Describe \code{Tadj} here~~
}
  \item{cutoff}{
  a logical value or a numeric value.
  The default value is FALSE, which means do not remove any P value.
  If the value is TRUE, P value > 0.05 will remove and showed as NA in result.
  If the value is numeric, P value > the number will remove and showed as NA in result.
}
  \item{\dots}{
Other arguments.
}
}

\value{
a data.frame of protein ID and Statistics information.

}

\author{
Kefu Liu
}

\examples{
data(imputedData)
group <- gsub("[0-9]+","", colnames(imputedData$intensity))
data <- imputedData
data$inf <- data$inf[1:100,]
data$intensity <- data$intensity[1:100,]
stat <- dataStatInf(data, group, meanmethod = "median",
                    T.test = "pairwise", Aadj = "fdr",
                    Tadj = "fdr", cutoff = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

