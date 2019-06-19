\name{Module_Enrich}
\alias{Module_Enrich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Module_Enrich
}
\description{
Enrichment analysis of a sets of proteins in all modules. The function offered two enrichment methods:ORA and FCS.
}
\usage{
Module_Enrich(module, classifiedID, enrichtype = "FCS",
              coln = "new.ID", datainf = NULL, p.adj.method = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{module}{
module information which is getted in \code{Module_inf} function.
}
  \item{classifiedID}{
a sets of protein IDs which is ordered by change value/ p value and so on.
}
  \item{enrichtype}{
enrichment method.
This must be (an abbreviation of) one of the strings "\code{FCS}","\code{ORA}".
"\code{FCS}" means  analyzes step-by-step a proteins list which is ordered by change ratio/ p value and so on.
"\code{ORA}" means  analyzes a proteins list by Fisher's extract test.
}
  \item{coln}{
column name of module contains protein IDs. it could be matched with "\code{classifiedID}"
}
  \item{datainf}{
proteomic data protein ID information.
The default value is "\code{NULL}". which is means that the "\code{classifiedID}" come from proteomic information is the same as the module construction proteomic information.
If they are different, proteomic data information should be given.
}
  \item{p.adj.method}{
p adjust methods of multiple comparisons.
it can seen in \code{p.adjust.methods}.
}
}

\value{
a list contains classifiedID enrichment information.
  \item{Counts }{the counts of classifiedID in module.}
  \item{module.size }{the number of module ID}
  \item{module.name }{module name}
  \item{precent }{counts divided module.size }
  \item{p }{enrichment p value in each module}
  \item{p.adj }{enrichment p.adj value in each module}
  \item{Z.score }{Z score is -log2 P value.}
%% ...
}

\author{
Kefu Liu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(net)
data(imputedData)
data <- imputedData
logD <- data$log2_value
rownames(logD) <- data$inf$ori.ID
group <- gsub("[0-9]+","", colnames(logD))
Module <- Module_inf(net, data$inf)
up <- changedID(logD, group, vs.set2 = "ad",vs.set1 = "ctl",
              rank = "foldchange",anova = FALSE, Padj = "none",cutoff = 1,
              datatype = "log2",fctype = "up")
FCSenrich <- Module_Enrich(Module, up, coln="ori.ID")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

