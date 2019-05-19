\name{getmoduleHub}
\alias{getmoduleHub}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getmoduleHub
}
\description{
extract PMFG information and get Module hub proteins.
}
\usage{
getmoduleHub(data, module, mod_num, coln = "new.ID",
             cor.sig = 0.05, cor.r = 0,
             adjustp = TRUE, hub.p = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
proteomic quantification data.
}
  \item{module}{
module information which is getted in \code{Module_inf} function.
}
  \item{mod_num}{
the module name which module will be calculate.
}
  \item{coln}{
column name of module contains protein IDs. it could be matched with "\code{classifiedID}"
}
  \item{cor.sig}{
a numberic value indicated that correlation p value less than cor.sig will be picked.
}
  \item{cor.r}{
a numberic value indicated that correlation r value larger than cor.r will be picked.
}
  \item{adjustp}{
a logical value indicating whether pick hub protein by FDR methods.
}
  \item{hub.p}{
a numberic value indicated that hub proteins are p value less than hub.p.
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
logD <- data$log2_value
rownames(logD) <- data$inf$ori.ID
group <- gsub("[0-9]+","", colnames(logD))
Module <- Module_inf(net, data$inf)
Mod10 <- getmoduleHub(logD, Module, 10, coln = "ori.ID", adjustp = FALSE)
if (requireNamespace("MEGENA", quietly = TRUE)) {
 library(MEGENA)
 plot_subgraph(module = Mod10$degreeStat$gene,
              hub = Mod10$hub,PFN = Mod10$PMFG,
              node.default.color = "black",
              gene.set = NULL,color.code = c("grey"),show.legend = TRUE,
              label.hubs.only = TRUE,hubLabel.col = "red",hubLabel.sizeProp = 0.5,
              show.topn.hubs = 10,node.sizeProp = 13,label.sizeProp = 13,
              label.scaleFactor = 10,layout = "kamada.kawai")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

