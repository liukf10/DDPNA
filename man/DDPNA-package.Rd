\name{DDPNA-package}
\alias{DDPNA-package}
\alias{DDPNA}
\docType{package}
\title{
Disease-drived Differential Proteins And Proteomic Co-expression Network Associated Analysis
}
\description{
disease drived proteins associated network in different species crosstalk.
The package is used to analysis differential proteomics consensus network in two or more datasets.
}
\details{
\tabular{ll}{
Package: \tab DHMDA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2019-03-18\cr
License: \tab GPL (>= 2)\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Kefu Liu

Maintainer: Kefu Liu <lkf1013@gmail.com>
}


\keyword{ package }

\examples{
rm(list = ls())
library(DDPNA)
#database file should be downloaded and extracted the rar file in C Disk.
#blast+ software should be downloaded and installed in C Disk.
#if you want to do the test, just make docompletetest is TRUE.
docompletetest <- FALSE;
if(docompletetest){
 data(Dforimpute)
data <- Data_impute(Dforimpute, miss.value=0, distmethod="manhattan")
logD <- data$log2_value
rownames(logD) <- data$inf$ori.ID
sft <- SoftThresholdScaleGraph(t(logD))
if(requireNamespace("WGCNA", quietly = TRUE)){
  library(WGCNA)
  net <- blockwiseModules(t(logD), power = sft$powerEstimate, maxBlockSize = 6000,
                       TOMType = "unsigned", deepSplit = 4,minModuleSize = 17,
                       reassignThreshold = 0.05, mergeCutHeight = 0.07,
                       numericLabels = TRUE,pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,loadTOMs = FALSE,
                       saveTOMFileBase = "blockwisetom",
                       verbose = 3)
} else {
data(net)
}
MEs <- ME_inf(net$MEs, logD)
Module <- Module_inf(net, data$inf)
group <- gsub("[0-9]+","", colnames(logD))
datamean <- groupmean(logD, group, name = FALSE)
up <- changedID(logD, group, vs.set2 = "ad",vs.set1 = "ctl",
              rank = "foldchange",anova = FALSE, Padj = "none",cutoff = 1,
              datatype = "log2",fctype = "up")
m5enrich <- single_mod_enrichplot(Module, 5, up, coln="ori.ID")
ORAenrich <- Module_Enrich(Module, up, coln="ori.ID", enrichtype = "ORA")
rowname <- ORAenrich$module.name;
ORAenrich <- data.frame(Counts=ORAenrich$Counts,
                        mosule.size=ORAenrich$module.size,
                        precent=ORAenrich$precent,
                        p=ORAenrich$p,p.adj=ORAenrich$p.adj,
                        Z.score=ORAenrich$Z.score,
                        stringsAsFactors = F)
rownames(ORAenrich) <- rowname;
rm(rowname)
FCSenrich <- Module_Enrich(Module, up, coln="ori.ID")
FCSenrichplot(FCSenrich)

Mod10<-getmoduleHub(logD,Module,10,coln="ori.ID",adjustp=FALSE)
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

}
