\name{DEP_Mod_net_plot}
\alias{DEP_Mod_net_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DEP_Mod_net_plot
}
\description{
remove hubs which is not in the IDsets and replot the PFG network
}
\usage{
DEP_Mod_net_plot(ModNet, IDsets = NULL, data = NULL, module = NULL,
                 plot = TRUE, filename = NULL, filetype = "pdf",
                 OnlyPlotLast = TRUE, BranchCut = TRUE,
                 reconstructNet = TRUE,
                 iteration = Inf, label.hubs.only = TRUE,
                 node.default.color = "grey",
                 hubLabel.col = "black", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ModNet}{
  data contains network information which get from \code{getmoduleHub}

}
  \item{IDsets}{
  ID sets information which get from \code{DEPsets}
%%     ~~Describe \code{IDsets} here~~
}
  \item{data}{
  the value should be defined only when \code{reconstructNet} is TRUE.
  The value is proteomic quantification data, which is same as the input in \code{getmoduleHub}.
}
  \item{module}{
  the value should be defined only when \code{reconstructNet} is TRUE.
  The value is module information which is same as the input in \code{getmoduleHub}.
}
  \item{plot}{
a logical value whether plot a picture.
}
  \item{filename}{
the filename of plot. The default value is \code{NULL} which means no file saving. The function is use \code{ggsave} to achieve.
}
  \item{filetype}{
the file type of plot.  the type should be one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).
}
  \item{OnlyPlotLast}{
a logical value whether plot the final network.
}
\item{BranchCut}{
a logical value whether remove unhub proteins which have no connection to DEPs.
}
  \item{reconstructNet}{
a logical value whether reconstruct network.
}
  \item{iteration}{
iteration times when reconstruct network.
}
  \item{label.hubs.only}{
a logical value whether show labels for hubs only.
}
  \item{node.default.color}{
Default node colors for those that do not intersect with signatures in gene.set.
}
  \item{hubLabel.col}{
Label color for hubs.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}

\value{
%%  ~Describe the value returned
a list contains network information
  \item{netgene }{all IDs in network.}
  \item{hub }{hub IDs}
  \item{PMFG }{PMFG graph data frame information}
%% ...
}

\author{
Kefu Liu%%  ~~who you are~~
}

\examples{
\donttest{
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
stat1 <- stat$ori.ID[stat$ad > 1]
stat2 <- stat$ori.ID[stat$asym > 1]
datalist <- list(stat1 = stat1, stat2 = stat2)
sets <- DEPsets(datalist)

logD <- imputedData$log2_value
rownames(logD) <- imputedData$inf$ori.ID
Mod3 <- getmoduleHub(logD, Module, 3, coln = "ori.ID", adjustp = FALSE)

newnet <- DEP_Mod_net_plot(Mod3, sets,
                           data = logD, module = Module,
                           plot = FALSE, filename = NULL, filetype = "pdf",
                           OnlyPlotLast = FALSE,reconstructNet = FALSE)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the

