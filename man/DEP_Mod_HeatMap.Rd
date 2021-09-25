\name{DEP_Mod_HeatMap}
\alias{DEP_Mod_HeatMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DEP_Mod_HeatMap
}
\description{
 get the DEP enrich fold in Module and plot a HeatMap
}
\usage{
DEP_Mod_HeatMap(DEP_Mod, xlab = "DEP", filter = c("p","p.adj"),
                cutoff = 0.05, filename = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DEP_Mod}{
  a list of DEP_Mod enrichment information. data.frame in list is get from \code{Module_Enrich} function.
}
\item{xlab}{
  it indicate x value in heatmap. it must be a value between "\code{DEP}" and "\code{MOD}".
}
  \item{filter}{
  p value or p.adjust value used to filter the enrich significant module.
%%     ~~Describe \code{filter} here~~
}
  \item{cutoff}{
  a numeric value is the cutoff of p value. Larger than the value will remove to show in plot.
%%     ~~Describe \code{cutoff} here~~
}
  \item{filename}{
  plot filename. If filename is null, it will print the plot.
%%     ~~Describe \code{filename} here~~
}
  \item{\dots}{
other argument.
}
}

\value{
a list of enrich fold heatmap information.
  \item{enrichFold}{enrichFold of DEP in Modules.}
  \item{textMatrix}{siginificant enrichment module information.}

%% ...
}

\author{
Kefu Liu
}



\examples{
\donttest{
data(net)
data(imputedData)
data <- imputedData
logD <- data$log2_value
rownames(logD) <- data$inf$ori.ID
group <- gsub("[0-9]+","", colnames(logD))
Module <- Module_inf(net, data$inf)
# define 2 DEP ID data: a and b
a <- Module$ori.ID[1:100]
b <- Module$ori.ID[50:100]
a <- Module_Enrich(Module, a, coln="ori.ID", enrichtype = "ORA")
b <- Module_Enrich(Module, b, coln="ori.ID", enrichtype = "ORA")
rowname <- a$module.name;
a <- data.frame(Counts = a$Counts, module.size = a$module.size,
                precent = a$precent, p = a$p, p.adj = a$p.adj,
                Z.score = a$Z.score, stringsAsFactors = FALSE)
rownames(a) <- rowname;
rowname <- b$module.name;
b <- data.frame(Counts = b$Counts, module.size = b$module.size,
                precent = b$precent, p = b$p, p.adj = b$p.adj,
                Z.score = b$Z.score, stringsAsFactors = FALSE)
rownames(b) <- rowname;
DEP_Mod <- list(a = a , b = b)
heatMapInf <- DEP_Mod_HeatMap(DEP_Mod)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

