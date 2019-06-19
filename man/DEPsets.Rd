\name{DEPsets}
\alias{DEPsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DEPsets
}
\description{
extract two or more IDsets interesection set and complementary set and define the colors.
}
\usage{
DEPsets(datalist, colors = c("red", "green", "blue"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datalist}{
a list contains more than two ID sets.
}
  \item{colors}{
 define each ID sets color.
}
}

\value{
a list contains interesection set and complementary set information and colors.
  \item{gene.set }{a list of each set ID information.}
  \item{color.code }{the colors of each set}
%% ...
}

\author{
Kefu Liu%%  ~~who you are~~
}

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
stat1 <- stat$new.ID[stat$ad > 1]
stat2 <- stat$new.ID[stat$asym > 1]
datalist <- list(stat1 = stat1, stat2 = stat2)
sets <- DEPsets(datalist)

}

