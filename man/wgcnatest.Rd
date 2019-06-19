\name{wgcnatest}
\alias{wgcnatest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 wgcnatest
}
\description{
The major parameter optimization in function \code{blockwiseModules} in WGCNA package. The function will do a series of network construction by change various parameter in \code{blockwiseModules} and record the result. (it will take a long time)
}
\usage{
wgcnatest(data, power = NULL, TOMType = "unsigned",
          detectCutHeight = NULL, maxBlockSize = 5000,
          deepSplit = TRUE, minModSize = TRUE,
          pamRespectsDendro = FALSE,
          minKMEtoStay = TRUE,
          minCoreKME = FALSE,
          reassignThreshold = FALSE,
          mergeCutHeight = FALSE,
          maxModNum = 30,
          minModNum = 8,
          MaxMod0ratio = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  protein quantification data used in network construction. Row is sample. Column is protein ID.
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{power}{
  Soft-thresholding power for network construction. The default value is NULL. it will run \code{pickSoftThreshold} function in WGCNA package to pick the lowest appropriate power.
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{TOMType}{
  one of "none", "unsigned", "signed".
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{detectCutHeight}{
  dendrogram cut height for module detection.
  The default value is NULL, which means it will calculate the cutheight through correlation r when p value is 0.05. When the value is larger than 0.995, it will set to detectCutHeight or 0.995.
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{maxBlockSize}{
  integer giving maximum block size for module detection.
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{deepSplit}{
  The default value is TRUE, which means the function will test deepSplit from 0 to 4. If the value is FALSE, deepSplit is 2. You also can setting integer value between 0 and 4 by yourself.
  integer value between 0 and 4.
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{minModSize}{
  minimum module size for module detection.
  The default value is TRUE, which means the function will test 15, 20, 30, 50. If the value is FALSE, minModSize is 20. You also can setting integer value by yourself.
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{pamRespectsDendro}{
  a logical value indicated that whether do pamStage or not.
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{minKMEtoStay}{
  The default value is TRUE, which means the function will test 0.1, 0.2, 0.3. If the value is FALSE, minKMEtoStay is 0.3. You also can setting value by yourself.
  Value between 0 to 1.
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{minCoreKME}{
  The default value is FALSE, minCoreKME is 0.5. If the value is TRUE, which means the function will test 0.4 and 0.5. You also can setting value by yourself.
  Value between 0 to 1.
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{reassignThreshold}{
  p-value ratio threshold for reassigning genes between modules.
  The default value is FALSE, reassignThreshold is 1e-6. If the value is TRUE, which means the function will test 0.01 and 0.05. You also can setting value by yourself.
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{mergeCutHeight}{
  dendrogram cut height for module merging.
  The default value is FALSE, mergeCutHeight is 0.15. If the value is TRUE, which means the function will test 0.15, 0.3 and 0.45. You also can setting value by yourself.
  More information can get from \code{blockwiseModules} in WGCNA package.
  }
  \item{maxModNum}{
  The maximum module number. If network construction make more than maxModnum of modules. The result will not record.
  }
  \item{minModNum}{
The mininum module number. If network construction make less than minModNum of modules. The result will not record.
}
  \item{MaxMod0ratio}{
  The maximum Mod0 protein numbers ratio in total proteins. If network construction make more than MaxMod0ratio in module 0. The result will not record.
}
  }
\details{
More information can get from \code{blockwiseModules} in WGCNA package.
}
\value{
a data.frame contains protein number in each module and the parameter information.

}

\author{
Kefu Liu
}

\examples{
\donttest{
data(imputedData)
wgcnadata <- t(imputedData$intensity)
sft <- SoftThresholdScaleGraph(wgcnadata)
# It will take a lot of time
if (requireNamespace("WGCNA", quietly = TRUE)){
require("WGCNA")
WGCNAadjust <- wgcnatest(wgcnadata, power = sft$powerEstimate)
}
}
}
