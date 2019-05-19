#library("WGCNA")


#prodata<-BLSAFC.LFQ.imp$log2_value
#data = t(prodata)
#draw graph of soft Threshold

SoftThresholdScaleGraph <- function(data,
                                    xlab = "Soft Threshold (power)",
                                    ylab = "Scale Free Topology Model Fit, signed R^2",
                                    main = "Scale independence",
                                    filename = NULL){
  sft <- pickSoftThreshold(data);
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2));
  if (is.character(filename) & length(filename) == 1){
    if (!dir.exists("plot")) dir.create("plot");
    pdf(paste0("plot/", "SoftThershold ", filename, ".pdf"))
  }
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[ ,2],
       xlab = xlab, ylab = ylab, type = "n",
       main = main);
  text(sft$fitIndices[ ,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=0.9,col="red");
  abline(h=0.90,col="red");
  if (is.character(filename) & length(filename) == 1) dev.off()
  sft;
}

# Modules_Build<-WGCNA::blockwiseModules
# net = blockwiseModules(data, power = sft$powerEstimate, maxBlockSize = 6000,
#                       TOMType = "unsigned", deepSplit = 4,minModuleSize = 17,
#                       reassignThreshold = 0.05, mergeCutHeight = 0.07,
#                       numericLabels = TRUE,pamRespectsDendro = FALSE,
#                       saveTOMs=F,loadTOMs=F,
#                       saveTOMFileBase = "blockwisetom",
#                       verbose = 3)

#table(net$colors)


#MEs = net$MEs


ME_inf <- function(MEs, data, intensity.type = "LFQ", rowname = NULL){
  intensity.type <- match.arg(intensity.type,
                              c("LFQ","iBAQ","intensity","none"));
  n <- pmatch(intensity.type,
              c("LFQ", "iBAQ", "intensity", "none"));
  name <- c("LFQ.intensity.", "iBAQ.", "Intensity.");
  if (length(rowname)!=0 & length(rowname) != ncol(data)) {
    rowname <- NULL;
    warning ("wrong rowname input, rowname will be set to NULL.");
  }
  if (n < 4) row.names(MEs) <- gsub(paste0(name[n], "(.*)"), "\\1", colnames(data))
  else if(n == 4 & is.null(rowname)) row.names(MEs) <- colnames(data)
  else row.names(MEs) <- rowname;
  MEs
}


Module_inf <- function(net, inf, inftype = "Convert", IDname = NULL, ...){
  inftype <- match.arg(inftype, c("Convert","MaxQ","none"));
  moduleColors <- labels2colors(net$colors);
  if (inftype == "Convert") inf <- data.frame(ori.ID = inf$ori.ID,
                                              new.ID = inf$new.ID,
                                              stringsAsFactors = FALSE)
  else if (any(colnames(inf) == IDname)) inf <- inf[ , IDname]
  else stop("IDname is uncorrect.");
  if (inftype == "MaxQ" & length(IDname) == 1) {
    inf <- P.G.extract(inf,...);
    names(inf) <- c("db.type","ID","ENTRY.NAME");
  }
  if (class(inf) == "data.frame") inf <- data.frame(inf,
                                                    moduleNum = net$colors,
                                                    moduleColors,
                                                    stringsAsFactors = FALSE)
  else inf <- data.frame(ProteinID = inf,
                         moduleNum = net$colors,
                         moduleColors,
                         stringsAsFactors = FALSE)
}



