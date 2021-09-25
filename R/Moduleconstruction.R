#draw graph of soft Threshold

SoftThresholdScaleGraph <- function(data,
                                    xlab = "Soft Threshold (power)",
                                    ylab = "Scale Free Topology Model Fit, signed R^2",
                                    main = "Scale independence",
                                    filename = NULL){
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop ("WGCNA package is needed. Please install it.",
          call. = FALSE)}
  sft <- WGCNA::pickSoftThreshold(data);
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2));
  if (is.character(filename) & length(filename) == 1){
    #if (!dir.exists("plot")) dir.create("plot");
    #pdf(paste0("plot/", "SoftThershold ", filename, ".pdf"))
    pdf(paste0( "SoftThershold ", filename, ".pdf")) #200703
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





#190607 add the question whether load WGCNA package
wgcnatest <- function(data, power = NULL, TOMType = "unsigned",
                      detectCutHeight = NULL, maxBlockSize = 5000,
                      deepSplit = TRUE, minModSize = TRUE,
                      pamRespectsDendro = FALSE,
                      minKMEtoStay = TRUE,
                      minCoreKME = FALSE,
                      reassignThreshold = FALSE,
                      mergeCutHeight = FALSE,
                      maxModNum = 30,
                      minModNum = 8,
                      MaxMod0ratio = 0.3) {

  if( ncol(data) < 50) stop("The number of proteins is too less to fit the function.")
  if( maxBlockSize > 45000) stop("maxBlockSize is too larger to fit the function.")
  if(!is.numeric(power) & length(power) > 1)
    stop("power is not a correct value.")
  if(!is.numeric(detectCutHeight) & length(detectCutHeight) > 1)
    stop("detectCutHeight is not a correct value.")
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop ("WGCNA package is needed. Please install it.",
          call. = FALSE)}
#  question <- function (...)#190607
#  {
#    yes <- c("Yes", "Definitely", "Sure",
#             "I agree", "Absolutely")
#    no <- c("No", "Nope", "Not Sure")
#    cat(paste0(..., collapse = ""))
#    qs <- c(sample(yes, 1), sample(no, 1))
#    menu(qs) == 2
#  }
  #cat("The function need load WGCNA package first.\n"); #190607
#  if (question("Did you already loaded WGCNA package?")) #190607
#    stop("Please install WGCNA package and run require(\"WGCNA\") first.")#190607

  cor <- WGCNA::cor;
  block <- ceiling(ncol(data) / maxBlockSize);
  .parameter <- function(x, num1, num2, name){
    if(is.null(x) || all(is.na(x)))
      x = num1 else if (isFALSE(x))
        x = num1 else if (isTRUE(x))
          x <- num2 else if (!is.numeric(x))
            stop(paste("wrong", name));
        x
  }

  deepSplit <- .parameter(deepSplit, 2, c(0,1,2,3,4), "deepSplit");
  minModSize <- .parameter(minModSize, 20, c(15, 20, 30, 50), "minModSize");
  minKMEtoStay <- .parameter(minKMEtoStay, 0.3, c(0.1, 0.2, 0.3), "minKMEtoStay");
  reassignThreshold <- .parameter(reassignThreshold, 1e-6, c(0.01, 0.05), "reassignThreshold");
  mergeCutHeight <- .parameter(mergeCutHeight, 0.15, c(0.15, 0.3, 0.45), "mergeCutHeight");
  minCoreKME <- .parameter(minCoreKME, 0.5, c(0.4, 0.5), "minCoreKME");

  if(is.null(power)) {
    sft <- WGCNA::pickSoftThreshold(data);
    power <- sft$powerEstimate}
  #r when p<0.05    detectcutheight
  if(is.null(detectCutHeight)){
    num <- nrow(data); p = 0; r <- 0.6
    while (p < 0.05) {
      r <- r-0.01
      p <- pt(-r*sqrt((num-2)/(1-r^2)),num-2)*2
    }
    while (p > 0.05) {
      r <- r+0.001
      p <- pt(-r*sqrt((num-2)/(1-r^2)),num-2)*2
    }
    while (p < 0.05) {
      r <- r-0.0001
      p <- pt(-r*sqrt((num-2)/(1-r^2)),num-2)*2
    }
    while (p > 0.05) {
      r <- r+0.00001
      p <- pt(-r*sqrt((num-2)/(1-r^2)),num-2)*2
    }
    detectCutHeight = 1-r^power
    if(detectCutHeight < 0.995) detectCutHeight = 0.995
    rm(r,p,num)
  }

  net <- try(WGCNA::blockwiseModules(data, power = power, maxBlockSize = maxBlockSize,
                                     TOMType = TOMType, deepSplit = 2,minModuleSize = 20,
                                     reassignThreshold = 0.05, mergeCutHeight = 0.15,
                                     numericLabels = TRUE, pamRespectsDendro = FALSE,
                                     saveTOMs = TRUE, saveTOMFileBase = "blockwisetom",
                                     verbose = 0), silent = TRUE)
  if (class(net) == "try-error") #190607
    stop("Run require(\"WGCNA\") first, if error again, please check the data.")

  module <- NULL; unmergedmodule <- NULL; ModNum <- NULL;
  Modnum <- rep(NA, maxModNum + 1); iter = 1;
  deepsplit <- NULL; size <- NULL;  KMEtostay <- NULL;
  reass <- NULL; ModH <- NULL; KMECore <- NULL;
  for (i1 in deepSplit) {
    for (i2 in minModSize) {
      for (i3 in minKMEtoStay) {
        for(i4 in reassignThreshold){
          for (i5 in mergeCutHeight) {
            for (i6 in minCoreKME) {
              net <- WGCNA::blockwiseModules(data, loadTOM = TRUE, power = power,
                                             maxBlockSize = maxBlockSize, TOMType = TOMType,
                                             deepSplit = i1, minModuleSize = i2,
                                             detectCutHeight = detectCutHeight,
                                             minKMEtoStay = i3, minCoreKME = i6,
                                             reassignThreshold = i4, mergeCutHeight = i5,
                                             numericLabels = TRUE,
                                             pamRespectsDendro = pamRespectsDendro,
                                             saveTOMs = FALSE,
                                             verbose = 0);
              if(table(net$colors)[1]/ncol(data) < MaxMod0ratio &
                 length(table(net$colors)) <= length(Modnum) &
                 length(table(net$colors)) >= minModNum) {
                if (is.null(module)) module <- net$colors else
                  module <- data.frame(module, net$colors);
                if (is.null(unmergedmodule)) unmergedmodule <- net$unmergedColors else
                  unmergedmodule <- data.frame(unmergedmodule, net$unmergedColors);
                if (is.null(ModNum)) ModNum <- data.frame(Modnum) else
                  ModNum <- data.frame(ModNum, Modnum);
                #if(length(table(net$colors)) > length(Modnum))
                # warning(paste0("deepSplit",i,",minModSize",j,",minKMEtoStay ",k,",h ",l, " have more than 100 mod"));
                mn <- min(length(table(net$colors)), length(Modnum));
                ModNum[1:mn, iter] <- table(net$colors)[1 : mn];
                if(length(deepSplit) != 1) deepsplit <- c(deepsplit, i1);
                if(length(minModSize) != 1) size <- c(size, i2);
                if(length(minKMEtoStay) != 1) KMEtostay <- c(KMEtostay, i3);
                if(length(reassignThreshold) != 1) reass <- c(reass, i4);
                if(length(mergeCutHeight) != 1) ModH <- c(ModH, i5);
                if(length(minCoreKME) != 1) KMECore <- c(KMECore, i6);
                #colnames(ModNum)[iter] <- paste0("split",i,",size",j,",KME ",k,",h ",l)
                iter = iter + 1
              }
            }
          }
        }
      }
    }
  }
  rownames(ModNum) <- 0:maxModNum
  ModNum <- rbind(deepSplit = if(!is.null(deepsplit)) deepsplit,
                  minModuleSize = if(!is.null(size)) size,
                  minKMEtoStay = if(!is.null(KMEtostay)) KMEtostay,
                  reassignThreshold = if(!is.null(reass)) reass,
                  mergeCutHeight = if(!is.null(ModH)) ModH,
                  minCoreKME = if(!is.null(KMECore)) KMECore,
                  ModNum)
  filetoremove <- paste("blockwisetom-block.", 1:block, ".RData", sep="")
  file.remove(filetoremove)
  ModNum
}






ME_inf <- function(MEs, data, intensity.type = "LFQ", rowname = NULL){
  intensity.type <- match.arg(intensity.type,
                              c("LFQ","iBAQ","intensity","none"));
  n <- pmatch(intensity.type,
              c("LFQ", "iBAQ", "intensity", "none"));
  name <- c("LFQ.intensity.", "iBAQ.", "Intensity.");
  if (length(rowname)!= 0 & length(rowname) != ncol(data)) {
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
                                                    stringsAsFactors = FALSE)
  else inf <- data.frame(ProteinID = inf,
                         moduleNum = net$colors,
                         stringsAsFactors = FALSE)
}



