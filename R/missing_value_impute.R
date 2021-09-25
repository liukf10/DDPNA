#missing value knn impute and outlier remove
Data_impute <- function(data, inf = "inf", intensity = "LFQ", miss.value = NA,
                        splNExt = TRUE, maxNAratio = 0.5,
                        removeOutlier = TRUE,
                        outlierdata = "intensity", iteration = NA, sdout = 2,
                        distmethod = "manhattan", A.IAC = FALSE,
                        dohclust = FALSE, treelabels = NA,
                        plot = TRUE, filename = NULL,
                        text.cex = 0.7, text.col = "red", text.pos = 1,
                        text.labels = NA, abline.col = "red", abline.lwd = 2,
                        impute = TRUE, verbose = 1, ...){
  if (is.null(data[[inf]]) | is.null(data[[intensity]])) stop("data, inf or intensity is wrong.")
  if (class(data[[intensity]]) != "data.frame") stop("data[[intensity]] must be a data.frame.")
  if (is.character(as.matrix(data[[intensity]]))) stop("data[[intensity]] contains character.")
  #extract sample name
  if (splNExt) {
    inten <- data[[intensity]];
    if (intensity == "LFQ") colnames(inten) <- gsub("LFQ.intensity.(.*)", "\\1",
                                                    colnames(inten));
    if (intensity == "iBAQ") colnames(inten) <- gsub("iBAQ.(.*)", "\\1",
                                                     colnames(inten));
    if (intensity == "intensity") colnames(inten)<-gsub("Intensity.(.*)", "\\1",
                                                        colnames(inten));
    data[[intensity]] <- inten;
    rm(inten)
  }
  data <- list(inf = data[[inf]], intensity = data[[intensity]]);
  data <- .NAnum.proteomic_data(data, miss.value = miss.value, verbose = verbose);
  #table(data$inf$NA_num)
  data_imp <- .deNA2.proteomic_data(data, trunc(maxNAratio * ncol(data$intensity)));
  #remove Outlier Sample
  if (removeOutlier) {
    data_deoutlier <- try(.OutlierDetect(data_imp, iteration = iteration,
                                         intensity = "intensity", #outlierdata
                                         distmethod = distmethod, A.IAC = A.IAC,
                                         maxNAratio = maxNAratio,
                                         miss.value = miss.value, sdout = sdout,
                                         dohclust = dohclust, treelabels = treelabels,
                                         plot = plot, filename = filename,
                                         text.cex = text.cex, text.col = text.col,
                                         text.pos = text.pos, text.labels = text.labels,
                                         abline.col = abline.col, abline.lwd = abline.lwd,
                                         ...)); #outlierdata change to intensity
    if (class(data_deoutlier) != "try-error")
      data_imp <- data_deoutlier$data else
        warning("removeOutlier is failed.");
  }
  if (impute)
  {#knn impute
  dimp <- try(.knn_impute3.proteomic_data(data_imp, miss.value = miss.value,
                                          maxp = nrow(data_imp$intensity), ...));
  if (class(dimp) != "try-error")
    data_imp <- dimp else
      warning("impute is failed.")
  }
  data_imp
}



#geomean, relative_value and log2 value
.proteomic_data <- function(inf, intensity){
  geo_mean <- function(data) {
    data[data < 0] <- NA;
    log_data <- log(data);
    gm <- NA;
    for (i in 1:length(log_data[ ,1])){
      num <- as.numeric(log_data[i, ]);
      gm[i] <- exp(mean(num[is.finite(num)]));
    }
    gm
  }
  proteomic_data <- list();
  proteomic_data$inf <- inf;
  proteomic_data$intensity <- intensity;
  proteomic_data$geo_mean <- geo_mean(intensity);
  proteomic_data$relative_value <- intensity/proteomic_data$geo_mean;
  proteomic_data$log2_value <- log2(proteomic_data$relative_value);
  class(proteomic_data) <- "proteomic_data";
  proteomic_data
}

#NA number count
#190716 fix the problem when use the function more than one time will appear many NA_num column
.NAnum.proteomic_data <- function(data, miss.value = NA, colmax = 0.8, verbose = 1){
  #calculate NA number every column
  if (is.numeric(miss.value)){
    NA_num <- as.numeric(drop(rep(1, nrow(data$intensity))
                              %*% (data$intensity == miss.value)));
    pos <- which(NA_num > trunc(colmax * nrow(data$intensity)));
    if (length(pos) > 0) {
      if(verbose > 0) warning (paste("The sample", colnames(data$intensity)[pos],
                                     "have been removed"));
      data$intensity <- data$intensity[ ,-pos];
      data$relative_value <- data$relative_value[ ,-pos];
      data$log2_value <- data$log2_value[ ,-pos];
    }
    #calculate NA number every row
    NA_num <- as.numeric(drop((data$intensity==miss.value)
                              %*% rep(1, ncol(data$intensity))));
  }
  else if (is.na(miss.value)) {
    NA_num <- as.numeric(drop(rep(1, nrow(data$intensity))
                              %*% is.na(data$intensity)));
    pos <- which(NA_num > trunc(colmax * nrow(data$intensity)));
    if (length(pos) > 0) {
      if(verbose > 0) warning (paste("The sample",colnames(data$intensity)[pos],
                                     "have been removed"));
      data$intensity <- data$intensity[ ,-pos];
      data$relative_value <- data$relative_value[ ,-pos];
      data$log2_value <- data$log2_value[ ,-pos];
    }
    NA_num <- as.numeric(drop(is.na(data$intensity)
                              %*% rep(1, ncol(data$intensity))));
  }
  else { stop ("wrong missing value setting"); }
  if(any(colnames(data$inf) == "NA_num")) data$inf <- data$inf[ ,-which(colnames(data$inf) == "NA_num")]; #190716
  data$inf <- data.frame(data$inf, NA_num, stringsAsFactors = FALSE);
  data
}

#remove NA number large than n
.deNA.proteomic_data <- function(data, n, miss.value = NA){
  #NA_num <- NA;
  #for (i in 1:nrow(data$inf)) {  #length(data$inf$X#
  #  NA_num[i] <- sum(is.na(data$log2_value[i,]));}
  #calculate NA number every row
  if (!is.na(miss.value) & miss.value == 0) {
    NA_num <- as.numeric(drop((data$intensity == 0)
                              %*% rep(1, ncol(data$intensity))));
  }
  else if(is.na(miss.value)){
    NA_num <- as.numeric(drop(is.na(data$intensity)
                              %*% rep(1, ncol(data$intensity))));
  }
  else { stop ("wrong missing value setting");}
  pos <- which(NA_num <= n);
  row_name <- row.names(data$inf)[pos];
  extract2.proteomic_data <- function(x, row_num){
    #serach postion row_num in x $inf, ordered by row_num
    row <- match(row_num, row.names(x$inf));
    inf <- x$inf[row, ];
    intensity <- x$intensity[row, ];
    relative_value <- x$relative_value[row,] ;
    log2_value <- x$log2_value[row, ];
    #list(inf=inf,intensity=intensity,relative_value=relative_value,log2_value=log2_value)
    proteomic_data <- list();
    proteomic_data$inf <- inf;
    proteomic_data$intensity <- intensity;
    proteomic_data$relative_value <- relative_value;
    proteomic_data$log2_value <- log2_value;
    class(proteomic_data) <- "proteomic_data";
    proteomic_data
  }
  extract2.proteomic_data(data,row_name)
}
.deNA2.proteomic_data <- function(data, n){
  pos <- which(data$inf$NA_num <= n);
  inf <- data$inf[pos,];
  intensity <- data$intensity[pos,];
  proteomic_data <- list();
  proteomic_data$inf <- inf;
  proteomic_data$intensity <- intensity;
  class(proteomic_data) <- "proteomic_data";
  proteomic_data
}

#missing value knn impute
.knn_impute.proteomic_data <- function(data){
  geo_mean <- function(data) {
    data[data < 0] <- NA;
    log_data <- log(data);
    gm <- NA;
    for (i in 1:length(log_data[ ,1])) {
      num <- as.numeric(log_data[i, ]);
      gm[i] <- exp(mean(num[is.finite(num)]));
    }
    gm
  }
  if (!requireNamespace("impute", quietly = TRUE)) {
    stop ("impute in Bioconductor needed for this function to work. Please install it.",
         call. = FALSE)
  }
  intensity <- impute::impute.knn(as.matrix(data$intensity));
  intensity <- as.data.frame(intensity$data);
  geo_mean <- geo_mean(intensity);
  relative_value <- intensity/geo_mean;
  #relative_value <- impute::impute.knn(as.matrix(data$relative_value));
  #relative_value <- as.data.frame(relative_value$data);
  log2_value <- log2(relative_value);
  #log2_value <- impute::impute.knn(as.matrix(data$log2_value));
  #log2_value<as.data.frame(log2_value<$data);
  proteomic_data <- list();
  proteomic_data$inf <- data$inf;
  proteomic_data$intensity <- intensity;
  proteomic_data$relative_value <- relative_value;
  proteomic_data$log2_value <- log2_value;
  class(proteomic_data) <- "proteomic_data";
  proteomic_data
}
.knn_impute2.proteomic_data <- function(data){
  geo_mean <- function(data) {
    data[data < 0] <- NA;
    log_data <- log(data);
    gm <- NA;
    for (i in 1:length(log_data[,1])){
      num <- as.numeric(log_data[i, ]);
      gm[i] <- exp(mean(num[is.finite(num)]));
    }
    gm
  }
  if (!requireNamespace("impute", quietly = TRUE)) {
    stop ("impute in Bioconductor needed for this function to work. Please install it.",
         call. = FALSE)
  }
  intensity <- impute::impute.knn(as.matrix(data$intensity));
  intensity <- as.data.frame(intensity$data);
  geo_mean <- geo_mean(intensity);
  #relative_value <- intensity/geo_mean;
  relative_value <- impute::impute.knn(as.matrix(data$relative_value));
  relative_value <- as.data.frame(relative_value$data);
  #log2_value <- log2(relative_value);
  log2_value <- impute::impute.knn(as.matrix(data$log2_value));
  log2_value <- as.data.frame(log2_value$data);
  proteomic_data <- list();
  proteomic_data$inf <- data$inf;
  proteomic_data$intensity <- intensity;
  proteomic_data$relative_value <- relative_value;
  proteomic_data$log2_value <- log2_value;
  class(proteomic_data) <- "proteomic_data";
  proteomic_data
}
.knn_impute3.proteomic_data <- function(data,miss.value = NA,...){
  geo_mean <- function(data) {
    data[data < 0] <- NA;
    log_data <- log(data);
    gm <- NA;
    for(i in 1:length(log_data[,1])){
      num <- as.numeric(log_data[i, ]);
      gm[i] <- exp(mean(num[is.finite(num)]));
    }
    gm
  }
  if (!requireNamespace("impute", quietly = TRUE)) {
    stop ("impute in Bioconductor needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (length(miss.value) != 1 ) stop("wrong miss value");
  if (!miss.value %in% c(0, 1, NA)) stop("wrong miss value");
  data$intensity[data$intensity == miss.value] <- NA;
  intensity <- impute::impute.knn(as.matrix(data$intensity),...);
  intensity <- as.data.frame(intensity$data);
  geo_mean <- geo_mean(intensity);
  relative_value <- intensity/geo_mean;
  #relative_value <- impute::impute.knn(as.matrix(data$relative_value));
  #relative_value <- as.data.frame(relative_value$data);
  log2_value <- log2(relative_value);
  #log2_value <- impute::impute.knn(as.matrix(data$log2_value));
  #log2_value<as.data.frame(log2_value<$data);
  proteomic_data <- list();
  proteomic_data$inf <- data$inf;
  proteomic_data$intensity <- intensity;
  proteomic_data$relative_value <- relative_value;
  proteomic_data$log2_value <- log2_value;
  class(proteomic_data) <- "proteomic_data";
  proteomic_data
}



#190615 fix a bug when doclust, it will error use the previous function
#190716 fix outlier detect problem when miss.value is not NA, the outlier detect will different.
#210728 the param:outilerdata change to param:intensity
#210728 add a new choose in distmethod: bicor
#210728 fix the plot and doclust plot
.OutlierDetect <- function(data, iteration = NA, intensity = "intensity",
                           distmethod = "manhattan", A.IAC = FALSE,
                           maxNAratio = 0.5,
                           miss.value = NA, sdout = 2,
                           dohclust = FALSE, treelabels = NA,
                           plot = TRUE, filename = NULL,
                           text.cex = 0.7, text.col = "red", text.pos = 1,
                           text.labels = NA, abline.col = "red", abline.lwd = 2,
                           ...){
  distmethod <- match.arg(distmethod,
                          c("manhattan","euclidean", "canberra","correlation","bicor"));
  #210728 del; outlierdata <- match.arg(outlierdata, c("intensity","relative_value","log2_value"));
  #210728 change intensity dataframe name to "intensity"
  if (intensity != "intensity") #210728 del; data <- .proteomic_data(data$inf, data$intensity);
    names(data)[names(data) == intensity] <- "intensity";
  outliersample = 1; outlier = NULL; iter = 0;
  if (!is.na(miss.value)) #210728 del; data[[outlierdata]][data[[outlierdata]] == miss.value] <- NA; #190716
    data[["intensity"]][data[["intensity"]] == miss.value] = NA;#210728
  while (outliersample > 0){
    iter <- iter + 1;
    if (distmethod %in% c("correlation","bicor"))
    {
      # intra-assay correlation
      #210728 del; IAC <- cor(data[[outlierdata]], use = "p");
      if(distmethod == "correlation") IAC <- cor(data[["intensity"]], use = "pairwise.complete.obs"); #210728
      #210728 add bicor function
      if(distmethod == "bicor") {
        if (!requireNamespace("WGCNA", quietly = TRUE)) {
          stop ("WGCNA in Bioconductor needed when using bicor method. Please install it.",
                call. = FALSE)}
        IAC <- WGCNA::bicor(data[["intensity"]], use='pairwise.complete.obs')}
      if (A.IAC) IAC <- ((1+IAC)/2)^2;
      #IAC histogram
      #hist(IAC,breaks=50,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
      # 1-IAC as distance
      dist <- as.dist(1 - IAC);
    } else {
      #dist <- dist(t(data[[outlierdata]]), method = distmethod);
      dist <- dist(t(data[["intensity"]]), method = distmethod);
      IAC <- as.matrix(dist);
    }
    # Another way to visualize outliers is to calculate the mean IAC for each array and examine this distribution
    meanIAC <- apply(IAC, 2, mean);#hist(meanIAC,breaks=10)
    sdCorr <- sd(meanIAC);#sd
    numbersd <- abs(meanIAC - mean(meanIAC)) / sdCorr;
    #if meanIAC is normal distribution,+-2 means 95%
    all_name <- dimnames(data[["intensity"]])[[2]];
    #210728 del; out_name <- dimnames(data[[outlierdata]])[[2]][numbersd > sdout];
    out_name <- all_name[numbersd > sdout];#210728
    out_pos <- as.numeric(which(numbersd > sdout));
    if(dohclust) {
      sampleTree <- hclust(dist, method = "average");
      #210728
      if (length(treelabels)!= dim(data$intensity)[[2]] & !all(is.na(treelabels))) {
        treelabels = NA;
        warning("the number of treelabels is differ with sample number and use sample name instead of.")
      }
      if (all(is.na(treelabels))) treelabels <- dimnames(data$intensity)[[2]];
      if (is.character(filename) & length(filename) == 1){
        if(!dir.exists("plot")) dir.create("plot"); #210615fix exist plot dir
        pdf(paste0("plot/",filename, " sampleTree ",iter,".pdf")) #210728 save file in plot dir
      }
      plot(sampleTree, main = "Sample clustering to detect outliers",
           sub="", labels=treelabels,...);
      if (is.character(filename) & length(filename) == 1) dev.off()
      #cluster <- hclust(as.dist(1-IAC),method="average")
      #plot(cluster,cex=0.7,labels=dimnames(data$intensity)[[2]])
    }
    #outliers
    if (plot) {
      #210728
      if (length(text.labels)!= dim(data$intensity)[[2]] & !all(is.na(text.labels))) {
        text.labels = NA;
        warning("the number of text.labels is differ with sample number and use sample name instead of.")
      }
      #210728 use textlabels instead of text.labels
      if (all(is.na(text.labels))) textlabels <- out_name;
      #210728
      if (length(text.labels) == dim(data$intensity)[[2]]) {
        textlabels <- text.labels[out_pos];
        text.labels <- text.labels[-out_pos];
      }
      col <- rep("black", length(numbersd));
      col[out_pos] <- "red";
      pch <- rep(20, length(numbersd));
      pch[out_pos] <- 21;
      if (is.character(filename) & length(filename) == 1){
        if(!dir.exists("plot")) dir.create("plot"); #210615fix exist plot dir
        pdf(paste0("plot/",filename, " outliersample ",iter,".pdf")) #210728 save file in plot dir
      }
      plot(numbersd, pch = pch, bg = "red", col = col, xlab = "");
      if (length(out_pos) > 0){
        text(out_pos, numbersd[out_pos], labels = textlabels,
             pos = text.pos, cex = text.cex, col = text.col);
        abline(h = sdout, col = abline.col, lwd = abline.lwd);
      }
      if (is.character(filename) & length(filename) == 1) dev.off()
      #text.labels <- NA;
    }
    outliersample <- length(out_pos);
    outlier <- c(outlier, out_name);
    if (outliersample > 0) {
      data[["intensity"]] <- data[["intensity"]][ ,-out_pos];
      treelabels <- treelabels[-out_pos];#190615
      data$inf <- data$inf[ ,-which(colnames(data$inf) == "NA_num")];
      #data <- .NAnum.proteomic_data(data, miss.value = miss.value);
      data <- .NAnum.proteomic_data(data, miss.value = NA);#190716
      data <- .deNA2.proteomic_data(data, trunc(maxNAratio*ncol(data$intensity)));
      #if (outlierdata != "intensity")
      #  data <- .proteomic_data(data$inf, data$intensity);#210728 del
    }
    if (!is.na(iteration) & iter >= iteration) outliersample <- 0; #190716
  }
  #if (!is.na(miss.value)) data[[outlierdata]][is.na(data[[outlierdata]])] <- miss.value; #190716
  if (!is.na(miss.value)) data[["intensity"]][is.na(data[["intensity"]])] <- miss.value; #210728
  list(outlier = outlier, data = data)
}
