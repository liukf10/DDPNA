#extract pca compont1,2
modpcomp <- function(data, colors, nPC = 2,
                     plot = FALSE, filename = NULL, group = NULL){
  modnum <- as.numeric(rownames(table(colors)));
  pc <- list();
  length(pc) <- length(modnum);
  names(pc) <- paste0("mod",modnum);
  if(plot & (is.null(group) | length(group) != ncol(data)))
    warning ("Don't plotting because group information is error.");
  for(i in modnum){
    modulegene <- as.matrix(data[colors == i, ]);
    modulepca <- prcomp(t(scale(t(modulegene))), center = FALSE, scale = FALSE);
    pc[[paste0("mod",i)]] <- modulepca$rotation[ , 1:nPC];
    #deprecated
    if(plot & length(group) == ncol(data)) {
      #20190614
      dfcomp <- t(modulegene);
      dfcomp <- data.frame(dfcomp, group = group)
      pic <- try(autoplot(prcomp(scale(t(modulegene))),data = dfcomp, colour = "group"),
                 silent = TRUE)
      if ( all(class(pic)=="try-error"))
        warning ("No plot. Please install ggfority package and run require(\"ggfority\") first.");
      if (FALSE){
        #install.packages("devtools", repo="http://cran.us.r-project.org")
        #library(devtools)
        #install_github("vqv/ggbiplot")
        #library(ggbiplot)
        #      if (!requireNamespace("ggbiplot", quietly = TRUE)) {
        #        stop("ggbiplot needed for this function to work. Please install it.", call. = FALSE)}
        modulepca2 <- prcomp(scale(t(modulegene)), center = FALSE, scale = FALSE);
        #pc[[paste0("mod",i)]] <- modulepca2$x[ , 1:nPC]
        pic <- .ggbiplot(modulepca2, obs.scale = 1, var.scale = 1,
                         groups = group, var.axes = FALSE) +
          scale_color_discrete(name = '') +
          theme_bw();}
      if ( all(class(pic)!="try-error") ) {
        if (is.character(filename) & length(filename) == 1){
          if (!dir.exists("plot")) dir.create("plot");
          pdf(paste0("plot/", filename, " PCA plot mod", i, ".pdf"))
          dev.off()
        } else print(pic)
      }
    }
  }
  pc
}


groupmean <- function(data, group, method = c("mean","median"), name = TRUE){
  if (length(group) != ncol(data))
    stop ("Group length and data col length is not the same.");
  group <- as.factor(group);
  mean.matrix <- NULL;
  method <- match.arg(method, c("mean","median"));
  met <- match.fun(method);
  for(i in 1:nlevels(group)){
    pos <- which(group == levels(group)[i]);
    mean <- apply(data[ , pos], 1, met, na.rm = T); #20190514
    if (is.null(mean.matrix)) mean.matrix <- data.frame(mean)
    else mean.matrix <- data.frame(mean.matrix, mean);
    if(name) colnames(mean.matrix)[i] <- paste0(levels(group)[i], method)
    else colnames(mean.matrix)[i] <- paste0(levels(group)[i]);
  }
  mean.matrix
}


##
#20190509 fix error when any group in groups have no replication
multi.t.test <- function(data, group,
                         sig = 0.05, Adj.sig = TRUE,
                         grpAdj = "bonferroni",
                         geneAdj="fdr",...)
{
  tcutoff <- function(MEs_t, sig = 0.05){
    MEs_t[MEs_t > sig] <- NA;
    NA_num <- as.numeric(drop(rep(1, nrow(MEs_t)) %*% is.na(MEs_t)));
    pos1 <- which(NA_num != nrow(MEs_t));
    NA_num <- as.numeric(is.na(MEs_t) %*%  drop(rep(1, ncol(MEs_t))));
    pos2 <- which(NA_num != ncol(MEs_t));
    list(pos2 = pos2, pos1 = pos1)
  }
  MEs_col <- t(data);
  if (any(table(group) == 1)) {
    delgroup <- names(table(group))[table(group) == 1];
    delpos <- which(group %in% delgroup);
    warning(paste("Group", delgroup, "is deleted because it have no replication. "));
    group <- group[-delpos];
    MEs_col <- MEs_col[-delpos, ];
  }
  if (nrow(MEs_col) != length(group))
    stop ("Group length and data col length is not the same.");
  if (length(table(group)) < 2)
    stop ("It should be at least two groups to do t.test.");
  if (is.null(colnames(MEs_col))) stop("wrong data rownames");
  MEs.t <- NULL; vsname <- NULL; Nnum <- NULL;
  if(!is.null(Nnum))
  {
    grpnum <- length(table(group)); k = 1;
    while((grpnum - k) > 1 ){
      Nnum <- c(Nnum, (grpnum - 1) * (grpnum - k - 1) + (1: (grpnum - k - 1)) )
      k = k + 1
    }
  }
  for (i in 1:ncol(MEs_col)) {
    t <- pairwise.t.test(MEs_col[ , i], group, p.adjust.method = grpAdj, ...);
    t.list <- as.vector(as.matrix(t$p.value))
    if (is.null(Nnum)) Nnum <- which(is.na(t.list) & !is.nan(t.list)); #20190514
    if (is.null(MEs.t)){
      MEs.t <- data.frame(M = t.list[-Nnum]);
      for (j in 1:ncol(t$p.value)) {
        for (k in 1:nrow(t$p.value)) {
          vsname <- c(vsname,
                      paste(colnames(t$p.value)[j], "vs", rownames(t$p.value)[k]));
        }
      }
      rownames(MEs.t) <- vsname[-Nnum];
    } else MEs.t <- data.frame(MEs.t, M = t.list[-Nnum]);
    colnames(MEs.t)[ncol(MEs.t)] <- colnames(MEs_col)[i];
  }
  geneAdj <- match.arg(geneAdj, p.adjust.methods);
  if (geneAdj != "none") {
    MEs.q <- MEs.t;
    rownames(MEs.q) <- paste(rownames(MEs.q), geneAdj);
    for (i in 1:nrow(MEs.t)) {
      MEs.q[i,] <- p.adjust(MEs.t[i, ], method = geneAdj);
    }
    if (Adj.sig) {
      pos <- tcutoff(MEs.q, sig = sig);
      #MEs.q[MEs.q > sig] <- NA;
      #NA_num <- as.numeric(drop(rep(1, nrow(MEs.q)) %*% is.na(MEs.q)));
      #pos <- which(NA_num != nrow(MEs.q));
      #NA_num <- as.numeric(is.na(MEs.q) %*%  drop(rep(1, ncol(MEs.q))));
      #pos2 <- which(NA_num != ncol(MEs.q));
      MEs.q <- MEs.q[pos$pos2, pos$pos1];
      MEs.t <- MEs.t[pos$pos2, pos$pos1];
    } else {
      pos <- tcutoff(MEs.t, sig = sig);
      #MEs.t[MEs.t > sig] <- NA;
      #NA_num <- as.numeric(drop(rep(1, nrow(MEs.t)) %*% is.na(MEs.t)));
      #pos <- which(NA_num != nrow(MEs.t));
      #NA_num <- as.numeric(is.na(MEs.t) %*%  drop(rep(1, ncol(MEs.t))));
      #pos2 <- which(NA_num != ncol(MEs.t));
      MEs.t <- MEs.t[pos$pos2, pos$pos1];
      MEs.q <- MEs.q[pos$pos2, pos$pos1];
    }
    MEs.all <- rbind(MEs.t, MEs.q);
  }
  if (geneAdj == "none") {
    pos <- tcutoff(MEs.t, sig = sig);
    #MEs.t[MEs.t > sig] <- NA;
    #NA_num <- as.numeric(drop(rep(1, nrow(MEs.t)) %*% is.na(MEs.t)));
    #pos <- which(NA_num!=nrow(MEs.t));
    #NA_num <- as.numeric(is.na(MEs.t) %*%  drop(rep(1, ncol(MEs.t))));
    #pos2 <- which(NA_num!=ncol(MEs.t));
    MEs.t <- MEs.t[pos$pos2, pos$pos1];
    MEs.all <- MEs.t;
  }
  t(MEs.all)
}



#anova analysis
##
#20190509 fix error when any group in groups have no replication
#20190514 fix some value is NA
anova_p <- function(data, group)
{
  if (ncol(data) != length(group))
    stop ("Group length and data col length is not the same.");
  data <- t(data);
  if (any(table(group) == 1)) {
    delgroup <- names(table(group))[table(group) == 1];
    delpos <- which(group %in% delgroup);
    warning(paste("Group", delgroup, "is deleted because it have no replication. "));
    group <- group[-delpos];
    data <- data[-delpos, ];
  }
  if (length(table(group)) < 2)
    stop ("It should be at least two groups to do anova analysis.");
  if (is.null(colnames(data))) stop("wrong data rownames");
  aov.p <- NULL;
  for (i in 1:ncol(data)) { #20190514
    gfactor<- levels(factor(group))
    barttest <- NULL; samevalue <- NULL;
    for(j in 1:length(gfactor)) {
      bart <- data[which(group == gfactor[j]), i]
      samevalue <- c( samevalue, length(table(bart)) )
      barttest <- c(barttest, sum(!is.na(bart)))
    }
    if (all(barttest > 1) & all(samevalue > 1) ) {
      Btlt <- bartlett.test(data[,i] ~ group);
      if(Btlt$p.value < 0.05){
        anova <- oneway.test(data[,i]~group);
        aov.p <- c(aov.p, anova$p.value);
      }
      if (Btlt$p.value > 0.05) {
        anova2 <- aov(data[,i]~group);
        anova3<-summary(anova2)[[1]];
        aov.p<-c(aov.p, anova3$`Pr(>F)`[1]);
      }
    } else aov.p<-c(aov.p, 2)
  }
  aov.p
}




#proteins which foldchange larger than cutoff   vs.set2 vs vs.set1
#20190531 fix the order when datatype is none
fc.pos <- function(fc, vs.set2, vs.set1 = "WT", cutoff = 1,
                   datatype = c("none", "log2"), fctype = "all",
                   order = TRUE){
  datatype <- match.arg(datatype, c("none", "log2"));
  fctype <- match.arg(fctype, c("all", "up", "down"));
  if (FALSE){
    if (datatype == "none") {
      fc <- log2(fc);
    }
    cutoff <- log2(cutoff);
    fc_1vs2 <- fc[ ,colnames(fc) == vs.set2] - fc[ , colnames(fc) == vs.set1];
    if (fctype == "up") fc.pos <- which(fc_1vs2 > cutoff);
    if (fctype == "down") fc.pos <- which(fc_1vs2 < -cutoff);
    if (fctype == "all")  fc.pos <- which(fc_1vs2 > cutoff | fc_1vs2 < -cutoff);
  }
  if (datatype == "log2") {
    cutoff <- log2(cutoff);
    fc_1vs2 <- fc[ ,colnames(fc) == vs.set2] - fc[ , colnames(fc) == vs.set1];
    if (fctype == "up") fc.pos <- which(fc_1vs2 > cutoff);
    if (fctype == "down") fc.pos <- which(fc_1vs2 < -cutoff);
    if (fctype == "all")  fc.pos <- which(fc_1vs2 > cutoff | fc_1vs2 < -cutoff);
  } else {
    fc_1vs2 <- fc[ , colnames(fc) == vs.set2] / fc[ , colnames(fc) == vs.set1];
    if (fctype == "up") fc.pos <- which(fc_1vs2 > cutoff);
    if (fctype == "down") fc.pos <- which(fc_1vs2 < 1/cutoff);
    if (fctype == "all")  fc.pos <- which(fc_1vs2 > cutoff | fc_1vs2 < 1/cutoff);
    fc_1vs2 <- log2(fc_1vs2);#190531
  }
  if(order) fc.pos <- fc.pos[order(abs(fc_1vs2[fc.pos]), decreasing = TRUE)]
  fc.pos
}



#siginificant proteins  vs.set2 vs vs.set1
##
#20190509 fix error when any group in groups have no replication
#20190513 fix t.test and pairwise t.test. do pairwise t.test when anova is TRUE.
#20190531 add padjust in anovaP and student t.test use Padj as p.adjust methods.
changedID <- function(relative_value, group, vs.set2, vs.set1 = "WT",
                      rank = "none",
                      anova = TRUE, anova.cutoff = 0.05,
                      T.cutoff = 0.05, Padj = "fdr",
                      cutoff = 1.5, datatype = c("none","log2"), fctype = "all",
                      ...){
  rank <- match.arg(rank, c("none", "foldchange", "anova", "t.test"));
  Padj <- match.arg(Padj, p.adjust.methods);
  if (rank == "foldchange") order = TRUE else order = FALSE;
  vs.name <- c(paste(vs.set1, "vs", vs.set2), paste(vs.set2, "vs", vs.set1));
  if (Padj != "none") vs.name <- paste(vs.name, Padj);

  if (ncol(relative_value) != length(group))
    stop ("Group length and data col length is not the same.");
  if (any(table(group) == 1)) {
    delgroup <- names(table(group))[table(group) == 1];
    delpos <- which(group %in% delgroup);
    warning(paste("Group", delgroup, "is deleted because it have no replication. "));
    group <- group[-delpos];
    relative_value <- relative_value[,-delpos];
  }
  if (length(table(group)) < 2)
    stop ("It should be at least two groups to do anova analysis.");
  vsgroup <- any(group == vs.set2, na.rm = TRUE) & any(group == vs.set1, na.rm = TRUE)
  if (!vsgroup) stop ("no vs.set1 or vs.set2 in  group.")
  sig.pos = NULL; pos_1vs2 = NULL;
  if (anova) {
    Anova <- anova_p(relative_value, group);
    anovapadj <- p.adjust(Anova[Anova != 2], method = Padj); #190531
    Anova[Anova != 2] <- anovapadj#190531
    sig.pos <- which(Anova < anova.cutoff);
    if (rank == "anova") sig.pos <- sig.pos[order(Anova[sig.pos])];
    if (is.numeric(T.cutoff)) {
      t.test <- multi.t.test(relative_value, group,
                             sig = 1, geneAdj = Padj, ...);
      vs.name <- vs.name[vs.name %in% colnames(t.test)];
      if (length(vs.name) != 1)
        stop ("vs.set1 and vs.set2 name is uncorrect");
      p_1vs2 <- as.numeric(t.test[ , vs.name]);
      ID_1vs2 <- rownames(t.test)[which(p_1vs2 < T.cutoff)]; #190531
      pos_1vs2 <- which(rownames(relative_value) %in% ID_1vs2)#190531
      if (rank == "t.test") pos_1vs2 <- pos_1vs2[order(p_1vs2[pos_1vs2])];
    }
  }
  if (!anova & is.numeric(T.cutoff)) {
    Ttestx <- which(group == vs.set1)
    Ttesty <- which(group == vs.set2)
    t.test <- NULL;
    for (i in 1:nrow(relative_value)) {
      TX <- as.numeric(relative_value[i,Ttestx]);
      TY <- as.numeric(relative_value[i,Ttesty]);
      samevalue <- c(length(table(TX)) > 1, length(table(TY)) > 1)
      if (sum(!is.na(TX)) > 2 & sum(!is.na(TY)) > 2 & any(samevalue)){
        Ttest <- t.test(TX, TY);
        t.test <- c(t.test, Ttest$p.value)
      } else t.test <- c(t.test, 2)
    }
    p_1vs2 <- p.adjust(t.test[t.test != 2], method = Padj);#190531
    t.test[t.test != 2] <- p_1vs2;#190531
    p_1vs2 <- t.test;#190531
    pos_1vs2 <- which(p_1vs2 < T.cutoff);
    if (rank == "t.test") pos_1vs2 <- pos_1vs2[order(p_1vs2[pos_1vs2])];
  }
  fc <- groupmean(relative_value, group, name = FALSE, ...);
  fc_1vs2 <- fc.pos(fc, vs.set2 = vs.set2, vs.set1 = vs.set1,
                    cutoff = cutoff, datatype = datatype,
                    fctype = fctype, order = order);
  if (is.null(sig.pos) & is.null(pos_1vs2))
    pos <- fc_1vs2 else if (is.null(sig.pos) & rank == "foldchange")
      pos <- fc_1vs2[fc_1vs2 %in% pos_1vs2] else if (is.null(sig.pos))
        pos <- pos_1vs2[pos_1vs2 %in% fc_1vs2] else if (is.null(pos_1vs2) &
                                                        rank == "foldchange")
          pos <- fc_1vs2[fc_1vs2 %in% sig.pos] else if (is.null(pos_1vs2))
            pos <- sig.pos[sig.pos %in% fc_1vs2] else if (rank == "foldchange") {
              pos <- fc_1vs2[fc_1vs2 %in% pos_1vs2];
              pos <- pos[pos %in% sig.pos];
            } else if (rank == "anova") {
              pos <- sig.pos[sig.pos %in% pos_1vs2];
              pos <- pos[pos %in% fc_1vs2];
  } else {
    pos <- pos_1vs2[pos_1vs2 %in% sig.pos];
    pos <- pos[pos %in% fc_1vs2];
  }
  rownames(relative_value)[pos]
}


#v0.2.1
#190605
#extract stat inf
dataStatInf <- function(prodata,  group, intensity = "intensity", Egrp = NULL, Cgrp = "ctl",
                        meanmethod = "mean", datatype = c("none", "log2"),
                        anova = TRUE, T.test = c("pairwise","student", "none"),
                        Aadj = "none", Tadj = "none", cutoff = FALSE,
                        ... ) {
  if (!intensity %in% names(prodata))
    stop(paste("intenstity is error. No",
               intensity, "data.frame in prodata."))
  data <- prodata[[intensity]];
  rownames(data) <- prodata$inf$ori.ID
  StudentT <- function(relative_value, group, Padj = p.adjust.methods) {
    Ttestx <- which(group == levels(as.factor(group))[1])
    Ttesty <- which(group == levels(as.factor(group))[2])
    ttest <- NULL;
    for (i in 1:nrow(relative_value)) {
      TX <- as.numeric(relative_value[i,Ttestx]);
      TY <- as.numeric(relative_value[i,Ttesty]);
      samevalue <- c(length(table(TX)) > 1, length(table(TY)) > 1)
      if (sum(!is.na(TX)) > 2 & sum(!is.na(TY)) > 2 & any(samevalue)){
        Ttest <- t.test(TX, TY);
        ttest <- c(ttest, Ttest$p.value)
      } else ttest <- c(ttest, 2)
    }
    ttest[ttest == 2] <- NA;
    p <- p.adjust(ttest, method = Padj);
    p
  }
  datatype <- match.arg(datatype, c("none", "log2"));
  meanmethod <- match.arg(meanmethod, c("mean", "median"));
  T.test <- match.arg(T.test, c("pairwise","student", "none"));
  Aadj <- match.arg(Aadj, p.adjust.methods);
  Tadj <- match.arg(Tadj, p.adjust.methods);
  if (isTRUE(cutoff)) cutoff = 0.05;
  if (length(cutoff) != 1 || (!is.numeric(cutoff) && !isFALSE(cutoff)) )
    stop ("wrong cutoff");
  groups <- levels(as.factor(group));
  if (!Cgrp %in% groups)
    stop ("wrong Cgrp or group");
  if (T.test == "student") {
    if ( is.null(Egrp))
      Egrp <- groups[groups != Cgrp] else if (!Egrp %in% groups)
        stop ("wrong Egrp");
    if (length(Egrp) != 1)
      stop("Egrp should be assigned when T.test is student.");
  }

  gmean <- groupmean(data, group, method = meanmethod, name = FALSE);
  if (datatype == "none")
    gmean <- gmean/gmean[, Cgrp];
  if (datatype == "log2")
    gmean <- gmean - gmean[, Cgrp];
  gmean <- gmean[ ,!colnames(gmean) %in% Cgrp];
  if (is.data.frame(gmean)) gmean <- as.matrix(gmean);
  if (anova) {
    Anova <- anova_p(data, group);
    Anova[Anova == 2] <- NA;
    anovapadj <- p.adjust(Anova, method = Aadj);
    Anova <- cbind(anovaP =Anova,
                   anovaPadj = if(Aadj != "none") anovapadj);
    if (cutoff) Anova[Anova > cutoff] <- NA;
  }
  if (T.test == "pairwise") {
    ttestP <- multi.t.test(data, group, sig = 1, geneAdj = Tadj, ...);
    P <- matrix(data = NA, nrow = nrow(data), ncol = ncol(ttestP),
                dimnames = list(rownames(data), colnames(ttestP)));
    P[rownames(P) %in% rownames(ttestP), ] <- ttestP;
    ttestP <- P;
    if (cutoff) ttestP[ttestP > cutoff] <- NA;
  }
  if (T.test == "student") {
    Ttestx <- which(group == Cgrp);
    Ttesty <- which(group == Egrp);
    dataforTtest <- data[, c(Ttestx, Ttesty)];
    ttestP <- StudentT(dataforTtest, group, Padj = "none");
    ttestPadj <- p.adjust(ttestP, method = Tadj);
    ttestP <- cbind(ttestP =ttestP,
                    ttestPadj = if(Tadj != "none") ttestPadj);
    if (cutoff) ttestP[ttestP > cutoff] <- NA;
  }
  final <- cbind(mean = gmean,
                 anovaP = if(anova) Anova,
                 ttestP = if(T.test != "none") ttestP);
  final <- data.frame(prodata$inf ,final, stringsAsFactors = FALSE);
  final
}
