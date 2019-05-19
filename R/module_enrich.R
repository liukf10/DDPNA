#Module enrichment analysis
#20190512 fix the datainf duplicationed id
Module_Enrich <- function(module, classifiedID, enrichtype = "FCS",
                          coln = "new.ID", datainf = NULL, p.adj.method = "BH")
  {
  enrichtype <- match.arg(enrichtype, c("FCS","ORA"));
  mod.name = NULL; pvalue = NULL; precentage = NULL;
  fcnum = NULL; modnum = NULL; Z = NULL;
  im <- length(table(module$moduleNum));
  classifiedID <- classifiedID[!duplicated(classifiedID)] #20190512

  for ( i in 0 : (im - 1)) {
    if (is.null(datainf)) {
      m <- module[module$moduleNum == i, coln];
      Total <- nrow(module);
      } else {
        datainf <- datainf[!duplicated(datainf)]#20190512
        m <- moduleID(datainf, module, i, coln = coln);
        Total <- length(datainf);
      }
    if (enrichtype == "FCS") {
      mod_enrich <- single_mod_enrichplot(module, i, classifiedID, coln = coln,
                            datainf = datainf,plot = FALSE);
      q <- mod_enrich$hits;
      p <- mod_enrich$p;
      z <- mod_enrich$Z.score;
      modnum <- c(modnum, length(m));
      mod.name <- c(mod.name, i);
      if(is.null(fcnum)) fcnum <- q
      else fcnum <- data.frame(fcnum, q, stringsAsFactors = FALSE);
      if(is.null(precentage)) precentage <- q/length(m)*100
      else precentage <- data.frame(precentage, q/length(m)*100,
                                    stringsAsFactors = FALSE);
      if(is.null(pvalue)) pvalue <- p
      else pvalue <- data.frame(pvalue, p, stringsAsFactors = FALSE);
      if(is.null(Z)) Z <- z
      else Z <- data.frame(Z, z, stringsAsFactors = FALSE);
    }
    if (enrichtype == "ORA") {
      q <- sum(classifiedID %in% m);
      n <- Total - length(m);
      k <- length(classifiedID);
      p <- phyper(q-1, length(m), n, k, lower.tail = FALSE);
      modnum <- c(modnum, length(m));
      fcnum <- c(fcnum, q);
      precentage <- c(precentage, q/length(m)*100);
      mod.name <- c(mod.name, i);
      pvalue <- c(pvalue, p);
      Z <- c(Z, -log10(p));
    }
  }
  if(enrichtype == "ORA") p.adj <- p.adjust(pvalue, method = p.adj.method)
  else{
    colnames(fcnum) <- mod.name;
    colnames(precentage) <- mod.name;
    colnames(pvalue) <- mod.name;
    colnames(Z) <- mod.name;
    p.adj <- pvalue;
    for (i in 1:length(classifiedID)) {
      p.adj[i, ] <- p.adjust(pvalue[i, ], method = p.adj.method);
    }
  }
  list(Counts = fcnum, module.size = modnum, module.name = mod.name,
       precent = precentage, p = pvalue, p.adj = p.adj, Z.score = Z)
}

#a set of proteins which is enrichment in single module.
#classfiedID is ordered ID, which is gotten from changedID function.
single_mod_enrichplot <- function(module, Mod_Nam, classifiedID,
                                  coln = "new.ID", datainf = NULL,
                                  plot = TRUE, filename = NULL, ...){
  pvalue = NULL; hits = NULL; Z = NULL;
  if (is.null(datainf)) {
    m <- module[module$moduleNum == Mod_Nam, coln];
    Total <- nrow(module);}
  else {
    m <- moduleID(datainf, module, Mod_Nam, coln = coln);
    Total <- length(datainf);}
  for (i in 1:length(classifiedID)){
    q <- sum(classifiedID[1:i] %in% m);
    n <- Total-length(m);
    k <- i;
    p <- phyper(q-1, length(m), n, k, lower.tail = FALSE)
    hits <- c(hits, q);
    pvalue <- c(pvalue, p);
    Z <- c(Z, -log10(p));
  }
  x <- 1:length(classifiedID);
  if (plot) {
    #spline interpolation
    sp1 = spline(x, Z, n=3*length(classifiedID), ...);
    if (is.character(filename) & length(filename) == 1){
      if (!dir.exists("plot")) dir.create("plot");
      pdf(paste0("plot/ module", Mod_Nam, " enrichment ", filename, ".pdf"))
    }
    plot(x, Z, type = "p", pch=20, ylab="Z score", xlab = "",
         main = paste("module", Mod_Nam, "enrichment analysis"));
    lines(sp1);
    abline(h = 2, col = "red");
    abline(h = -log10(0.05), col = "blue");
    if (is.character(filename) & length(filename) == 1) dev.off()
    }
  list(hits = hits, p = pvalue, Z.score = Z)
}




#extract intersection ID between dataset and module
moduleID <- function(inf, module, num, coln = "new.ID"){
  m <- module[module$moduleNum == num, "new.ID"]
  #m <- gsub("(tr|sp)\\|(.*)\\|.*","\\2",m)
  ori <- length(m)
  m <- inf[inf %in% m]
  #if (ori != length(m)) {warning (paste("old# is ",ori, ",now is",length(m)))}
  m
  }
