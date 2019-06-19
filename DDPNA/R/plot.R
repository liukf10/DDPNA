
#need ggplot2 and ggalt function
FCSenrichplot <- function(FCSenrich, count = 1, p = 0.05, filter = "p",
                          plot = TRUE, filename = NULL, filetype = "pdf", ...){
  filter <- match.arg(filter, c("p", "p.adj", "none"));
  if(filter == "p") pvalue <- FCSenrich$p;
  if(filter == "p.adj") pvalue <- FCSenrich$p.adj;
  Counts <- FCSenrich$Counts;
  Z <- FCSenrich$Z.score;
  pickmod <- 1:ncol(Z);
  pickmod <- pickmod[Counts[nrow(Counts), ] > count];
  X<-as.numeric(rownames(Z));
  x = NULL; y = NULL; module = NULL;
  for(i in pickmod){
    if(filter != "none"){
      if(any(pvalue[ ,i] < p)){
        y <- c(y, Z[ ,i]);
        module <- c(module, rep(colnames(Z)[i], nrow(Z))); #190603
        x <- c(x, X);
      }
    } else {
      y <- c(y, Z[,i]);
      module <- c(module, rep(colnames(Z)[i], nrow(Z))); #190603
      x <- c(x, X);
    }
  }
  Zscore <- data.frame(number = x, Zs = y, module, stringsAsFactors = FALSE);
  if(length(Zscore) == 0) {
    warning("No module satisfy the rule.")
    plot = FALSE}
  #library(ggplot2)
  #library(ggalt)
  if(plot){
    number <- Zscore$number;
    Z <- Zscore$Zs
    pic <- ggplot(Zscore, aes(number, Z, group = module, color = module)) +
      geom_point(size = 1, shape = 16) +
      geom_xspline(spline_shape = -0.4, size = 1) +
      geom_hline(yintercept = -log10(p), col = "red") +
      theme_classic() +
      theme(text = element_text(size = 17));
    if (is.character(filename) & length(filename) == 1) {
      if (!dir.exists("plot")) dir.create("plot");
      filename = paste0("plot/module enrichment ", filename, ".", filetype);
      ggsave(pic, filename = filename, ...);
      unlink(pic)
    } else print(pic)
#    if (is.character(filename) & length(filename) == 1) {
#      if (!dir.exists("plot")) dir.create("plot");
#      pdf(paste0("plot/module enrichment ", filename, ".pdf"))}
#    print(pic)
#    if (is.character(filename) & length(filename) == 1) dev.off()

  }
  Zscore
}



DEP_Mod_HeatMap <- function(DEP_Mod, filter = c("p","p.adj"),
                            cutoff = 0.05, filename = NULL, ...) {
  filter <- match.arg(filter, c("p","p.adj"));
  if (!is.list(DEP_Mod)) stop("DEP_Mod is not the list class.")
  x <- matrix(data = NA, nrow = length(DEP_Mod[[1]][[1]]), ncol = length(DEP_Mod))
  x <- as.data.frame(x)
  rownames(x) <- rownames(DEP_Mod[[1]])
  if (is.null(names(DEP_Mod))) names(DEP_Mod) <- 1:length(DEP_Mod);
  colnames(x) <- names(DEP_Mod)
  x1 = x2 = x3 = x4 = x5 = x
  for (i in 1:length(DEP_Mod) ) {
    if (all(c("precent","Counts","module.size","module.size","p","p.adj") %in%
            names(DEP_Mod[[i]]) ) ) {
      if(length(DEP_Mod[[i]][["precent"]]) != nrow(x1))
        stop ( paste(names(DEP_Mod)[i], "is not the same module number with others.") );
      x1[ , i ] = DEP_Mod[[i]][["precent"]];
      x2[ , i ] = DEP_Mod[[i]][["Counts"]];
      x3[ , i ] = DEP_Mod[[i]][["module.size"]];
      x4[ , i ] = DEP_Mod[[i]][["p"]];
      x5[ , i ] = DEP_Mod[[i]][["p.adj"]];
    } else stop( paste(names(DEP_Mod)[i], "is not the correct data."))

  }
  if (filter == "p")
    connect <- list(precent = x1, Counts =x2,
                    module.size = x3, p = x4)
  if (filter == "p.adj")
    connect <- list(precent = x1, Counts =x2,
                    module.size = x3, p = x5)

  p <- connect$p
  p <- signif(p, digits = 2)
  p <- format(p, scientific = TRUE)
  p[connect$p > cutoff] <- ""
  Counts <- connect$Counts
  Counts[connect$p > cutoff] <- ""
  #colnames(p) <- gsub("(.*).p","\\1",colnames(p))
  #rownames(p) <- paste0("M",rownames(p))
  textMatrix = paste(as.matrix(Counts),"\n(",as.matrix(p),")",sep="")
  dim(textMatrix) = dim(Counts)
  textMatrix[textMatrix=="\n()"]<-""
  precent <- connect$precent
  #precent <- round(precent, digits = 3)
  ratio <- colSums(connect$Counts) / colSums(connect$module.size)
  #precent2 <- apply(precent, 1, function(x) x/ratio/100)
  enrichFold <- t(t(precent)/ratio/100)
  if (!is.null(filename)) pdf(paste0("plot/",filename,".pdf"))
  WGCNA::labeledHeatmap(Matrix = enrichFold, xLabels = colnames(p),
                        yLabels = rownames(p),
                        cex.lab = 1, colorLabels = TRUE,
                        colors = WGCNA::blueWhiteRed(100)[51:100],
                        textMatrix = textMatrix,
                        setStdMargins = FALSE,
                        cex.text = 1, ...)
  if (!is.null(filename)) dev.off()
  list(enrichFold = enrichFold,
       textMatrix = textMatrix)
}

#creat Rd file need devtools, roxygen2
#library(devtools)
#devtools::install_github("klutometis/roxygen")
#library(roxygen2)

#disease drived human-mouse difference protein associated network analysis
#setwd("C:/Rpackage")

