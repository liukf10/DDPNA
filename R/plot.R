
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
        module <- c(module, rep(paste0("mod", colnames(Z)[i]), nrow(Z)));
        x <- c(x, X);
      }
    } else {
      y <- c(y, Z[,i]);
      module <- c(module, rep(paste0("mod", colnames(Z)[i]), nrow(Z)));
      x <- c(x, X);
    }
  }
  Zscore <- data.frame(number = x, Zs = y, module, stringsAsFactors = FALSE);
  if(length(Zscore) == 0) {
    warning("no module satisfy the rule.")
    plot = FALSE}
  #library(ggplot2)
  #library(ggalt)
  if(plot){
    pic <- ggplot(Zscore, aes(Zscore$number, Zscore$Zs, group = module, color = module)) +
      geom_point(size = 1, shape = 16) +
      geom_xspline(spline_shape = -0.4, size = 1) +
      geom_hline(yintercept = -log10(p), col = "red") +
      theme_classic() +
      theme(text = element_text(size = 17));
    if (is.character(filename) & length(filename) == 1) {
      if (!dir.exists("plot")) dir.create("plot");
      filename = paste0("plot/module enrichment ", filename, ".", filetype);
      ggsave(pic, filename = filename, ...);
    } else print(pic)
#    if (is.character(filename) & length(filename) == 1) {
#      if (!dir.exists("plot")) dir.create("plot");
#      pdf(paste0("plot/module enrichment ", filename, ".pdf"))}
#    print(pic)
#    if (is.character(filename) & length(filename) == 1) dev.off()

  }
  Zscore
}





#creat Rd file need devtools, roxygen2
#library(devtools)
#devtools::install_github("klutometis/roxygen")
#library(roxygen2)

#disease drived human-mouse difference protein associated network analysis
#setwd("C:/Rpackage")

