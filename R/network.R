#PMFG
getmoduleHub <- function(data, module, mod_num, coln = "new.ID",
                         cor.sig = 0.05, cor.r = 0,
                         adjustp = TRUE, hub.p = 0.05){
  mod.name <- module[module$moduleNum == mod_num, coln];
  #mod <- data[ , colnames(data) %in% mod.name]
  mod <- data[rownames(data) %in% mod.name, ];
  mod <- t(mod);

  #library("Hmisc")
  cor_matrix <- Hmisc::rcorr(mod, type = "pearson");
  #
  cor_matrix_sig <- .sig_cor(cor_matrix, p = cor.sig, r = cor.r, type = 1);
  cor_order <- .orderCor(abs(cor_matrix_sig$r), cor_matrix_sig$P);

  # make igraph objects
  #library("igraph")
  g <- graph.data.frame(cor_order[ ,1:3], directed = FALSE);
  #plot(g)

  #library(MEGENA)
  el <- calculate.PFN(cor_order[,1:3], doPar = FALSE, num.cores = 1);
  gg <- graph.data.frame(el, directed = FALSE);
  #plot(gg)
  subnet <- induced.subgraph(gg, mod.name);
  Stat <- get.DegreeHubStatistic(subnet, n.perm = 100);
  #Stat2 <- .get_HubStatistic(subnet, n.perm = 100)
  if(adjustp) hubgene <- Stat$gene[Stat$pvalue.BH < hub.p]
  else hubgene <- Stat$gene[Stat$pvalue < hub.p];
  hubgene <- as.character(hubgene);
  Stat$gene <- as.character(Stat$gene);
  list(hub = hubgene, degreeStat = Stat, graph = g, PMFG = gg)
}










##extract correlation
#P 0.05,0.01,0.001  r  0.1,0.3,0.5,0.8
.sig_cor <- function(cor_matrix, p = 0.05, r = NA, type = 1){
  #lower.tri
  cor_matrix$P[upper.tri(cor_matrix$P,diag = TRUE)] <- NA;
  cor_matrix$r[upper.tri(cor_matrix$P,diag = TRUE)] <- NA;
  #order
  #rank_r<-matrix(rank(1-abs(cor_matrix$r)),nrow=nrow(cor_matrix$r));
  #rank_p<-matrix(rank(abs(cor_matrix$P)),nrow=nrow(cor_matrix$P));
  #rank_r[upper.tri(rank_r,diag = TRUE)]<-NA;
  #rank_p[upper.tri(rank_p,diag = TRUE)]<-NA;
  if (!is.na(p)) {
    #p value less than p
    logi <- cor_matrix$P > p;
    cor_matrix$P[logi] <- NA;
    cor_matrix$r[logi] <- NA;
  }
  if (!is.na(r)) {
    logi2 <- abs(cor_matrix$r) < r;
    #r value >r
    cor_matrix$P[logi2] <- NA;
    cor_matrix$r[logi2] <- NA;
  }
  if (!is.na(p) | !is.na(r)) {
    #
    posit<-which(!is.na(cor_matrix$P), arr.ind=TRUE);
    if(type==2){
      #
      column<-sort(unique(c(posit[,1],posit[,2])));
      cor_matrix$P<-cor_matrix$P[column,column];
      cor_matrix$r<-cor_matrix$r[column,column];
    }
    else if (type ==1 ) {
      #
      roworder <- sort(unique(posit[ ,1]));
      column <- sort(unique(posit[ ,2]));
      cor_matrix$P <- cor_matrix$P[roworder, column];
      cor_matrix$r <- cor_matrix$r[roworder, column];
    }
  }
  cor_matrix
}

# orderCorr
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
.orderCor <- function(cormat, pmat) {
  postion <- which(!is.na(cormat), arr.ind = TRUE);
  pos <- which(!is.na(cormat));
  ro <- attr(cormat, "dimnames")[[1]][postion[ ,1]];
  co <- attr(cormat, "dimnames")[[2]][postion[ ,2]];
  r <- cormat[pos];
  p <- pmat[pos];
  flat <- data.frame(row = as.character(ro),
                     column = as.character(co),
                     cor = r,
                     p = p,
                     stringsAsFactors = FALSE);
  orderCor <- flat[order(abs(flat$cor), decreasing = TRUE), ];
  rownames(orderCor) <- order(abs(orderCor$cor), decreasing = TRUE);
  orderCor
}


