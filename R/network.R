#PMFG
#190604
#allow extract 2 or more Modules to as a big net
#the previous code will show the wrong result when extract 2 or more Modules
#190612 fix the error when PFG is not the whole mod.name

getmoduleHub <- function(data, module, mod_num, coln = "new.ID",
                         cor.sig = 0.05, cor.r = 0,
                         adjustp = TRUE, hub.p = 0.05){
  #mod.name <- module[module$moduleNum == mod_num, coln];
  mod.name <- module[module$moduleNum %in% mod_num, coln]; #190604
  if(length(mod.name) == 0) stop(paste("module information have no", mod_num))
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
  mod.name.net <- c(el[,1], el[,2]); #190612
  mod.name.net <- mod.name.net[!duplicated(mod.name.net)]; #190612
  subnet <- induced.subgraph(gg, mod.name.net); #190612
  #subnet <- induced.subgraph(gg, mod.name);
  Stat <- get.DegreeHubStatistic(subnet, n.perm = 100);
  #Stat2 <- .get_HubStatistic(subnet, n.perm = 100)
  if(adjustp) hubgene <- Stat$gene[Stat$pvalue.BH < hub.p]
  else hubgene <- Stat$gene[Stat$pvalue < hub.p];
  hubgene <- as.character(hubgene);
  Stat$gene <- as.character(Stat$gene);
  list(hub = hubgene, degreeStat = Stat, graph = g, PMFG = gg)
}



#v0.2.1 add
#rename the duplicated newID in moduleinf and renew the ID in DEPstat
rename_dupnewID <- function(DEPstat, moduleinf, DEPfromMod = FALSE ) {
  #find duplicated ID and add _1 after the ID name.
  #If the ID duplicated more than one time, then add _n. n means times.
  if (sum(colnames(DEPstat) %in%
          "new.ID") * sum(colnames(DEPstat) %in%
                          "ori.ID") != 1)
    stop("DEPstat have no column: ori.ID or new.ID.")
  if (sum(colnames(moduleinf) %in%
          "new.ID") * sum(colnames(DEPstat) %in%
                          "ori.ID") != 1)
    stop("moduleinf have no column: ori.ID or new.ID.")
  rename_newID <- function(data)
  {
    newID <- data$new.ID; dup <- which(duplicated(newID));
    dupID <- newID[dup]; i = 1;
    while (length(dup) != 0) {
      newID[dup] <- paste0(newID[dup], "_", i);
      dup <- which(duplicated(newID));
      i = i+1;
    }
    newID
  }
#  question <- function (...)
#  {
#    yes <- c("Yes", "Definitely", "Sure",
#             "I agree", "Absolutely")
#    no <- c("No", "Nope", "Not Sure")
#    cat(paste0(..., collapse = ""))
#    qs <- c(sample(yes, 1), sample(no, 1))
#    menu(qs) == 1
#  }
  dupID <- moduleinf$new.ID[which(duplicated(moduleinf$new.ID))];
  newID <- rename_newID(moduleinf);
  #DEPstat and module construction is the same datasets.
  if (DEPfromMod) {
    posLogical <- match(DEPstat$ori.ID, moduleinf$ori.ID);
    if ( !any(is.na(posLogical)) )
      DEPstat$new.ID <- newID[posLogical] else
        stop("Some DEPstat ID are not in the moduleinf ID list.")
        #{
        #cat("Some DEPstat ID are not in the moduleinf ID list.\n")
        #if (question("Do you want to reset DEPfromMod to FALSE?"))
        #  DEPfromMod = FALSE else stop("Please check the input and rerun!")
      #}
  }
  #DEPstat and module construction is the different datasets.
  if (!DEPfromMod) {
    DEPstat_intersect <- DEPstat[DEPstat$new.ID %in% moduleinf$new.ID, ];
    if (sum(dupID %in% DEPstat_intersect$new.ID) != 0) {
      pos <- which(DEPstat_intersect$new.ID %in% dupID)
      for (i in pos) {
        id <- newID[moduleinf$new.ID %in% DEPstat_intersect$new.ID[i]]
        idd <- rbind(DEPstat_intersect[i, ], DEPstat_intersect[i, ])
        while (length(id) > nrow(idd)) {
          idd <- rbind(idd, DEPstat_intersect[i, ])
        }
        idd$new.ID <- id;
        if(i == pos[1]) {
          if( i != 1) intersect <- rbind(DEPstat_intersect[1:(i-1), ],
                                         idd)
          if( i == 1) intersect <- idd
        } else {
          pre <- j+1; now <- i-1;
          if( pre <= now )
            intersect <- rbind(intersect,
                               DEPstat_intersect[(pre:now), ],
                               idd)
          if( pre > now )
            intersect <- rbind(intersect, idd)
        }
        j = i;
      }
      if (pos[length(pos)] != nrow(DEPstat_intersect))
        intersect <- rbind(intersect,
                           DEPstat_intersect[(pos[length(pos)]+1):nrow(DEPstat_intersect), ])
      DEPstat <- intersect
    } else DEPstat <- DEPstat_intersect
  }
  DEPstat
}



#v0.2.1 add
# get two or more IDsets interesection set and complementary set
# give each IDsets a color and use overlaped color as intersection set color
# The algorithm is Photoshop filter algorithm
DEPsets <- function(datalist,colors = c("red","green","blue")){
  if (length(datalist) < 2)
    stop("datalist should have at least 2 sets of DEPs.")
  VennInf <- VennDiagram::get.venn.partitions(datalist)
  RGB <- col2rgb(colors)
  if (length(datalist) > length(colors) )
    stop("colors number is less than datalist number.")

  gene.set <- vector("list", nrow(VennInf))
  Colr <- NULL;
  for (i in 1:nrow(VennInf)) {
    gene.set[[i]] <- VennInf$..values..[[i]];
    pos<- as.logical(VennInf[i,1:length(datalist)]);
    name <- colnames(VennInf)[1:length(datalist)][pos];
    name <- paste0(name,collapse = "_");
    names(gene.set)[i] <- name;
    if (sum(pos) == 1) {
      color = colors[1:length(datalist)][pos];
      #color = rgb(col2rgb(color)[1], col2rgb(color)[2], col2rgb(color)[3],
      #            maxColorValue = 255);
    }
    if (sum(pos) != 1) {
      oricol3 = RGB[, 1:length(datalist)][ ,pos];
      while (is.matrix(oricol3)) {
        oricol3[,1] = round(oricol3[,1] + oricol3[,2] - (oricol3[,1]*oricol3[,2])/255);
        oricol3 <- oricol3[,-2];
      }
      color = rgb(oricol3[1], oricol3[2], oricol3[3], maxColorValue = 255);
    }
    Colr <- c(Colr, color);
  }
  list(gene.set = gene.set,
       color.code = Colr)
}


#v0.2.1 add
# remove hubs which not in the IDsets and replot the PFG network
#190611 add value BranchCut, a logical value indicating whether remove proteins
#which have no connection to DEPs.
# 190611 fix a bug that the proteins only choose el$from when reconstructNet.
DEP_Mod_net_plot <- function(ModNet, IDsets = NULL,
                             data = NULL, module = NULL,
                             plot = TRUE,
                             filename = NULL, filetype = "pdf",
                             OnlyPlotLast = TRUE,
                             BranchCut = TRUE,
                             reconstructNet = TRUE, iteration = Inf,
                             label.hubs.only = TRUE,
                             node.default.color = "grey",
                             hubLabel.col = "black",...)
{
  continue = TRUE;
  if (OnlyPlotLast) plot = FALSE;
  if (is.null(IDsets)) continue = FALSE;
  if (!continue) {
    test <- try(plot_subgraph(module = ModNet$degreeStat$gene,
                              hub = ModNet$hub, PFN = ModNet$PMFG,
                              node.default.color = node.default.color,
                              show.legend = TRUE,
                              label.hubs.only = label.hubs.only,
                              hubLabel.col = hubLabel.col,
                              hubLabel.sizeProp = 1,
                              show.topn.hubs = 10,
                              node.sizeProp = 13, label.sizeProp = 13,
                              label.scaleFactor = 10,
                              layout = "kamada.kawai"), silent = TRUE)
    if (class(test) == "try-error")
      stop("plot_subgraph is error. Please check ModNet.")
    netgene = ModNet$degreeStat$gene;
    hub = ModNet$hub;
    PMFG = ModNet$PMFG;
    #cat("No IDsets, it will only plot the module network.");
    warning("No IDsets, it will only plot the module network.")
    if (plot | OnlyPlotLast) {
      pic <- plot_subgraph(module = ModNet$degreeStat$gene,
                           hub = ModNet$hub, PFN = ModNet$PMFG,
                           node.default.color = node.default.color,
                           show.legend = TRUE,
                           label.hubs.only = label.hubs.only,
                           hubLabel.col = hubLabel.col, hubLabel.sizeProp = 1,
                           show.topn.hubs = 10,
                           node.sizeProp = 13, label.sizeProp = 13,
                           label.scaleFactor = 10,
                           layout = "kamada.kawai")
      if (is.character(filename) & length(filename) == 1) {
        if (!dir.exists("plot"))
          dir.create("plot")
        filename2 = paste0("plot/DEP_Mod_net", filename, "_ori",
                           ".", filetype)
        ggsave(pic$pnet, filename = filename2, ...)
        unlink(pic$pnet)
      }
      else print(pic$pnet)
    }
  }
  if (continue) {
    if (reconstructNet & (is.null(data) | is.null(module)))
      stop ("data and module should be defined when reconstructNet is TRUE.");
    test <- try(plot_subgraph(module = ModNet$degreeStat$gene,
                              hub = ModNet$hub, PFN = ModNet$PMFG,
                              node.default.color = node.default.color,
                              gene.set = IDsets$gene.set,
                              color.code = IDsets$color.code,
                              show.legend = TRUE,
                              label.hubs.only = label.hubs.only,
                              hubLabel.col = hubLabel.col,
                              hubLabel.sizeProp = 1,
                              show.topn.hubs = 10,
                              node.sizeProp = 13, label.sizeProp = 13,
                              label.scaleFactor = 10,
                              layout = "kamada.kawai"), silent = TRUE)
    if (class(test) == "try-error")
      stop("plot_subgraph is error. Please check ModNet and IDsets.")
    if (plot) {
      pic <- plot_subgraph(module = ModNet$degreeStat$gene,
                           hub = ModNet$hub, PFN = ModNet$PMFG,
                           node.default.color = node.default.color,
                           gene.set = IDsets$gene.set,
                           color.code = IDsets$color.code, show.legend = TRUE,
                           label.hubs.only = label.hubs.only,
                           hubLabel.col = hubLabel.col, hubLabel.sizeProp = 1,
                           show.topn.hubs = 10,
                           node.sizeProp = 13, label.sizeProp = 13,
                           label.scaleFactor = 10,
                           layout = "kamada.kawai")
      if (is.character(filename) & length(filename) == 1) {
        if (!dir.exists("plot"))
          dir.create("plot")
        filename2 = paste0("plot/DEP_Mod_net", filename, "_ori",
                           ".", filetype)
        ggsave(pic$pnet, filename = filename2, ...)
        unlink(pic$pnet)
      }
      else print(pic$pnet)
    }
    #not reconstrunction PFG and replot
    #IDs in DEPsets which have no any connection will not showed in picture
    if (!reconstructNet)
    {
      node.features<- plot_subgraph(module = ModNet$degreeStat$gene,
                                    hub = ModNet$hub, PFN = ModNet$PMFG,
                                    node.default.color = "grey",
                                    gene.set = IDsets$gene.set,
                                    color.code = IDsets$color.code,
                                    show.legend = TRUE,
                                    label.hubs.only = TRUE,
                                    hubLabel.col = "black", hubLabel.sizeProp = 1,
                                    show.topn.hubs = 10,
                                    node.sizeProp = 13, label.sizeProp = 13,
                                    label.scaleFactor = 10,
                                    layout = "kamada.kawai")$node.features
      removedHub <- as.character(node.features$id[node.features$node.shape == "hub" &
                                                    node.features$node.stat == "NA"])

      el <- as_data_frame(ModNet$PMFG);
      el2 <- el[!(el[,1] %in% removedHub | el[,2]  %in% removedHub), ];
      if (BranchCut) {
        markedID<- as.character(node.features$id[node.features$node.stat !="NA"])
        el2 <- el2[el2$from %in% markedID | el2$to %in% markedID, ]
      }
      PMFG <- graph.data.frame(el2, directed = FALSE);
      hubs <- ModNet$hub[!ModNet$hub %in% removedHub];
      modulegene <- c(as.character(el2[,1]),as.character(el2[,2]))
      modulegene <- as.factor(modulegene);
      modulegene <- levels(modulegene);
      if (plot | OnlyPlotLast) {
        pic <- plot_subgraph(module = modulegene,
                             hub = hubs, PFN = PMFG,
                             node.default.color = node.default.color,
                             gene.set = IDsets$gene.set,
                             color.code = IDsets$color.code, show.legend = TRUE,
                             label.hubs.only = label.hubs.only,
                             hubLabel.col = hubLabel.col,
                             hubLabel.sizeProp = 1,
                             show.topn.hubs = 10,
                             node.sizeProp = 13, label.sizeProp = 13,
                             label.scaleFactor = 10,
                             layout = "kamada.kawai")
        if (is.character(filename) & length(filename) == 1) {
          if (!dir.exists("plot"))
            dir.create("plot")
          filename2 = paste0("plot/DEP_Mod_net", filename, "_new",
                             ".", filetype)
          ggsave(pic$pnet, filename = filename2, ...)
          unlink(pic$pnet)
        }
        else print(pic$pnet)
      }
      if (length(removedHub) == 0 )
        warning("all Hubs ID belong IDsets.")

    }
    #reconstrunction PFG and replot
    #remove unDEPsets hubs ID and related unDEPsets ID.
    #the function can run a cycle until no unDEPsets hubs.
    if (reconstructNet)
    {
      node.features<- plot_subgraph(module = ModNet$degreeStat$gene,
                                    hub = ModNet$hub, PFN = ModNet$PMFG,
                                    node.default.color = node.default.color,
                                    gene.set = IDsets$gene.set,
                                    color.code = IDsets$color.code, show.legend = TRUE,
                                    label.hubs.only = label.hubs.only,
                                    hubLabel.col = hubLabel.col,
                                    hubLabel.sizeProp = 1,
                                    show.topn.hubs = 10,
                                    node.sizeProp = 13, label.sizeProp = 13,
                                    label.scaleFactor = 10,
                                    layout = "kamada.kawai")$node.features;
      if (sum(node.features$node.shape == "hub" & node.features$node.stat != "NA") == 0)
        stop("No hub is in the IDsets.")
      markedID <- as.character(node.features$id[node.features$node.stat != "NA"])
      removedHub <- as.character(node.features$id[node.features$node.shape == "hub" &
                                                    node.features$node.stat == "NA"]);
      iter = 1;
      if (FALSE) {
        while (length(removedHub) > 0) {
          if (iter > iteration ) break;
          el <- as_data_frame(ModNet$PMFG);
          #el2 <- el[el[,1] %in% markedID | el[,2]  %in% markedID, ]
          el <- el[!(el[,1] %in% removedHub | el[,2]  %in% removedHub), ];
          if (BranchCut) el <- el[el$from %in% markedID | el$to %in% markedID, ];
          modulegene <- c(as.character(el$from),as.character(el$to), markedID);
          modulegene <- as.factor(modulegene);
          modulegene <- levels(modulegene);
          moduleinf <- module;
          moduleinf$moduleNum[moduleinf$new.ID %in% modulegene] <- "a";
          ModNet <- try(getmoduleHub(data, moduleinf, "a",
                                     coln = "new.ID",
                                     adjustp = FALSE),silent = TRUE);
          if (class(ModNet) == "try-error")
            stop ("data or module have some problem and cannot run getmoduleHub.")
          node.features<- plot_subgraph(module = ModNet$degreeStat$gene,
                                        hub = ModNet$hub, PFN = ModNet$PMFG,
                                        node.default.color = "grey",
                                        gene.set = IDsets$gene.set,
                                        color.code = IDsets$color.code, show.legend = TRUE,
                                        label.hubs.only = TRUE,
                                        hubLabel.col = "black", hubLabel.sizeProp = 1,
                                        show.topn.hubs = 10,
                                        node.sizeProp = 13, label.sizeProp = 13,
                                        label.scaleFactor = 10,
                                        layout = "kamada.kawai")$node.features;
          removedHub <- as.character(node.features$id[node.features$node.shape == "hub" &
                                                        node.features$node.stat == "NA"]);
          if (plot) {
            pic <- plot_subgraph(module = ModNet$degreeStat$gene,
                                 hub = ModNet$hub, PFN = ModNet$PMFG,
                                 node.default.color = node.default.color,
                                 gene.set = IDsets$gene.set,
                                 color.code = IDsets$color.code, show.legend = TRUE,
                                 label.hubs.only = label.hubs.only,
                                 hubLabel.col = hubLabel.col,
                                 hubLabel.sizeProp = 1,
                                 show.topn.hubs = 10,
                                 node.sizeProp = 13, label.sizeProp = 13,
                                 label.scaleFactor = 10,
                                 layout = "kamada.kawai");
            if (is.character(filename) & length(filename) == 1) {
              if (!dir.exists("plot"))
                dir.create("plot")
              filename2 = paste0("plot/DEP_Mod_net", filename, "_",iter,
                                 ".", filetype)
              ggsave(pic$pnet, filename = filename2, ...)
              unlink(pic$pnet)
            }
            else print(pic$pnet)
          }
          iter <- iter + 1;
        }
      }
      if (TRUE) {
        while (iter < 1000) {
          if (iter > iteration ) break;
          el <- as_data_frame(ModNet$PMFG);
          previous <- nrow(el);
          if (length(removedHub) > 0) {
            el <- el[!(el[,1] %in% removedHub | el[,2]  %in% removedHub), ];
            after_rmhub <- nrow(el);
          }
          if (BranchCut) el <- el[el$from %in% markedID | el$to %in% markedID, ];
          after_all <- nrow(el);
          if (previous == after_all) break;
          if (iter == 1) modulegene2 = ModNet$degreeStat$gene;
          if (iter != 1) modulegene2 = modulegene;
          modulegene <- c(as.character(el$from),as.character(el$to), markedID);
          modulegene <- as.factor(modulegene);
          modulegene <- levels(modulegene);
          if (length(modulegene) == length(modulegene2)) break;
          moduleinf <- module;
          moduleinf$moduleNum[moduleinf$new.ID %in% modulegene] <- "a";
          ModNet <- try(getmoduleHub(data, moduleinf, "a",
                                     coln = "new.ID",
                                     adjustp = FALSE),silent = TRUE);
          if (class(ModNet) == "try-error")
            stop ("data or module have some problem and cannot run getmoduleHub.");
          if (FALSE) {
            node.features<- plot_subgraph(module = ModNet$degreeStat$gene,
                                          hub = ModNet$hub, PFN = ModNet$PMFG,
                                          node.default.color = "grey",
                                          gene.set = IDsets$gene.set,
                                          color.code = IDsets$color.code, show.legend = TRUE,
                                          label.hubs.only = TRUE,
                                          hubLabel.col = "black", hubLabel.sizeProp = 1,
                                          show.topn.hubs = 10,
                                          node.sizeProp = 13, label.sizeProp = 13,
                                          label.scaleFactor = 10,
                                          layout = "kamada.kawai")$node.features;
            removedHub <- as.character(node.features$id[node.features$node.shape == "hub" &
                                                          node.features$node.stat == "NA"]);
          }
          #cat(paste(" - # of genes:", length(modulegene), "\n"))
          #cat(paste(" - # of hubs:", length(ModNet$hub), "\n"))
          removedHub <- ModNet$hub[!ModNet$hub %in% markedID];
          if (plot) {
            pic <- plot_subgraph(module = ModNet$degreeStat$gene,
                                 hub = ModNet$hub, PFN = ModNet$PMFG,
                                 node.default.color = node.default.color,
                                 gene.set = IDsets$gene.set,
                                 color.code = IDsets$color.code, show.legend = TRUE,
                                 label.hubs.only = label.hubs.only,
                                 hubLabel.col = hubLabel.col,
                                 hubLabel.sizeProp = 1,
                                 show.topn.hubs = 10,
                                 node.sizeProp = 13, label.sizeProp = 13,
                                 label.scaleFactor = 10,
                                 layout = "kamada.kawai");
            if (is.character(filename) & length(filename) == 1) {
              if (!dir.exists("plot"))
                dir.create("plot")
              filename2 = paste0("plot/DEP_Mod_net", filename, "_",iter,
                                 ".", filetype)
              ggsave(pic$pnet, filename = filename2, ...)
              unlink(pic$pnet)
            }
            else print(pic$pnet)
          }
          iter <- iter + 1;
          if ( iter == 1000) warning("It's already run 1000 cycle, but no final result. Please email the author.")
        }
      }
      if (OnlyPlotLast) {
        pic <- plot_subgraph(module = ModNet$degreeStat$gene,
                             hub = ModNet$hub, PFN = ModNet$PMFG,
                             node.default.color = node.default.color,
                             gene.set = IDsets$gene.set,
                             color.code = IDsets$color.code, show.legend = TRUE,
                             label.hubs.only = label.hubs.only,
                             hubLabel.col = hubLabel.col,
                             hubLabel.sizeProp = 1,
                             show.topn.hubs = 10,
                             node.sizeProp = 13, label.sizeProp = 13,
                             label.scaleFactor = 10,
                             layout = "kamada.kawai")
        if (is.character(filename) & length(filename) == 1) {
          if (!dir.exists("plot"))
            dir.create("plot")
          filename2 = paste0("plot/DEP_Mod_net", filename, "_",iter,
                             ".", filetype)
          ggsave(pic$pnet, filename = filename2, ...)
          unlink(pic$pnet)
        }
        else print(pic$pnet)
      }
      modulegene = ModNet$degreeStat$gene;
      hubs = ModNet$hub;
      PMFG = ModNet$PMFG;
      if (iteration == 0)
        warning("iteration is 0, it's not remove any Hubs.")
      if (iter == 1 & iteration != 0 )
        warning("all Hubs ID belong IDsets.")
    }
  }

  output <- list(netgene = modulegene, hub = hubs, PMFG = PMFG)
  return(output)
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


