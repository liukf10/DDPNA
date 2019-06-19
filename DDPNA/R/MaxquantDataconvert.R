#one-step to extract Maxquant quantification information and convert
#20190418version
#20190517
MaxQdataconvert <- function(pgfilename, IDname = "Majority.protein.IDs",
                            IDtype = c("MaxQ","none"), CONremove = TRUE,
                            justID = TRUE, status1 = TRUE, ENTRY1 = TRUE,
                            db1.path = NULL, db2.path = NULL,
                            out.folder = NULL, blast.path = NULL,
                            savecsvpath = NULL, csvfilename = NULL,
                            verbose = 1, ...)
{
  if (!file.exists(pgfilename)) stop("pgfilename is not exist.")
  my_fasta <- read.delim(pgfilename, stringsAsFactors = FALSE)
  IDtype <- match.arg(IDtype, c("MaxQ","none"))
  data <- MaxQprotein(my_fasta, IDname = IDname,
                      IDtype = IDtype, remove = CONremove,
                      verbose = verbose, ...);
  NEXTtoIDmatch = FALSE;
  if (IDtype == "MaxQ"){
    inf <- try(P.G.extract(data$protein_IDs, justID = justID,
                           status1 = status1, ENTRY1 = ENTRY1,
                           verbose = verbose-1));
    if (class(inf) != "try-error") {
      inf <- data.frame(inf[ ,1:3], stringsAsFactors = FALSE);
      colnames(inf) <- c("db.type", "ori.ID", "ENTRY.NAME");
      NEXTtoIDmatch = TRUE;
      }
  }


  if (IDtype == "none") {
    #inf <- my_fasta[ ,IDname];
    inf <- data$protein_IDs;
    inf <- gsub("(.*?);.*","\\1",inf);
    inf <- gsub(".*(tr|sp)\\|(.*)\\|.*", "\\2", inf);
    if (!is.null(db2.path)) {
      if(file.exists(db2.path)){
        my_sequences <- try(.uniprot_database(input_file = db2.path))
        if (class(data_match) != "try-error") {
          pos <- match(inf, my_sequences$pro_info$Uniprot.ID);
          inf <- my_sequences$pro_info[pos, 1:3];
          con <- which(is.na(inf[ ,2]));
          #if (length(con) != 0) inf[con, 2] <- my_fasta[con, IDname];
          if (length(con) != 0) inf[con, 2] <- data$protein_IDs[con];
          if (CONremove & length(con)!= 0) {
            inf <- inf[-con, ];
            data$protein_IDs <- data$protein_IDs[-con]
            data$intensity <- data$intensity[-con, ];
          }
          rownames(inf) <- 1:nrow(inf);
          colnames(inf) <- c("db.type", "ori.ID", "ENTRY.NAME");
          NEXTtoIDmatch = TRUE;
          if (length(con) == length(inf)) {
            inf <- data$protein_IDs;
            if(verbose > 0) warning("protein ID is no match for db2.");
            NEXTtoIDmatch = FALSE;
          }
        }
      }
    }
  }

  if(1!=1)
  {if (is.null(db1.path) | is.null(db2.path) | is.null(out.folder) | is.null(blast.path)) {
    warning("Break short. ID match is not execute.");
    NEXTtoIDmatch = FALSE;
  } else if (file.exists(db1.path) & file.exists(db2.path) &  dir.exists(blast.path)) {
    warning("Break short. ID match is not execute.");
    NEXTtoIDmatch = FALSE;
  }}

  if (NEXTtoIDmatch) {
    data_match <- try(ID_match(data = inf,
                               db1.path = db1.path,
                               db2.path = db2.path,
                               out.folder=out.folder,
                               blast.path=blast.path,
                               verbose = verbose));
    if (class(data_match) != "try-error") {
      inf <- data.frame(data_match, db.type = inf$db.type,
                        protein_IDs = data$protein_IDs,
                        stringsAsFactors = FALSE);
    } else warning("ID match is not execute.");
  } else {
    data_match <- try(ID_match(data = inf), silent = TRUE);
    warning("ID match is not execute.")
  }

  if(IDtype == "none")
    data <- list(inf = inf, intensity = data$intensity) else
      data <- list(inf = inf, intensity = data$intensity,
                   iBAQ = data$iBAQ, LFQ = data$LFQ);
  if(NEXTtoIDmatch & class(data_match) != "try-error") {
    if (is.character(savecsvpath) & is.character(csvfilename)){
      my_fasta_match <- data.frame(inf, data$LFQ);
      if (!dir.exists(savecsvpath)) dir.create(savecsvpath);
      csvfilename <- paste0(savecsvpath,"/",csvfilename);
      write.csv(my_fasta_match, file = csvfilename);
    } else if (is.character(savecsvpath)) {
      if (verbose > 0) warning("wrong csvfilename");
    } else if (verbose > 0) warning("Don't save csv file.")
  }
  data
}


#read Maxquant proteingroups data, remove contaminant and reverse protein
#20190517 add verbose
MaxQprotein <- function(proteinGroups, IDname = "Majority.protein.IDs",
                        IDtype = "MaxQ", remove = TRUE,
                        QuanCol = NULL, verbose = 1){
  removeConRev <- function(infile){
    ##  remove contaminant, reverse proteinID
    if (is.element("Potential.contaminant", colnames(infile)) &
        is.element("+", unique(infile$Potential.contaminant))) {
      infile <- infile[-which(infile$Potential.contaminant %in% "+"), ]
    }
    if (is.element("Reverse", colnames(infile)) &
        is.element("+", unique(infile$Reverse))) {
      infile <- infile[-which(infile$Reverse %in% "+"), ]
    }
    infile
  }
  if (length(IDname) != 1 & is.vector(IDname) )
    stop ("IDname only allow one column name.")
  if (is.na(IDname)) stop ("IDname is NA.")
  IDtype <- match.arg(IDtype, c("MaxQ","none"));
  if (remove & IDtype == "MaxQ") {
    proteinGroups <- removeConRev(proteinGroups);
    if(verbose > 0) message('* + Contaminant, + Reverse, proteins are removed.')
    }
  colname <- colnames(proteinGroups);
  if (is.element(IDname, colname))
    protein_IDs <- proteinGroups[, IDname] else
      stop(paste("no", IDname));
  if (IDtype == "MaxQ"){
    intensity <- grep("Intensity..*", colname);
    if (length(intensity) == 0)
      intensity <- NULL else
        intensity <- proteinGroups[, intensity];
    iBAQ <- grep("iBAQ..*", colname);
    if (length(iBAQ) == 0)
      iBAQ <- NULL else
        iBAQ <- proteinGroups[, iBAQ];
    LFQ <- grep("LFQ.intensity..*", colname);
    if (length(LFQ) == 0)
      LFQ <- NULL else
        LFQ <- proteinGroups[, LFQ];
    data <- list(protein_IDs = protein_IDs, intensity = intensity, iBAQ = iBAQ, LFQ = LFQ)
  }
  if (IDtype == "none"){
    if (is.null(QuanCol)) {
      IDpos <- which(colname == IDname);
      intensity <- proteinGroups[, -IDpos];
    }
    if (is.numeric(QuanCol)) {
      if (max(QuanCol) <= ncol(proteinGroups) & min(QuanCol) >= 1)
        intensity <- proteinGroups[, QuanCol] else
          stop ("QuanCol is wrong.
Some of the numbers are larger than proteinGroups column number or less than 1.")
    }
    if (is.character(QuanCol)) {
      if (all(is.element(QuanCol, colname)))
        intensity <- proteinGroups[, QuanCol] else
          stop ("QuanCol is wrong.")
    }
    data <- list(protein_IDs = protein_IDs, intensity = intensity)
  }
  data
}


#extract protein groups information
#20190105version
#20190517 add verbose
P.G.extract <- function(inf, ncol = 4,
                        justID = FALSE, status1 = FALSE,
                        ENTRY1 = FALSE, verbose = 0) {
  n.gene <- length(inf)
  if (n.gene == 0) stop("wrong inf.")
  inf_matrix <- matrix(rep("", ncol*n.gene), nrow = n.gene, ncol = ncol)
  for (i in 1:n.gene){
    gene <- unlist(strsplit(inf[i], split = ";"))
    status <- NULL; ENTRY <- NULL;
    seprow1 <- unlist(strsplit(gene[1], split = "\\|"));
    if (status1) status <- seprow1[1];
    #status<-gsub("(tr|sp)\\|.*","\\1",gene[1])
    if (ENTRY1) ENTRY <- seprow1[3];
    #ENTRY<-gsub("tr|sp\\|.*\\|(.*)","\\1",gene[1])
    if (justID) gene <- gsub(".*(tr|sp)\\|(.*)\\|.*", "\\2", gene);
    gene <- c(status, gene[1], ENTRY, gene[-1]);
    if (ncol >= length(gene))
      inf_matrix[i, 1:length(gene)] <- gene else {
        inf_matrix[i,] <- gene[1:ncol];
        if(verbose > 0) warning (paste("The", i, "proteins have more than",
                                       ncol, "protein ID"));
        }
  }
  inf_matrix
}

################################
#need install blast+   blast.path="e:/blast/ncbi-blast-2.7.1+/bin/"
#out.folder and db1.path should be in the same folder path
#path should have no special character
#db1.path,db2.path,out.folder are both need the complete path
#input should have colname: ori.ID,ENTRY.NAME  output have 4 cols:ori.ID,ENTRY.NAME,new.ID,match.type
#db1.path is transfered species fasta file, db2.path is original species fasta file
#blast+ download website:ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#blast+ version need 2.7.1
#evalue lower means more rigorous;verbose:0 no verbose:1 have
#190601 reset the default working directory
#190605 add alternative solution without use blast+

#db1.path="c:/workingdirectory/human-uniprot-20180108.fasta"
#db2.path="c:/workingdirectory/mouse-uniprot-20180108.fasta"
#out.folder="c:/workingdirectory"
#blast.path="e:/blast/ncbi-blast-2.7.1+/bin/"
ID_match <- function(data, db1.path = NULL, db2.path = NULL,
                     out.folder = NULL, blast.path = NULL,
                     evalue = 0.1, verbose = 1)
{

  question <- function (...)#190605
  {
    yes <- c("Yes", "Definitely", "Sure",
               "I agree", "Absolutely")
    no <- c("No", "Nope", "Not Sure")
    cat(paste0(..., collapse = ""))
    qs <- c(sample(yes, 1), sample(no, 1))
    menu(qs) == 1
  }
  #extract local database
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings in Bioconductor needed for this function to work. Please install it.",
         call. = FALSE)
  }
  old.wd <- getwd(); on.exit(setwd(old.wd)); #190601
  if (is.null(db1.path)) stop ("db1.path is null.");
  if (is.null(db2.path)) stop ("db2.path is null.");
  if (is.null(out.folder)) stop ("out.folder path is null.");
  if (!dir.exists(out.folder)) dir.create(out.folder);
  if (!file.exists(db2.path)) stop ("db2.path is wrong.");
  if (!file.exists(db1.path)) stop ("db1.path is wrong.");
  db1.folder <- gsub("(.*/+).*","\\1", db1.path);
  db1.folder <- c(db1.folder, gsub("(.*)/","\\1", db1.folder));
  #db1.folder <- gsub(paste0("(",out.folder,").*"), "\\1", db1.path);
  if (all(out.folder != db1.folder)) stop ("db1.path is not in the out.folder");
  Alignment = FALSE; #190605
  #190605
  if (is.null(blast.path)) {
    message("blast.path is null.\n");
    #cat("blast.path is null.\n");
    if (question("Do you want to do ID_match without using blast+ software?"))
      {if(question("It will take large amouts of time.  Are you sure?"))
        Alignment <- TRUE else stop("Please give the blast.path and re-run the function.");}
    else stop("Please give the blast.path and re-run the function.");
  } else if (!dir.exists(blast.path)) {
    message("blast.path is not exists.\n");
    #cat("blast.path is not exists.\n");
    if (question("Do you want to do ID_match without using blast+ software?"))
      {if(question("It will take large amouts of time. Are you sure?"))
        Alignment <- TRUE else stop("Please give the correct blast.path and re-run the function.");
    } else stop("Please give the correct blast.path and re-run the function.");
  }#190605
  my_sequences <- .uniprot_database(input_file = db1.path);
  if(class(my_sequences) == "try-error") #190605
    stop("db1.path file is not the correct format.") else{ #190605
      database1 <- my_sequences$pro_info;
      database1_seq <- my_sequences$pro_seq;
      if (verbose > 0) message("db1 extraction is completed")
    }
  my_sequences <- .uniprot_database(input_file = db2.path);
  if(class(my_sequences) == "try-error") #190605
    stop("db1.path file is not the correct format.") else{ #190605
      database2 <- my_sequences$pro_info;
      database2_seq <- my_sequences$pro_seq;
      if (verbose > 0) message("db2 extraction is completed")
    }
  #190605 add error information
  if (sum(colnames(data) %in% c("ori.ID")) * sum(colnames(data) %in% c("ENTRY.NAME")) != 1 )
    stop("data should have column: ori.ID and ENTRY.NAME.")


  #1.ENTRY.NAME match between data and database1, no matched IDs saved in datano1,match.type is 1
  #delete species information in ENTRY.NAME.
  protein_name_despecies <- function(data){
    pro_name <- gsub("(.*)_.*", "\\1", as.character(data$ENTRY.NAME), perl = TRUE);
    data$ENTRY.NAME <- pro_name;
    data;
  }
  #data database1 despecies
  data_despecies <- protein_name_despecies(data);
  database1_despecies <- protein_name_despecies(database1);
  #order is ENTRY.NAME,ori.ID,new.ID
  data_IDmatch <- merge(data_despecies[ , c('ori.ID', 'ENTRY.NAME')],
                        database1_despecies[ , c('Uniprot.ID', 'ENTRY.NAME')],
                      by='ENTRY.NAME',all.x=TRUE,sort = FALSE);
  names(data_IDmatch)[names(data_IDmatch) == 'Uniprot.ID'] <- 'new.ID';
  #match.type is 1
  match.type <- "1";
  data_IDmatch <- data.frame(data_IDmatch, match.type, stringsAsFactors = FALSE);
  #data_IDmatch$match.type[is.na(data_IDmatch$new.ID)]<-NA;
  data_IDmatch <- merge(data[ , c('ori.ID', 'ENTRY.NAME')],
                      data_IDmatch[ ,-1],
                      by = 'ori.ID', sort = FALSE);
  #data_IDmatch$new.ID<-as.character(data_IDmatch$new.ID);#factor to character
  datano1 <- data_IDmatch[is.na(data_IDmatch$new.ID), ];
  if (verbose > 0) message("ID match is completed, match.type is 1")
  #2.gn match, no matched IDs in datano2, match.type is 2
  if (nrow(datano1) != 0){
    #add GN information
    data_withgn <- merge(datano1, database2[ ,c('Uniprot.ID', 'GN')],
                       by.x = 'ori.ID',by.y = 'Uniprot.ID',sort = FALSE);
    nogn <- which(data_withgn$GN == "") #181210
    data_withgn <- data_withgn[-nogn, ] #181210
    # extract GN information
    data_GN <- gsub("(.*)","^\\1", as.character(data_withgn$GN), perl = TRUE);
    #limitation maxium matched protein ID number is 14
    GN <- matrix(nrow = length(data_GN), ncol = 14);
    ID <- NULL;
    #GN match: if no match, delete the last character and match again
    #delete at most 4 characters, if match more than 1 ID, do pairwise
    for (j in 1:length(data_GN))
    {
      gn <- NULL;
      ##maxium delete character:4
      max <- nchar(data_GN[j]);
      for (i in 1:5)
      {
        ##match first, if no then delete character
        if (length(gn) == 0){
          data_GN[j] <- substr(data_GN[j], start = 1, stop = max+1-i);
          gn <- grep(data_GN[j], database1$GN,
                     ignore.case = TRUE, perl = TRUE, value = TRUE);
          ##if delete 4 character and no match,then gn<-NA
          if (i == 5 && length(gn) == 0) gn <- NA;
          if (i == 5 && length(gn) > 14) gn <- NA;
        }
        ##if match, how many ID have been matched.
        else{
          if (length(gn) > 14) gn<-NA;
          break;
        }
      }
      GN[j, 1:length(gn)] <- gn;
      id <- as.character(database1$Uniprot.ID[database1$GN %in% gn]);
      ##if ID# more than 1, pairwise
      if (length(id) != 1 && length(id) != 0){
        similar_seq <- as.character(database1_seq[names(database1_seq) %in% id]);
        name <- data_withgn$ori.ID[j];
        ori_seq <- as.character(database2_seq[names(database2_seq) == name]);
        #remove the unregular amino acid
        aa <- unlist(Biostrings::strsplit(unname(ori_seq),"")); #190605
        aa <- aa[aa %in% c("A","C","D","E","F","G","H","I",
                           "K","L","M","N","P","Q","R","S",
                           "T","V","W","Y")]; #190605
        aa <- paste0(aa,collapse = "");   #190605
        ori_seq <- aa;  #190605
        scorem <- 0;
        type2_results <- list();
        for (k in 1:length(similar_seq)) {#190605
          pwalign <- try(Biostrings::pairwiseAlignment(as.character(ori_seq),
                                                       as.character(similar_seq[k]),
                                                       type = "overlap",
                                                       substitutionMatrix = "BLOSUM62",
                                                       gapOpening = 9.5,
                                                       gapExtension = 0.5,
                                                       scoreOnly = TRUE), silent = TRUE);
          if (class(pwalign) == "try-error") {
            aa <- unlist(Biostrings::strsplit(unname( as.character(similar_seq[k]) ),""));
            aa <- aa[aa %in% c("A","C","D","E","F","G","H","I",
                               "K","L","M","N","P","Q","R","S",
                               "T","V","W","Y")];
            aa <- paste0(aa,collapse = "")
            align_seq <- aa;
            pwalign <- try(Biostrings::pairwiseAlignment(as.character(ori_seq),
                                                         align_seq,
                                                         type = "overlap",
                                                         substitutionMatrix = "BLOSUM62",
                                                         gapOpening = 9.5,
                                                         gapExtension = 0.5,
                                                         scoreOnly = TRUE), silent = TRUE);
            if (class(pwalign) == "try-error")
            {pwalign = 0;
            warning("Please email the warning to the author. Thank you!")
            }
          }
          scorem[k] <- pwalign;
        }#190605
        type2_results[[name]] <- data.frame(rep(name, length(similar_seq)),
                                            similar_seq,
                                            scorem);
        ID[j] <- id[which.max(scorem)];
      }
      ##if ID# is 1 then record
      else if(length(id) != 0) ID[j] <- id
      else {id <- NA; ID[j] <- id;}
    }
    #no matched ID in datano2
    new.ID <- ID;
    match.type <- rep('2', length(new.ID));#181210
    data_withgn <- data.frame(data_withgn[,c(-3,-4,-5)],
                              new.ID,
                              match.type,
                              stringsAsFactors = FALSE);#181210
    datano1[datano1$ori.ID %in% data_withgn$ori.ID, ] <- data_withgn;#181210
    #match.type<-rep('2',length(datano1$ori.ID)); #181210
    #datano1<-data.frame(datano1[,c(-3,-4)],new.ID,match.type,stringsAsFactors = FALSE);#181210
    datano1$match.type[is.na(datano1$new.ID)] <- NA;
    datano2 <- datano1[is.na(datano1$new.ID), ];
    if(verbose > 0) message("GN similar match is completed, match.type is 2");
    #
    data_IDmatch[is.na(data_IDmatch$new.ID), ] <- datano1;
    #3.blast match.type is 3
    #ceshi#if(length(datano2$ori.ID)==length(data$ori.ID)){#ceshi#
    if (length(datano2$ori.ID) != 0) {
      ID <- NULL;
      data_ori.ID <- as.character(datano2$ori.ID);
      #if("micropan" %in% rownames(installed.packages()) == FALSE)
      #{install.packages("micropan")}
      #suppressMessages(library("micropan"))
      if (!Alignment) { #190605
        log.fil <- file.path(out.folder, "log.txt");
        db.fil <- file.path(out.folder, "blastDB");
        command <- paste("makeblastdb -logfile", log.fil,
                         "-dbtype prot -out", db.fil, "-in", db1.path);
        setwd(blast.path)
        system(command)
        Continue <- TRUE; #190605
        #190605
        if ( !file.exists(log.fil) ) {
          warning("blast.path have no blast+ software.\n");
          Continue <- FALSE;
          if (question("Do you want to do ID_match without using blast+ software?"))
            if(question("It will take large amouts of time. Are you sure?"))
              Alignment <- TRUE;
        } else if (!file.exists(paste0(db.fil,".psq"))) {
          file.remove( log.fil, paste( log.fil, ".perf", sep=""));
          warning("blast+ software cannot extract db1.bath information.\n")
          #cat("Please email the error to author. Thank you!\n");
          #cat("The author Email is liukefu19@163.com\n")
          stop("Please email the error to author. Thank you!\n The author Email is liukefu19@163.com");
        }
        if (Continue) {
          type3_results <- list();
          for (j in 1:length(data_ori.ID)) {
            name <- data_ori.ID[j];
            ori_seq <- as.character(database2_seq[names(database2_seq) == name]);
            write(ori_seq, file = file.path(out.folder, "input.fasta"));
            input <- paste( "-query ", file.path(out.folder, "input.fasta"), sep = "" );
            dbase <- paste( "-db ", db.fil, sep = "" );
            output <- paste( "-out ",
                             file.path( out.folder, "blast_result.txt" ),
                             sep = "" );
            # command <- paste( "blastp -matrix BLOSUM45 -evalue", 0.001, "-num_threads", 1,
            #                  "-outfmt 6", input, dbase, output )
            command <- paste( "blastp -evalue", evalue,
                              "-outfmt 6", input, dbase, output );
            system(command)
            if(length(scan(file.path(out.folder, "blast_result.txt"),
                           what = character(), quiet = TRUE))== 0)
              ID[j] <- NA
            else{
              similar_seq <- read.table(file.path( out.folder, "blast_result.txt" ));
              similar_seq[,1] <- data_ori.ID[j];
              type3_results[[name]] <- similar_seq;
              similar_ID <- as.character(similar_seq[1, 2]);
              similar_ID <- unlist(strsplit(similar_ID, split = "\\|"));
              ID[j] <- similar_ID[2];
            }
          }
          file.remove( paste( db.fil, ".pin", sep="" ) );
          file.remove( paste( db.fil, ".phr", sep="" ) );
          file.remove( paste( db.fil, ".psq", sep="" ) );
          file.remove( log.fil, paste( log.fil, ".perf", sep=""));
          file.remove(file.path( out.folder, "input.fasta"),
                      file.path( out.folder, "blast_result.txt"));
        }
      }#190605
      #190605
      if (Alignment) {
        type3_results <- list();
        for (j in 1:length(data_ori.ID)) {
          name <- data_ori.ID[j];
          ori_seq <- as.character(database2_seq[names(database2_seq) == name]);
          aa <- unlist(Biostrings::strsplit(unname(ori_seq),""));
          aa <- aa[aa %in% c("A","C","D","E","F","G","H","I",
                             "K","L","M","N","P","Q","R","S",
                             "T","V","W","Y")];
          aa <- paste0(aa,collapse = "");
          ori_seq <- aa;
          scorem <- 0;
          for (k in 1:length(database1_seq)) {
            pwalign <- try(Biostrings::pairwiseAlignment(as.character(ori_seq),
                                                         as.character(database1_seq[k]),
                                                         type = "overlap",
                                                         substitutionMatrix = "BLOSUM62",
                                                         gapOpening = 9.5,
                                                         gapExtension = 0.5,
                                                         scoreOnly = TRUE), silent = TRUE);
            if (class(pwalign) == "try-error") {
              aa <- unlist(Biostrings::strsplit(unname(as.character(database1_seq[k])),""));
              aa <- aa[aa %in% c("A","C","D","E","F","G","H","I",
                                 "K","L","M","N","P","Q","R","S",
                                 "T","V","W","Y")];
              aa <- paste0(aa,collapse = "")
              align_seq <- aa;
              pwalign <- try(Biostrings::pairwiseAlignment(as.character(ori_seq),
                                                           align_seq,
                                                           type = "overlap",
                                                           substitutionMatrix = "BLOSUM62",
                                                           gapOpening = 9.5,
                                                           gapExtension = 0.5,
                                                           scoreOnly = TRUE), silent = TRUE);
              if (class(pwalign) == "try-error")
              {pwalign = 0;
              warning("Please email the warning to the author. Thank you!")
              }
            }
            scorem[k] <- pwalign;
          }
          ID[j] <- names(database1_seq)[which.max(scorem)];
        }
      }
      if (Alignment | Continue) { #190605
        new.ID <- ID;
        match.type <- "3";
        datano2 <- data.frame(datano2[,c(-3,-4)],
                              new.ID,
                              match.type,
                              stringsAsFactors = FALSE);
        #
        data_IDmatch[is.na(data_IDmatch$new.ID), ] <- datano2;
        if (verbose > 0) message("seq blast is completed, match.type is 3")
      } else message("match.type 3 is not run.")#190605
    }
  }
  data_IDmatch
}



#Uniprot database fasta file(type2 suited for human, mouse, Monkey )
##version2.0 no factor format
#190605 add error when the file is not correct.
.uniprot_database <- function(input_file, type = 1) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings in Bioconductor needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # read fasta file
  my_fasta <- try(Biostrings::readAAStringSet(input_file), silent = TRUE); #190605
  if(class(my_fasta) != "try-error") { #190605
    protein_inf <- names(my_fasta);
    #seperate information based by "|"
    if (type == 1) {
      #protein_inf<-unlist(strsplit(protein_inf,split="_HUMAN.*|_MOUSE.*|_MACMU.*|_MACFA.*"))
      protein_inf <- unlist( strsplit( protein_inf, split = "_.*"))
      protein_inf <- data.frame(matrix( unlist( strsplit( protein_inf, split="\\|")),
                                        ncol = 3,byrow = TRUE),
                                stringsAsFactors = FALSE)
      names(protein_inf)<- c("status", "Uniprot.ID", "ENTRY-NAME");
      GN <- gsub(".*GN=(.*) PE=.*", "\\1",
                 as.character(names(my_fasta)), perl = TRUE)
      GN <- gsub(".*\\|.*", "", GN, perl = TRUE)
      protein_inf<-data.frame(protein_inf,GN,stringsAsFactors = FALSE);
    }
    if(type == 2){
      protein_inf <- data.frame(matrix(unlist(strsplit(protein_inf, split = "\\|")),
                                       ncol = 3,byrow = T),
                                stringsAsFactors = FALSE);
      names(protein_inf) <- c("status", "Uniprot.ID", "ENTRY-NAME");
      #extract entryname and GN
      entryname <- gsub("(.*_HUMAN).*|(.*_MOUSE).*|(.*_MACMU).*|(.*_MACFA).*",
                        "\\1\\2", as.character(protein_inf[,3]), perl = TRUE);
      GN <- gsub(".*GN=(.*) PE=.*", "\\1",
                 as.character(protein_inf[ ,3]), perl = TRUE);
      protein_inf[ ,3] <- entryname;
      protein_inf <- data.frame(protein_inf, GN, stringsAsFactors = FALSE);
    }
    my_id_sequence <- Biostrings::readAAStringSet(input_file);
    names(my_id_sequence) <- as.character(protein_inf$Uniprot.ID);
    list(pro_info = protein_inf, pro_seq = my_id_sequence);
  } else my_fasta; #190605
}
