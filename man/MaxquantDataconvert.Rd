\name{MaxQdataconvert}
\alias{MaxQdataconvert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
one-step to extract Maxquant quantification data and convert
}
\description{
Maxquant quantification data extract and homolog protein Uniprot ID match.
}
\usage{
MaxQdataconvert(pgfilename, IDname = "Majority.protein.IDs",
                IDtype = c("MaxQ","none"), CONremove = TRUE,
                justID = TRUE, status1 = TRUE, ENTRY1 = TRUE,
                db1.path = NULL, db2.path = NULL,
                out.folder = NULL, blast.path = NULL,
                savecsvpath = NULL, csvfilename = NULL,
                verbose = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pgfilename}{
  Maxquant quantification file "protein groups.txt"
  }
  \item{IDname}{
  The column name of uniprot ID. The default value is "\code{Majority.protein.IDs}"" which is the column name in MaxQuant data.
  }
  \item{IDtype}{
  "\code{MaxQ}" means \code{proteinGroups} is Maxquant quantification data, "\code{none}" means other type data. This must be (an abbreviation of) one of the strings: "\code{MaxQ}","\code{none}".
  }
  \item{CONremove}{
  a logical value indicated whether remove contaminant IDs. When IDtype is "none", it will remove unmatch ID compared with database2.
  }
  \item{justID}{
  a logical value indicated whether only extract ID when IDtype is "MaxQ".
  }
  \item{status1}{
  a logical value indicated whether extract the first ID status when IDtype is "MaxQ".
  }
  \item{ENTRY1}{
  a logical value indicated whether extract the first ID ENTRY NAME when IDtype is "MaxQ".
  }
  \item{db1.path}{
  fasta file, database of transfered species
  }
  \item{db2.path}{
  fasta file, database of original species
  }
  \item{out.folder}{
  blast result output folder, the folder path should be the same with db1.path
  }
  \item{blast.path}{
  blast+ software install path
  }
  \item{savecsvpath}{
  the information of csv file name output path. The default value means don't save csv file.
  }
  \item{csvfilename}{
  the name of csv file which the data are to be output. The default value means don't save csv file.
  }
  \item{verbose}{
integer level of verbosity. Zero means silent, higher values make the output progressively more and more verbose.
}
  \item{\dots}{
Other arguments.
}
  }
\details{
one-step to extract MaxQuant or other quantification data and convert.
The function contain \link{ID_match} function.
}
\value{
a list of proteomic information.
  \item{protein_IDs}{Portein IDs which is \code{IDname} column information.}
  \item{intensity}{Quantification intensity informaton. When  \code{IDtype} is "\code{none}", it is the  \code{QuanCol} columns information.}
  \item{iBAQ}{Quantification iBAQ intensity informaton.(only for \code{IDtype} is "\code{MaxQ}")}
  \item{LFQ}{Quantification LFQ intensity informaton.(only for \code{IDtype} is "\code{MaxQ}")}
%% ...
}

\author{
Kefu Liu
}
\note{
The function should install blast+ software, Version 2.7.1.
blast+ download website:ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
db1.path, db2.path, out.folder are both need the complete path. Out.folder and db1.path should be in the same folder. Path should have no special character.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{ID_match}
}

\examples{
rm(list = ls())
library(DDPNA)
#database file should be downloaded and extracted the rar file in C Disk.
#blast+ software should be downloaded and installed in C Disk.
#the example file should be downloaded and decompress in C Disk.
out.folder = "C:/DDPNA_test/";
blast.path = "C:/blast/ncbi-blast-2.7.1+/bin/";
pgfilename = "C:/DDPNA_test/Sample_protein_groups.txt"
if (file.exists(pgfilename) & file.exists(out.folder) & dir.exists(blast.path))
  Maxdata <- MaxQdataconvert(pgfilename, IDtype = "MaxQ",
                             db1.path = "C:/DDPNA_test/new.fasta",
                             db2.path = "C:/DDPNA_test/ori.fasta",
                             out.folder = out.folder,
                             blast.path = blast.path)
}
