\name{ID_match}
\alias{ID_match}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
homolog protein Uniprot ID transformation
}
\description{
homolog protein Uniprot ID match
}
\usage{
ID_match(data, db1.path = NULL, db2.path = NULL,out.folder = NULL,
         blast.path = NULL,evalue = 0.1, verbose = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  dataset of protein information.Column Names should contain "ori.ID" and "ENTRY.NAME". "ori.ID" is Uniprot ID
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
  \item{evalue}{
  blast threshold, the lower means more rigorous
  }
  \item{verbose}{
  integer level of verbosity. Zero means silent, 1 means have Diagnostic Messages.
  }
  }
\details{
homolog protein Uniprot ID match is based on the ENTRY.NAME, gene name and sequence homophyly in two different species or different version of database.
}
\value{
a data.frame included 4 columns: ori.ID, ENTRY.NAME, new.ID, match.type.
}

\author{
Kefu Liu
}
\note{
This function should install blast+ software, Version 2.7.1.
blast+ download website:ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
db1.path, db2.path, out.folder are both need the complete path. Out.folder and db1.path should be in the same folder. Path should have no special character.
data should have colname: ori.ID, ENTRY.NAME.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
rm(list = ls())

library(DDPNA)
#database file should be downloaded and extracted the rar file in C Disk.
#blast+ software should be downloaded and installed in C Disk.
#the example file should be downloaded and decompress in C Disk.
data(Sample_ID_data)
# need to install blast+ software
out.folder = "C:/DDPNA_test/";
blast.path = "C:/blast/ncbi-blast-2.7.1+/bin/";
if(file.exists(out.folder) & dir.exists(blast.path))
  data <- ID_match(ID_data,
                   db1.path = "C:/DDPNA_test/new.fasta",
                   db2.path = "C:/DDPNA_test/ori.fasta",
                   out.folder = out.folder,
                   blast.path = blast.path,
                   evalue = 0.1, verbose = 1)
}
