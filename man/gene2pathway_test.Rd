\name{gene2pathway_test}
\alias{gene2pathway_test}

\title{
A wrapper function to perform gene2pathway test.
}
\description{
The function includes two part, one runs the classical Fisher's exact test, the other runs novel gene2pathway test.
}
\usage{
gene2pathway_test(dat, DataBase="GOterm", FisherTest=TRUE, EmpiricalTest=FALSE, 
                 method=c("FAIME", "KS-rank", "cumulative-rank"), 
                 genome=c("hg38","hg19","mm10","mm9"), alpha=5, 
                 logCheck=FALSE, na.rm=FALSE, B=100, min_Intersect_Count=5)
}

\arguments{
  \item{dat}{
  A data frame of gene expression or a matrix of sequencing-derived gene-level measurements. The rows of dat correspond to genes, and the columns correspond to sample profile (eg. Chip-seq peak scores, somatic mutation p-values, RNS-seq or micro-array gene expression values). Note that the rows must be labeled by official gene symbol. The values contained in dat should be either finite or NA.
}
  \item{DataBase}{
  A character string assigns an R GSA.genesets object to define gene-set. User can call GSA.read.gmt to load customized gene-sets with a .gmt format. If not specified, GO defined gene sets (BP,MF,CC) will be used.
  }
  \item{FisherTest}{
  A Boolean value. By default is TRUE to excute the function of the Fisher's exact test. Otherwise, only excutes the function of gene2pathway test.
  }
  \item{EmpiricalTest}{
  A Boolean value. By default is FALSE for multiple-sample dat. When true, gene2pathway_test calculates empirical p-values for gene-sets.
  }
  \item{method}{
  A character string determines the method to calculate the pathway scores. Currently, "FAIME" (default), "KS-rank", and "cumulative-rank" are supported.
  }
  \item{genome}{
  A character specifies the genome type. Currently, choice of "hg38", "hg19", "mm10", and "mm9" is supported.
  }
  \item{alpha}{
  A positive integer, 5 by default. This is a FAIME-specific parameter. A higher value puts more weights on the most highly-expressed ranks than the lower expressed ranks.
  }
  \item{logCheck}{
  A Boolean value. By default is FALSE. When true, the function takes the log-transformed values of gene if the maximum value of sample profile is larger than 20.
  }
  \item{na.rm}{
  A Boolean value indicates whether to keep missing values or not when method="FAIME". By default is FALSE.
  }
  \item{B}{
  A positive integer assigns the total number of random sampling trials to calculate the empirical pvalues. By default is 100.
  }
  \item{min_Intersect_Count}{
  A number decides the cutoff of the minimum number of intersected genes when reporting Fisher's exact tested results.
  }
}

\value{
A list or data frame. If the parameter "FisherTest" is true, the result is a list including both reports for Fisher's exact test and the gene2pathway test. Otherwise, only reports the gen2pathway tested results.
}

\author{
Xinan Yang
}
\examples{
data(dat_chip)
data(MsigDB_C5,package="seq2pathway.data")
res<-gene2pathway_test(dat=head(dat_chip), DataBase=MsigDB_C5, FisherTest=FALSE, EmpiricalTest=FALSE, 
                 method="FAIME", genome="hg19", alpha=5, 
                 logCheck=FALSE, na.rm=FALSE, B=100, min_Intersect_Count=1)
}
\keyword{methods}