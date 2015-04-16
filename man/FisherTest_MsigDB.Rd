\name{FisherTest_MsigDB}
\alias{FisherTest_MsigDB}

\title{
A wrapper function to perform conditional Fisher's exact test, using custom-defined genesets.
}
\description{
A wrapper function to perform conditional Fisher's exact test, using custom-defined genesets.
}
\usage{
FisherTest_MsigDB(gsmap, gs, genome=c("hg38","hg19","mm10","mm9"), 
                  min_Intersect_Count=5)
}
\arguments{
  \item{gsmap}{
  An R GSA.genesets object defined by the package "GSA" for functional gene-set (or termed as pathway for simplification). User can call the GSA.read.gmt function in R GSA package to load customized gene-sets with a .gmt format. 
}
  \item{gs}{A characteristic vector of gene symbols, the input genelist.}
  \item{genome}{
  A character specifies the genome type. Currently, choice of "hg38", "hg19", "mm10", and "mm9" is supported.
  }
  \item{min_Intersect_Count}{
  A number decides the cutoff of the minimum number of intersected genes when reporting Fisher's exact tested results.
  }


}

\value{
A data frame of Fisher's exact tested result with the following columns:
  \item{GeneSet}{MSigDB gene-set names (IDs)}
  \item{Description}{MSigDB definition and description for gene-sets}
  \item{Fisher_Pvalue}{the raw Pvalues} 
  \item{Fisher_odds}{estimate of the odds ratios} 
  \item{FDR}{the multi-test adjusted Pvalues using the Benjamini and Hochberg method} 
  \item{Intersect_Count}{the sizes of the overlap between the MSigDB gene-set genes and the input genelist} 
  \item{MsigDB_gene_inBackground}{the counts of genes among each MSigDB gene-set that are also within genome background}
  \item{MsigDB_gene_raw_Count}{the original counts of genes in each MSigDB gene-set}  
  \item{Intersect_gene}{the intersecting genes' symbols}
}

\author{
Bin Wang
}





\examples{
  data(dat_chip)
  data(MsigDB_C5,package="seq2pathway.data")
  FS_test<-FisherTest_MsigDB(gsmap=MsigDB_C5, gs=rownames(dat_chip)[1:5], genome="hg19", min_Intersect_Count=1)

}

\keyword{methods}

