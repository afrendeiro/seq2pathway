\name{addDescription}
\alias{addDescription}

\title{
Retrieve "gene description" attributes for gene symbol.
}
\description{
A function wrappered from Rpackage "biomaRt". Get gene description information from gene symbol information.
}
\usage{
addDescription(genome=c("hg19","mm10","mm9"), genevector)
}

\arguments{
  \item{genome}{
A character specifies the genome type. Currently, choice of "hg19", "mm10", and "mm9" is supported.
}
  \item{genevector}{
A characteristic vector of gene symbols.
}
}

\value{
A characteristic matrix of gene symbols and descriptions.
}

\references{
  Durinck S, Spellman P, Birney E and Huber W (2009)
  Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt.
  \emph{Nature Protocols}, \bold{4}, 1184--1191.

  Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A and Huber W (2005)
  BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis.
  \emph{Bioinformatics}, \bold{21}, 3439--3440.
}

\author{
Bin Wang
}


\examples{
  data(dat_chip)
  gene_description<-addDescription(genome="hg19",genevector=rownames(dat_chip))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

