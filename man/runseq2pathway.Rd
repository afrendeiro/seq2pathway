\name{runseq2pathway}
\alias{runseq2pathway}

\title{
An function to perform the runseq2pathway algorithm(s).
}
\description{
A wrapper function to perform seq2gene and gene2pathway in series.
}
\usage{
runseq2pathway(inputfile,
               search_radius=150000, promoter_radius=200, promoter_radius2=100, 
               genome=c("hg38","hg19","mm10","mm9"), adjacent=FALSE, SNP= FALSE, 
               PromoterStop=FALSE, NearestTwoDirection=TRUE,
               DataBase=c("GOterm"), FAIMETest=FALSE, FisherTest=TRUE,
               collapsemethod=c("MaxMean","function","ME",
               "maxRowVariance","MinMean","absMinMean","absMaxMean","Average"), 
               alpha=5, logCheck=FALSE, B=100, na.rm=FALSE, min_Intersect_Count=5)
}

\arguments{
  \item{inputfile}{
  An R object input file that records genomic region information (coordinates). The file format could be data frame defined as:  
  \enumerate{
    \item{column 1}{
    the unique IDs of genomic regions of interest (peaks, mutations, or SNPs)
    }
    \item{column 2}{
    the chromosome IDs (eg. chr5 or 5)
    }
    \item{column 3}{
    the start of genomic regions
    }
    \item{column 4}{
    the end of genomic regions (for SNP and point mutations, the difference of start and end is 1bp)
    }
    \item{column 5...}{
    Other custom defined information (option)
    }
  }
  Or, the input format should be RangedData object(from R package IRanges) with value column.
    \enumerate{
    \item{column 1: space}{
    the chromosome IDs (eg. chr5 or 5)
    }
    \item{column 2: ranges}{
    the ranges of genomic regions
    }
    \item{column 3: name}{
    the unique IDs of genomic regions of interest (peaks, mutations, or SNPs)
    }
    \item{more columns:}{
    Other custom defined information (optional)
    }
  }
  }

  \item{search_radius}{
  A non-negative integer, with which the input genomic regions can be assigned not only to the matched or nearest gene, but also with all genes within a search radius for some genomic region type. This parameter works only when the parameter "SNP" is FALSE. Default is 150000. 
  }
  \item{promoter_radius}{
  A non-negative integer. Default is 200. Promoters are here defined as upstream regions of the transcription start sites (TSS). User can assign the promoter radius, a suggested value is between 200 to 2000. 
  }
  \item{promoter_radius2}{
  A non-negative integer. Default is 100. Promoters are here defined as downstream regions after the transcription start sites (TSS).
  }
  \item{genome}{
  A character specifies the genome type. Currently, choice of "hg38", "hg19", "mm10", and "mm9" is supported.
  }
  \item{adjacent}{
  A Boolean. Default is FALSE to search all genes within the search_radius. Using "TRUE" to find the adjacent genes only and ignore the parameters "SNP" and "search_radius".
  }
  \item{SNP}{
  A Boolean specifies the input object type. FALSE by default to keep on searching for intron and neighboring genes. Otherwise, runseq2gene stops searching when the input genomic region is residing on exon of a coding gene.
  }
  \item{PromoterStop}{
  A Boolean, "FALSE" by default to keep on searching neighboring genes using the parameter "search_radius". Otherwise, runseq2gene stops searching neighboring genes. This parameter has function only if an input genomic region maps to promoter of coding gene(s).
  }
  \item{NearestTwoDirection}{
  A boolean, "TRUE" by default to output the closest left and closest right coding genes with directions. Otherwise, output only the nearest coding gene regardless of direction.
  }
  \item{DataBase}{
  A character string assigns an R GSA.genesets object to define gene-set. User can call GSA.read.gmt to load customized gene-sets with a .gmt format. If not specified, GO defined gene sets (BP,MF,CC) will be used.
  }
  \item{FAIMETest}{
  A boolean values. By default is FALSE. When true, executes function of gene2pathway test using the FAIME method, which only functions when the fifth column of input file exsists and is a vector of scores or values.
  }
  \item{FisherTest}{
  A Boolean value. By default is TRUE to excute the function of the Fisher's exact test. Otherwise, only excutes the function of gene2pathway test.
  }
  \item{collapsemethod}{
  A character for determining which method to use when call the function collapseRows in package WGCNA. The function "collapsemethod" uses this paramter to call the collapseRows() function in package "WGCNA".
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
An R list of several data frames. The results of function seq2gene, Fisher's exact test and gene2pathway test results are included.
}

\references{
  Langfelder P, Horvath S (2008)
  WGCNA: an R package for weighted correlation network analysis.
  \emph{BMC Bioinformatics},  9:559.

  Miller JA, Cai C, Langfelder P, Geschwind DH, Kurian SM, Salomon DR, Horvath S (2011)
  Strategies for aggregating gene expression data: The collapseRows R function.
  \emph{BMC Bioinformatics}, 12:322.
  
  Lawrence M, Huber W, Pages H, Aboyoun P, Carlson M, Gentleman R, Morgan M and Carey V (2013) 
  "Software for Computing and Annotating Genomic Ranges.".
  \emph{PLoS Computational Biology}, \bold{9}.
}


\author{
Bin Wang
}

\keyword{methods}

