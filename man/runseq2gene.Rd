\name{runseq2gene}
\alias{runseq2gene}

\title{
R wrapped python function to map genomic regions on the sequence-level to genes.
}
\description{
Annotate genome regions of interest to either the nearest TSS or a broader range of neighboring genes.
}
\usage{
runseq2gene(inputfile, 
            search_radius=150000, promoter_radius=200, promoter_radius2=100, 
            genome=c("hg38","hg19","mm10","mm9"), adjacent=FALSE, SNP=FALSE, 
            PromoterStop=FALSE,NearestTwoDirection=TRUE)
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
}


\value{
A matrix with multiple columns.
  \item{Columns 1 to 4}{The same as the first four columns in the input file. }
  \item{PeakLength}{An integer gives the length of the input genomic region. It is the number of base pairs between the start and end of the region. }
  \item{PeakMtoStart_Overlap}{
  An integer gives the distance from the TSS of mapped gene to the middle of genomic region. A negative value indicates that TSS of the mapped gene is at the right of the peak. Otherwise, PeakMtoStart_Overlap reports a numeric range showing the location of overlapped coordinates (exon, intron, CDS, or UTR). }
  \item{type}{
  A character specifies the relationship between the genomic region and the mapped gene.
  \enumerate{
    \item{"Exon"}{
    any part of a genomic region overlaps the exon region of the mapped gene
    }
    \item{"Intron"}{
    any part of a genomic region overlaps an intron region of the mapped gene
    }  
    \item{"cds"}{
    any part of a genomic region overlaps the CDS region
    }
    \item{"utr"}{
    any part of a genomic region overlaps a UTR region
    }
    \item{"promoter"}{
    any part of a genomic region overlaps the promoter region of the mapped gene when an intergenic region of mapped gene covers the input genomic region
    }
    \item{"promoter_internal"}{
    any part of a genomic region overlaps the promoter region of the mapped gene when an adjacent TTS region of mapped gene covers the input genomic region
    }
    \item{"Nearest"}{
    the mapped gene is the nearest gene if the genomic region is located in an intergenic region
    }
    \item{"L" and "R"}{
    show the relative location of mapped genes when the input genomic region resides within a bidirectional region
    }
    \item{"Neighbor"}{
    any mapped gene within the search radius but belongs to none of the prior types
    }
  }
  }
  \item{BidirectionalRegion}{
  A Boolean indicates whether or not the input genomic region is in bidirectional region. \cite{"A 'bidirectional gene pair' refers to two adjacent genes coded on opposite strands, with their 5' UTRs oriented toward one another." (from wiki http://en.wikipedia.org/wiki/Promoter_(genetics) )}. NA means the genomic region is at exon or intron region. }
  \item{Chr}{An integer gives chromosome number of mapped gene. }
  \item{TSS}{An integer indicates transcription start site of mapped gene regardless of strand.  }
  \item{TTS}{An integer indicates transcription termination site of mapped gene regardless of strand. }
  \item{strand}{A character indicates whether mapped gene is in forward (+) or reverse (-) direction on chromosome.}
  \item{gene_name}{A character gives official gene symbol of mapped genes.   }
  \item{source}{A character gives gene source (Ensembl classification) of mapped genes.  }
  \item{transID}{A character gives Ensemble transcript ID of mapped genes. }
}

\references{
  Lawrence M, Huber W, Pages H, Aboyoun P, Carlson M, Gentleman R, Morgan M and Carey V (2013) 
  "Software for Computing and Annotating Genomic Ranges.".
  \emph{PLoS Computational Biology}, \bold{9}.
}

\author{
Bin Wang
}



\examples{
  \dontrun{res=runseq2gene(inputfile=Chipseq_Peak_demo)    }
}

\keyword{methods}

