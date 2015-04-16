\name{Chipseq_Peak_demo}
\alias{Chipseq_Peak_demo}
\docType{data}
\title{
chip seq loci data example
}
\description{
chip seq loci data example
}
\usage{data("Chipseq_Peak_demo")}
\format{
  A data frame with 5 observations on the following 5 variables.
  \describe{
    \item{\code{peakID}}{unique chip peak name information}
    \item{\code{chrom}}{chromosome information}
    \item{\code{start}}{loci start}
    \item{\code{end}}{loci end}
    \item{\code{signalvalue}}{a numeric vector}
  }
}
\value{
a data frame of chip sequence peak information
}

\examples{
data(Chipseq_Peak_demo)
head(Chipseq_Peak_demo)
}
\keyword{datasets}
