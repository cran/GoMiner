% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GoMiner.R
\name{GOheatmap}
\alias{GOheatmap}
\title{GOheatmap}
\usage{
GOheatmap(sampleList, x, thresh, fdrThresh = 0.105, verbose)
}
\arguments{
\item{sampleList}{character list of gene names}

\item{x}{DB component of return value of GOtable3()}

\item{thresh}{output of GOthresh()}

\item{fdrThresh}{numeric value of FDR acceptance threshold}

\item{verbose}{integer vector representing classes}
}
\value{
returns a matrix to be used as input to a heat map
}
\description{
generate a matrix to be used as input to a heat map
}
\examples{
\dontrun{
# GOGOA3.RData is too large to include in the R package
# so I need to load it from a file that is not in the package.
# Since this is in a file in my own file system, I could not
# include this as a regular example in the package.
# you can generate it using the package 'minimalistGODB'
# or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
heatmap<-GOheatmap(cluster52,GOGOA3$ontologies[["biological_process"]],x_thresh,verbose=1)
}

}
