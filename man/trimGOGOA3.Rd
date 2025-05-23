% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GoMiner.R
\name{trimGOGOA3}
\alias{trimGOGOA3}
\title{trimGOGOA3}
\usage{
trimGOGOA3(GOGOA3, mn, mx, verbose)
}
\arguments{
\item{GOGOA3}{return value of subsetGOGOA()}

\item{mn}{integer min size threshold for a category}

\item{mx}{integer max size threshold for a category}

\item{verbose}{integer vector representing classes}
}
\value{
returns trimmed version of GOGOA3
}
\description{
remove categories from GOGOA3 that are too small or too large
}
\details{
If a category is too small, it is unreliable for statistical evaluation
Also, in the extreme case of size = 1, then that category is essentially
equivalent to a gene rather than a category. Same is partially true for size = 2.
If a category is too large, it is too generic to be useful for categorization.
Finally, by trimming the database, analyses will run faster.
}
\examples{
\dontrun{
# GOGOA3.RData is too large to include in the R package
# so I need to load it from a file that is not in the package.
# Since this is in a file in my own file system, I could not
# include this as a regular example in the package.
# This example is given in full detail in the package vignette.
# You can generate GOGOA3.RData using the package 'minimalistGODB'
# or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases

GOGO3tr<-trimGOGOA3(GOGOA3,mn=2,mx=200,1)
}

}
