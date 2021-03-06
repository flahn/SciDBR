% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbops.R
\name{sort.scidb}
\alias{sort.scidb}
\title{Sort a SciDB array}
\usage{
\method{sort}{scidb}(x, decreasing = FALSE, ...)
}
\arguments{
\item{x}{a SciDB array}

\item{decreasing}{set to \code{TRUE} to sort in decreasing order}

\item{...}{optional SciDB-specific character vector of SciDB array attribute names to sort by}
}
\value{
a SciDB array
}
\description{
Sort a SciDB array
}
\examples{
\dontrun{
# Create a copy of the iris data frame in a 1-d SciDB array named "iris."
# Note that SciDB attribute names will be changed to conform to SciDB
# naming convention.
x <- as.scidb(iris,name="iris")

# Sort x by Petal_Width and Petal_Length:
a <- sort(x, attributes=c("Petal_Width","Petal_Length"))
}
}

