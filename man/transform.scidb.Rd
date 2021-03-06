% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbops.R
\name{transform.scidb}
\alias{transform.scidb}
\title{Transform SciDB array values}
\usage{
\method{transform}{scidb}(`_data`, ...)
}
\arguments{
\item{_data}{SciDB array}

\item{...}{named transformations}
}
\value{
a SciDB array
}
\description{
Use \code{transform} to add new derived attributes to a SciDB array, or to
replace an existing attribute. New attribute names must not conflict with array
dimension names.
}
\note{
Expressions that can't be evaluated in R are passed to SciDB as is. Explicitly
quote expressions to guarantee that they will be evaluated only by SciDB.
}
\examples{
\dontrun{
x <- scidb("build(<v:double>[i=1:5,5,0], i)")
transform(x, a="2 * v")
# Note replacement in this example:
transform(x, v="3 * v")
# Illustration of quoting expressions to force them to evaluate in SciDB:
v <- pi  # local R assignment of variable 'v'
transform(x, b=sin(v), c="sin(v)")
}
}

