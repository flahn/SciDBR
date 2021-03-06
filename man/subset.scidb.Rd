% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbops.R
\name{subset.scidb}
\alias{subset.scidb}
\title{Filter SciDB array values or dimensions}
\usage{
\method{subset}{scidb}(x, ...)
}
\arguments{
\item{x}{SciDB array}

\item{...}{filter expression (see notes)}
}
\value{
a SciDB array object
}
\description{
Filter SciDB array values or dimensions
}
\note{
Perform a SciDB \code{filter} operation on a SciDB array.  The \code{subset}
argument can be an R expression or a character string representing an explicit
SciDB filter operation.  The R expression form can include R scalar values and
can generate more efficient SciDB queries in some cases as shown in the
examples.

When \code{subset} is an R expression, conditions involving array dimensions
will be translated to SciDB \code{between} statements when possible.  The R
expression it must use valid R syntax, although no distinction are made between
scalar and vector forms of logical operators.  For instance, \code{|} and
\code{||} are both translated to SciDB \code{or}.
 
Simple R scalars and constants may be used in R expressions and they will
be translated appropriately in the generated SciDB query. More complex
R objects like functions can't be used, however, because the logical
expressions are ultimately evaluated by SciDB. Dimension values are
treated as integer values. Values are evaulated using R dynamic scoping/
nonstandard evaluation (NSE). Values are evaluated in the enclosing R environments
first, then among the names of SciDB attributes and dimensions. Quote the entire
expression to avoid NSE
and force the expression to be evaluated verbatim in SciDB (see examples).

Explicit grouping by parenthesis may be required to generate most
optimal queries when attribute and dimension conditions are mixed together
in an expression.
}
\examples{
\dontrun{
# Create a copy of the iris data frame in a 1-d SciDB array named "iris."
# Variable names are changed to conform to SciDB attribute naming convention.
x <- as.scidb(iris)
# Filter the array explicitly using SciDB filter syntax
y <- subset(x, "Species = 'setosa'")
# Using an R expression form is equivalent in this example
z <- subset(x, Species == "setosa")

# The R expression form can sometimes generate better-optimized SciDB
# expressions than the explicit form.
# Compare a filter involving the 'row' dimension and
# an attribute. Note the difference in the generated queries:

y <- subset(x, "Species = 'setosa' and row > 40")
# [1] "filter(R_array5494563bc4e1101849601199,Species = 'setosa' and row > 40)"

i <- 40
z <- subset(x, Species == 'setosa' & row > i)
# [1] "filter(between(R_array5494563bc4e1101849601199,41,null),Species = 'setosa' )"

# Important things to note:
# 1. The R expression form uses R syntax.
# 2. The R expression form generates a SciDB query using between on
#    the dimensions when possible.
# 3. Simple R scalars may be used in the R expression form.
}
}

