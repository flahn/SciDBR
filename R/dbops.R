# The functions and methods defined below are based closely on native SciDB
# functions, some of which have weak or limited analogs in R.

# SciDB rename wrapper
# Note! that the default garbage collection option here is to *not* remove.
#' Rename a SciDB array
#' @param x \code{scidb} array object
#' @param name optional new name
#' @param gc set \code{TRUE} to connect to result R's garbage collector (default is \code{FALSE})
#' @return a \code{scidb} array object
#' @export
rename = function(x, name=x@name, gc=FALSE)
{
  if(!(inherits(x, "scidb"))) stop("`x` must be a scidb object.")
  query = sprintf("rename(%s,%s)", x@name, name)
  iquery(query)
  scidb(name, gc=gc)
}

#' Remove old SciDB array versions
#' @param x \code{scidb} array object
#' @return \code{NULL}
#' @export
remove_old_versions = function(x)
{
  if(!(inherits(x, "scidb"))) stop("`x` must be a scidb object.")
  versions_query = sprintf("aggregate(versions(%s), max(version_id) as max_version)", x@name)
  versions = iqdf(versions_query, n=Inf)
  max_version = versions$max_version
  if (is.na(max_version))
  {
    stop("The array has no versions to remove")
  }
  remove_query = sprintf("remove_versions(%s, %i)", x@name, max_version)
  iquery(remove_query)
}

#' Unpack a SciDB array into a data frame-like shape
#' @param x \code{scidb} array object
#' @return a new unpacked \code{scidb} array object
#' @export
unpack = function(x)
{
  dimname = make.unique_(c(dimensions(x), scidb_attributes(x)), "i")
  query = sprintf("unpack(%s, %s)", x@name, dimname)
  ans = .scidbeval(query, FALSE, depend=list(x))
  class(ans) = "scidb"
  ans
}

#' Rename SciDB array attributes
#' @param x \code{scidb} array object
#' @param old Either a numeric vector of attribute positions, or a character vector of attribute names
#' @param new Character vector of new attribute names of same length as \code{old}
#' @return a new \code{scidb} array object
#' @export
attribute_rename = function(x, old, `new`)
{
  if(!(is.scidb(x))) stop("Requires a scidb object")
  atr = scidb_attributes(x)
  if(missing(old)) old=scidb_attributes(x)
# Positional attributes
  if(!is.numeric(old))
  {
    old = which(atr %in% old)
  }
  idx = old
  atr[idx] = `new`
  query = sprintf("cast(%s, %s%s)", x@name, build_attr_schema(x,newnames=atr),
             build_dim_schema(x))
  .scidbeval(query, FALSE, depend=list(x))
}

#' Rename SciDB array dimensions
#' @param x \code{scidb} array object
#' @param old Either a numeric vector of dimensions positions, or a character vector of dimension names
#' @param new Character vector of new dimension names of same length as \code{old}
#' @return a new \code{scidb} array object
#' @export
dimension_rename = function(x, old, `new`)
{
  if(!(is.scidb(x))) stop("Requires a scidb object")
  if(missing(old)) old = dimensions(x)
  dnames = dimensions(x)
  if(!is.numeric(old))
  {
    old = which(dnames %in% old)
  }
  idx = old
  if(length(idx)!=length(new)) stop("Invalid old dimension name specified")
  dnames[idx] = `new`
  query = sprintf("cast(%s, %s%s)", x@name, build_attr_schema(x),
             build_dim_schema(x, newnames=dnames))
  .scidbeval(query, FALSE, depend=list(x))
}

#' SciDB 'slice' operator, slice out a particular dimension value
#' @param x \code{scidb} array object
#' @param d name or numeric position of dimension to slice
#' @param n coordinate value to slice on
#' @return a new \code{scidb} array object
#' @export
slice = function(x, d, n)
{
  if(!(is.scidb(x))) stop("Requires a scidb object")
  N = length(dimensions(x))
  i = d
  if(is.character(d))
  {
    i = which(dimensions(x) %in% d)
  }
  if(length(i)==0 || i > N)
  {
    stop("Invalid dimension specified")
  }
  if(missing(n)) n = scidb_coordinate_bounds(x)$start[i]
  query = sprintf("slice(%s, %s)",x@name,paste(paste(dimensions(x)[i], noE(n), sep=","), collapse=","))
  .scidbeval(query, FALSE, depend=list(x))
}

#' SciDB 'substitute' operator, substitiute SciDB NULL values with new values
#' @param x \code{scidb} array object
#' @param value value with which to replace \code{NULL} values
#' @param attribute Either a numeric vector of attribute positions or character vector of attribute names. Substitution is limited to the specified attributes.
#' @return a new \code{scidb} array object
#' @export
replaceNA = function(x, value, `attribute`)
{
  if(!(is.scidb(x))) stop("Requires a scidb object")
  if(!any(scidb_nullable(x))) return(x)
  if(missing(attribute))
  {
    attribute = 1:length(scidb_attributes(x))
  }
  if(!is.numeric(attribute))
  {
    attribute = which(scidb_attributes(x) %in% attribute)
  }
  if(length(attribute)<1) stop("Invalid attribute(s)")
  if(missing(value))
  {
    ba = paste("build(",lapply(attribute,function(a) build_attr_schema(x,I=a,nullable=FALSE,newnames="____")),"[i=0:0,1,0]",",",.scidb_default_subst[scidb_types(x)],")",sep="")
  } else
  {
    ba = paste("build(",lapply(attribute,function(a) build_attr_schema(x,I=a,nullable=FALSE),newnames="____"),"[i=0:0,1,0]",",",value,")",sep="")
  }
  query = x@name
  i = which(scidb_types(x) == "binary")
  if(length(i) > 0)
  {
    ba[i] = "build(<____:binary>[i=0:0,1,0],'{0}[()]',true)"
  } 
  for(j in 1:length(ba))
  {
    query = sprintf("substitute(%s,%s,%s)",query,ba[j],scidb_attributes(x)[j])
  }
  .scidbeval(query, FALSE, depend=list(x))
}

#' SciDB 'subarray/between' operator, select subsets of arrays
#' @param x \code{scidb} array object
#' @param limits vector of coordinate ranges to select
#' @param between use SciDB \code{between} operator instead of \code{subarray} when \code{TRUE}
#' @return a new \code{scidb} array object
#' @export
subarray = function(x, limits, between=FALSE)
{
  if(!is.scidb(x))stop("Invalid SciDB object")
  if(missing(limits)) limits=paste(rep("null",2*length(dimensions(x))),collapse=",")
  else if(is.character(limits))
  {
# Assume user has supplied a schema string
    limits = paste(between_coordinate_bounds(limits),collapse=",")
  } else if(is.scidb(limits))
  {
# User has supplied an array
    limits = paste(between_coordinate_bounds(schema(limits)),collapse=",")
  } else
  {
# Assume a vector of limits
    if(length(limits)!=2*length(dimensions(x))) stop("Mismatched bounds length")
    limits = paste(noE(limits),collapse=",")
  }
  limits = gsub("\\*",.scidb_DIM_MAX,limits)
  if(between)
    query = sprintf("between(%s,%s)",x@name,limits)
  else
    query = sprintf("subarray(%s,%s)",x@name,limits)
  .scidbeval(query, FALSE, depend=list(x))
}

#' SciDB cast operator
#' @param x \code{scidb} array object
#' @param schema new SciDB array schema to use
#' @return a new \code{scidb} array object
#' @export
cast = function(x, schema)
{
  if(!(class(x) %in% c("scidb"))) stop("Invalid SciDB object")
  if(missing(schema)) stop("Missing cast schema")
  if(is.scidb(schema)) schema = schema(schema) # wow!
  query = sprintf("cast(%s,%s)",x@name,schema)
  .scidbeval(query, FALSE, depend=list(x))
}

#' Count the number of non-empty cells in a SciDB array
#' @param x \code{scidb} array object
#' @return count of non-empty cells
#' @export
count = function(x)
{
  if(!(class(x) %in% c("scidb"))) stop("Invalid SciDB object")
  iquery(sprintf("aggregate(%s, count(*) as count)", x@name), return=TRUE)$count
}

#' Project onto a subset of SciDB array attributes
#'
#' @param x \code{scidb} array object
#' @param attributes a character vector of array attribute names or a numeric vector of attribute positions
#' @return a new \code{scidb} array object
#' @export
project = function(x, attributes)
{
  xname = x
  if(is.logical(attributes))
    attributes = scidb_attributes(x)[which(attributes)]
  if(is.numeric(attributes))
    attributes = scidb_attributes(x)[attributes]
  if(class(x) %in% c("scidb")) xname = x@name
  query = sprintf("project(%s,%s)", xname,paste(attributes,collapse=","))
  .scidbeval(query, depend=list(x))
}

#' Flexible SciDB array subsetting
#' @param x \code{scidb} array object
#' @param expr Either a quoted SciDB filter expression, or an R expression involving array attributes and dimensions
#' @note The \code{expr} value can include scalar R values, but not more complicated expressions since the expression
#' is evaluated on the server inside SciDB (not R). Scalar R values are translated to constants in the SciDB expression
#' using R dynamic scoping/nonstandard evaluation (NSE). Quote full expressions to avoid NSE and force evaluation of
#' the quoted expression in SciDB (see examples).
#' @return a new \code{scidb} array object
#' @keywords internal
filter_scidb = function(x, expr)
{
  if(!inherits(x,"scidb")) stop("x must be a scidb object")
  xname = x@name
  isdf = "scidb" %in% class(x)
  ischar = tryCatch( is.character(expr), error=function(e) FALSE)
  if(ischar)
  {
# Check for special filter cases and adjust expr
    if(length(scidb_attributes(x)) == 2 && nchar(expr) == 1)
    {
      expr = paste(scidb_attributes(x), collapse=expr)
    }
    query = sprintf("filter(%s,%s)", xname,expr)
  } else
  {
    query = rewrite_subset_expression(substitute(expr), x, parent.frame(2))
  }
  .scidbeval(query, depend=list(x))
}


#' SciDB index lookup operator
#'
#' The \code{index_lookup} function is a wrapper to the SciDB `index_lookup` operator.
#' It produces a new SciDB dimension array that joins the unqiue indices defined in the array
#' \code{I} with values looked up in array \code{X} for attribute \code{attr}. Use
#' the \code{index_lookup} with the \code{unique} and \code{sort} functions.
#' @param x SciDB array
#' @param I single-attribute SciDB array
#' @param attr character string attribute name from the \code{x} array
#' @param new_attr optional name for the new attribute
#' @return a SciDB array
#' @note If \code{attr} is missing the first listed attribute in \code{x} will be used. If \code{I} has more than
#' one attribute, only its first listed attribute is used.
#' @examples
#' \dontrun{
#' data("iris")
#' x <- as.scidb(iris)
#' 
#' # Create a unique list of elements of the "Species" attribute.
#' # Note that we choose to defer evaluation of this expression.
#' y <- unique(sort(project(x,"Species")))
#'
#' # Append a new attribute to the array x called "Species_index" that
#' # enumerates the unique values of the "Species" attribute:
#' z <- index_lookup(x, y, "Species")
#'
#' print(head(z))
#' }
#' @export
`index_lookup` = function(x, I, attr, new_attr)
{
  if(missing(attr)) attr = scidb_attributes(x)[[1]]
  if(missing(new_attr)) new_attr=paste(attr,"index",sep="_")
  al = scidb_alias(x,I)
  xname = x
  if(class(x) %in% c("scidb")) xname=x@name
  iname = I
  if(class(I) %in% c("scidb"))
  {
    if(length(scidb_attributes(I))>1) I = project(I,1)
    if(scidb_nullable(I)) I = replaceNA(I)
    iname=I@name
  }
  query = sprintf("index_lookup(%s as %s, %s as %s, %s.%s, %s)",xname, al[1], iname, al[2], al[1], attr, new_attr)
  .scidbeval(query, depend=list(x,I))
}

# supports transform/within
bind = function(x, name, FUN)
{
  oldname = name
  aname = x
  if(class(x) %in% c("scidb")) aname=x@name
# Auto-generate names like x_n:
  if(missing(name))
  {
    name = rep("x", length(FUN))
  }
  name = make.unique_(c(scidb_attributes(x)), name)
  if(length(name) != length(FUN)) stop("name and FUN must be character vectors of identical length")
  expr = paste(paste(name, FUN, sep=","), collapse=",")
  query = sprintf("apply(%s, %s)", aname, expr)
  ans = .scidbeval(query, depend=list(x))
  if(!all(oldname == name))
  {
    drop = setdiff(oldname, name)
    a    = scidb_attributes(ans)
    ans  = project(ans, setdiff(a, drop))
    idx  = oldname != name
    a1   = name[idx]
    a2   = oldname[idx]
    ans  = attribute_rename(ans, a1, a2)
  }
  ans
}

unique_scidb = function(x, incomparables=FALSE, sort=TRUE, ...)
{
  mc = list(...)
  if(incomparables!=FALSE) warning("The incomparables option is not available yet.")
  if(any(scidb_attributes(x) %in% "i"))
  {
    new_attrs = scidb_attributes(x)
    new_attrs = new_attrs[scidb_attributes(x) %in% "i"] = make.unique_(scidb_attributes(x),"i")
    x = attribute_rename(x,scidb_attributes(x),new_attrs)
  }
  got_cu = any(grepl("^cu$", .scidbenv$ops$name))
  if(sort)
  {
    dimname = make.unique_(scidb_attributes(x),"n")
    if(length(scidb_attributes(x))>1)
    {
      if(got_cu && scidb_types(x)[[1]] == "string")
        query = sprintf("uniq(sort(cu(project(%s,%s))))",x@name,scidb_attributes(x)[[1]])
      else
        query = sprintf("uniq(sort(project(%s,%s)))",x@name,scidb_attributes(x)[[1]])
    }
    else
    {
      if(got_cu && scidb_types(x)[[1]] == "string")
        query = sprintf("uniq(sort(cu(%s)))",x@name)
      else
        query = sprintf("uniq(sort(%s))",x@name)
    }
  } else
  {
    query = sprintf("uniq(%s)", x@name)
  }
  .scidbeval(query, depend=list(x))
}

sort_scidb = function(x, decreasing=FALSE, ...)
{
  mc = list(...)
  if(!is.null(mc$na.last))
    warning("na.last option not supported by SciDB sort. Missing values are treated as less than other values by SciDB sort.")
  dflag = ifelse(decreasing, 'desc', 'asc')
# Check for ridiculous SciDB name conflict problem
  if(any(scidb_attributes(x) %in% "n"))
  {
    new_attrs = scidb_attributes(x)
    new_attrs = new_attrs[scidb_attributes(x) %in% "n"] = make.unique_(scidb_attributes(x), "n")
    x = attribute_rename(x, scidb_attributes(x), new_attrs)
  }
  if(is.null(mc$attributes))
  {
    if(length(scidb_attributes(x)) > 1)
      warning("Array contains more than one attribute, sorting on all of them.\nUse the attributes= option to restrict the sort.")
    mc$attributes = scidb_attributes(x)
  }
  `eval` = ifelse(is.null(mc$eval), FALSE, mc$eval)
  a = paste(paste(mc$attributes, dflag, sep=" "),collapse=",")
  if(!is.null(mc$chunk_size)) a = paste(a, mc$chunk_size, sep=",")

#  rs = sprintf("%s[n=0:%s,%s,0]",build_attr_schema(x,I=1),.scidb_DIM_MAX,noE(min(1e6,prod(dim(x)))))
#  query = sprintf("redimension(sort(%s,%s),%s)", x@name,a,rs)
  query = sprintf("sort(%s,%s)", x@name,a)
  .scidbeval(query,eval,depend=list(x))
}

#' SciDB \code{merge}, \code{cross_join}, and \code{join} operations.
#'
#' Only one of either \code{by} or both \code{by.x} and \code{by.y} may be
#' specified.  If none of the \code{by.x},\code{by.y} arguments are specified, and
#' \code{by=NULL} the result is the Cartesian cross product of \code{x} and
#' \code{y}.  The default value of \code{by} performs a \code{cross_join} or
#' \code{join} along common array dimensions. The \code{by} arguments may
#' be specified by name or 1-based integer dimension index.
#'
#' If only \code{by} is specified, the dimension names or attribute name in
#' \code{by} are assumed to be common across \code{x} and \code{y}.  Otherwise
#' dimension names or attribute names are matched across the names listed in
#' \code{by.x} and \code{by.y}, respectively.
#'
#' If dimension names are specified and \code{by} contains all the dimensions
#' in each array, then the SciDB \code{join} operator is used, otherwise SciDB's
#' \code{cross_join} operator is used. In each either case, the output is a cross
#' product set of the two arrays along the specified dimensions.
#'
#' If \code{by} or each of \code{by.x} and \code{by.y} list a single attribute
#' name, the indicated attributes will be lexicographically ordered as categorical
#' variables and SciDB will redimension each array along new coordinate systems
#' defined by the attributes, and then those redimensioned arrays will be joined.
#' This method limits joins along attributes to a single attribute from
#' each array. The output array will contain additional columns showing the
#' attribute factor levels used to join the arrays.
#'
#' Specify \code{merge=TRUE} to perform a SciDB merge operation instead
#' of a SciDB join.
#'
#' If \code{all=FALSE} (the default), then a SQL-like `natural join`
#' (an inner join) is performed. If \code{all=TRUE} then SQL-like `outer join`
#' is performed, but this case has some limitiations; in particular the
#' outer join is not available yet for the \code{merge=TRUE} case, for
#' joining on SciDB attributes, or for joining on subsets of dimensions.
#'
#' The various SciDB \code{join} operators generally require that the arrays have
#' identical partitioning (coordinate system bounds, chunk size, etc.) in the
#' common dimensions.  The \code{merge} method attempts to rectify SciDB
#' arrays along the specified dimensions as required before joining. Those
#' dimensions must at least have common lower index bounds.
#'
#' The merge function may rename SciDB attributes and dimensions as required
#' to avoid name conflicts in SciDB. See the last example for an
#' illustration.
#' @param x A \code{scidb} array
#' @param y A \code{scidb} array
#' @param by optional vector of common dimension or attribute names or dimension indices
#'            to join on. See details below.
#' @param ... optional additional agruments:
#' \code{fillin} is an optional argument specifying a value used to fill
#'            attributes as required by merge, it defaults to null;
#' \code{all} is an optional argument that, if TRUE, indicates outer join. It only
#'            applies in limited settings (the default is inner join);
#' \code{merge} if \code{TRUE}, perform a SciDB merge operation instead of join.
#' \code{by.x} optional vector of dimension or attribute names or dimension indices
#'                of array \code{x} to join on;
#' \code{by.y} optional vector of dimension or attribute names or dimension indices
#'                   of array \code{y} to join on;
#' @return a \code{scidb} object
#' @examples
#' \dontrun{
#' x <- as.scidb(iris,name="iris")
#'
#' a <- x$Species
#' b <- x$Petal_Length
#'
#' c <- merge(a, b, by="row")
#' merge(b, b, by="row", merge=TRUE)
#'
#'
#' # Here is an example that joins on SciDB array attributes instead of
#' # dimensions. It works by enumerating the attribute values and
#' # redimensioning along those.
#' set.seed(1)
#' a <- as.scidb(data.frame(a=sample(10, 5), b=rnorm(5)))
#' b <- as.scidb(data.frame(u=sample(10, 5), v=rnorm(5)))
#' merge(x=a, y=b, by.x="a", by.y="u")[]
#'
#'
#' # The following example joins on a subset of coordinate axes:
#' x <- build(5.5, c(3, 3));                   print(schema(x))
#' y <- build(1.1, c(3, 3), chunksize=c(2,1)); print(schema(y))
#' z <- merge(x, y, by="i")
#' print(schema(z))
#' }
#' @export
`merge.scidb` = function(x, y, by=intersect(dimensions(x), dimensions(y)), ...) merge_scidb(x, y, by, ...)

#' Sort a SciDB array
#'
#' @param x a SciDB array
#' @param decreasing set to \code{TRUE} to sort in decreasing order
#' @param ... optional SciDB-specific character vector of SciDB array attribute names to sort by
#' @return a SciDB array
#' @examples
#' \dontrun{
#' # Create a copy of the iris data frame in a 1-d SciDB array named "iris."
#' # Note that SciDB attribute names will be changed to conform to SciDB
#' # naming convention.
#' x <- as.scidb(iris,name="iris")
#'
#' # Sort x by Petal_Width and Petal_Length:
#' a <- sort(x, attributes=c("Petal_Width","Petal_Length"))
#' }
#' @export
`sort.scidb` = function(x, decreasing=FALSE, ...) sort_scidb(x, decreasing, ...)

#' Filter unique elements from a SciDB array
#' @param x SciDB array
#' @param incomparables ignored
#' @param ... optional arguments: specify \code{sort=TRUE} to first sort values
#' @return a SciDB array
#' @export
`unique.scidb` = function(x, incomparables=FALSE, ...) unique_scidb(x, incomparables, ...)

#' Filter SciDB array values or dimensions
#' @param x SciDB array
#' @param ... filter expression (see notes)
#' @note
#' Perform a SciDB \code{filter} operation on a SciDB array.  The \code{subset}
#' argument can be an R expression or a character string representing an explicit
#' SciDB filter operation.  The R expression form can include R scalar values and
#' can generate more efficient SciDB queries in some cases as shown in the
#' examples.
#'
#' When \code{subset} is an R expression, conditions involving array dimensions
#' will be translated to SciDB \code{between} statements when possible.  The R
#' expression it must use valid R syntax, although no distinction are made between
#' scalar and vector forms of logical operators.  For instance, \code{|} and
#' \code{||} are both translated to SciDB \code{or}.
#'  
#' Simple R scalars and constants may be used in R expressions and they will
#' be translated appropriately in the generated SciDB query. More complex
#' R objects like functions can't be used, however, because the logical
#' expressions are ultimately evaluated by SciDB. Dimension values are
#' treated as integer values. Values are evaulated using R dynamic scoping/
#' nonstandard evaluation (NSE). Values are evaluated in the enclosing R environments
#' first, then among the names of SciDB attributes and dimensions. Quote the entire
#' expression to avoid NSE
#' and force the expression to be evaluated verbatim in SciDB (see examples).
#'
#' Explicit grouping by parenthesis may be required to generate most
#' optimal queries when attribute and dimension conditions are mixed together
#' in an expression.
#'
#' @export
#' @return a SciDB array object
#' @examples
#' \dontrun{
#' # Create a copy of the iris data frame in a 1-d SciDB array named "iris."
#' # Variable names are changed to conform to SciDB attribute naming convention.
#' x <- as.scidb(iris)
#' # Filter the array explicitly using SciDB filter syntax
#' y <- subset(x, "Species = 'setosa'")
#' # Using an R expression form is equivalent in this example
#' z <- subset(x, Species == "setosa")
#'
#' # The R expression form can sometimes generate better-optimized SciDB
#' # expressions than the explicit form.
#' # Compare a filter involving the 'row' dimension and
#' # an attribute. Note the difference in the generated queries:
#'
#' y <- subset(x, "Species = 'setosa' and row > 40")
#' # [1] "filter(R_array5494563bc4e1101849601199,Species = 'setosa' and row > 40)"
#'
#' i <- 40
#' z <- subset(x, Species == 'setosa' & row > i)
#' # [1] "filter(between(R_array5494563bc4e1101849601199,41,null),Species = 'setosa' )"
#'
#' # Important things to note:
#' # 1. The R expression form uses R syntax.
#' # 2. The R expression form generates a SciDB query using between on
#' #    the dimensions when possible.
#' # 3. Simple R scalars may be used in the R expression form.
#' }
`subset.scidb` = function(x, ...) filter_scidb(x, ...)


#' Transform SciDB array values
#'
#' Use \code{transform} to add new derived attributes to a SciDB array, or to
#' replace an existing attribute. New attribute names must not conflict with array
#' dimension names.
#' @param _data SciDB array
#' @param ... named transformations
#' @note
#' Expressions that can't be evaluated in R are passed to SciDB as is. Explicitly
#' quote expressions to guarantee that they will be evaluated only by SciDB.
#' @return a SciDB array
#' @examples
#' \dontrun{
#' x <- scidb("build(<v:double>[i=1:5,5,0], i)")
#' transform(x, a="2 * v")
#' # Note replacement in this example:
#' transform(x, v="3 * v")
#' # Illustration of quoting expressions to force them to evaluate in SciDB:
#' v <- pi  # local R assignment of variable 'v'
#' transform(x, b=sin(v), c="sin(v)")
#' }
#' @export
`transform.scidb` = function(`_data`, ...)
{
  `_val` = as.list(match.call())[-(1:2)]
  if(length(`_val`) == 0) return()
  n = names(`_val`)
  v = unlist(Map(function(x) tryCatch(eval(x), error=function(e) x), `_val`))
  names(v) = c()
  bind(`_data`, n, v)
}
