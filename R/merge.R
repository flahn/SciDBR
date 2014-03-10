#
#    _____      _ ____  ____
#   / ___/_____(_) __ \/ __ )
#   \__ \/ ___/ / / / / __  |
#  ___/ / /__/ / /_/ / /_/ / 
# /____/\___/_/_____/_____/  
#
#
#
# BEGIN_COPYRIGHT
#
# This file is part of SciDB.
# Copyright (C) 2008-2014 SciDB, Inc.
#
# SciDB is free software: you can redistribute it and/or modify
# it under the terms of the AFFERO GNU General Public License as published by
# the Free Software Foundation.
#
# SciDB is distributed "AS-IS" AND WITHOUT ANY WARRANTY OF ANY KIND,
# INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY,
# NON-INFRINGEMENT, OR FITNESS FOR A PARTICULAR PURPOSE. See
# the AFFERO GNU General Public License for the complete license terms.
#
# You should have received a copy of the AFFERO GNU General Public License
# along with SciDB.  If not, see <http://www.gnu.org/licenses/agpl-3.0.html>
#
# END_COPYRIGHT
#

# SciDB cross_join wrapper internal function to support merge on various
# classes (scidb, scidbdf). This is an internal function to support
# R's merge on various SciDB objects.
#
# X and Y are SciDB array references of any kind (scidb, scidbdf)
# by is either a single character indicating a dimension name common to both
# arrays to join on, or a two-element list of character vectors of array
# dimensions to join on.
# eval=TRUE means run the query and return a scidb object.
# eval=FALSE means return a promise object representing the query.
#
# Examples:
# merge(x,y)        # Cross-product of x and y
# merge(x,y,by='i') # Natural join on common dimension i
# merge(x,y,by.x='i',by.y='i') # equiv. to the last expression
`merge_scidb` = function(x,y,...)
{
  mc = list(...)
  `by` = by.x = by.y = NULL
  `all` = FALSE
  scidbmerge = FALSE
  if(!is.null(mc$by)) `by` = mc$by
  if(!is.null(mc$all)) `all` = mc$all
  if(!is.null(mc$by.x)) by.x = mc$by.x
  if(!is.null(mc$by.y)) by.y = mc$by.y
  if(!is.null(mc$merge)) scidbmerge = mc$merge
  `eval` = ifelse(is.null(mc$eval), FALSE, mc$eval)
  xname = x@name
  yname = y@name

# Check input
  if(!is.null(`by`) && (!is.null(by.x) || !is.null(by.y)))
  {
    stop("`by` may not also be specified with `by.x` or `by.y`")
  }

# Check for simple cross case
  if(is.null(`by`) && is.null(by.x) && is.null(by.y))
  {
    query = sprintf("cross(%s, %s)",xname,yname)
    return(.scidbeval(query,eval,depend=list(x,y)))
  }

# Convert identically specified by int separate by.x by.y
  if(length(by)>0)
  {
    by.x = `by`
    by.y = `by`
  }
# Check for special join on attributes case (limited applicability)
# In particular:
# - join on only one attribute per array
# - only inner join
  if(all(by.x %in% x@attributes) && all(by.y %in% y@attributes))
  {
    by.x = by.x[[1]]  # Limitation: only one attribute for now
    by.y = by.y[[1]]  # Ditto
    lkup = unique(project(x,by.x),attributes=by.x)
    XI = index_lookup(x,lkup,by.x,`eval`=FALSE)
    YI = index_lookup(y,lkup,by.y,`eval`=FALSE)

    new_dim_name = make.unique_(c(x@D$name,y@D$name),"row")
    a = XI@attributes %in% paste(by.x,"index",sep="_")
    n = XI@attributes[a]
    redim = paste(paste(n,"=-1:",.scidb_DIM_MAX,",100000,0",sep=""), collapse=",")
    S = build_attr_schema(x, I=!(x@attributes %in% by.x))
    D = sprintf("[%s,%s]",redim,build_dim_schema(x,bracket=FALSE))
    q1 = sprintf("redimension(substitute(%s,build(<_i_:int64>[_j_=0:0,1,0],-1),%s),%s%s)",XI@name,n,S,D)

    a = YI@attributes %in% paste(by.y,"index",sep="_")
    n = YI@attributes[a]
    redim = paste(paste(n,"=-1:",.scidb_DIM_MAX,",100000,0",sep=""), collapse=",")
    S = build_attr_schema(y)
    D = sprintf("[%s,%s]",redim,build_dim_schema(y,bracket=FALSE))
    D2 = sprintf("[%s,_%s]",redim,build_dim_schema(y,bracket=FALSE))
    q2 = sprintf("cast(redimension(substitute(%s,build(<_i_:int64>[_j_=0:0,1,0],-1),%s),%s%s),%s%s)",YI@name,n,S,D,S,D2)
    query = sprintf("unpack(cross_join(%s as _cazart, %s as _yikes, _cazart.%s_index, _yikes.%s_index),%s)",q1,q2,by.x,by.y,new_dim_name)
    return(.scidbeval(query,eval,depend=list(x,y)))
  }

# New attribute schema for y that won't conflict with x:
  newas = build_attr_schema(y,newnames=make.unique_(x@attributes,y@attributes))
# Check for join case (easy case)
  if((length(by.x) == length(by.y)) && all(x@D$name %in% by.x) && all(y@D$name %in% by.y))
  {
    newds = build_dim_schema(y,newnames=x@D$name)
    castschema = sprintf("%s%s", newas, newds)
    reschema = sprintf("%s%s", newas,build_dim_schema(x))
# Cast and redimension y conformably with x:
    z = redimension(cast(y,castschema),reschema)
    if(all)
    {
# Experimental outer join XXX XXX
      x = make_nullable(x)
      z = make_nullable(z)
# Form a null-valued version of each array in the alternate array coordinate system
      xnames = make.unique_(c(z@D$name,z@attributes),x@attributes)
      vals = paste(x@types, rep("(null)",length(x@types)))
      xnull = attribute_rename(project(bind(z,xnames,vals),xnames),xnames,x@attributes)
      znames = make.unique_(c(x@D$name,x@attributes),z@attributes)
      vals = paste(z@types, rep("(null)",length(z@types)))
      znull = attribute_rename(project(bind(x,znames,vals),znames),znames,z@attributes)
# Merge each array with its nullified counterpart, then join:
      query = sprintf("join(merge(%s,%s),merge(%s,%s))",x@name,xnull@name,z@name,znull@name)
    }
    else
      if(scidbmerge)
      {
        query = sprintf("merge(%s,%s)",x@name,z@name)
      } else
      {
        query = sprintf("join(%s,%s)",x@name,z@name)
      }
    return(.scidbeval(query,eval,depend=list(x,y)))
  }

# Cross-join case (trickiest)
# Cast and redimension y conformably with x along join dimensions:
  idx.x = which(x@D$name %in% by.x)
  msk.y = y@D$name %in% by.y
  newds = lapply(1:length(y@D$name),
    function(j) {
      if(!msk.y[j])
      {
        d = build_dim_schema(y,I=j,bracket=FALSE)
      } else
      {
        k = sum(msk.y[1:j])
        d = build_dim_schema(x,I=idx.x[k],newnames=y@D$name[j],bracket=FALSE)
      }
    })
  newds = sprintf("[%s]",paste(newds,collapse=","))
  reschema = sprintf("%s%s", build_attr_schema(y),newds)
  castschema = sprintf("%s%s",newas,newds)
  z = cast(redimension(subarray(y,schema=reschema,between=TRUE),reschema),castschema)

# Join on dimensions.
  query = sprintf("cross_join(%s as __X, %s as __Y", xname, z@name)
  cterms = paste(c("__X","__Y"), as.vector(rbind(by.x,by.y)), sep=".")
  cterms = paste(cterms,collapse=",")
  query  = paste(query,",",cterms,")",sep="")
  ans = .scidbeval(query,eval,depend=list(x,y))
# Fix up bogus dimension names ... todo XXX
  ans
}