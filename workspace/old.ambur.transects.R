ambur.transects <-
function(userinput1=50, userinput2=500, userinput3=5, userinput4=5) {

require(shapefiles)
require(locfit)
require(spatstat)
require(tcltk)

# Establish the inputs
transpace <- userinput1
tranlength <- userinput2
innersample <- userinput3
outersample <- userinput4

#transpace <- 50         
#tranlength <- 500
#innersample <- 5
#outersample <- 5
####################################################################################################
###Old Spatstat functions for PPP and PSP
###########################################

#
#	ppp.R
#
#	A class 'ppp' to define point patterns
#	observed in arbitrary windows in two dimensions.
#
#	$Revision: 4.79 $	$Date: 2011/04/17 05:50:13 $
#
#	A point pattern contains the following entries:	
#
#		$window:	an object of class 'owin'
#				defining the observation window
#
#		$n:	the number of points (for efficiency)
#	
#		$x:	
#		$y:	vectors of length n giving the Cartesian
#			coordinates of the points.
#
#	It may also contain the entry:	
#
#		$marks:	a vector of length n
#			whose entries are interpreted as the
#			'marks' attached to the corresponding points.	
#	
#--------------------------------------------------------------------------
ambur.ppp <- function(x, y, ..., window, marks, check=TRUE ) {
  # Constructs an object of class 'ppp'
  #
  if(!missing(window))
    verifyclass(window, "owin")
  else
    window <- owin(...)

  if(missing(x) && missing(y)) 
    x <- y <- numeric(0)
  
  n <- length(x)
  if(length(y) != n)
    stop("coordinate vectors x and y are not of equal length")
  
  # validate x, y coordinates
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  if(any(!ok)) {
    nbg <- is.na(x) | is.na(y)
    if(any(nbg)) {
      howmany <- if(all(nbg)) "all" else paste(sum(nbg),  "out of", length(nbg))
      stop(paste(howmany, "coordinate values are NA or NaN"))
    }
    howmany <- if(!any(ok)) "all" else paste(sum(!ok),  "out of", length(ok))
    stop(paste(howmany, "coordinate values are infinite"))
  }

  names(x) <- NULL
  names(y) <- NULL
  
  # check (x,y) points lie inside window
  if(check && n > 0) {
    ok <- inside.owin(x, y, window)
    nout <- sum(!ok)
    if(nout > 0) {
      warning(paste(nout,
                    ngettext(nout, "point was", "points were"),
                    "rejected as lying outside the specified window"))
      rr <- ripras(x,y)
      bb <- bounding.box.xy(x,y)
      bb <- bounding.box(rr, bb, window)
      rejectwindow <-
        if(!is.null(rr)) rebound.owin(rr, bb) else bb
      rejects <- ppp(x[!ok], y[!ok], window=rejectwindow, check=FALSE)
      # discard illegal points
      x <- x[ok]
      y <- y[ok]
      n <- length(x)
    }
  } else nout <- 0
  # initialise ppp object
  pp <- list(window=window, n=n, x=x, y=y)
  # coerce marks to appropriate forma
  if(missing(marks))
    marks <- NULL
  if(is.hyperframe(marks)) 
    stop("Hyperframes of marks are not implemented for ppp objects; use ppx")
  if(is.matrix(marks)) 
    marks <- as.data.frame(marks)
  if(is.data.frame(marks)) {
    nc <- ncol(marks)
    if(nc == 0)
      marks <- NULL
    else if(nc == 1)
      marks <- marks[,,drop=TRUE]
  }
  # attach marks 
  if(is.null(marks)) {
    # no marks
    pp$markformat <- "none"
  } else if(is.data.frame(marks)) {
    # data frame of marks
    pp$markformat <- "dataframe"
    if(nout > 0) {
      marks <- marks[ok, ]
      marks(rejects) <- marks[!ok,]
    }
    if(nrow(marks) != n)
      stop("number of rows of marks != length of x and y")
    pp$marks <- marks
  } else {
    # should be a vector or factor
    # To recognise vector, strip attributes
    if(!is.factor(marks))
      attributes(marks) <- NULL
    if(!(is.vector(marks) || is.factor(marks)))
      stop("Format of marks not understood")
    # OK, it's a vector or factor
    pp$markformat <- "vector"
    if(nout > 0) {
      marks(rejects) <- marks[!ok]
      marks <- marks[ok]
    }
    if(length(marks) != n)
      stop("length of marks vector != length of x and y")
    names(marks) <- NULL
    pp$marks <- marks
  }
  class(pp) <- "ppp"
  if(check && any(duplicated(pp)))
    warning("data contain duplicated points")
  if(nout > 0) 
    attr(pp, "rejects") <- rejects
  pp
}

#
#--------------------------------------------------------------------------
#

is.ppp <- function(x) { inherits(x, "ppp") }

#
#--------------------------------------------------------------------------
#

as.ppp <- function(X, ..., fatal=TRUE) {
  UseMethod("as.ppp")
}

as.ppp.ppp <- function(X, ..., fatal=TRUE) {
  check <- resolve.defaults(list(...), list(check=FALSE))$check
  return(ppp(X$x, X$y, window=X$window, marks=X$marks, check=check))
}

as.ppp.quad <- function(X, ..., fatal=TRUE) {
  return(union.quad(X))
}

as.ppp.data.frame <- function(X, W = NULL, ..., fatal=TRUE) {
  check <- resolve.defaults(list(...), list(check=TRUE))$check
  if(ncol(X) < 2) 
    return(complaining("X must have at least two columns",
                       fatal, value=NULL))

  if(is.null(W))
    return(complaining("x,y coords given but no window specified",
                       fatal, value=NULL))
    
  if(is.function(W))
    Z <- cobble.xy(X[,1], X[,2], W, fatal)
  else {
    win <- as.owin(W)
    Z <- ppp(X[,1], X[,2], window = win, check=check)
  }

  # add marks from other columns
  if(ncol(X) > 2)
    marks(Z) <- X[, -(1:2)]

  return(Z)
}
    
as.ppp.matrix <- function(X, W = NULL, ..., fatal=TRUE) {
  check <- resolve.defaults(list(...), list(check=TRUE))$check
  if(!verifyclass(X, "matrix", fatal=fatal)
     || !is.numeric(X))
    return(complaining("X must be a numeric matrix",
                       fatal, value=NULL))

  if(ncol(X) < 2)
    return(complaining("X must have at least two columns",
                       fatal, value=NULL))

  if(is.null(W))
    return(complaining("x,y coords given but no window specified",
                       fatal, value=NULL))
    
  if(is.function(W))
    Z <- cobble.xy(X[,1], X[,2], W, fatal)
  else {
    win <- as.owin(W)
    Z <- ppp(X[,1], X[,2], window = win, check=check)
  }

  # add marks from other columns
  if(ncol(X) > 2)
    marks(Z) <- X[, -(1:2)]

  return(Z)
}
    
as.ppp.default <- function(X, W=NULL, ..., fatal=TRUE) {
	# tries to coerce data X to a point pattern
	# X may be:
	#	1. a structure with entries x, y, xl, xu, yl, yu
	#	2. a structure with entries x, y, area where
        #                    'area' has entries xl, xu, yl, yu
	#	3. a structure with entries x, y
        #       4. a vector of length 2, interpreted as a single point.
	# The second argument W is coerced to an object of class 'owin' by the 
	# function "as.owin" in window.S
        # If X also has an entry X$marks
        # then this will be interpreted as the marks vector for the pattern.
	#
  check <- resolve.defaults(list(...), list(check=TRUE))$check
  if(checkfields(X, c("x", "y", "xl", "xu", "yl", "yu"))) {
		xrange <- c(X$xl, X$xu)
		yrange <- c(X$yl, X$yu)
		if(is.null(X$marks))
			Z <- ppp(X$x, X$y, xrange, yrange, check=check)
		else
			Z <- ppp(X$x, X$y, xrange, yrange, 
				marks=X$marks, check=check)
		return(Z)
        } else if(checkfields(X, c("x", "y", "area"))
                  && checkfields(X$area, c("xl", "xu", "yl", "yu"))) {
                win <- as.owin(X$area)
                if (is.null(X$marks))
                  Z <- ppp(X$x, X$y, window=win, check=check)
                else
                  Z <- ppp(X$x, X$y, window=win, marks = X$marks, check=check)
                return(Z)
	} else if(checkfields(X, c("x", "y"))) {
                if(is.function(W))
                  return(cobble.xy(X$x, X$y, W, fatal))
		if(is.null(W)) {
                  if(fatal)
                    stop("x,y coords given but no window specified")
                  else
                    return(NULL)
                }
		win <- as.owin(W)
		if(is.null(X$marks))
                  Z <- ppp(X$x, X$y, window=win, check=check)
                else
                  Z <- ppp(X$x, X$y, window=win, marks=X$marks, check=check)
                return(Z)
        } else if(is.vector(X) && length(X) == 2) {
                win <- as.owin(W)
                Z <- ppp(X[1], X[2], window=win, check=check)
                return(Z)
	} else {
          if(fatal)
            stop("Can't interpret X as a point pattern")
          else
            return(NULL)
        }
}

cobble.xy <- function(x, y, f=ripras, fatal=TRUE) {
  if(!is.function(f))
    stop("f is not a function")
  w <- f(x,y)
  if(!is.owin(w)) {
    gripe <- "Supplied function f did not return an owin object"
    if(fatal)
      stop(gripe)
    else {
      warning(gripe)
      return(NULL)
    }
  }
  return(ppp(x, y, window=w))
}
  

# --------------------------------------------------------------

"[.ppp" <-
  function(x, i, j, drop, ...) {

        verifyclass(x, "ppp")
        
        if(missing(i) && missing(j))
          return(x)

        if(!missing(i)) {
          if(inherits(i, "owin")) {
            # i is a window
            window <- i
            ok <- inside.owin(x$x, x$y, window)
            x <- ppp(x$x[ok], x$y[ok], window=window, #SIC
                     marks=marksubset(x$marks, ok),
                     check=FALSE)
          } else if(inherits(i, "im")) {
            # i is an image
            if(i$type != "logical")
              stop(paste("Subset operator X[i] undefined",
                         "when i is a pixel image",
                         "unless it has logical values"), call.=FALSE)
            # convert logical image to window
            e <- sys.frame(sys.nframe())
            window <- solutionset(i, e)
            ok <- inside.owin(x$x, x$y, window)
            x <- ppp(x$x[ok], x$y[ok], window=window, #SIC
                     marks=marksubset(x$marks, ok),
                     check=FALSE)
          } else {
            # assume i is a subset index
            if(x$n == 0)
              return(x)
            subset <- i
            x <- ppp(x$x[subset], x$y[subset], window=x$window,
                     marks=marksubset(x$marks, subset),
                     check=FALSE)
          } 
        }

        if(!missing(j))
          x <- x[j]   # invokes code above

        return(x)
}


# ------------------------------------------------------------------
#
#
scanpp <- function(filename, window, header=TRUE, dir="", multitype=FALSE) {
  filename <- if(dir=="") filename else
              paste(dir, filename, sep=.Platform$file.sep)
  df <- read.table(filename, header=header)
  if(header) {
    # check whether there are columns named 'x' and 'y'
    colnames <- dimnames(df)[[2]]
    xycolumns <- match(c("x", "y"), colnames, 0)
    named <- all(xycolumns > 0)
  } else {
    named <- FALSE
  }
  if(named) {
    x <- df$x
    y <- df$y
  } else {
    # assume x, y given in columns 1, 2 respectively
    x <- df[,1]
    y <- df[,2]
    xycolumns <- c(1,2)
  }
  if(ncol(df) == 2) 
      X <- ppp(x, y, window=window)
  else {
    marks <- df[ , -xycolumns]
    if(multitype) 
      marks <- factor(marks)
    X <- ppp(x, y, window=window, marks = marks)
  }
  X
}

#-------------------------------------------------------------------

"markspace.integral" <-
  function(X) {
  verifyclass(X, "ppp")
  if(!is.marked(X, dfok=TRUE))
    return(1)
  if(is.multitype(X))
    return(length(levels(marks(X))))
  else
    stop("Don't know how to compute total mass of mark space")
}

#-------------------------------------------------------------------

print.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  ism <- is.marked(x, dfok=TRUE)
  cat(paste(if(ism) "marked" else NULL,
            "planar point pattern:",
            x$n,
            ngettext(x$n, "point", "points"),
            "\n"))
  if(ism) {
    mks <- marks(x, dfok=TRUE)
    if(is.data.frame(mks)) {
      # data frame of marks
      cat(paste("Mark variables: ", paste(names(mks), collapse=", "), "\n"))
    } else {
      # vector of marks
      if(is.factor(mks)) {
        cat("multitype, with ")
        cat(paste("levels =", paste(levels(mks), collapse="\t"),"\n"))
      } else {
        cat(paste("marks are",
                  if(is.numeric(mks)) "numeric,",
                  "of type", sQuote(typeof(mks)), "\n"))
      }
    }
  }
  print(x$window)
  if(!is.null(rejects <- attr(x, "rejects"))) {
    nrejects <- rejects$n
    cat(paste("\n***",
              nrejects,
              ngettext(nrejects, "illegal point", "illegal points"),
              "stored in",
              paste("attr(,", dQuote("rejects"), ")", sep=""),
              "***\n"))
  }
  return(invisible(NULL))
}

summary.ppp <- function(object, ..., checkdup=TRUE) {
  verifyclass(object, "ppp")
  result <- list()
  result$is.marked <- is.marked(object, dfok=TRUE)
  result$n <- object$n
  result$window <- summary(object$window)
  result$intensity <- result$n/result$window$area
  if(checkdup)
    result$nduplicated <- sum(duplicated(object))
  if(result$is.marked) {
    mks <- marks(object, dfok=TRUE)
    if(result$multiple.marks <- is.data.frame(mks)) {
      result$marknames <- names(mks)
      result$is.numeric <- FALSE
      result$marktype <- "dataframe"
      result$is.multitype <- FALSE
    } else {
      result$is.numeric <- is.numeric(mks)
      result$marknames <- "marks"
      result$marktype <- typeof(mks)
      result$is.multitype <- is.multitype(object)
    }
    if(result$is.multitype) {
      tm <- as.vector(table(mks))
      tfp <- data.frame(frequency=tm,
                        proportion=tm/sum(tm),
                        intensity=tm/result$window$area,
                        row.names=levels(mks))
      result$marks <- tfp
    } else 
      result$marks <- summary(mks)
  }
  class(result) <- "summary.ppp"
  if(!is.null(rejects <- attr(object, "rejects"))) 
    result$rejects <- rejects$n
  return(result)
}

print.summary.ppp <- function(x, ..., dp=3) {
  verifyclass(x, "summary.ppp")
  cat(paste(if(x$is.marked) "Marked planar " else "Planar ",
            "point pattern: ",
            x$n,
            " points\n",
            sep=""))
  oneline <- resolve.defaults(list(...), list(oneline=FALSE))$oneline
  if(oneline) return(invisible(NULL))
  unitinfo <- summary(x$window$units)
  cat(paste("Average intensity",
            signif(x$intensity,dp),
            "points per square",
            unitinfo$singular,
            unitinfo$explain,
            "\n"))
  ndup <- x$nduplicated
  if((!is.null(ndup)) && (ndup > 0))
    cat("\n*Pattern contains duplicated points*\n")
  if(x$is.marked) {
    if(x$multiple.marks) {
      cat(paste("Mark variables: ", paste(x$marknames, collapse=", "), "\n"))
      cat("Summary:\n")
      print(x$marks)
    } else if(x$is.multitype) {
      cat("Multitype:\n")
      print(signif(x$marks,dp))
    } else {
      cat(paste("marks are ",
                if(x$is.numeric) "numeric, ",
                "of type", sQuote(x$marktype), "\n"))
      cat("Summary:\n")
      print(x$marks)
    }
  }
  cat("\n")
  print(x$window)
  if(!is.null(nrejects <- x$rejects)) 
    cat(paste("\n***",
              nrejects,
              ngettext(nrejects, "illegal point", "illegal points"),
              "stored in",
              paste("attr(,", dQuote("rejects"), ")", sep=""),
              "***\n"))
  return(invisible(x))
}

# ---------------------------------------------------------------

identify.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  if(!is.marked(x) || "labels" %in% names(list(...)))
    identify(x$x, x$y, ...)
  else {
    marx <- marks(x, dfok=FALSE)
    marques <- if(is.numeric(marx)) paste(signif(marx, 3)) else paste(marx)
    id <- identify(x$x, x$y, labels=marques, ...)
    mk <- marx[id]
    if(is.factor(marx)) mk <- levels(marx)[mk]
    cbind(id=id, marks=mk)
  }
}

rebound <- function(x, rect) {
  UseMethod("rebound")
}

rebound.ppp <- function(x, rect) {
  verifyclass(x, "ppp")
  x$window <- rebound.owin(x$window, rect)
  return(x)
}

as.data.frame.ppp <- function(x, row.names=NULL, ...) {
  df <- data.frame(x=x$x, y=x$y, row.names=row.names)
  marx <- marks(x, dfok=TRUE)
  if(is.null(marx))
    return(df)
  if(is.data.frame(marx))
    df <- cbind(df, marx)
  else
    df <- data.frame(df, marks=marx)
  return(df)
}

is.empty.ppp <- function(x) { return(x$n == 0) }

npoints <- function(x) {
  UseMethod("npoints")
}

nobjects <- function(x) {
  UseMethod("nobjects")
}

nobjects.ppp <- npoints.ppp <- function(x) { x$n }
##########################################################################################
#
#  psp.R
#
#  $Revision: 1.61 $ $Date: 2011/05/19 06:14:18 $
#
# Class "psp" of planar line segment patterns
#
#
#################################################
# creator
#################################################
ambur.psp <- function(x0, y0, x1, y1, window, marks=NULL,
                check=spatstat.options("checksegments")) {
  stopifnot(is.numeric(x0))
  stopifnot(is.numeric(y0))
  stopifnot(is.numeric(x1))
  stopifnot(is.numeric(y1))
  stopifnot(is.vector(x0))
  stopifnot(is.vector(y0))
  stopifnot(is.vector(x1))
  stopifnot(is.vector(y1))
  stopifnot(length(x0) == length(y0))
  stopifnot(length(x1) == length(y1))
  stopifnot(length(x0) == length(x1))
  ends <- data.frame(x0=x0,y0=y0,x1=x1,y1=y1)
  if(!missing(window))
    verifyclass(window,"owin")
  if(check) {
    ok <- inside.owin(x0,y0, window) & inside.owin(x1,y1,window)
    if((nerr <- sum(!ok)) > 0)
      stop(paste(nerr, ngettext(nerr, "segment does not", "segments do not"),
                 "lie entirely inside the window.\n"), call.=FALSE)
  }
  out <- list(ends=ends,
              window=window,
              n = nrow(ends))

# add marks if any
  if(!is.null(marks)) {
    if(is.matrix(marks))
      marks <- as.data.frame(marks)
    if(is.data.frame(marks)) {
      omf <- "dataframe"
      nmarks <- nrow(marks)
      rownames(marks) <- seq_len(nmarks)
      whinge <- "The number of rows of marks"
    } else {
      omf <- "vector"
      names(marks) <- NULL
      nmarks <- length(marks)
      whinge <- "The length of the marks vector"
    }
    if(nmarks != out$n) stop(paste(whinge, "!= length of x and y.\n"))
    out$marks <- marks
    out$markformat <- omf
  } else {
    out$markformat <- "none"
  }

  class(out) <- c("psp", class(out))
  return(out)
}

######################################################
#  conversion
######################################################

is.psp <- function(x) { inherits(x, "psp") }

as.psp <- function(x, ..., from=NULL, to=NULL) {
  # special case: two point patterns
  if(is.null(from) != is.null(to))
    stop(paste("If one of", sQuote("from"), "and", sQuote("to"),
               "is specified, then both must be specified.\n"))
  if(!is.null(from) && !is.null(to)) {
    verifyclass(from, "ppp")
    verifyclass(to, "ppp")
    if(from$n != to$n)
      stop(paste("The point patterns", sQuote("from"), "and", sQuote("to"),
                 "have different numbers of points.\n"))
    uni <- union.owin(from$window, to$window)
    Y <- do.call("psp",
                 resolve.defaults(list(from$x, from$y, to$x, to$y),
                                  list(...),
                                  list(window=uni)))
    return(Y)
  }
  UseMethod("as.psp")
}

as.psp.psp <- function(x, ..., check=FALSE, fatal=TRUE) {
  if(!verifyclass(x, "psp", fatal=fatal))
    return(NULL)
  ends <- x$ends
  psp(ends$x0, ends$y0, ends$x1, ends$y1, window=x$window,
      marks=x$marks, check=check)
}

as.psp.data.frame <- function(x, ..., window=NULL, marks=NULL,
                              check=spatstat.options("checksegments"), fatal=TRUE) {
  window <- suppressWarnings(as.owin(window,fatal=FALSE))
  if(!is.owin(window)) {
    if(fatal) stop("Cannot interpret \"window\" as an object of class owin.\n")
    return(NULL)
  }

  if(checkfields(x,"marks")) {
    if(is.null(marks)) marks <- x$marks
    else warning(paste("Column named \"marks\" ignored;\n",
                       "argument named \"marks\" has precedence.\n",sep=""))
    x$marks <- NULL
  }

  if(checkfields(x, c("x0", "y0", "x1", "y1"))) {
    out <- psp(x$x0, x$y0, x$x1, x$y1, window=window,
               check=check)
    x <- x[-match(c("x0","y0","x1","y1"),names(x))]
  }
  else if(checkfields(x, c("xmid", "ymid", "length", "angle"))) {
    rr <- x$length/2
    dx <- cos(x$angle) * rr
    dy <- sin(x$angle) * rr
    bb <- bounding.box(window)
    rmax <- max(rr)
    bigbox <- owin(bb$xrange + c(-1,1) * rmax, bb$yrange + c(-1,1) * rmax)
    pattern <- psp(x$x - dx, x$y - dy, x$x + dx, x$y + dy,
                   window=bigbox,check=FALSE)
    out <- pattern[window]
    x <- x[-match(c("xmid","ymid","length","angle"),names(x))]
  }
  else if(ncol(x) >= 4) {
    out <- psp(x[,1], x[,2], x[,3], x[,4], window=window,
               check=check)
    x <- x[-(1:4)]
  }
  else if(fatal)
    stop("Unable to interpret x as a line segment pattern.\n")
  else out <- NULL

  if(!is.null(out)) {
    if(is.null(marks) & ncol(x) > 0) marks <- x
    if(is.null(marks)) {
       out$markformat <- "none"
    } else {
       out$marks <- marks
       out$markformat <- if(is.data.frame(marks)) "dataframe" else "vector"
       out <- as.psp(out,check=FALSE)
    }
  }
  return(out)
}

as.psp.matrix <- function(x, ..., window=NULL, marks=NULL,
                          check=spatstat.options("checksegments"), fatal=TRUE) {
   x <- as.data.frame(x)
   as.psp(x,...,window=window,marks=marks,check=check,fatal=fatal)
}

as.psp.default <- function(x, ..., window=NULL, marks=NULL,
                           check=spatstat.options("checksegments"), fatal=TRUE) {
  if(checkfields(x,"marks")) {
	if(is.null(marks)) marks <- x$marks
	else warning(paste("Component of \"x\" named \"marks\" ignored;\n",
                             "argument named \"marks\" has precedence.\n",sep=""))
  }
  if(checkfields(x, c("x0", "y0", "x1", "y1")))
    return(psp(x$x0, x$y0, x$x1, x$y1, window=window, marks=marks,
               check=check))
  else if(checkfields(x, c("xmid", "ymid", "length", "angle"))) {
    rr <- x$length/2
    dx <- cos(x$angle) * rr
    dy <- sin(x$angle) * rr
    window <- as.owin(window)
    bb <- bounding.box(window)
    rmax <- max(rr)
    bigbox <- owin(bb$xrange + c(-1,1) * rmax, bb$yrange + c(-1,1) * rmax)
    pattern <- psp(x$x - dx, x$y - dy, x$x + dx, x$y + dy,
                   window=bigbox, marks=marks, check=FALSE)
    clipped <- pattern[window]
    return(clipped)
  }
  else if(fatal)
    stop("Unable to interpret x as a line segment pattern")
  return(NULL)
}

as.psp.owin <- function(x, ..., check=spatstat.options("checksegments"), fatal=TRUE) {
  verifyclass(x, "owin")
  # can't use as.rectangle here; still testing validity
  xframe <- owin(x$xrange, x$yrange)
  switch(x$type,
         rectangle = {
           xx <- x$xrange[c(1,2,2,1)]
           yy <- x$yrange[c(1,1,2,2)]
           nxt <- c(2,3,4,1)
           out <- psp(xx, yy, xx[nxt], yy[nxt], window=x, check=check)
           return(out)
         },
         polygonal = {
           x0 <- y0 <- x1 <- y1 <- numeric(0)
           bdry <- x$bdry
           for(i in seq_along(bdry)) {
             po <- bdry[[i]]
             ni <- length(po$x)
             nxt <- c(2:ni, 1)
             x0 <- c(x0, po$x)
             y0 <- c(y0, po$y)
             x1 <- c(x1, po$x[nxt])
             y1 <- c(y1, po$y[nxt])
           }
           out <- psp(x0, y0, x1, y1,  window=xframe, check=check)
           return(out)
         },
         mask = {
           if(fatal) stop("x is a mask")
           else warning("x is a mask - no line segments returned")
           return(psp(numeric(0), numeric(0), numeric(0), numeric(0),
                      window=xframe, check=FALSE))
         })
  return(NULL)
}


#################

as.data.frame.psp <- function(x, row.names=NULL, ...) {
  df <- as.data.frame(x$ends, row.names=row.names)
  if(is.marked(x))
    df <- cbind(df, if(x$markformat=="dataframe") marks(x)
                    else data.frame(marks=marks(x)))
  return(df)
}

#######  manipulation ##########################

append.psp <- function(A,B) {
  verifyclass(A, "psp")
  verifyclass(B, "psp")
  stopifnot(identical(A$window, B$window))
  marks <- marks(A) %mapp% marks(B)
  ends <- rbind(A$ends, B$ends)
  out  <- as.psp(ends,window=A$window,marks=marks,check=FALSE)
  return(out)
}

rebound.psp <- function(x, rect) {
  verifyclass(x, "psp")
  x$window <- rebound.owin(x$window, rect)
  return(x)
}


#################################################
#  marks
#################################################

is.marked.psp <- function(X, ...) {
  marx <- marks(X, ...)
  return(!is.null(marx))
}

marks.psp <- function(x, ..., dfok = TRUE) {
  # data frames of marks are as of 19/March 2011 implemented for psp
    ma <- x$marks
    if ((is.data.frame(ma) || is.matrix(ma)) && !dfok) 
        stop("Sorry, not implemented when the marks are a data frame.\n")
    return(ma)
}

"marks<-.psp" <- function(x, ..., value) {
  stopifnot(is.psp(x))
  if(is.null(value)) {
    return(unmark(x))
  }
  m <- value
  if(!(is.vector(m) || is.factor(m) || is.data.frame(m) || is.matrix(m)))
    stop("Incorrect format for marks")

    if (is.hyperframe(m)) 
        stop("Hyperframes of marks are not supported in psp objects.\n")
    nseg <- nsegments(x)
    if (!is.data.frame(m) && !is.matrix(m)) {
        if (length(m) == 1) 
            m <- rep(m, nseg)
        else if (nseg == 0) 
            m <- rep(m, 0)
        else if (length(m) != nseg) 
            stop("Number of marks != number of line segments.\n")
        marx <- m
    }
    else {
        m <- as.data.frame(m)
        if (ncol(m) == 0) {
            marx <- NULL
        }
        else {
            if (nrow(m) == nseg) {
                marx <- m
            }
            else {
                if (nrow(m) == 1 || nseg == 0) {
                  marx <- as.data.frame(lapply(as.list(m),function(x,k) {
                    rep(x, k)}, k = nseg))
                }
                else stop("Number of rows of data frame != number of points.\n")
            }
        }
    }
    Y <- as.psp(x$ends, window = x$window, marks = marx, check = FALSE)
    return(Y)
}

markformat.psp <- function(x) {
    mf <- x$markformat
    if(is.null(mf)) 
      mf <- markformat(marks(x))
    return(mf)
}

unmark.psp <- function(X) {
  X$marks <- NULL
  X$markformat <- "none"
  return(X)
}

#################################################
#  plot and print methods
#################################################

plot.psp <- function(x, ..., add=FALSE, which.marks=1,
                     ribbon=TRUE, ribsep=0.15, ribwid=0.05, ribn=1024) {
  main <- deparse(substitute(x))
  verifyclass(x, "psp")
  #
  n <- nsegments(x)
  marx <- marks(x)
  #
  use.colour <- !is.null(marx) && (n != 0)
  do.ribbon <- identical(ribbon, TRUE) && use.colour && !add
  #
  if(!add) {
    # create plot region
    if(!do.ribbon) {
      # window of x
      do.call.matched("plot.owin", 
                      resolve.defaults(list(x=x$window),
                                       list(...),
                                       list(main=deparse(substitute(x)))))
    } else {
      # enlarged window with room for colour ribbon
      # x at left, ribbon at right
      bb <- as.rectangle(as.owin(x))
      xwidth <- diff(bb$xrange)
      xheight <- diff(bb$yrange)
      xsize <- max(xwidth, xheight)
      bb.rib <- owin(bb$xrange[2] + c(ribsep, ribsep+ribwid) * xsize,
                     bb$yrange)
      bb.all <- bounding.box(bb.rib, bb)
      # establish coordinate system
      do.call.matched("plot.default",
                      resolve.defaults(list(x=0, y=0, type="n",
                                            axes=FALSE, asp=1,
                                            xlim=bb.all$xrange,
                                            ylim=bb.all$yrange),
                                       list(...),
                                       list(main=main, xlab="", ylab="")))
      # now plot window of x
      do.call.matched("plot.owin", 
                      resolve.defaults(list(x=x$window, add=TRUE),
                                       list(...)))
    }
  }

  # plot segments
  if(n == 0)
    return(invisible(NULL))
  
  # determine colours if any
  if(!use.colour) {
    # black
    col <- colmap <- NULL
  } else {
    # multicoloured 
    marx <- as.data.frame(marx)[, which.marks]
    if(is.character(marx) || length(unique(marx)) == 1)
      marx <- factor(marx)
    if(is.factor(marx)) {
      lev <- levels(marx)
      colmap <- colourmap(col=rainbow(length(lev)), inputs=factor(lev))
    } else {
      colmap <- colourmap(col=rainbow(ribn), range=range(marx))
    }
    col <- colmap(marx)
  }
  # plot segments
  do.call("segments",
          resolve.defaults(as.list(x$ends),
                           list(...),
                           list(col=col),
                           .StripNull=TRUE))
  # plot ribbon
  if(do.ribbon) 
    plot(colmap, vertical=TRUE, add=TRUE,
         xlim=bb.rib$xrange, ylim=bb.rib$yrange)
  
  # return colour map
  return(invisible(colmap))
}

print.psp <- function(x, ...) {
   verifyclass(x, "psp")
    ism <- is.marked(x, dfok = TRUE)
    cat(paste(if(ism)
        "marked"
    else NULL, "planar line segment pattern:", x$n, "line", ngettext(x$n, "segment",
        "segments"), "\n"))
    if (ism) {
        mks <- marks(x, dfok = TRUE)
        if (is.data.frame(mks)) {
            cat(paste("Mark variables: ", paste(names(mks), collapse = ", "),
                "\n"))
        }
        else {
            if (is.factor(mks)) {
                cat("multitype, with ")
                cat(paste("levels =", paste(levels(mks), collapse = "\t"), 
                  "\n"))
            }
            else {
                cat(paste("marks are", if (is.numeric(mks)) 
                  "numeric,", "of type", sQuote(typeof(mks)), 
                  "\n"))
            }
        }
    }
    print(x$window)
    return(invisible(NULL))
}

unitname.psp <- function(x) {
  return(unitname(x$window))
}

"unitname<-.psp" <- function(x, value) {
  w <- x$window
  unitname(w) <- value
  x$window <- w
  return(x)
}

####################################################
#    summary information
####################################################

endpoints.psp <- function(x, which="both") {
  verifyclass(x, "psp")
  ends <- x$ends
  n <- x$n
  switch(which,
         both={
           first <- second <- rep(TRUE, n)
         },
         first={
           first <- rep(TRUE, n)
           second <- rep(FALSE, n)
         },
         second={
           first <- rep(FALSE, n)
           second <- rep(TRUE, n)
         },
         left={
           first <- (ends$x0 < ends$x1)
           second <- !first
         },
         right={
           first <- (ends$x0 > ends$x1)
           second <- !first
         },
         lower={
           first <- (ends$y0 < ends$y1)
           second <- !first
         },
         upper={
           first <- (ends$y0 > ends$y1)
           second <- !first
         },
         stop(paste("Unrecognised option: which=", sQuote(which)))
         )
  ok <- rbind(first, second)
  xmat <- rbind(ends$x0, ends$x1)
  ymat <- rbind(ends$y0, ends$y1)
  idmat <- col(ok)
  xx <- as.vector(xmat[ok])
  yy <- as.vector(ymat[ok])
  id <- as.vector(idmat[ok])
  result <- ppp(xx, yy, window=x$window, check=FALSE)
  attr(result, "id") <- id
  return(result)
}

midpoints.psp <- function(x) {
  verifyclass(x, "psp")
  xm <- eval(expression((x0+x1)/2), envir=x$ends)
  ym <- eval(expression((y0+y1)/2), envir=x$ends)
  win <- x$window
  ok <- inside.owin(xm, ym, win)
  if(any(!ok)) {
    warning(paste("Some segment midpoints lie outside the original window;",
                  "window replaced by bounding box"))
    win <- bounding.box(win)
  }
  ppp(x=xm, y=ym, window=win, check=FALSE)
}

lengths.psp <- function(x) {
  verifyclass(x, "psp")
  eval(expression(sqrt((x1-x0)^2 + (y1-y0)^2)), envir=x$ends)
}

angles.psp <- function(x, directed=FALSE) {
  verifyclass(x, "psp")
  a <- eval(expression(atan2(y1-y0, x1-x0)), envir=x$ends)
  if(!directed) 
    a <- a %% pi
  return(a)
}

summary.psp <- function(object, ...) {
  verifyclass(object, "psp")
  len <- lengths.psp(object)
  out <- list(n = object$n,
              len = summary(len),
              totlen = sum(len),
              ang= summary(angles.psp(object)),
              w = summary.owin(object$window),
              marks=if(is.null(object$marks)) NULL else summary(object$marks),
              unitinfo=summary(unitname(object)))
  class(out) <- c("summary.psp", class(out))
  return(out)
}

print.summary.psp <- function(x, ...) {
  cat(paste(x$n, "line segments\n"))
  cat("Lengths:\n")
  print(x$len)
  unitblurb <- paste(x$unitinfo$plural, x$unitinfo$explain)
  cat(paste("Total length:", x$totlen, unitblurb, "\n"))
  cat(paste("Length per unit area:", x$totlen/x$w$area, "\n"))
  cat("Angles (radians):\n")
  print(x$ang)
  print(x$w)
  if(!is.null(x$marks)) {
    cat("Marks:\n")
    print(x$marks)
  }
  return(invisible(NULL))
}

  
########################################################
#  subsets
########################################################

"[.psp" <-
  function(x, i, j, drop, ...) {

    verifyclass(x, "psp")
    
    if(missing(i) && missing(j))
      return(x)
        
    if(!missing(i)) {
      style <- if(inherits(i, "owin")) "window" else "index"
      switch(style,
             window={
               x <- clip.psp(x, window=i)
             },
             index={
               subset <- i
               markformat <- markformat(x)
               x <- as.psp(x$ends[subset, ],
                           window=x$window,
                           marks=switch(markformat,
                             none=NULL,
                             vector=x$marks[subset],
                             dataframe=x$marks[subset,]))
             })
    }

    if(!missing(j))
      x <- x[j] # invokes code above
    
    return(x)
 }
  


####################################################
# affine transformations
####################################################

affine.psp <- function(X,  mat=diag(c(1,1)), vec=c(0,0), ...) {
  verifyclass(X, "psp")
  W <- affine.owin(X$window, mat=mat, vec=vec)
  E <- X$ends
  ends0 <- affinexy(list(x=E$x0,y=E$y0), mat=mat, vec=vec)
  ends1 <- affinexy(list(x=E$x1,y=E$y1), mat=mat, vec=vec)
  psp(ends0$x, ends0$y, ends1$x, ends1$y, window=W, marks=marks(X, dfok=TRUE),
      check=FALSE)
}

shift.psp <- function(X, vec=c(0,0), ..., origin=NULL) {
  verifyclass(X, "psp")
  if(!is.null(origin)) {
    stopifnot(is.character(origin))
    if(!missing(vec))
      warning("Argument vec ignored; argument origin has precedence.\n")
    origin <- pickoption("origin", origin, c(centroid="centroid",
                                             midpoint="midpoint",
                                             bottomleft="bottomleft"))
    W <- as.owin(X)
    locn <- switch(origin,
                   centroid={ unlist(centroid.owin(W)) },
                   midpoint={ c(mean(W$xrange), mean(W$yrange)) },
                   bottomleft={ c(W$xrange[1], W$yrange[1]) })
    return(shift(X, -locn))
  }
  W <- shift.owin(X$window, vec=vec, ...)
  E <- X$ends
  ends0 <- shiftxy(list(x=E$x0,y=E$y0), vec=vec, ...)
  ends1 <- shiftxy(list(x=E$x1,y=E$y1), vec=vec, ...)
  psp(ends0$x, ends0$y, ends1$x, ends1$y, window=W, marks=marks(X, dfok=TRUE),
      check=FALSE)
}

rotate.psp <- function(X, angle=pi/2, ...) {
  verifyclass(X, "psp")
  W <- rotate.owin(X$window, angle=angle, ...)
  E <- X$ends
  ends0 <- rotxy(list(x=E$x0,y=E$y0), angle=angle, ...)
  ends1 <- rotxy(list(x=E$x1,y=E$y1), angle=angle, ...)
  psp(ends0$x, ends0$y, ends1$x, ends1$y, window=W, marks=marks(X, dfok=TRUE),
      check=FALSE)
}

is.empty.psp <- function(x) { return(x$n == 0) } 

identify.psp <- function(x, ..., labels=seq_len(nsegments(x)), n=nsegments(x), plot=TRUE) {
  Y <- x
  W <- as.owin(Y)
  mids <- midpoints.psp(Y)
  if(!(is.numeric(n) && (length(n) == 1) && (n %% 1 == 0) && (n >= 0)))
    stop("n should be a single integer")
  out <- integer(0)
  while(length(out) < n) {
    xy <- locator(1)
    # check for interrupt exit
    if(length(xy$x) == 0)
      return(out)
    # find nearest segment
    X <- ppp(xy$x, xy$y, window=W)
    ident <- project2segment(X, Y)$mapXY
    # add to list
    if(ident %in% out) {
      cat(paste("Segment", ident, "already selected\n"))
    } else {
      if(plot) {
        # Display
        mi <- mids[ident]
        li <- labels[ident]
        text(mi$x, mi$y, labels=li)
      }
      out <- c(out, ident)
    }
  }
  # exit if max n reached
  return(out)
}

nsegments <- function(x) {
	UseMethod("nsegments")
}

nobjects.psp <- nsegments.psp <- function(x) {
   x$n
}

as.ppp.psp <- function (X, ..., fatal=TRUE) 
{
  Y <- endpoints.psp(X, which="both")
  m  <- marks(X)
  marks(Y) <- markappend(m, m)
  return(Y)
}


 #######################
 ### BEGIN ACTUAL FUNCTION
#######################################################################################################################################
tkmessageBox(message = "Please select the inner baseline...")

path1 <- tk_choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path1)

path1a <- paste(substr(path1,1,max(del.ext)-4),".shp",sep="")

my.shapefile <- read.shp(path1a)

mydata_inner <- convert.to.simple(my.shapefile)

path2 <- dirname(path1)
setwd(path2)

tkmessageBox(message = "Please select the outter baseline...")
path3 <- tk_choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path3)

path3a <- paste(substr(path3,1,max(del.ext)-4),".shp",sep="")

my.shapefile2 <- read.shp(path3a)

mydata_outer <- convert.to.simple(my.shapefile2)

baseline.dbf <- data.frame(read.dbf(paste(substr(path3,1,max(del.ext)-4),".dbf",sep="")))

colnames(baseline.dbf) <- gsub("ID", "Id", colnames(baseline.dbf))

baseline.dbf$dbf.Id <- seq(1,length(baseline.dbf$dbf.Id),by=1)


time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_transects", showWarnings=FALSE)
setwd("AMBUR_transects")

dir.create(paste(time.stamp2," ","transects",sep=""))
setwd(paste(time.stamp2," ","transects",sep="")) 

###########build the status bar
pb <- tkProgressBar("AMBUR: progress bar", "Some information in %", 0, 100, 50)

    info <- sprintf("%d%% done", 0)
    setTkProgressBar(pb, 0, sprintf("AMBUR: Transect casting (%s)", info), info)

#####################



#build function for points along a line (modified to add the start and end points)

pts.along <- function(x,y,pspace) {

#VARIABLES: x = x coord, y = y coord, pspace = spacing between points

Cx <- x
Cy <- y

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Segment.Length <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

Dx <- Cx2 - Cx
Dy <- Cy2 - Cy

Segment.Azimuth <- ifelse(Dx >= 0, 90 -(180/pi) * atan(Dy/Dx),270 -(180/pi) * atan(Dy/Dx))

Cumltiv.Sum <- ave(Segment.Length, FUN=cumsum)

pts.int <- pspace  ##defines the point spacing along the line

pts.bin <- seq(from=pts.int,to=max(Cumltiv.Sum),by=pts.int)

Cumltiv.Sum2 <- c(0,Cumltiv.Sum[-max(length(Cumltiv.Sum))])


pts.bin.up <- length(pts.bin)
pts.bin.upx <- length(pts.bin)
pts.bin.upy <- length(pts.bin)
pts.bin.upx2 <- length(pts.bin)
pts.bin.upy2 <- length(pts.bin)
pts.bin.diff <- length(pts.bin)

for (i in 1:length(pts.bin)) {

pts.bin.up[i] <-  which.max(Cumltiv.Sum2[Cumltiv.Sum2  <=  pts.bin[i]])
pts.bin.upx[i] <- x[pts.bin.up[i]]
pts.bin.upy2[i] <- y[pts.bin.up[i]+1]
pts.bin.upx2[i] <- x[pts.bin.up[i]+1]
pts.bin.upy[i] <- y[pts.bin.up[i]]
pts.bin.diff[i] <- pts.bin[i] - Cumltiv.Sum2[pts.bin.up[i]]



}

Dx2 <- pts.bin.upx2 - pts.bin.upx
Dy2 <- pts.bin.upy2 - pts.bin.upy


pts.bin.az  <- ifelse(Dx2 >= 0, 90 -(180/pi) * atan(Dy2/Dx2),270 -(180/pi) * atan(Dy2/Dx2))

t.azimuth <- c(pts.bin.az[1],pts.bin.az,pts.bin.az[max(length(pts.bin.az))])

t.startx <- c(Cx[1],sin((pts.bin.az * pi/180)) * pts.bin.diff + pts.bin.upx,Cx[max(length(Cx))])
t.starty <- c(Cy[1],cos((pts.bin.az * pi/180)) * pts.bin.diff + pts.bin.upy,Cy[max(length(Cy))])


cbind(t.startx,t.starty,t.azimuth)



}



 
##############################################################################################################
#set up the variables for the analyses


for (m in 1:length(mydata_outer$Id)) {

test1 <- which(mydata_outer$Id[m] == baseline.dbf$dbf.Id)

mydata_outer$Id[m] <- baseline.dbf$dbf.BaseOrder[test1]


}

#######status checkpoint (5%)
    info <- sprintf("%d%% done", 5)
    setTkProgressBar(pb, 5, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################    




#setup table to hold the results
finaltable2 <- matrix(ncol=16)
colnames(finaltable2) <- c("startx","starty","perpx","perpy","perpaz","trimx","trimy","trimdist","nearx","neary","nearaz","neardist","basemaxbnum","baseorder","baseoffshore","basecastdir")

texttable2 <- matrix(ncol=3)
colnames(texttable2) <- c("startx","baseloc1","baseloc2")


nbaselines <- sort(unique(c(mydata_outer$Id)))

transtart <- 0

Cx <- 0
Cy <- 0
Cid <- 0

Cx2 <- 0
Cy2 <- 0
 
 pts.id <- 0 
  int.ptsX <- 0
    int.ptsY <- 0
      trans.id <- 0

#######status checkpoint (10%)
    info <- sprintf("%d%% done", 10)
    setTkProgressBar(pb, 10, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################      


for (b in 1:length(nbaselines)) {


Ax <- mydata_inner$X
Ay <- mydata_inner$Y

Bx <- mydata_outer$X[mydata_outer$Id == nbaselines[b]]
By <- mydata_outer$Y[mydata_outer$Id == nbaselines[b]]

baseloc1x <- as.character(baseline.dbf$dbf.Location[baseline.dbf$dbf.BaseOrder == nbaselines[b]])
basemaxbnumx <-  baseline.dbf$dbf.MaxBNum[baseline.dbf$dbf.BaseOrder == nbaselines[b]]
baseorderx <-  baseline.dbf$dbf.BaseOrder[baseline.dbf$dbf.BaseOrder == nbaselines[b]]
baseoffshorex <- baseline.dbf$dbf.OFFshore[baseline.dbf$dbf.BaseOrder == nbaselines[b]]
basecastdirx <-   baseline.dbf$dbf.CastDir[baseline.dbf$dbf.BaseOrder == nbaselines[b]]
baseloc2x <- as.character(baseline.dbf$dbf.BASE_LOC[baseline.dbf$dbf.BaseOrder == nbaselines[b]])

##############################################################################################################



Cx <- Bx
Cy <- By

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Segment.Length <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

Dx <- Cx2 - Cx
Dy <- Cy2 - Cy

Segment.Azimuth <- ifelse(Dx >= 0, 90 -(180/pi) * atan(Dy/Dx),270 -(180/pi) * atan(Dy/Dx))



####cast transects

outer.basepts <- pts.along(Bx,By,transpace)  #for transects



#cast near transects
Cx <- mydata_inner$X
Cy <- mydata_inner$Y

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Cx3 <- c(Cx[1],Cx[-length(Cx)])
Cy3 <- c(Cy[1],Cy[-length(Cy)])

Cxo <- outer.basepts[,1]
Cyo <- outer.basepts[,2]

Cx2o <- c(Cxo[-1],Cxo[length(Cxo)])
Cy2o <- c(Cyo[-1],Cyo[length(Cyo)])

Cx3o <- c(Cxo[1],Cxo[-length(Cxo)])
Cy3o <- c(Cyo[1],Cyo[-length(Cyo)])


TY.w <- owin()
TY.w <- owin(c(min(Cx-100000),max(Cx+100000)), c(min(Cy-100000),max(Cy+100000)))
TY <- ambur.psp(Cx,Cy,Cx2,Cy2,window=TY.w)

TX.w <- owin()
TX.w <- owin(c(min(outer.basepts[,1]-100000),max(outer.basepts[,1]+100000)), c(min(outer.basepts[,2]-100000),max(outer.basepts[,2]+100000)))
TX <- ambur.ppp(outer.basepts[,1],outer.basepts[,2],window=TX.w)

################################### adjusted project2segment function to handle NA values
adj.project2segment <-       function (X, Y, action = "project", check = FALSE) 
{
    stopifnot(is.ppp(X))
    stopifnot(is.psp(Y))
    stopifnot(action %in% c("distance", "identify", "project"))
    if (Y$n == 0) 
        stop("Segment pattern Y contains 0 segments; projection undefined")
    if (X$n == 0) {
        nowt <- numeric(0)
        none <- integer(0)
        switch(action, identify = return(none), distance = return(list(dist = nowt, 
            which = none)), project = return(list(Xproj = X, 
            mapXY = none, d = nowt, tp = nowt)))
    }
    XX <- as.matrix(as.data.frame(unmark(X)))
    YY <- as.matrix(as.data.frame(unmark(Y)))
    d <- distppllmin(XX, YY)
    mapXY <- d$min.which
    if (action == "identify") 
        return(mapXY)
    else if (action == "distance") 
        return(data.frame(dist = d$min.d, which = mapXY))
    alldata <- as.data.frame(cbind(XX, YY[mapXY, , drop = FALSE]))
    colnames(alldata) <- c("x", "y", "x0", "y0", "x1", "y1")
    dx <- with(alldata, x1 - x0)
    dy <- with(alldata, y1 - y0)
    leng <- sqrt(dx^2 + dy^2)
    
    co <- dx/leng
    co[is.na(co)==TRUE] <- 0.000001  #added 6/27/2011 to remove NA values from crashing the function  
    
    si <- dy/leng
    si[is.na(si)==TRUE] <- 0.000001   #added 6/27/2011 to remove NA values from crashing the function  
    
    xv <- with(alldata, x - x0)
    yv <- with(alldata, y - y0)
    xpr <- xv * co + yv * si
    ypr <- -xv * si + yv * co
    left <- (xpr <= 0)
    right <- (xpr >= leng)
    xr <- with(alldata, ifelse(left, 0, ifelse(right, leng, xpr)))
    xproj <- with(alldata, x0 + xr * co)
    yproj <- with(alldata, y0 + xr * si)
    Xproj <- ambur.ppp(xproj, yproj, window = X$window, marks = X$marks, 
        check = check)
    tp <- xr/leng
    tp[!is.finite(tp)] <- 0
    return(list(Xproj = Xproj, mapXY = mapXY, d = d$min.d, tp = tp))
}


####################################end adj.project2segment function
v <- adj.project2segment(TX,TY)
  Xproj <- v$Xproj

Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #2nd iteration:this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #3nd iteration:this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)

#transect filter using average window of azimuths
tfilter.tab <- cbind(TX$x,TX$y,Xproj$x,Xproj$y,v$d,v$mapXY,v$tp)
colnames(tfilter.tab) <- c("T_x","T_y","In_x","In_y","Dist","Near_Seg","Pos")


Dx2_in <- tfilter.tab[,3] - tfilter.tab[,1]
Dy2_in <- tfilter.tab[,4] - tfilter.tab[,2]



in.az  <- ifelse(Dx2_in >= 0, 90 -(180/pi) * atan(Dy2_in/Dx2_in),270 -(180/pi) * atan(Dy2_in/Dx2_in))

in.length <- ((tfilter.tab[,3]- tfilter.tab[,1])^2 +  (tfilter.tab[,4] - tfilter.tab[,2])^2)^(1/2)


Inx2 <- c(tfilter.tab[-1,3],tfilter.tab[length(tfilter.tab[,3]),3])
Iny2 <- c(tfilter.tab[-1,4],tfilter.tab[length(tfilter.tab[,4]),4])

Inx3 <- c(tfilter.tab[1,3],tfilter.tab[-length(tfilter.tab[,3]),3])
Iny3 <- c(tfilter.tab[1,4],tfilter.tab[-length(tfilter.tab[,4]),4])

end.tspace <-  ((Inx2- tfilter.tab[,3])^2 +  (Iny2 - tfilter.tab[,4])^2)^(1/2)
end.tspace2 <-  ((Inx3- tfilter.tab[,3])^2 +  (Iny3 - tfilter.tab[,4])^2)^(1/2)

#######status checkpoint (25%)
    info <- sprintf("%d%% done", 25)
    setTkProgressBar(pb, 25, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################

#cast perpendicular transects along outer baseline
t.startx <- outer.basepts[,1]
t.starty <- outer.basepts[,2]
t.azimuth <- ifelse(outer.basepts[,3] + 90 >= 360, outer.basepts[,3] + 90 - 360, outer.basepts[,3] + 90) 
t.length <- tranlength * basecastdirx

t.endx <- sin((t.azimuth * pi/180)) * t.length + t.startx
t.endy <- cos((t.azimuth * pi/180)) * t.length + t.starty


#######status checkpoint (50%)
    info <- sprintf("%d%% done", 50)
    setTkProgressBar(pb, 50, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################


####################################test to trim transects

test.wx <- c(t.startx,t.endx,Cx,Cx2,mydata_outer$X)
test.wy <- c(t.starty,t.endy,Cy,Cy2,mydata_outer$Y)

Test.w <- owin()
Test.w <- owin(c(min(test.wx-100000),max(test.wx+100000)), c(min(test.wy-100000),max(test.wy+100000)))

TY.w <- owin()
TY.w <- owin(c(min(Cx-100000),max(Cx+100000)), c(min(Cy-100000),max(Cy+100000)))
TY <- ambur.psp(Cx,Cy,Cx2,Cy2,window=Test.w)

trim.x <- numeric(length(tfilter.tab[,1]))
trim.y  <- numeric(length(tfilter.tab[,1]))
trim.length  <- numeric(length(tfilter.tab[,1]))

for (i in 1:length(tfilter.tab[,1])) {
 
b.x <- ambur.psp(t.startx[i], t.starty[i], t.endx[i], t.endy[i], window=Test.w)

 inner.trim <- crossing.psp(b.x,TY)

pts.dists <- ((inner.trim$x - tfilter.tab[,1][i])^2 +  (inner.trim$y - tfilter.tab[,2][i])^2)^(1/2)

 trim.x[i] <- ifelse(inner.trim$n == 0, t.endx[i],inner.trim$x[which.min(pts.dists)] )
 trim.y[i] <- ifelse(inner.trim$n == 0, t.endy[i],inner.trim$y[which.min(pts.dists)] )
 trim.length[i] <- ifelse(inner.trim$n == 0,t.length,min(pts.dists))
 } 

filter.length6 <- trim.length
filter.x6 <- trim.x
filter.y6 <- trim.y
#####################################################################


#construct master data table

finaltable1 <- cbind(t.startx,t.starty,t.endx,t.endy,t.azimuth,trim.x,trim.y,trim.length,Xproj$x,Xproj$y,in.az,in.length,basemaxbnumx,baseorderx,baseoffshorex,basecastdirx)
colnames(finaltable1) <- c("startx","starty","perpx","perpy","perpaz","trimx","trimy","trimdist","nearx","neary","nearaz","neardist","basemaxbnum","baseorder","baseoffshore","basecastdir")

finaltable2 <- rbind(finaltable2,finaltable1)


texttable1 <- cbind(t.startx,baseloc1x,baseloc2x)

colnames(texttable1) <- c("ttx","baseloc1","baseloc2")

texttable2 <- rbind(texttable2,texttable1)

transtart <- max(length(t.startx))

}
finaltable3 <- finaltable2[-1,]

texttable3 <- texttable2[-1,]

attach(data.frame(finaltable3))

#######status checkpoint (75%)
    info <- sprintf("%d%% done", 75)
    setTkProgressBar(pb, 75, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################

#plots for fun
par(mfrow=(c(3,1)))
par(pty= "m")
plot(c(startx,perpx),c(starty,perpy), type="l",asp=1,col="white",main="Perpendicular Transects",xlab="X",ylab="Y")
lines(mydata_inner$X,mydata_inner$Y, type="l",col="gray")
segments(startx,starty,perpx,perpy,col="blue")

plot(c(startx,perpx),c(starty,perpy), type="l",asp=1,col="white",main="Perpendicular Trimmed Transects",xlab="X",ylab="Y")
lines(mydata_inner$X,mydata_inner$Y, type="l",col="gray")
segments(startx,starty,trimx,trimy,col="blue")

plot(c(startx,perpx),c(starty,perpy), type="l",asp=1,col="white",main="Near Transects",xlab="X",ylab="Y")
lines(mydata_inner$X,mydata_inner$Y, type="l",col="gray")
segments(startx,starty,nearx,neary,col="blue")


#write shapefiles of the transects
library(shapefiles)
id.field <- seq(1,length(startx),by=1)
dd <- data.frame(Id=c(id.field,id.field),X=c(startx,perpx),Y=c(starty,perpy))
ddTable <- data.frame(Id=id.field,Transect=id.field,TranSpace=transpace,TranDist=tranlength,Location=texttable3[,"baseloc1"],MaxBNum=basemaxbnum,BaseOrder=baseorder,OFFshore=baseoffshore,CastDir=basecastdir,BASE_LOC=texttable3[,"baseloc2"],StartX=startx,StartY=starty,EndX=perpx,EndY=perpy,Azimuth=perpaz,Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("b",b,"transects_perp",sep=""), arcgis=T)

id.field <- seq(1,length(startx),by=1)
dd <- data.frame(Id=c(id.field,id.field),X=c(startx,trimx),Y=c(starty,trimy))
ddTable <- data.frame(Id=id.field,Transect=id.field,TranSpace=transpace,TranDist=trimdist,Location=texttable3[,"baseloc1"],MaxBNum=basemaxbnum,BaseOrder=baseorder,OFFshore=baseoffshore,CastDir=basecastdir,BASE_LOC=texttable3[,"baseloc2"],StartX=startx,StartY=starty,EndX=trimx,EndY=trimy,Azimuth=perpaz,Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("b",b,"transects_perp_trim",sep=""), arcgis=T)


id.field <- seq(1,length(startx),by=1)
dd <- data.frame(Id=c(id.field,id.field),X=c(startx,nearx),Y=c(starty,neary))
ddTable <- data.frame(Id=id.field,Transect=id.field,TranSpace=transpace,TranDist=neardist,Location=texttable3[,"baseloc1"],MaxBNum=basemaxbnum,BaseOrder=baseorder,OFFshore=baseoffshore,CastDir=basecastdir,BASE_LOC=texttable3[,"baseloc2"],StartX=startx,StartY=starty,EndX=nearx,EndY=neary,Azimuth=nearaz,Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("b",b,"transects_near_inner",sep=""), arcgis=T)




#detach("package:shapefiles")
detach(data.frame(finaltable3))

#######status checkpoint (100%)
    info <- sprintf("%d%% done", 100)
    setTkProgressBar(pb, 100, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################


}

