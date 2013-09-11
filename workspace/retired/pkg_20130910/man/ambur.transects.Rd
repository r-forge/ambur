\name{ambur.transects}
\alias{ambur.transects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Construct transects }
\description{
  Construct transects using multiple methods
}
\usage{
ambur.transects(userinput1 = 50, userinput2 = 500, userinput3 = 5, userinput4 = 5, userinput5 = 90)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{userinput1}{ Transect spacing }
  \item{userinput2}{ Transect length }
  \item{userinput3}{ Not implemented. Inner baseline sampling distance for new vertices }
  \item{userinput4}{ Not implemented. Outer baseline sampling distance for new vertices }
  \item{userinput5}{ Beta test. Angle(s) to cast transects }
}
\details{
  Construct transects using the perpendicular, perpendicular-trimmed, and near methods.
}
\value{
Shapefiles:
transect_perp.shp
transect_perp_trim.shp
transect_near_inner.shp
}
\references{ }
\note{ }
\seealso{ }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ transect }
\keyword{ near }% __ONLY ONE__ keyword per line