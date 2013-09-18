\name{ambur.critical}
\alias{ambur.critical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate critical points for a polyline shapefile }
\description{
 Determines the critical points of a single polyline using a linear spline function.
}
\usage{
ambur.critical(ncritpts = 50, sampledist = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ncritpts}{ Number of critical points for output}
  \item{sampledist}{Sampling distance along the polyline for spline function}
}
\details{
A linear spline function is used to approximate the shape of a single polyline in shapefile.  The linear spline is constructed using the number of critical points assigned to 'ncritpts'.
}
\value{
A shapefile containing critical points.  The shapefile is named "criticalpts.shp" and is stored in a time-stamped folder in side a folder named "AMBUR_critical" where the original polyine shapefile is located.
}

\author{ Chester W. Jackson Jr. }


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ critical }
\keyword{ point }% __ONLY ONE__ keyword per line