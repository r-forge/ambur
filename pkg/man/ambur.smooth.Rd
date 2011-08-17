\name{ambur.smooth}
\alias{ambur.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smooth a single polyline }
\description{
  Smooth a single polyine in a shapefile using a local polynomial fitting algorithm.
}
\usage{
ambur.smooth(alpha = 0.7, degree = 2, sampledist = 5, presame = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{ Alpha value for smoothing }
  \item{degree}{ Degree of polynomial }
  \item{sampledist}{ Sampling distance along the original polyline }
  \item{sampledist}{ Perform sampling along the original polyline. 1 = yes and 0 = no }
}
\details{
  A polyline is smoothed using a local polynomial fit based on a modified 'locfit' algorithm.  See the 'locfit' package for details on the original 'locfit' function.
}
\value{
 
 A shapefile with a smoothed polyline.
 
}
\references{ }
\note{ ~~ This function is a beta version ~~ }
\seealso{ }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
\keyword{ generalize }% __ONLY ONE__ keyword per line