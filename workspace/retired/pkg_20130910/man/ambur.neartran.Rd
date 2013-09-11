\name{ambur.neartran}
\alias{ambur.neartran}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Construct transects }
\description{
  Construct transects using the 'near' method between to polylines
}
\usage{
ambur.neartran(sampledist = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampledist}

}
\details{
  Construct transects using the the near method.  User must supply two different polyline shapefiles that represent an inner and outer baseline to cast transects in between.  Baseline shapefiles must be properly prepared with the correct fields in the attribute table.  See ambur.addfields function.
}
\value{
Shapefiles:
near_transects.shp
}
\references{ }
\note{ }
\seealso{ }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ transect }
\keyword{ near }% __ONLY ONE__ keyword per line