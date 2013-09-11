\name{ambur.tranrepair}
\alias{ambur.tranrepair}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Repair transects }
\description{
  Construct transects using the 'near' method between to polylines
}
\usage{
ambur.tranrepair(test = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test}{ test }

}
\details{
  Repairs transects orientations (azimuths), distances, and starting/ending coordinates for transect positions manually adjusted by the user within GIS software. 'test' argument is unimplemented.
}
\value{
Shapefiles:
shapefile with shapefile name appended with '_repaired.shp'
}
\references{ }
\note{ }
\seealso{ }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ transect }
\keyword{ repair }% __ONLY ONE__ keyword per line