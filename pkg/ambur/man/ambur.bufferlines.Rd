\name{ambur.bufferlines}
\alias{ambur.bufferlines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create buffer lines }
\description{
  Create buffer lines from a polyline shapefile
}
\usage{
ambur.bufferlines(buffdist=10,buffnum=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{buffdist}{distance to buffer polylines}
  \item{buffnum}{not implemented }

}
\details{
  Creates a new polyline shapefile that contains buffered lines given a user defined distance.  The process involves buffering lines initially resulting in multiple (potentially overlapping) polygons.  The polygons are dissolve where they overlap and are converted into polylines. In case the buffers will be used as baselines (to make transects from), the attribute table is updated with the necessary fields for baselines.  The values in these fields will need to be edited based on the users needs.
}
\value{
Shapefiles:
filename_buffers_dist.shp
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ buffer }
\keyword{ polylines }% __ONLY ONE__ keyword per line