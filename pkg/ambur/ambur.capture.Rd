\name{ambur.capture}
\alias{ambur.capture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Capture intersection points between transects and historical boundaries }
\description{
Generates a point shapefile that contains the locations of where a shapefile with transects intersect another shapefile with boundaries (historical shorelines)  }
\usage{
ambur.capture(userinput1 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{userinput1}{ Not currently used by the function. Default = 1 }
}
\details{
A shapefile with transects is intersected with another shapefile containing historical boundaries (shorelines) to create a point shapefile.  The transect and boundary shapefiles must have the proper fields in their attribute tables.  See 'ambur.check' and 'ambur.addfields' command for determining if proper fields are present in shapefiles. The resulting point shapefile will be used by 'ambur.analysis' to generate statistics/graphics of changes and rates.
}
\value{
A shapefile containing intersection points.  The shapefile is named "shorepts.shp" and is stored in a time-stamped folder in side a folder named "AMBUR_capture" where the boundary shapefile is located.

}
\references{ }
\author{ Chester W. Jackson Jr. }
\note{ }
\seealso{  }
\examples{
Enter the following at the command prompt to start the function:

ambur.capture()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ intersection }
\keyword{ capture }% __ONLY ONE__ keyword per line
