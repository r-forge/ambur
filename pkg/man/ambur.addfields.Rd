\name{ambur.addfields}
\alias{ambur.addfields}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Check shapefiles for missing attribute fields and add missing fields}
\description{
Determines if necessary fields are present in boundary, transect, or capture points shapefiles.  
}
\usage{
ambur.addfields(checktype = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{checktype}{ Check for missing fields.  Enter '1' for boundaries (shorelines) or '2' for baselines }
}
\details{
  Adds all necessary fields must be present in order to run various AMBUR functions.  The following fields are required for:

boundaries (shorelines):  "Id","DATE_","ACCURACY","SHORE_LOC","CLASS_1","CLASS_2","CLASS_3","GROUP"

baselines:  "Id","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC"
 


}
\value{
Adds missing fields to the shapefile's data table.
}

\author{ Chester W. Jackson Jr. }
\note{ The shapefile must contain polyline data.}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ field }
\keyword{ missing }% __ONLY ONE__ keyword per line