\name{ambur.check}
\alias{ambur.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Check shapefiles for missing attribute fields }
\description{
Determines if necessary fields are present in boundary, transect, or capture points shapefiles.  
}
\usage{
ambur.check(checktype = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{checktype}{ Check for missing fields.  Enter: '1' for "boundaries; '2' for transects; '3' for capture points }
}
\details{
  All necessary fields must be present in order to run various AMBUR functions.  The following fields are required for:

boundaries:  
"Id","DATE_","ACCURACY","SHORE_LOC","CLASS_1","CLASS_2","CLASS_3","GROUP"

transects:  
"Id","Transect","TranSpace","TranDist","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC","StartX","StartY","EndX","EndY","Azimuth"
 
capture points:
"Id","DATE_","ACCURACY","SHORE_LOC","CLASS_1","CLASS_2","CLASS_3","GROUP","Transect","TranSpace","TranDist","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC","StartX","StartY","EndX","EndY","Azimuth","Point_X","Point_Y","Distance"


}
\value{
Returns a list of missing field names.
}
\references{  }
\author{ Chester W. Jackson Jr. }
\note{ }
\seealso{  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ field }
\keyword{ missing }% __ONLY ONE__ keyword per line