\name{ambur.forecast}
\alias{ambur.forecast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Forecast shoreline position }
\description{
 Estimates the location of future or past shoreline postions.
}
\usage{
ambur.forecast(years = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{years}{ Number of years into the future }
}
\details{
  Provides a forecast of shoreline position for n years in to the future for the EPR, LRR, and WLR rate models.  The method of forecasting is a simple linear extrapolation and should be used with caution.  
}
\value{
Shapefiles for each rate model: 
EPR_forecast.shp 
LRR_forecast.shp
WLR_forecast.shp)
}

\author{ Chester W. Jackson Jr. }
\note{ ~~ This function is a beta version ~~ }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ setback }
\keyword{ extrapolation }% __ONLY ONE__ keyword per line