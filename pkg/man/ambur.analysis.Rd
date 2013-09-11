\name{ambur.analysis}
\alias{ambur.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Analyze changes and rates of boundary movement }
\description{
  The function requires an ESRI Shapefile with point data representing historical boundary positions along a transect.  The shapefile must be properly formatted and contain the proper fields. See links below for additional information.  The function will generate data tables and graphics containing the results of analysis statics given user specified parameters.
}
\usage{
ambur.analysis(userinput1 = "first", userinput2 = 95, 
userinput3 = "m", userinput4 = "", userinput5 = 1, 
userinput6 = "all", userinput7 = "basic",  userinput8 = "yr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{userinput1}{ Use first or last intersections along the transect.  DEFAULT = "first" }
  \item{userinput2}{ Confidence interval. DEFAULT = 95  }
  \item{userinput3}{ Map units. DEFAULT = "m" for meters (mainly used for labeling graphic axes) }
  \item{userinput4}{ Name of the analysis.  DEFAULT = nothing. (will add a name to the time-stamped folder where results are stored }
  \item{userinput5}{ Start transect. DEFAULT = 1 }
  \item{userinput6}{ End transect. DEFAULT = all }
  \item{userinput7}{ Basic or advanced analysis.  DEFAULT = basic }
  \item{userinput8}{ Time units for rates.  DEFAULT = "yr" (options = "day", "hr", "min", "sec")}
}
\details{
  If an "advanced" analysis is specified, more statics are generated and included in the output data tables.  However, this will increase the processing time.
}
\value{
  Creates a date and time stamped folder inside a parent folder named "AMBUR_results", located in the same folder as the shapefile.  The folder contains all of the output data and graphics.
}

\author{ Chester W. Jackson Jr. }


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ambur }
\keyword{ shoreline }% __ONLY ONE__ keyword per line
