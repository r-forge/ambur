\name{ambur.filtertran}
\alias{ambur.filtertran}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Advanced baseline and transect generation (Experimental) }
\description{
 Generate smoothed baselines and transects using experimental techniques. 
}
\usage{
ambur.filtertran(winsize=5, indv=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{winsize}{ Transect filter window size: Value must be an odd number }
  \item{indv}{ Filter individual baseline transects: Value of 1 means 'yes' and 0 is 'no'}

}
\details{
A new transect azimuth orienation is assigned to the middle transect of and odd numbered moving (sliding) window of transects.  The new transect azimuth takes into account the azimuths of adjacent transects within the window.
}
\value{
  A shapefile with adjusted transect azimuths.
}
\references{  }
\author{ Chester W. Jackson Jr. }
\note{ ~~ This function is a beta version ~~ }
\seealso{  }
\examples{
Filter transects created by the near method using a 5-transect window and do not filter transects across baselines:

ambur.filtertran(winsize=5, indv=1)

Filter transects created by the near method using a 3-transect window and filter transects across baselines:

ambur.filtertran(winsize=3, indv=0)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ filter }
\keyword{ transect }% __ONLY ONE__ keyword per line