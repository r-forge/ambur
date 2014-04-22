\name{ambur.hva}
\alias{ambur.hva}

\title{
Coastal Hazards Vulnerability Assessment 
}
\description{
  The function requires four ESRI Shapefiles with properly formatted field names and attributes:   \cr
 
 1. post-AMBUR analysis envelope transect shapefile  \cr
 2. FEMA Q3/FIRM shapefile  \cr
 3. NOAA/NWS MOM-AGL or USACOE SLOSH shapefile  \cr
 4. NOAA SoVI shapefile  \cr 
 
  The function will generate shapefiles for five vulnerability analyses:   \cr
  
  shoreline change   \cr
  inundation   \cr
  SoVI scaled   \cr
  inundation/SoVI   \cr
  shoreline change/inundation/SoVI   \cr
   
  The function will generate a text file during the analysis:   \cr
  
  parameters (a text file with analysis paramaters/variables than can be   \cr
  'sourced' in R or cut & pasted in R at the R prompt)   \cr
}
\usage{
ambur.hva(loc = "location", bufferdist = 150, bufferdist2 = 10, removeanthro = "no", 
 file1 = "path", file2 = "path", file3 = "path", file4 = "path", dir1 = "path", 
 scr = c(-1, -0.2, 0.2, 1), tscr = c(0.5, 0.2), sscr = c(0.16, 0.074), 
 sovi = c(6.77, 4.06, -1.37, -4.09))
}

\arguments{
  \item{loc}{
Defines the geographic location of the study area that will be used to name the output folder tag. (optional input)
}
  \item{bufferdist}{
Buffer radius in map units that is used to search adjacent transects to collect spatial shoreline change variability.
}
  \item{bufferdist2}{
Buffer distance to convert transect lines into polygons for intersecting with inundation and SoVi data.
}
  \item{removeanthro}{
Remove transects that are classified as "anthropogenic" from the analyses.
}
  \item{file1}{
File path (in quotes) to the post-AMBUR analysis envelope transect shapefile.  If left as "path" then a file will be chosen interactively.
}
  \item{file2}{
File path (in quotes) to the FEMA Q3/FIRM shapefile.  If left as "path" then a file will be chosen interactively.
}
  \item{file3}{
File path (in quotes) to the NOAA/NWS MOM-AGL or USACOE SLOSH shapefile.  If left as "path" then a file will be chosen interactively.
}
  \item{file4}{
File path (in quotes) to the NOAA SoVI shapefile.  If left as "path" then a file will be chosen interactively.
}
  \item{dir1}{
Directory/folder path (in quotes) to store the results in.  If left as "path" then a directory/folder will be chosen interactively.
}
  \item{scr}{
Shoreline change rate breakpoints in map units/year (high to low risk)
}
  \item{tscr}{
Temporal shoreline change rate breakpoints (high to low risk)
}
  \item{sscr}{
Spatial shoreline change rate breakpoints (high to low risk)
}
  \item{sovi}{
Sovi score breakpoints (high to low risk)
}
}
\details{
none
}
\value{
 date and time stamped folder containing the following shapefiles:

 ambur_hva_shoreline_change 
 ambur_hva_inundation
 ambur_hva_inundation_sovi
 ambur_hva_sc_inundation_sovi


}
\references{
http://ambur.r-forge.r-project.org/ \cr

 https://msc.fema.gov/ \cr

 http://www.fema.gov/floodplain-management/flood-zones \cr

 http://slosh.nws.noaa.gov/ \cr
 (Use a user name of Gustav2000 and a password of Ike2008) \cr

 http://webra.cas.sc.edu/hvri/products/sovi.aspx  \cr


}
\author{
Chester W. Jackson Jr., Ph.D.
 Georgia Southern University - Applied Coastal Research Lab
}
\note{
 Some data are available at NOAA's Digital Coastal website.
}


\examples{

A <- system.file("external/Jekyll_envelope_transects_original.shp", package="amburhva")
B <- system.file("external/Jekyll_Q3_DFIRM_GEMA.shp", package="amburhva")
C <- system.file("external/Jekyll_SLOSH_USACOE_GEMA.shp", package="amburhva")
D <- system.file("external/Jekyll_SoVI.shp", package="amburhva")
E <- dirname(system.file("external/Jekyll_envelope_transects_original.shp", package="amburhva"))

ambur.hva(loc = "location", bufferdist = 150, bufferdist2 = 10, removeanthro = "no", 
 file1 = A, file2 = B, file3 = C, file4 = D, dir1 = E, scr = c(-1, -0.2, 0.2, 1), 
 tscr = c(0.5, 0.2), sscr = c(0.16, 0.074), sovi = c(6.77, 4.06, -1.37, -4.09))

}

\keyword{ vulnerability }
\keyword{ ambur }
