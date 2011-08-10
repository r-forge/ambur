require(rgdal) # also loads sp and lattice packages
 require(tcltk)
tkmessageBox(message = "Please select the shapefile...")

getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)

 
   shapedata <- readOGR(getdata,layer=shapename)
   
   attrtable <- data.frame(shapedata)
attrtable
  
  
  
# plot the three geo-objects with a single, simple plot

   plot(shapedata)
  
   title("Map")  
 # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata) # contains projection info  
 
# finally, write a shape file (with .prj component)
  outputname <- paste("ambur_",shapename,sep="")
  
   writeOGR(shapedata, ".", outputname, driver="ESRI Shapefile") 
   message("done with cleaning shapefile") 