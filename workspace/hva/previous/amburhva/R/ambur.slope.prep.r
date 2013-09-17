ambur.slope.prep <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1

require(raster)
require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select the percent slope raster image...")
filters <- matrix(c("Imagine", ".img"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files(filter = filters,multi = FALSE)
shapename <- gsub(".img", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)

shapedata <- raster(getdata)

#spplot(shapedata)

#shapedata <- readGDAL(getdata)
#attrtable <- data.frame(shapedata)


tkmessageBox(message = "Please select a post-AMBUR analysis transect shapefile generated from ambur.statshape function...")
filters2 <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata2 <- tk_choose.files(filter = filters2,multi = FALSE)
shapename2 <- gsub(".shp", "", basename(getdata2))

shapedata2 <- readOGR(getdata2,layer=shapename2)
attrtable2 <- data.frame(shapedata2)

 ########
time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_HVA_data_prep", showWarnings=FALSE)
setwd("AMBUR_HVA_data_prep")

#shapedata2A <- as(shapedata2, "SpatialLines")

 raster.tran <-  extract(shapedata, shapedata2,fun=mean, na.rm=TRUE, sp=TRUE, progress='text')
  #extract(x, y, fun=NULL, na.rm=FALSE, cellnumbers=FALSE, df=FALSE, layer, nl, factors=FALSE, along=FALSE, sp=FALSE, ...)
######setup HVA calcs

shapetable <-  data.frame(raster.tran)

colnames(shapetable)[length(colnames(shapetable))] <- "mn_slp_pct"

hva_cat <- numeric(length(shapetable[,"mn_slp_pct"] ))

hva_cat[shapetable[,"mn_slp_pct"] > 70] <- 1
hva_cat[shapetable[,"mn_slp_pct"] > 40 & shapetable[,"mn_slp_pct"] < 70] <- 2
hva_cat[shapetable[,"mn_slp_pct"] > 20 & shapetable[,"mn_slp_pct"] < 40] <- 3
hva_cat[shapetable[,"mn_slp_pct"] > 10 & shapetable[,"mn_slp_pct"] < 20] <- 4
hva_cat[shapetable[,"mn_slp_pct"] < 10] <- 5

shapetable$hva_cat <- hva_cat          

slope_hva_table <- cbind(shapetable$Transect, hva_cat)

colnames(slope_hva_table) <- c("Transect","hva_cat")

write.table(slope_hva_table, file = "ambur_hva_slope_data.csv", sep = ",", row.names = FALSE)

###########################################

lns.output <- SpatialLinesDataFrame(raster.tran,shapetable)


 writeOGR(lns.output, ".", "ambur_hva_slope_lns", driver="ESRI Shapefile")
 


 




}