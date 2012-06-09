require(raster)
require(rgdal) # also loads sp and lattice packages
 require(tcltk)
tkmessageBox(message = "Please select the raster file...")

getdata <- tk_choose.files(default = " ",multi = TRUE)
shapename <- gsub(".", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)


# read the individual rasters into a list of RasterLayer objects
#input.rasters <- lapply(list.files(pattern="*[.]asc$"), raster)

input.rasters <- lapply(getdata, raster)



# create an empty output raster that spans the full extent of all input
# rasters, and uses the same coordinate reference system; in this case
# we know that all input rasters share the same CRS, so we can
# arbitrarily extract CRS information from the first one
#full.extent <- unionExtent(input.rasters)
#bounding.raster <- raster(full.extent,
 #   crs=projection(input.rasters[[1]]))

# set the output resolution to match the center tile (somewhat
# arbitrarily); this can also be specified manually if preferred
#res(bounding.raster) <- res(input.rasters[[1]])

# for each input raster, extract the corresponding sub-extent of the
# empty output raster, and use this as the basis for resampling
# !! note that this may take several minutes to run !!
#resampled.rasters <- lapply(input.rasters, function(input.raster) {
 #   target.raster <- crop(bounding.raster, input.raster)
  #  resample(input.raster, target.raster, method="bilinear")
#})

# mosaic all 9 resampled rasters, taking the maximum pixel value in
# cases where there is overlap
#raster.mosaic <- mosaic(resampled.rasters, fun=max)

raster.mosaic <- mosaic(input.rasters, fun=min)

#  generate output map

plot(raster.mosaic, col=grey((0:256)/256))


# write the result out to a GeoTIFF, rounding pixel values to the
# nearest integer to preserve the original 1-byte integer format
writeRaster(raster.mosaic, filename="raster_mosaic.img", datatype="FLT4S")





 
 #spplot(raster.mosaic)