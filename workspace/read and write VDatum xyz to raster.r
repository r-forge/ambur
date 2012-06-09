require(raster)
require(rgdal) # also loads sp and lattice packages
require(tcltk)
require(sqldf)


tkmessageBox(message = "Please select the '*.txt' file...")


data.path <- tk_choose.files(default = "*.txt",multi = FALSE)



f <- file(data.path)
bigdf <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = T, row.names = F, sep=" "))

blah <- rasterFromXYZ(bigdf[,1:3])

spplot(blah)

writeRaster(blah, filename="test.img", format="HFA", overwrite=TRUE)