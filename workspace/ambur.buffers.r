ambur.buffer <-
function(buffdist=5,buffnum=1) {

require(tcltk)
require(rgdal)
require(rgeos)

#buffdist <- 5 #for testing
#buffnum <- 10 #for testing

totbuffers <- seq(buffdist,buffdist*buffnum,buffdist)

tkmessageBox(message = "Please select the polyline shapefile...")
getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
shapedata <- as(shapedata, "SpatialLinesDataFrame")
attrtable <- data.frame(shapedata)

workingdir <- dirname(getdata)
setwd(workingdir)
time.stamp1 <- as.character(Sys.time())
time.stamp2 <- gsub("[:]", "_", time.stamp1)


dir.create("AMBUR_buffer", showWarnings=FALSE)
setwd("AMBUR_buffer")

dir.create(paste(time.stamp2," ","buffers",sep=""))
setwd(paste(time.stamp2," ","buffers",sep=""))



 ###### break down outer baseline into simply points with IDs
crdl0 <- coordinates(shapedata)
crd.1 <- sapply(crdl0, function(x) do.call("rbind", x),simplify = FALSE)
crd.2    <- do.call("rbind", crd.1)
crd.3 <- as.numeric(sapply(crd.1, function(x) max(row(x)) ,simplify = TRUE))
crd.len.test <- as.numeric(length(crd.3))
if(crd.len.test <= 1) crd.rep <-  1 else crd.rep <- seq(1, length(crd.3),1)
basepointIDs <- rep(crd.rep,crd.3)
baseshapeIDs <- basepointIDs - 1
sortshapeIDs <- seq(1,length(basepointIDs),1)
basex <- crd.2[,1]
basey <- crd.2[,2]

outerbase.tab <- data.frame(sortshapeIDs,baseshapeIDs,basepointIDs,basex,basey)
colnames(outerbase.tab) <- c("sortshapeID","shapeID","baseID","baseX", "baseY")


 Baseline.Factor <- unique(outerbase.tab$baseID)
 Buffer.Factor <- totbuffers


width.factor <- rep(Buffer.Factor, nrow(shapedata))
line.id <- sort(rep(1:nrow(shapedata),buffnum))

l <- vector("list", nrow(shapedata)*buffnum)

for (i in seq_len(nrow(shapedata)*buffnum)) {

        l[[i]] <- gBuffer(shapedata[line.id[i],], width = width.factor[i], byid=TRUE, capStyle="FLAT")
        l[[i]] <- spChFIDs(l[[i]], as.character(i))

}

vec2 <- unlist(l)


out <- do.call("rbind", vec2)
rn <- row.names(out)
nrn <- do.call("rbind", strsplit(rn, " "))

final_buffers <- as(out, "SpatialLines")

plot(final_buffers)


buff.tab <- data.frame(line.id,width.factor)
colnames(buff.tab) <- c("baseID","distance")

out.buffers <- SpatialLinesDataFrame(final_buffers, buff.tab)


writeOGR(out.buffers, ".", "polyline_buffers", driver="ESRI Shapefile")

}