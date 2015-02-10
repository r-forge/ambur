ambur.bufferlines <-
function(buffdist=10,buffnum=1) {

#require(tcltk)
#require(rgdal)
#require(rgeos)
#buffdist <- 10 #for testing
#buffnum <- 1

totbuffers <- seq(buffdist,buffdist*buffnum,buffdist)

tkmessageBox(message = "Please select the shoreline shapefile...")
filetype <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files("","Choose file",multi = FALSE,filetype,1)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
shapedata <- as(shapedata, "SpatialLinesDataFrame")
attrtable <- data.frame(shapedata)

workingdir <- dirname(getdata)
setwd(workingdir)
time.stamp1 <- as.character(Sys.time())
time.stamp2 <- gsub("[:]", "_", time.stamp1)


#dir.create("AMBUR_buffer", showWarnings=FALSE)
#setwd("AMBUR_buffer")

#dir.create(paste(time.stamp2," ","buffers",sep=""))
#setwd(paste(time.stamp2," ","buffers",sep=""))



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

        l[[i]] <- gBuffer(shapedata[line.id[i],], width = width.factor[i], byid=TRUE, capStyle="SQUARE")
        l[[i]] <- spChFIDs(l[[i]], as.character(i))

}

vec2 <- unlist(l)


out <- do.call("rbind", vec2)
rn <- row.names(out)
nrn <- do.call("rbind", strsplit(rn, " "))


dissolve <- gUnaryUnion(out)  #added to dissolve polygons

dissolve_buffers <- as(dissolve, "SpatialLines")

#  Extract the individual 'lines'
final_buffers <- disaggregate(dissolve_buffers)

plot(final_buffers)




line.ids <- length(final_buffers)
width.factors <- rep(Buffer.Factor, length(final_buffers))


buff.tab <- data.frame(line.ids,width.factors)
colnames(buff.tab) <- c("Id","buffdist")

######add field names that are required for baselines
reqfields <- c("Id","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC")
#check for missing fields
fieldcheck <- toupper(colnames(buff.tab))

 presentfields <-  which(toupper(reqfields) %in% fieldcheck)

   missingfields <- as.character(reqfields[-presentfields])

add.fields <- data.frame(matrix(data= NA,ncol=length(missingfields),nrow=length(buff.tab[,1])) )

colnames(add.fields) <- missingfields

## adjust the character field widths
fwidth <- rep(" ",30)
fwidth1 <- paste(fwidth, collapse = "")

if ("Location" %in% missingfields) add.fields$Location <- fwidth1
if ("MaxBNum" %in% missingfields) add.fields$MaxBNum <- 0
if ("SHORE_LOC" %in% missingfields) add.fields$SHORE_LOC <- fwidth1
if ("BaseOrder" %in% missingfields) add.fields$BaseOrder <- 1
if ("OFFshore" %in% missingfields) add.fields$OFFshore <- 1
if ("CastDir" %in% missingfields) add.fields$CastDir <- -1
if ("BASE_LOC" %in% missingfields) add.fields$BASE_LOC <- fwidth1
buff.tab$Id <- seq(1,length(buff.tab[,1]),1)


   new.table <-  data.frame(buff.tab,add.fields) 




out.buffers <- SpatialLinesDataFrame(final_buffers, new.table)


  locname <- shapename
      outputname <- paste(locname,"_buffers_",buffdist,sep="")

writeOGR(out.buffers , ".", outputname, driver="ESRI Shapefile")




#get.exp <- as.character(round(length(buff.tab$distance)/c(length(buff.tab$distance),8,4,2,1),0))  # try to get exponential values
#buff.tab2 <- buff.tab[rownames(buff.tab) %in% get.exp,]


#out.buffers2 <- SpatialLinesDataFrame(final_buffers, buff.tab)

#out.buffers3 <- out.buffers2[as.numeric(get.exp),]

#writeOGR(out.buffers3, ".", "polyline_buffers2", driver="ESRI Shapefile")

}

