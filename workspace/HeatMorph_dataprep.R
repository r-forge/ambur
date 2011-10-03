require(tcltk)
require(rgdal)


tkmessageBox(message = "Please select the outer baseline shapefile...")
getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

workingdir <- dirname(getdata)
setwd(workingdir)


tkmessageBox(message = "Please select the inner baseline shapefile...")
getdata2 <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename2 <- gsub(".shp", "", basename(getdata2))
shapedata2 <- readOGR(getdata2,layer=shapename2)
attrtable2 <- data.frame(shapedata2)


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


###### break down inner baseline into simply points with IDs
crdl0 <- coordinates(shapedata2)
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

innerbase.tab <- data.frame(sortshapeIDs,baseshapeIDs,basepointIDs,basex,basey)
colnames(innerbase.tab) <- c("sortshapeID","shapeID","baseID","baseX", "baseY")



inner_numpts <- length(innerbase.tab$baseX)
outer_numpts <- length(outerbase.tab$baseX)

inner_xy <- rep(0,inner_numpts *2)
inner_xy[1:length(inner_xy) %% 2 ==1] <- innerbase.tab$baseX
inner_xy[1:length(inner_xy) %% 2 ==0] <- innerbase.tab$baseY

outer_xy <- rep(0,outer_numpts *2)
outer_xy[1:length(outer_xy) %% 2 ==1] <- outerbase.tab$baseX
outer_xy[1:length(outer_xy) %% 2 ==0] <- outerbase.tab$baseY

#####give them the same start and end points by adding to the outer line
outer_xy_adj <- c(inner_xy[1:2],outer_xy,inner_xy[(length(inner_xy)-1):length(inner_xy)])
outer_numpts_adj <- outer_numpts + 2

output_data <- c(inner_numpts,inner_xy,outer_numpts_adj,outer_xy_adj)

write.table(output_data, "curves.txt", sep="\n", row.names = FALSE, col.names = FALSE)
