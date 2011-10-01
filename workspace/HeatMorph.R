require(tcltk)
require(rgdal)

tkmessageBox(message = "Please select the '*.txt' file...")

data.path <- tk_choose.files(default = "*.txt",multi = FALSE)

dir.path <- dirname(data.path)
setwd(dir.path)

mydata <- read.csv(data.path, header=FALSE, sep=",")


id <- numeric(length(mydata[,1]))
numptsremove <- numeric(length(mydata[,1]))

for (i in 1:length(mydata[,1])) {

id[1] <- 1

testspot <- ifelse(id[i] > 0 && id[i+1] == 0, 1,0)

numptsremove[i] <- ifelse(testspot == 1, i, 0)


idnum <- ifelse(testspot == 1, numptsremove[i], id*1)

range1 <- ifelse(testspot == 1, (mydata[,1][i]) * 2, 1) 

ifelse(testspot == 1, id[(i+1):(i+1+range1)] <- idnum, 0)

}

id <- id[-1]

ptsremove <- numptsremove[numptsremove > 0]

mydata2 <- mydata[,1][-ptsremove]
id1 <- id[-ptsremove]

###replace sequence values 1111888 to 111122222
recode <- function(var, from, to){
x <- as.vector(var)
for (i in 1:length(from)){
  x <- replace(x, x == from[i], to[i])}
if(is.factor(var)) factor(x) else x} 

id2 <- recode(id1, unique(id1), 1:length(unique(id1)))




id3 <- id2[1:length(id2) %% 2 ==1]
xx <- mydata2[1:length(mydata2) %% 2 ==1]
yy <- mydata2[1:length(mydata2) %% 2 ==0]
transectID <- rep(1:length(xx[id3==1]),max(id3))



plot(xx,yy,type="l",col=unique(id3))
points(xx,yy,col=unique(id3))


LineID.Factor <- factor(id3)
smooth.final <- sapply(levels(LineID.Factor), function(x)
list(Lines(list(Line(list(x=c(xx[id3 == x]), y=c(yy[id3 == x])))), ID=(as.numeric(x))))
,simplify = TRUE)
smooth.final2 <- SpatialLines(smooth.final)

smooth1.tab <- data.frame(ID1=1,test1=1)

smooth1.tab2 <-  smooth1.tab[rep(1, length(unique(LineID.Factor))),]
smooth1.tab2$baseID <- seq(1, length(unique(LineID.Factor)),1)
row.names(smooth1.tab2) <- seq(1, length(unique(LineID.Factor)),1)

smooth.dataframe <- smooth1.tab2

smooth.final3 <- SpatialLinesDataFrame(smooth.final2, smooth.dataframe)


#create shapefile and write it to the working directory
writeOGR(smooth.final3, ".", "ambur_HeatMorph", driver="ESRI Shapefile")

####################################################################################################################
LineID.Factor <- factor(transectID)
smooth.final <- sapply(levels(LineID.Factor), function(x)
list(Lines(list(Line(list(x=c(xx[transectID == x]), y=c(yy[transectID == x])))), ID=(as.numeric(x))))
,simplify = TRUE)
smooth.final2 <- SpatialLines(smooth.final)

smooth1.tab <- data.frame(ID1=1,test1=1)

smooth1.tab2 <-  smooth1.tab[rep(1, length(unique(LineID.Factor))),]
smooth1.tab2$Transect <- seq(1, length(unique(LineID.Factor)),1)
row.names(smooth1.tab2) <- seq(1, length(unique(LineID.Factor)),1)

smooth.dataframe <- smooth1.tab2

smooth.final3 <- SpatialLinesDataFrame(smooth.final2, smooth.dataframe)


#create shapefile and write it to the working directory
writeOGR(smooth.final3, ".", "ambur_HeatMorph_transects", driver="ESRI Shapefile")
