ambur.pts2line <-
function(npts="all") {

require(tcltk)
require(rgdal)

tkmessageBox(message = "Please select the *.csv file containing point data...")
getdata <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata, header=TRUE, sep=",")
attach(mydata)

dir.path <- dirname(getdata)
setwd(dir.path)

file.name <- substr(basename(getdata),1,nchar(basename(getdata))-4) 

  edit(mydata)
  
 



target.field <- select.list(colnames(mydata), multiple = FALSE, title = "Choose field with X (easting) data:")  

field.x <- mydata[,target.field]

target.field2 <- select.list(colnames(mydata), multiple = FALSE, title = "Choose field with Y (northing) data:")  

field.y <- mydata[,target.field2]

target.field3 <- select.list(colnames(mydata), multiple = FALSE, title = "Choose field with numeric ID to assign points to unique line(s):")  

mydata$Baseline.Location <- ifelse(is.na(mydata[,target.field3]) == T, 0, mydata[,target.field3])



Baseline.Factor <- factor(mydata$Baseline.Location)
LINE.final <- sapply(levels(Baseline.Factor), function(x)
list(Lines(list(Line(list(x=c(field.x[mydata$Baseline.Location == x]), y=c(field.y[mydata$Baseline.Location == x])))), ID=(as.numeric(x))))
,simplify = TRUE)
LINE.final2 <- SpatialLines(LINE.final)


LINE.tab <- data.frame(Meth="Points to Line",Source=dir.path,Creator="R - AMBUR")

LINE.tab2 <-  LINE.tab[rep(1, length(unique(Baseline.Factor))),]
row.names(LINE.tab2) <- getSLLinesIDSlots(LINE.final2)   ##fixes bug and gets IDs of spatial lines without having to manually coding them (4/2012)

LINE.final3 <- SpatialLinesDataFrame(LINE.final2, LINE.tab2)
#create shapefile and write it to the working directory
output.fname <- paste(file.name,"line", sep = "_", collapse = NULL)

writeOGR(LINE.final3, ".", output.fname, driver="ESRI Shapefile")

}