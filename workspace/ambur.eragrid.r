
ambur.eragrid <-
function(tstart=1, tend="all") {

#tstart = 25
#tend = 55

require(lattice)
require(grid)
require(tcltk)
tkmessageBox(message = "Please select the raw_positions.csv file...")
getdata <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata, header=TRUE, sep=",")
attach(mydata)


tkmessageBox(message = "Please select the raw_EPR_rates_eras.csv file...")
getdata2 <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata2 <- read.table(getdata2, header=TRUE, sep=",")
attach(mydata2)


dir.path <- dirname(getdata)
setwd(dir.path)




#cull the dataset based on the range of transects

#Extract (cull) the transects for analysis based on user's input range                                                               
Start.Transect <- tstart
End.Transect <- tend

End.Transect2 <- ifelse(End.Transect == "all", max(mydata[,1]), End.Transect)

mydata <- mydata[mydata[,1] >= Start.Transect & mydata[,1] <= End.Transect2,]
mydata2 <- data.frame(mydata2[mydata2[,1] >= Start.Transect & mydata2[,1] <= End.Transect2,],row.names=NULL)


chng.eras <- as.matrix(mydata2[,-1])
   
#get the date for the eras and labels
dates1 <- substr(colnames(mydata)[-1],2,11)
dates2 <- c(dates1[-1],"9999")
dates3 <- paste(dates1,dates2, sep = "-", collapse = NULL)[-(length(dates1))]

Window.Labels  <- gsub("[.]", "/", dates3)
Transect.Labels  <- pretty(as.numeric(mydata2[,1]))

#generate the grid plots
pdf("PDF/graph_eras_grid1.pdf", bg="white")


print(levelplot(chng.eras, pretty=TRUE,na.rm=TRUE,colorkey = list(space="right",tick.number=1,height=0.4,width=1,cex=1, at=c(-0.001,-0.0005,0,0.0005,0.001),labels=c("","Erosion","","Accretion","")),at= c(min(chng.eras,na.rm=TRUE)
,0,max(chng.eras,na.rm=TRUE)),col.regions=c("red","red","blue","blue")
,xlab=list(label="Transect", cex=0.8),ylab=list(label="Era", cex=0.8),scale = list(x = list(labels=Transect.Labels),
y = list(labels=Window.Labels,cex=0.5)),page = function(n) grid.text(paste(""), x = 0.5275, y = 0.00001, default.units = "npc", just = c("center", "top")),aspect=1))


  dev.off()
  
  
  pdf("PDF/graph_eras_grid2.pdf", bg="white")


print(levelplot(chng.eras, pretty=TRUE,na.rm=TRUE,colorkey = list(space="right",tick.number=1,height=0.4,width=1,cex=1, at=c(-0.001,-0.0005,0,0.0005,0.001),labels=c("","Erosion","","Accretion","")),at= c(min(chng.eras,na.rm=TRUE)
,0,max(chng.eras,na.rm=TRUE)),col.regions=c("red","red","blue","blue")
,xlab=list(label="Transect", cex=0.8),ylab=list(label="Era", cex=0.8),scale = list(x = list(labels=Transect.Labels),
y = list(labels=Window.Labels,cex=0.5)),page = function(n) grid.text(paste(""), x = 0.5275, y = 0.00001, default.units = "npc", just = c("center", "top")),aspect="iso"))


  dev.off()
  
  
    pdf("PDF/graph_eras_grid3.pdf", bg="white")

col.l <- colorRampPalette(c('red', 'blue'))(30)


print(levelplot(chng.eras, pretty=TRUE,na.rm=TRUE,col.regions=col.l
,xlab=list(label="Transect", cex=0.8),ylab=list(label="Era", cex=0.8),scale = list(x = list(labels=Transect.Labels),
y = list(labels=Window.Labels,cex=0.5)),page = function(n) grid.text(paste(""), x = 0.5275, y = 0.00001, default.units = "npc", just = c("center", "top")),aspect=1,level.colors(chng.eras, at = c(min(chng.eras,na.rm=TRUE),0,max(chng.eras,na.rm=TRUE)) )))


  dev.off()
  

    pdf("PDF/graph_eras_grid4.pdf", bg="white")

col.l <- colorRampPalette(c('red', 'blue'))(30)


print(levelplot(chng.eras, pretty=TRUE,na.rm=TRUE,col.regions=col.l
,xlab=list(label="Transect", cex=0.8),ylab=list(label="Era", cex=0.8),scale = list(x = list(labels=Transect.Labels),
y = list(labels=Window.Labels,cex=0.5)),page = function(n) grid.text(paste(""), x = 0.5275, y = 0.00001, default.units = "npc", just = c("center", "top")),aspect="iso",level.colors(chng.eras, at = c(min(chng.eras,na.rm=TRUE),0,max(chng.eras,na.rm=TRUE)) )))


  dev.off()
    
  
  
  
  }
  
  
  
