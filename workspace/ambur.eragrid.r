
ambur.eragrid <-
function(tstart=1, tend="all", downsample=0) {


#tstart = 15     #### starting transect
#tend = 305     #### ending transect
#downsample=25  ####sampling frequency to trim down the dataset for every n-th transect

require(lattice)
require(grid)
require(tcltk)

tkmessageBox(message = "This function might not work properly if 'rgeos' and 'sp' packages are loaded.")

tkmessageBox(message = "Please select the raw_positions.csv file...")
getdata <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata, header=TRUE, sep=",")    #this provides the list of raw dates for the EPR eras file
#attach(mydata)


dir.path <- dirname(getdata)
setwd(dir.path)


tkmessageBox(message = "Please select the raw_EPR_rates_eras.csv file...")
getdata2 <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata2 <- read.table(getdata2, header=TRUE, sep=",")        #actual EPR eras file with data
#attach(mydata2)






#cull the dataset based on the range of transects

#Extract (cull) the transects for analysis based on user's input range                                                               
Start.Transect <- tstart
End.Transect <- tend

End.Transect2 <- ifelse(End.Transect == "all", max(mydata[,1]), End.Transect)



mydata <- mydata[mydata[,1] >= Start.Transect & mydata[,1] <= End.Transect2,]
mydata <- if(downsample == 0) mydata[1:length(mydata[,1]),]  else mydata[c(1,seq(0,length(mydata[,1]),downsample)),]

mydata2 <- data.frame(mydata2[mydata2[,1] >= Start.Transect & mydata2[,1] <= End.Transect2,],row.names=NULL)
mydata2 <- if(downsample == 0) mydata2[1:length(mydata2[,1]),]  else mydata2[c(1,seq(0,length(mydata2[,1]),downsample)),]



chng.eras <- as.matrix(mydata2[,-1])
   
#get the date for the eras and labels
dates1 <- substr(colnames(mydata)[-1],2,11)
dates2 <- c(dates1[-1],"9999")
dates3 <- paste(dates1,dates2, sep = "-", collapse = NULL)[-(length(dates1))]

Window.Labels  <- gsub("[.]", "/", dates3)
#Transect.Labels  <- c(pretty(min((mydata2[,"TRANSECT"])):max(as.numeric(mydata2[,"TRANSECT"]))))  
#Transect.Labels  <- c(min((mydata2[,"TRANSECT"])),pretty(as.numeric(mydata2[,"TRANSECT"]))[-1]+min((mydata2[,"TRANSECT"] - 1)))  #changed to start at 1 instead of 0 

Transect.scale <-  pretty(1:length(mydata2[,1]),n=5)[-1]  #changed to start at 1 insteat of 0
Transect.Labels  <- mydata2[,"TRANSECT"][Transect.scale]
#Transect.scale  <- pretty(as.numeric(row.names(mydata2))) ####changed to label the n-th number of transect sample
#chng.eras[is.na(chng.eras) == T] <- 1    for test with NA data

#generate the grid plots
pdf("PDF/graph_eras_grid1.pdf", bg="white")



print(levelplot(chng.eras, pretty=TRUE,na.rm=TRUE,colorkey = list(space="right",tick.number=1,height=0.4,width=1,cex=1, at=c(-0.001,-0.0005,0,0.0005,0.001),labels=c("","Erosion","","Accretion","")),at= c(min(chng.eras,na.rm=TRUE)
,0,max(chng.eras,na.rm=TRUE)),col.regions=c("red","red","blue","blue")
,xlab=list(label="Transect", cex=0.8),ylab=list(label="Era", cex=0.8),scale = list(x = list(rot=0, at=Transect.scale,labels=Transect.Labels),
y = list(labels=Window.Labels,cex=0.5)),page = function(n) grid.text(paste(""), x = 0.5275, y = 0.00001, default.units = "npc", just = c("center", "top")),aspect=1))


#levelplot(chng.eras, pretty=TRUE,na.rm=TRUE,colorkey = list(space="right",tick.number=1,height=0.4,width=1,cex=1, at=c(-0.001,-0.0005,0,0.0005,0.001),labels=c("","Erosion","","Accretion","")),at= c(min(chng.eras,na.rm=TRUE)
#,0,max(chng.eras,na.rm=TRUE)),col.regions=c("red","red","blue","blue")
#,xlab=list(label="Transect", cex=0.8),ylab=list(label="Era", cex=0.8),scale = list(x = list(rot=0, at=Transect.scale,labels=Transect.Labels),
#y = list(labels=Window.Labels,cex=0.5)),page = function(n) grid.text(paste(""), x = 0.5275, y = 0.00001, default.units = "npc", just = c("center", "top")),aspect=1)


#levelplot(chng.eras, pretty=TRUE)


  dev.off()
  
  
  pdf("PDF/graph_eras_grid2.pdf", bg="white")


print(levelplot(chng.eras, pretty=TRUE,na.rm=TRUE,colorkey = list(space="right",tick.number=1,height=0.4,width=1,cex=1, at=c(-0.001,-0.0005,0,0.0005,0.001),labels=c("","Erosion","","Accretion","")),at= c(min(chng.eras,na.rm=TRUE)
,0,max(chng.eras,na.rm=TRUE)),col.regions=c("red","red","blue","blue")
,xlab=list(label="Transect", cex=0.8),ylab=list(label="Era", cex=0.8),scale = list(x = list(rot=0, at=Transect.scale,labels=Transect.Labels),
y = list(labels=Window.Labels,cex=0.5)),page = function(n) grid.text(paste(""), x = 0.5275, y = 0.00001, default.units = "npc", just = c("center", "top")),aspect="iso"))


  dev.off()
  
  
    pdf("PDF/graph_eras_grid3.pdf", bg="white")

col.l <- colorRampPalette(c('red', 'blue'))(30)


print(levelplot(chng.eras, pretty=TRUE,na.rm=TRUE,col.regions=col.l
,xlab=list(label="Transect", cex=0.8),ylab=list(label="Era", cex=0.8),scale = list(x = list(rot=0, at=Transect.scale,labels=Transect.Labels),
y = list(labels=Window.Labels,cex=0.5)),page = function(n) grid.text(paste(""), x = 0.5275, y = 0.00001, default.units = "npc", just = c("center", "top")),aspect=1,level.colors(chng.eras, at = c(min(chng.eras,na.rm=TRUE),0,max(chng.eras,na.rm=TRUE)) )))




  dev.off()
  

    pdf("PDF/graph_eras_grid4.pdf", bg="white")

col.l <- colorRampPalette(c('red', 'blue'))(30)


print(levelplot(chng.eras, pretty=TRUE,na.rm=TRUE,col.regions=col.l
,xlab=list(label="Transect", cex=0.8),ylab=list(label="Era", cex=0.8),scale = list(x = list(rot=0, at=Transect.scale,labels=Transect.Labels),
y = list(labels=Window.Labels,cex=0.5)),page = function(n) grid.text(paste(""), x = 0.5275, y = 0.00001, default.units = "npc", just = c("center", "top")),aspect="iso",level.colors(chng.eras, at = c(min(chng.eras,na.rm=TRUE),0,max(chng.eras,na.rm=TRUE)) )))


  dev.off()
    
  
 #tidy up and remove objects
rm(list = ls()) 
  
  }
  
  
  
