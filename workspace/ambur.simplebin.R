`ambur.simplebin` <-
function(  ) {

#open results stats file
winDialog("ok","Please select the 'results_stats.csv' file...")

data.path <- choose.files(default = "*.csv",multi = FALSE)

mydata <- read.csv(data.path, header=TRUE, sep=",")
attach(mydata)

##require the stats and lattice package
require(stats)
require(lattice)

#cull EPR values of "NaN'

mydata2 <- mydata[mydata$EPR != "NaN",]

EPR.cut <- cut(EPR,c(-Inf,0,+Inf),labels=FALSE)

EPR.cut.prep1 <- c(0,EPR.cut[-max(length(EPR.cut))])

EPR.cut.prep2 <- c(EPR.cut[-1],0)

Start.Zone <- which(EPR.cut != EPR.cut.prep1)
End.Zone <- which(EPR.cut != EPR.cut.prep2)


Zone.Num <- numeric(length(Start.Zone))
Zone.Mean <- numeric(length(Start.Zone))
Zone.Start.Tran <- numeric(length(Start.Zone))
Zone.End.Tran <- numeric(length(Start.Zone))

Zone.Num.Tran <- numeric(length(Transect))
Zone.Mean.Tran <- numeric(length(Transect))

for (i in 1:length(Start.Zone)) {

Zone.Num[i] <- i
Zone.Mean[i] <- mean(EPR[Start.Zone[i]:End.Zone[i]],na.rm=TRUE)
Zone.Start.Tran[i] <- Transect[Start.Zone[i]]
Zone.End.Tran[i] <- Transect[End.Zone[i]]


for (j in 1:length(Transect)) {

Zone.Num.Tran[j] <- ifelse(Transect[j] >= Start.Zone[i] & Transect[j] <= End.Zone[i],Zone.Num[i],Zone.Num.Tran[j]*1)

Zone.Mean.Tran[j] <- ifelse(Transect[j] >= Start.Zone[i] & Transect[j] <= End.Zone[i],Zone.Mean[i],Zone.Mean.Tran[j]*1)

}
}


Zone.Summary <- cbind(Zone.Num,Start.Zone,End.Zone,Zone.Start.Tran,Zone.End.Tran,Zone.Mean)



Transect.Zone.Summary <- cbind(Transect,EPR,EPR.cut,EPR.cut.prep1,EPR.cut.prep2,Zone.Num.Tran,Zone.Mean.Tran)




############plot with transect zones
par(mfrow=(c(1,2)))
par(pty= "m")
plot(c(Transect.StartX,rev(Transect.EndX)), c(Transect.StartY,rev(Transect.EndY)),col=0,asp=TRUE,xlab="X",ylab="Y",main="Zones (Transects)")

Transect.StartX.p <- numeric(length(Transect))
Transect.StartY.p <- numeric(length(Transect))
Transect.EndX.p <- numeric(length(Transect))
Transect.EndY.p <- numeric(length(Transect)) 


for (i in 1:length(Start.Zone)) {

Zone.Num[i] <- i

for (j in 1:length(Transect)) {

Transect.StartX.p[j] <- ifelse(Transect[j] >= Start.Zone[i] & Transect[j] <= End.Zone[i],Transect.StartX[j],Transect.StartX.p[j]*0)

Transect.StartY.p[j] <- ifelse(Transect[j] >= Start.Zone[i] & Transect[j] <= End.Zone[i],Transect.StartY[j],Transect.StartY.p[j]*0)

Transect.EndX.p[j] <- ifelse(Transect[j] >= Start.Zone[i] & Transect[j] <= End.Zone[i],Transect.EndX[j],Transect.EndX.p[j]*0)

Transect.EndY.p[j] <- ifelse(Transect[j] >= Start.Zone[i] & Transect[j] <= End.Zone[i],Transect.EndY[j],Transect.EndY.p[j]*0)

X1 <- Transect.StartX.p[Transect.StartX.p != 0]
X2 <- rev(Transect.EndX.p[Transect.StartX.p != 0])
Y1 <- Transect.StartY.p[Transect.StartY.p != 0]
Y2 <- rev(Transect.EndY.p[Transect.StartY.p != 0])


}
polygon(c(X1,X2), c(Y1,Y2),col=ifelse(Zone.Mean[i] <0,"red","blue"),border = NA)
}
lines(Max.Date.Xcoord, Max.Date.Ycoord,pch=20,col="yellow")



######plot with shoreline envelope
plot(c(Transect.Outer.Xcoord,rev(Transect.Inner.Xcoord)), c(Transect.Outer.Ycoord,rev(Transect.Inner.Ycoord)),asp=TRUE,xlab="X",ylab="Y",col="white",main="Zones (Envelope)")

Transect.Outer.Xcoord.p <- numeric(length(Transect))
Transect.Outer.Ycoord.p <- numeric(length(Transect))
Transect.Inner.Xcoord.p <- numeric(length(Transect))
Transect.Inner.Ycoord.p <- numeric(length(Transect)) 


for (i in 1:length(Start.Zone)) {

Zone.Num[i] <- i

for (j in 1:length(Transect)) {

Transect.Outer.Xcoord.p[j] <- ifelse(Transect[j] >= Start.Zone[i] & Transect[j] <= End.Zone[i],Transect.Outer.Xcoord[j],Transect.Outer.Xcoord.p[j]*0)

Transect.Outer.Ycoord.p[j] <- ifelse(Transect[j] >= Start.Zone[i] & Transect[j] <= End.Zone[i],Transect.Outer.Ycoord[j],Transect.Outer.Ycoord.p[j]*0)

Transect.Inner.Xcoord.p[j] <- ifelse(Transect[j] >= Start.Zone[i] & Transect[j] <= End.Zone[i],Transect.Inner.Xcoord[j],Transect.Inner.Xcoord.p[j]*0)

Transect.Inner.Ycoord.p[j] <- ifelse(Transect[j] >= Start.Zone[i] & Transect[j] <= End.Zone[i],Transect.Inner.Ycoord[j],Transect.Inner.Ycoord.p[j]*0)

X1 <- Transect.Outer.Xcoord.p[Transect.Outer.Xcoord.p != 0]
X2 <- rev(Transect.Inner.Xcoord.p[Transect.Outer.Xcoord.p != 0])
Y1 <- Transect.Outer.Ycoord.p[Transect.Outer.Ycoord.p != 0]
Y2 <- rev(Transect.Inner.Ycoord.p[Transect.Outer.Ycoord.p != 0])


}
polygon(c(X1,X2), c(Y1,Y2),col=ifelse(Zone.Mean[i] <0,"red","blue"),border = NA)
}
lines(Max.Date.Xcoord, Max.Date.Ycoord,col="green")

print(Zone.Summary)

#tidy up and remove all objects
detach(mydata)
rm(list = ls())

}

