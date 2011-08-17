ambur.smooth <-
function(alpha=0.7, degree=2, sampledist=5,presamp=1) {

 require(tcltk)
require(rgdal)
require(locfit)
require(spatstat)

# Establish the inputs

alpha=0.7
 degree=2 
 sampledist=5
presamp = 0
bufferdist <- 10
maxsimalpha <- 0.05
maxsimdeg <- 20
nsims <- 10 # number of simulations

outer.alpha.rec <- alpha
outer.deg.rec <- degree
outersample <- sampledist
presample  <- presamp




tkmessageBox(message = "Please select a baseline to smooth...")
getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

workingdir <- dirname(getdata)
setwd(workingdir)


time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_smoothing", showWarnings=FALSE)
setwd("AMBUR_smoothing")

dir.create(paste(time.stamp2," ","smooth",sep=""))
setwd(paste(time.stamp2," ","smooth",sep="")) 



#build function for points along a line (modified to add the start and end points)

pts.along <- function(x,y,pspace) {

#VARIABLES: x = x coord, y = y coord, pspace = spacing between points

Cx <- x
Cy <- y

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Segment.Length <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

Dx <- Cx2 - Cx
Dy <- Cy2 - Cy

Segment.Azimuth <- ifelse(Dx >= 0, 90 -(180/pi) * atan(Dy/Dx),270 -(180/pi) * atan(Dy/Dx))

Cumltiv.Sum <- ave(Segment.Length, FUN=cumsum)

pts.int <- pspace  ##defines the point spacing along the line

pts.bin <- seq(from=pts.int,to=max(Cumltiv.Sum),by=pts.int)

Cumltiv.Sum2 <- c(0,Cumltiv.Sum[-max(length(Cumltiv.Sum))])


pts.bin.up <- length(pts.bin)
pts.bin.upx <- length(pts.bin)
pts.bin.upy <- length(pts.bin)
pts.bin.upx2 <- length(pts.bin)
pts.bin.upy2 <- length(pts.bin)
pts.bin.diff <- length(pts.bin)

for (i in 1:length(pts.bin)) {

pts.bin.up[i] <-  which.max(Cumltiv.Sum2[Cumltiv.Sum2  <=  pts.bin[i]])
pts.bin.upx[i] <- x[pts.bin.up[i]]
pts.bin.upy2[i] <- y[pts.bin.up[i]+1]
pts.bin.upx2[i] <- x[pts.bin.up[i]+1]
pts.bin.upy[i] <- y[pts.bin.up[i]]
pts.bin.diff[i] <- pts.bin[i] - Cumltiv.Sum2[pts.bin.up[i]]



}

Dx2 <- pts.bin.upx2 - pts.bin.upx
Dy2 <- pts.bin.upy2 - pts.bin.upy


pts.bin.az  <- ifelse(Dx2 >= 0, 90 -(180/pi) * atan(Dy2/Dx2),270 -(180/pi) * atan(Dy2/Dx2))

t.azimuth <- c(pts.bin.az[1],pts.bin.az,pts.bin.az[max(length(pts.bin.az))])

t.startx <- c(Cx[1],sin((pts.bin.az * pi/180)) * pts.bin.diff + pts.bin.upx,Cx[max(length(Cx))])
t.starty <- c(Cy[1],cos((pts.bin.az * pi/180)) * pts.bin.diff + pts.bin.upy,Cy[max(length(Cy))])


cbind(t.startx,t.starty,t.azimuth)



}



  #build the locfit smoothing function
smlocfit.line <- function(x,y,n=20,alpha.z =0.7,kern.z="tricube",deg.z = 2,span.z = 2/3, fam.z = "gaussian") {
#VARIABLES: x = x coord, y = y coord, z = degrees of freedom, n  = densify number of points

Ox <- x
Oy <- y
n <- n #densification threshold
span.z <- span.z
deg.z <- deg.z
fam.z <- fam.z
alpha.z <- alpha.z
kern.z <- kern.z

   xy <- cbind(Ox,Oy)

   nP <- length(xy[,1])
   if((nP < 3) || (is.na(any(as.logical(xy[,1])))))
       return(print("need more than 3 points"))
       #return(list(x=numeric(0), y=numeric(0))) #old return call

  ## else :

   nP <- length(xy[,1])
   z <- n*(nP-1)

   i <- 1:nP
   Sx <- function(x) {fitted.values(locfit.raw(i, xy[,1],alpha=alpha.z,deg=deg.z,kern=kern.z,kt="sph"))}
   Sy <- function(x) {fitted.values(locfit.raw(i, xy[,2],alpha=alpha.z,deg=deg.z,kern=kern.z,kt="sph"))}
   ti <- seq(1, nP, length = z)
   opspl <- cbind(x = Sx(ti), y = Sy(ti))
   opspl
   }





#################################################################
###### break down outer polyline into simply points with IDs
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

attrtable$Id <- seq(0,max(length(attrtable[,1]))-1,1)  #repair baseline attr table for sequential IDs to match with bbb$shapeID
bbb <- data.frame(merge(outerbase.tab,attrtable,by.x= "shapeID" ,by.y = "Id", all.x = TRUE,sort=FALSE))

### cast transects for individual polylines based on unique IDs
Baseline.Factor <- factor(bbb$shapeID)

#to even out the points for smoothing
outer.sample <- sapply(levels(Baseline.Factor), function(x) data.frame(pts.along(bbb$baseX[bbb$shapeID == x],bbb$baseY[bbb$shapeID == x],outersample)) ,simplify = FALSE)
outer.basepts4sm <- data.frame(do.call("rbind", outer.sample))

#round down to get ID numbers
outer.basepts4sm$t.ID <- floor(as.numeric(row.names(outer.basepts4sm)))




sim.alpha <- seq(maxsimalpha/nsims,maxsimalpha,maxsimalpha/nsims)
sim.deg <- seq(2,maxsimdeg,1)

outer.baseptsSmooth <- numeric(length(sim.deg))
outersm.max <- numeric(length(sim.deg))
out.alpha  <- numeric(length(sim.deg))
out.deg  <- numeric(length(sim.deg))
outer.alpha <- 0
outer.deg <- 0
outer.maxdist <- 0

plotx <- rep(sim.alpha,length(sim.deg))
ploty <- rep(sim.deg,length(sim.alpha))
 plot(plotx,ploty,,xlab="Alpha",ylab="Degree",main="LocFit Curve Deviations (Outer Baseline)",col="light gray")


for (i in 1:length(sim.alpha)) {


for (j in 1:length(sim.deg)) {
#outer.baseptsSmooth <- smlocfit.line(outer.basepts4sm[,1], outer.basepts4sm[,2],alpha=outer.alpha.rec,kern="epan",n=1,deg=outer.deg.rec)
BaselineSM.Factor <- factor(outer.basepts4sm$t.ID)

outer.baseptsSmooth1 <- sapply(levels(BaselineSM.Factor), function(x) data.frame(smlocfit.line(outer.basepts4sm[,1][outer.basepts4sm$t.ID == x],outer.basepts4sm[,2][outer.basepts4sm$t.ID == x],alpha.z=sim.alpha[i],kern.z="epan",n=1,deg.z=sim.deg[j])) ,simplify = FALSE)


#outer.baseptsSmooth <- smlocfit.line(outer.basepts4sm[,1], outer.basepts4sm[,2],alpha=sim.alpha[i],kern="epan",n=1,deg=sim.deg[j])

outer.baseptsSmooth <- data.frame(do.call("rbind", outer.baseptsSmooth1))

outersm.distdiff <-  (((outer.baseptsSmooth[,1] - outer.basepts4sm[,1])^2 +  (outer.baseptsSmooth[,2] - outer.basepts4sm[,2])^2)^(1/2))

outersm.max[j] <- max(outersm.distdiff)

out.alpha[j] <- sim.alpha[i]

out.deg[j] <- sim.deg[j]



#plot(sim.alpha[i],sim.deg[j])


points(cbind(outer.alpha[-1],out.alpha),cbind(outer.deg[-1],out.deg), cex=1,pch=19,col=ifelse((cbind(outer.maxdist[-1],outersm.max))>=bufferdist,"red","green"))

}

outer.alpha <- c(outer.alpha,out.alpha)
outer.deg <- c(outer.deg,out.deg)
outer.maxdist <- c(outer.maxdist,outersm.max)



}


outer.tests <- cbind(outer.alpha,outer.deg,outer.maxdist)[-1,]


#experimental for ternary plots (convert each column to percentages) #variable needed
ternary.prep <- cbind(outer.alpha/max(outer.alpha),outer.deg/max(outer.deg),outer.maxdist/max(outer.maxdist))[-1,]


outer.alpha.rec <- max(outer.tests[,1][outer.tests[,3]<bufferdist])
outer.deg.rec <- min(outer.tests[,2][outer.tests[,3]<bufferdist & outer.tests[,1] == outer.alpha.rec])


plot(outer.tests[,1],outer.tests[,2],cex=(outer.tests[,3]/max(outer.tests[,3]))*3,pch=19,col=ifelse(outer.tests[,3]>=bufferdist,"red","green"),xlab="Alpha",ylab="Degree",main="LocFit Curve Deviations (Outer Baseline)")
points(outer.alpha.rec,outer.deg.rec,pch=1,cex=3,col="dark green")

dev.new()


plot(outer.basepts4sm[,1], outer.basepts4sm[,2],type="l",col="green" ,asp=1,xlab="X",ylab="Y")
lines(outer.baseptsSmooth[,1],outer.baseptsSmooth[,2],col="blue")







}

