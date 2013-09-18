ambur.analysis <-
function(userinput1="first", userinput2=95, userinput3="m", userinput4="", userinput5=1, userinput6="all",userinput7="basic",userinput8="yr") {


#require(tcltk)
#require(rgdal)
#require(rgeos)



#userinput1 <- "first"
#userinput2 <- 95
#userinput3 <- "m"
#userinput4 <- ""
#userinput5 <- 1
#userinput6 <- "all"
#userinput7 <- "basic"
#userinput8 <- "yr"


# Establish the inputs
if(userinput1 == "first") Intersect.Position <- min  else  Intersect.Position <- max
ConfInter <- userinput2 / 100
Map.Units <- userinput3
Start.Transect <- userinput5
End.Transect <- userinput6
Basic.Analysis <- ifelse(userinput7 == "basic",1,0)
Time.Units <- userinput8

  if(userinput8=="yr"){time.adj <- 31557600} else
    if(userinput8=="day"){time.adj <- 86400} else
     if(userinput8=="hr"){time.adj <- 3600} else
      if(userinput8=="min"){time.adj <- 60} else
       if(userinput8=="sec"){time.adj <- 1}



#kill any open windows/devices
graphics.off()

# NOW COMPATIBLE FOR TRANSECTS FROM AN OFFSHORE & ONSHORE BASELINE


#turn on the foreign, stats, and lattice package"
#require(foreign)
#require(stats)
#require(lattice)
#require(tcltk)
#require(nlme)
#require(MASS)
#require(grid)

#choose dbf file to import
tkmessageBox(message = "Please select the capture points shapefile...")
filetype <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files("","Choose file",multi = FALSE,filetype,1)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
mydata <- data.frame(shapedata)


workingdir <- dirname(getdata)
setwd(workingdir)
path <- getdata

#mydata <- foreign::read.dbf(path)

#status log and checkpoint
dev.new(width = 8, height = 4)
par("usr")[1]
par(mar=c(5,25,5,5))
plot(1, type="n", axes=F, xlab="", ylab="")
par("usr")
textalign <- par("usr")[1] - ((par("usr")[2] - par("usr")[1]) * 1.464759)
mtext("AMBUR: Shoreline Change Analysis",line= 3, cex= 0.75, at=textalign,adj=0)
mtext("-starting analyses...", cex= 0.75,line=2, at=textalign,adj=0)




#start timing the analysis
time0 <- Sys.time()

#set working directory and place results in a time stamped folder

time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

path2 <- dirname(path)
setwd(path2)

dir.create("AMBUR_results", showWarnings=FALSE)
setwd("AMBUR_results")

dir.create(paste(time.stamp2," ",userinput4,sep=""))
setwd(paste(time.stamp2," ",userinput4,sep=""))

dir.create("PDF")


results.path <- getwd()

#make sure all column names are upper case
colnames(mydata) <- toupper(colnames(mydata))

#Repair ArcGIS "DATE_", and XY fields

colnames(mydata) <- gsub("DATE_", "DATE", colnames(mydata))
colnames(mydata) <- gsub("POINT_X", "X_COORD", colnames(mydata))
colnames(mydata) <- gsub("POINT_Y", "Y_COORD", colnames(mydata))

#add unique ID to each point in the shoreline field  (added 6/2008)
mydata[,"ID"] <- seq(1:length(mydata[,"ID"]))

#Fix DATE field to include time
if(nchar(as.character(mydata[ ,"DATE"]), type = "chars", allowNA = FALSE) < 15){mydata[ ,"DATE"] <- paste(mydata[ ,"DATE"],"12:00:01 AM",sep=" ")}

#interactive Date Selection and removal of extaneous dates

DateList1 <- as.POSIXlt(mydata[ ,"DATE"],format="%m/%d/%Y %I:%M:%S %p")

DateList2 <- as.character(as.POSIXlt(sort(unique(DateList1)),format="%m/%d/%Y %I:%M:%S %p"))

DateList3 <- select.list(DateList2, multiple = TRUE, title = "AMBUR! Choose dates:")


mydata <- mydata[as.character(DateList1) %in% DateList3,]


#Extract (cull) the transects for analysis based on user's input range
End.Transect2 <- ifelse(End.Transect == "all", max(mydata$TRANSECT), End.Transect)

mydata <- mydata[mydata$TRANSECT >= Start.Transect & mydata$TRANSECT <= End.Transect2,]






#converts (overwrites) DATE column to true numerical dates readable by R
Setup.Date <- as.POSIXlt(mydata[ ,"DATE"], format="%m/%d/%Y %I:%M:%S %p")
Setup.Date2 <- as.numeric(Setup.Date - as.POSIXlt("01/01/1970 12:0:0 AM", format="%m/%d/%Y %I:%M:%S %p"),units="secs")



#repair the datatable by keep only the first or last intersections of shorelines (fixed 6/2008)

testsort <- cbind(mydata[ ,"TRANSECT"], Setup.Date2, mydata[ ,"DISTANCE"], mydata[ ,"ACCURACY"],mydata[ ,"ID"])[ order(mydata[ ,"TRANSECT"], Setup.Date2, mydata[ ,"DISTANCE"], mydata[ ,"ACCURACY"],mydata[ ,"ID"]) ,]

colnames(testsort) <- c("TRANSECT", "DATE", "DISTANCE", "ACCURACY","ID")
row.names(testsort) <- rownames(testsort, do.NULL = FALSE, prefix = "")




mini.pts <- aggregate(testsort, by=list(TR=testsort[,"TRANSECT"], DA=testsort[, "DATE"]), Intersect.Position)

sor.mini.pts <- (mini.pts)[order(mini.pts[,"TRANSECT"]),]

final.mini.pts <- sor.mini.pts[,-(1:2)]


#repair the dates in both tables before merging tables
final.mini.pts[,"DATE"] <- as.character(format(as.POSIXct(final.mini.pts[,"DATE"], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p"))

mydata[,"DATE"] <- as.character(format(as.POSIXct(Setup.Date2, origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p"))

#join the worktable by ID (added 6/2008)
worktable1 <- merge(final.mini.pts, mydata, by="ID")
colnames(worktable1) <- gsub(".x", "", colnames(worktable1))
colnames(worktable1) <- gsub(".y", "", colnames(worktable1))


#!#write.table(worktable1, file = "GIS_table_cleaned.csv", sep = ",", row.names = FALSE)



#tidy up
mydata <- worktable1

# correct for onshore and offshore transect positions
mydata$DISTANCE <- ifelse(mydata$OFFSHORE == 1, mydata$DISTANCE * -1, mydata$DISTANCE * 1)


remove(final.mini.pts, mini.pts, Setup.Date, Setup.Date2, sor.mini.pts, testsort, worktable1)




#status checkpoint
mtext("-building data tables...", side=3, line= 1, adj= 0, cex= 0.75, at=textalign)



# BEGIN ANALYSES

#make sure all column names are upper case
colnames(mydata) <- toupper(colnames(mydata))



#converts (overwrites) DATE column to true numerical dates readable by R
Setup.Date <- as.POSIXlt(mydata[ ,"DATE"], format="%m/%d/%Y %I:%M:%S %p")
Setup.Date2 <- as.numeric(Setup.Date - as.POSIXlt("01/01/1970 12:0:0 AM", format="%m/%d/%Y %I:%M:%S %p"),units="secs")



#################################################
#get all changes and rates

 #get unique dates and transects

unq.dates <- sort(unique(Setup.Date2))
unq.transects <- sort(unique(mydata$TRANSECT))

t1 <- unq.transects
t2 <- unq.dates
t3 <- merge(t1,t2)
colnames(t3) <- c("TRANSECT", "DATE")

mydata[,"DATE"] <- as.numeric(Setup.Date - as.POSIXlt("01/01/1970 12:0:0 AM", format="%m/%d/%Y %I:%M:%S %p"),units="secs")

mydataPrep <- mydata[,!duplicated(colnames(mydata))] ###added to remove duplicate column names 11/19/2012

t4 <- merge(t3, mydataPrep, all.x=T,by.y=c(2,3)) ###fixed to remove duplicate column names 11/19/2012

t5 <- t4[order(t4$TRANSECT, t4$DATE),]



#create matrix to fill the distance values
testmatrix <- matrix(data = t5$DISTANCE, ncol = length(unq.transects), nrow = length(unq.dates), byrow = FALSE, dimnames = NULL)

row.names(testmatrix) <- rownames(unq.dates, do.NULL = FALSE, prefix = "")

colnames(testmatrix) <- unq.transects

testmatrix2 <-  t(testmatrix)

row.names(testmatrix2) <- rownames(unq.transects, do.NULL = FALSE, prefix = "")

colnames(testmatrix2) <- unq.dates

testmatrix3 <- cbind(unq.transects, testmatrix2)

colnames(testmatrix3)[-1] <- as.character(format(as.POSIXct(unq.dates, origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p"))
colnames(testmatrix3)[1] <- "TRANSECT"


#this works to get all changes
tt3 <- as.numeric(colnames(testmatrix2))
tt1 <- sort(rep(tt3, length(tt3)))
tt2 <- rep(tt3, length(tt3))
subs.dates <- which(tt2 >= tt1)

subs.dates.eras <- which(tt2 > tt1)

###NEW stuff to get consecutive eras and rates only
tt4 <- c(tt3[-1],0)
subs.dates.eras <- which(tt4 > tt3 & tt4 != 0)

chng.eras <- as.matrix((testmatrix2[,as.character(tt4)[subs.dates.eras]] - testmatrix2[,as.character(tt3)[subs.dates.eras]]))

fname.eras <- paste(as.character(format(as.POSIXct((tt3)[subs.dates.eras], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")),"to",as.character(format(as.POSIXct((tt4)[subs.dates.eras], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")))




colnames(chng.eras) <- (fname.eras)

ellyrs.eras <- as.matrix(((tt4)[subs.dates.eras] - (tt3)[subs.dates.eras]) / time.adj)

ellyrs2.eras <- matrix(data = ellyrs.eras, ncol = length(ellyrs.eras), nrow = length(unq.transects), byrow = TRUE, dimnames = NULL)

rates.eras <- chng.eras / ellyrs2.eras

fname2.eras <- paste(as.character(format(as.POSIXct((tt3)[subs.dates.eras], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")),"to",as.character(format(as.POSIXct((tt4)[subs.dates.eras], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")), "EPR")

colnames(rates.eras) <- (fname2.eras)

#bind them to transects for output file
chng.eras2 <- cbind(unq.transects,chng.eras)
colnames(chng.eras2)[1] <- "TRANSECT"

rates.eras2 <- cbind(unq.transects,rates.eras)
colnames(rates.eras2)[1] <- "TRANSECT"

######### get all possible combinations of rates and changes
ellyrs <- as.matrix(((tt2)[subs.dates] - (tt1)[subs.dates]) / time.adj)

ellyrs2 <- matrix(data = ellyrs, ncol = length(ellyrs), nrow = length(unq.transects), byrow = TRUE, dimnames = NULL)


#delete chng multiply value of *-1 to adjust for onshore transects
chng <- as.matrix((testmatrix2[,as.character(tt2)[subs.dates]] - testmatrix2[,as.character(tt1)[subs.dates]]))

fname <- paste(as.character(format(as.POSIXct((tt1)[subs.dates], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")),"to",as.character(format(as.POSIXct((tt2)[subs.dates], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")))



colnames(chng) <- (fname)


rates <- chng / ellyrs2

fname2 <- paste(as.character(format(as.POSIXct((tt1)[subs.dates], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")),"to",as.character(format(as.POSIXct((tt2)[subs.dates], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")), "EPR")

colnames(rates) <- (fname2)

chng <- cbind(unq.transects,chng)
colnames(chng)[1] <- "TRANSECT"

rates <- cbind(unq.transects,rates)
colnames(rates)[1] <- "TRANSECT"


finalchngs <- cbind(testmatrix3, chng, rates)


#set up matrix for each class and changes
#Class 1:

Shore.Class1 <- matrix(data = t5$CLASS_1, ncol = length(unq.transects), nrow = length(unq.dates), byrow = FALSE, dimnames = NULL)

row.names(Shore.Class1) <- rownames(unq.dates, do.NULL = FALSE, prefix = "")

colnames(Shore.Class1) <- unq.transects

Shore.Class1.tran <-  t(Shore.Class1)

row.names(Shore.Class1.tran) <- rownames(unq.transects, do.NULL = FALSE, prefix = "")

colnames(Shore.Class1.tran) <- unq.dates

Shore.Class1.raw <- cbind(unq.transects, Shore.Class1.tran)

colnames(Shore.Class1.raw)[-1] <- as.character(format(as.POSIXct(unq.dates, origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p"))
colnames(Shore.Class1.raw)[1] <- "TRANSECT"

tt4 <- c(tt3[-1],0)
subs.dates.eras <- which(tt4 > tt3 & tt4 != 0)

Shore.Class1.chng <- as.matrix((Shore.Class1.tran[,as.character(tt4)[subs.dates.eras]] != Shore.Class1.tran[,as.character(tt3)[subs.dates.eras]]))

Shore.Class1.chng <- ifelse(Shore.Class1.chng == TRUE,1,0)

fname.eras <- paste(as.character(format(as.POSIXct((tt3)[subs.dates.eras], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")),"to",as.character(format(as.POSIXct((tt4)[subs.dates.eras], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")))

colnames(Shore.Class1.chng) <- (fname.eras)

Shore.Class1.chng.out <- cbind(unq.transects,Shore.Class1.chng)

colnames(Shore.Class1.chng.out)[1] <- "TRANSECT"

#Class 2:

Shore.Class2 <- matrix(data = t5$CLASS_2, ncol = length(unq.transects), nrow = length(unq.dates), byrow = FALSE, dimnames = NULL)

row.names(Shore.Class2) <- rownames(unq.dates, do.NULL = FALSE, prefix = "")

colnames(Shore.Class2) <- unq.transects

Shore.Class2.tran <-  t(Shore.Class2)

row.names(Shore.Class2.tran) <- rownames(unq.transects, do.NULL = FALSE, prefix = "")

colnames(Shore.Class2.tran) <- unq.dates

Shore.Class2.raw <- cbind(unq.transects, Shore.Class2.tran)

colnames(Shore.Class2.raw)[-1] <- as.character(format(as.POSIXct(unq.dates, origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p"))
colnames(Shore.Class2.raw)[1] <- "TRANSECT"

tt4 <- c(tt3[-1],0)
subs.dates.eras <- which(tt4 > tt3 & tt4 != 0)

Shore.Class2.chng <- as.matrix((Shore.Class2.tran[,as.character(tt4)[subs.dates.eras]] != Shore.Class2.tran[,as.character(tt3)[subs.dates.eras]]))

Shore.Class2.chng <- ifelse(Shore.Class2.chng == TRUE,1,0)

fname.eras <- paste(as.character(format(as.POSIXct((tt3)[subs.dates.eras], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")),"to",as.character(format(as.POSIXct((tt4)[subs.dates.eras], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")))

colnames(Shore.Class2.chng) <- (fname.eras)

Shore.Class2.chng.out <- cbind(unq.transects,Shore.Class2.chng)

colnames(Shore.Class2.chng.out)[1] <- "TRANSECT"

#Class 3:

Shore.Class3 <- matrix(data = t5$CLASS_3, ncol = length(unq.transects), nrow = length(unq.dates), byrow = FALSE, dimnames = NULL)

row.names(Shore.Class3) <- rownames(unq.dates, do.NULL = FALSE, prefix = "")

colnames(Shore.Class3) <- unq.transects

Shore.Class3.tran <-  t(Shore.Class3)

row.names(Shore.Class3.tran) <- rownames(unq.transects, do.NULL = FALSE, prefix = "")

colnames(Shore.Class3.tran) <- unq.dates

Shore.Class3.raw <- cbind(unq.transects, Shore.Class3.tran)

colnames(Shore.Class3.raw)[-1] <- as.character(format(as.POSIXct(unq.dates, origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p"))
colnames(Shore.Class3.raw)[1] <- "TRANSECT"

tt4 <- c(tt3[-1],0)
subs.dates.eras <- which(tt4 > tt3 & tt4 != 0)

Shore.Class3.chng <- as.matrix((Shore.Class3.tran[,as.character(tt4)[subs.dates.eras]] != Shore.Class3.tran[,as.character(tt3)[subs.dates.eras]]))

Shore.Class3.chng <- ifelse(Shore.Class3.chng == TRUE,1,0)

fname.eras <- paste(as.character(format(as.POSIXct((tt3)[subs.dates.eras], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")),"to",as.character(format(as.POSIXct((tt4)[subs.dates.eras], origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")))

colnames(Shore.Class3.chng) <- (fname.eras)

Shore.Class3.chng.out <- cbind(unq.transects,Shore.Class3.chng)

colnames(Shore.Class3.chng.out)[1] <- "TRANSECT"




# write all of the output tables

finalchngs.trans <- t(finalchngs)

write.table(testmatrix3, file = "raw_positions.csv", sep = ",", row.names = FALSE)

write.table(chng, file = "raw_changes.csv", sep = ",", row.names = FALSE)

write.table(rates, file = "raw_EPR_rates.csv", sep = ",", row.names = FALSE)

write.table(chng.eras2, file = "raw_changes_eras.csv", sep = ",", row.names = FALSE)

write.table(rates.eras2, file = "raw_EPR_rates_eras.csv", sep = ",", row.names = FALSE)

write.table(finalchngs, file = "all_combos.csv", sep = ",", row.names = FALSE)

write.table(finalchngs.trans, file = "all_combos_transposed.csv", sep = ",", row.names = TRUE)

write.table(Shore.Class1.raw, file = "Class_1_raw.csv", sep = ",", row.names = FALSE)

write.table(Shore.Class1.chng.out, file = "Class_1_change.csv", sep = ",", row.names = FALSE)

write.table(Shore.Class2.raw, file = "Class_2_raw.csv", sep = ",", row.names = FALSE)

write.table(Shore.Class2.chng.out, file = "Class_2_change.csv", sep = ",", row.names = FALSE)

write.table(Shore.Class3.raw, file = "Class_3_raw.csv", sep = ",", row.names = FALSE)

write.table(Shore.Class3.chng.out, file = "Class_3_change.csv", sep = ",", row.names = FALSE)


#repair the date field
mydata[,"DATE"] <- as.character(format(as.POSIXct(Setup.Date2, origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p"))



#################################################
# get min date and distance changes for each individual transect to merge with mydata table
Setup.Transect <- unique(mydata[ ,"TRANSECT"])
Setup.T.Min.Date <- numeric(length(Setup.Transect))
Setup.Min.Date.Position <- numeric(length(Setup.Transect))
Setup.Min.Date.Dist <- numeric(length(Setup.Transect))
Setup.Date <- as.numeric(Setup.Date - as.POSIXlt("01/01/1970 12:0:0 AM", format="%m/%d/%Y %I:%M:%S %p"),units="secs")


  for (i in 1:length(Setup.Transect)) {

if (i == 1) {screen(1)}

Setup.T.Min.Date[i] <- min((Setup.Date)[mydata[ ,"TRANSECT"] == Setup.Transect[i]])

Setup.Min.Date.Position[i] <- which.min((Setup.Date)[mydata[ ,"TRANSECT"] == Setup.Transect[i]])

Setup.Min.Date.Dist[i] <- mydata[ ,"DISTANCE"][mydata[ ,"TRANSECT"] == Setup.Transect[i]][Setup.Min.Date.Position[i]]


}

#create minimum change table for merge
ChangeTable <- cbind(Setup.Transect, Setup.T.Min.Date, Setup.Min.Date.Position, Setup.Min.Date.Dist)

#merge vectors
ChangeTable2 <- ChangeTable[,!duplicated(colnames(ChangeTable))]   #####added to remove duplicates  11/19/2012

WorkTable1 <- merge(mydata, ChangeTable2, by.x = 2, by.y = "Setup.Transect", all = FALSE, sort = TRUE)

Setup.Date_repair <- as.character(as.POSIXlt(WorkTable1[ ,"DATE"], origin= "01/01/1970 12:0:0 AM", format="%m/%d/%Y %I:%M:%S %p"))

WorkTable1 <- (WorkTable1)[order(WorkTable1[ ,"TRANSECT"], Setup.Date_repair) ,]





#now calculate the new Changes fields for distances and dates (correct for onshore * -1)

Changes.Distances <- (WorkTable1[ ,"DISTANCE"] - WorkTable1[ ,"Setup.Min.Date.Dist"])

chng.dt1 <- as.numeric(as.POSIXlt(WorkTable1[ ,"DATE"], format="%m/%d/%Y %I:%M:%S %p") - as.POSIXlt("01/01/1970 12:0:0 AM", format="%m/%d/%Y %I:%M:%S %p"),units="secs")


Changes.Years <- (chng.dt1 - WorkTable1[ ,"Setup.T.Min.Date"])/time.adj


Changes.Transects <- WorkTable1[ ,"TRANSECT"]

Changes.Dates <- WorkTable1[ ,"DATE"]

Changes.Dist.Eras.Raw <- c(0, Changes.Distances[-1] - Changes.Distances[1:(length(Changes.Distances)-1)])

Changes.Dist.Eras <- ifelse(Changes.Distances == 0, Changes.Dist.Eras.Raw * 0,  Changes.Dist.Eras.Raw * 1)

Changes.Years.Eras.Raw <- c(0, Changes.Years[-1] - Changes.Years[1:(length(Changes.Years)-1)])

Changes.Years.Eras <- ifelse(Changes.Years == 0, Changes.Years.Eras.Raw * 0,  Changes.Years.Eras.Raw * 1)

Rates.Consec.Eras.Raw <- Changes.Dist.Eras / Changes.Years.Eras

Rates.Consec.Eras <-  ifelse(Rates.Consec.Eras.Raw == "NaN", 0,  Rates.Consec.Eras.Raw * 1)


Changes.Dates2 <- WorkTable1[ ,"DATE"]


WorkTable1a <-cbind(Changes.Transects, Changes.Dates, Changes.Dates2, Changes.Years, Changes.Years.Eras, Changes.Distances, Changes.Dist.Eras, Rates.Consec.Eras)

######googleViz output table
Setup.Date <- as.POSIXlt(WorkTable1[ ,"DATE"], origin= "01/01/1970 12:0:0 AM", format="%m/%d/%Y %I:%M:%S %p")

Year.axis <- as.numeric(format(Setup.Date, "%Y"))


#fix transect id to get it to plot alphabetically

fill_zeros <- paste("%","0",max(nchar(WorkTable1[ ,"TRANSECT"])),"d",sep="")
Transect.axis.setup <- sprintf(fill_zeros, WorkTable1[ ,"TRANSECT"])
Transect.axis <- paste("T",Transect.axis.setup,sep="")

googleVizData <- data.frame(Transect.axis, Year.axis, WorkTable1[ ,"TRANSECT"],	as.Date(Setup.Date, format="%m/%d/%Y %I:%M:%S %p"),	WorkTable1[ ,"ACCURACY"],	WorkTable1[ ,"TRANSPACE"],	WorkTable1[ ,"TRANDIST"],	WorkTable1[ ,"LOCATION"],	WorkTable1[ ,"BASE_LOC"],	WorkTable1[ ,"STARTX"],	WorkTable1[ ,"STARTY"],	WorkTable1[ ,"ENDX"],	WorkTable1[ ,"ENDY"],	WorkTable1[ ,"AZIMUTH"],	WorkTable1[ ,"SHORE_LOC"],	WorkTable1[ ,"CLASS_1"],	WorkTable1[ ,"CLASS_2"],	WorkTable1[ ,"CLASS_3"],	WorkTable1[ ,"X_COORD"],	WorkTable1[ ,"Y_COORD"],	WorkTable1[ ,"DISTANCE"],	WorkTable1a[ ,"Changes.Years"],	WorkTable1a[ ,"Changes.Years.Eras"],	WorkTable1a[ ,"Changes.Distances"],	WorkTable1a[ ,"Changes.Dist.Eras"],	WorkTable1a[ ,"Rates.Consec.Eras"])


colnames(googleVizData) <- c("Transect",	"Year",	"Transect_number",	"Date",	"Accuracy",	"Transect_spacing",	"Transect_distance",	"Location",	"Baseline_location",	"Start_X",	"Start_Y",	"End_X",	"End_Y",	"Azimuth",	"Shoreline_location",	"Class_1",	"Class_2",	"Class_3",	"Position_X",	"Position_Y",	"Change_Distance_Raw",	"Change_Years",	"Change_Years_Eras",	"Change_Distance_Cumulative",	"Change_Distance_Eras",	"Rate_Consecutive_Eras")


d.Transect <- "00"
d.Year <- min(googleVizData$Year)
d.Transect_number <- 0
d.Date <- googleVizData$Date[1]
d.Accuracy <- 0
d.Transect_spacing <- 0
d.Transect_distance <- 0
d.Location <- "undefined"
d.Baseline_location <- "undefined"

testdiff <- ((max(googleVizData$Start_X) - min(googleVizData$Start_X)) - (max(googleVizData$Start_Y) - min(googleVizData$Start_Y)))
d.Start_X <- ifelse(testdiff > 0, max(googleVizData$Start_X), max(googleVizData$Start_X + abs(testdiff)))
d.Start_Y <- ifelse(testdiff < 0, max(googleVizData$Start_Y), max(googleVizData$Start_Y + abs(testdiff)))

testdiff <- ((max(googleVizData$End_X) - min(googleVizData$End_X)) - (max(googleVizData$End_Y) - min(googleVizData$End_Y)))
d.End_X <-  ifelse(testdiff > 0, max(googleVizData$End_X), max(googleVizData$End_X + abs(testdiff)))
d.End_Y <-  ifelse(testdiff < 0, max(googleVizData$End_Y), max(googleVizData$End_Y + abs(testdiff)))

d.Azimuth <- 0
d.Shoreline_location <- "undefined"
d.Class_1 <- "undefined"
d.Class_2 <- "undefined"
d.Class_3 <- "undefined"

testdiff <- ((max(googleVizData$Position_X) - min(googleVizData$Position_X)) - (max(googleVizData$Position_Y) - min(googleVizData$Position_Y)))
d.Position_X <- ifelse(testdiff > 0, max(googleVizData$Position_X), max(googleVizData$Position_X + abs(testdiff)))
d.Position_Y <- ifelse(testdiff < 0, max(googleVizData$Position_Y), max(googleVizData$Position_Y + abs(testdiff)))

d.Distance <- 0
d.Change_Years <- 0
d.Change_Years_Eras <- 0
d.Change_Distance <- 0
d.Change_Distance_Eras <- 0
d.Rate_Consecutive_Eras <- 0

d.googleVizData <- data.frame(d.Transect,	d.Year,	d.Transect_number,	d.Date,	d.Accuracy,	d.Transect_spacing,	d.Transect_distance,	d.Location,	d.Baseline_location,	d.Start_X,	d.Start_Y,	d.End_X,	d.End_Y,	d.Azimuth,	d.Shoreline_location,	d.Class_1,	d.Class_2,	d.Class_3,	d.Position_X,	d.Position_Y,	d.Distance,	d.Change_Years,	d.Change_Years_Eras,	d.Change_Distance,	d.Change_Distance_Eras,	d.Rate_Consecutive_Eras)


colnames(d.googleVizData) <- c("Transect",	"Year",	"Transect_number",	"Date",	"Accuracy",	"Transect_spacing",	"Transect_distance",	"Location",	"Baseline_location",	"Start_X",	"Start_Y",	"End_X",	"End_Y",	"Azimuth",	"Shoreline_location",	"Class_1",	"Class_2",	"Class_3",	"Position_X",	"Position_Y",	"Change_Distance_Raw",	"Change_Years",	"Change_Years_Eras",	"Change_Distance_Cumulative",	"Change_Distance_Eras",	"Rate_Consecutive_Eras")

googleVizData.tab <- rbind(googleVizData, d.googleVizData)

#adjust coordinates to display in googleViz


googleVizData.tab$Start_X <- googleVizData.tab$Start_X - min(googleVizData.tab$Start_X)
googleVizData.tab$Start_Y <- googleVizData.tab$Start_Y - min(googleVizData.tab$Start_Y)
googleVizData.tab$End_X <- googleVizData.tab$End_X - min(googleVizData.tab$End_X)
googleVizData.tab$End_Y <- googleVizData.tab$End_Y - min(googleVizData.tab$End_Y)
googleVizData.tab$Position_X <- googleVizData.tab$Position_X - min(googleVizData.tab$Position_X)
googleVizData.tab$Position_Y <- googleVizData.tab$Position_Y - min(googleVizData.tab$Position_Y)

googleVizData.tab$Size_1 <- 1
googleVizData.tab$Size_1[length(googleVizData.tab$Size_1)] <- 100

##fix NA values of text fields
googleVizData.tab$Baseline_location[is.na(googleVizData.tab$Baseline_location)] <- "undefined"
googleVizData.tab$Shoreline_location[is.na(googleVizData.tab$Shoreline_location)] <- "undefined"
googleVizData.tab$Location[is.na(googleVizData.tab$Location)] <- "undefined"
googleVizData.tab$Class_1[is.na(googleVizData.tab$Class_1)] <- "undefined"
googleVizData.tab$Class_2[is.na(googleVizData.tab$Class_2)] <- "undefined"
googleVizData.tab$Class_3[is.na(googleVizData.tab$Class_3)] <- "undefined"


write.table(googleVizData.tab, file = "ambur_motionchart_data.csv", sep = ",", row.names = TRUE)

#########

remove(Setup.T.Min.Date, Setup.Min.Date.Position, Setup.Min.Date.Dist, "i")

#attach(WorkTable1)

WorkTable1df <- as.data.frame(WorkTable1)

write.table(WorkTable1, file = "debugging3WorkTable1.csv", sep = ",", row.names = TRUE)
#1#write.table(WorkTable1a, file = "debugging3WorkTable1a.csv", sep = ",", row.names = TRUE)



#converts (overwrites) DATE column to true numerical dates readable by R
WorkTable1df$DATE <- as.POSIXlt(WorkTable1df$DATE, format="%m/%d/%Y %I:%M:%S %p")
WorkTable1df$DATE2 <- as.numeric(WorkTable1df$DATE - as.POSIXlt("01/01/1970 12:0:0 AM", format="%m/%d/%Y %I:%M:%S %p"),units="secs")





#status checkpoint
mtext("-constructing master data table...", side=3, line= 0, adj= 0, cex= 0.75, at=textalign)

pb <- tkProgressBar("AMBUR: progress bar", "Some information in %", 0, max(length(unique(WorkTable1df$TRANSECT))), 50)






#begin contructing the master data table



#setup the number of transect fields
Transect <- unique(WorkTable1df$TRANSECT)

# get the baseline location
Baseline.Offshore <- numeric(length(Transect))

#add start time stamps

Time.Stamp <- numeric(length(Transect))
Time.Stamp1 <- Sys.time()

# get stats for each individual transect
Transect.Spacing <- numeric(length(Transect))
Transect.Distance  <- numeric(length(Transect))
Transect.StartX <- numeric(length(Transect))
Transect.StartY <- numeric(length(Transect))
Transect.EndX <- numeric(length(Transect))
Transect.EndY <- numeric(length(Transect))
Transect.Min.Date <- numeric(length(Transect))
Transect.Max.Date <- numeric(length(Transect))
Elapsed.Years <- numeric(length(Transect))
Transect.Means <- numeric(length(Transect))
Min.Date.Position <-  numeric(length(Transect))
Max.Date.Position <-  numeric(length(Transect))
Min.Date.Dist <-  numeric(length(Transect))
Max.Date.Dist <-  numeric(length(Transect))
Net.Change <- numeric(length(Transect))
EPR <- numeric(length(Transect))
Mean.EPR.Eras <- numeric(length(Transect))
StDev.EPR.Eras <- numeric(length(Transect))
Mean.EPR.Eras.L <- numeric(length(Transect))
Mean.EPR.Eras.U <- numeric(length(Transect))
Min.Date.Acc <- numeric(length(Transect))
Max.Date.Acc <- numeric(length(Transect))
EPR.Error <- numeric(length(Transect))
Number.Dates <- numeric(length(Transect))
Range.Distance <- numeric(length(Transect))
Stdev.Change <- numeric(length(Transect))
Stdev.Eras.Distance <- numeric(length(Transect)) #new added 20130916
Mean.Eras.Distance <- numeric(length(Transect)) #new added 20130916
CoVar.Eras.Distance <- numeric(length(Transect)) #new added 20130916

#locate transects with multiple intersections of dates
Transect.Flag <- numeric(length(Transect))
Uniq.Dates <- numeric(length(Transect))
Trans.Dates <- numeric(length(Transect))

#for linear regression help, see page#25 in SimpleR - Using R for introductory Statistics
LRR.obj <- numeric(length(Transect))
LRR.slope <- numeric(length(Transect))
LRR.intercept <- numeric(length(Transect))
LRR.Rsquared <- numeric(length(Transect))
LRR.SECoef <- numeric(length(Transect))
LRR.SEResi <- numeric(length(Transect))
LRR.Pval <- numeric(length(Transect))
LRR.CI.L <- numeric(length(Transect))
LRR.CI.U <- numeric(length(Transect))

#for weighted linear regression
WLR.obj <- numeric(length(Transect))
WLR.slope <- numeric(length(Transect))
WLR.intercept <- numeric(length(Transect))
WLR.Rsquared <- numeric(length(Transect))
WLR.SECoef <- numeric(length(Transect))
WLR.SEResi <- numeric(length(Transect))
WLR.Pval <- numeric(length(Transect))
WLR.CI.L <- numeric(length(Transect))
WLR.CI.U <- numeric(length(Transect))

#for robust linear regression
RLR.slope <- numeric(length(Transect))
#RLR.intercept <- numeric(length(Transect))
#RLR.Rsquared <- numeric(length(Transect))
#RLR.SECoef <- numeric(length(Transect))
#RLR.SEResi <- numeric(length(Transect))
#RLR.tval <- numeric(length(Transect))
#RLR.CI <- numeric(length(Transect))

#for least median of squares
LMS.slope <- numeric(length(Transect))

#for jack knife
JK.avg <- numeric(length(Transect))
JK.min <- numeric(length(Transect))
JK.max <- numeric(length(Transect))

jackknife.lm.avg<-function (lmobj)
{
        n <- length(resid(lmobj))
        jval <- t(apply(as.matrix(1:n), 1, function(y) coef(update(lmobj, subset = -y))))
        mean(jval[,2])
}

jackknife.lm.min<-function (lmobj)
{
        n <- length(resid(lmobj))
        jval <- t(apply(as.matrix(1:n), 1, function(y) coef(update(lmobj, subset = -y))))
        min(jval[,2])
}

jackknife.lm.max<-function (lmobj)
{
        n <- length(resid(lmobj))
        jval <- t(apply(as.matrix(1:n), 1, function(y) coef(update(lmobj, subset = -y))))
        max(jval[,2])
}


#get attributes for old and new shoreline
Min.Date.Class1 <- numeric(length(Transect))
Max.Date.Class1 <- numeric(length(Transect))

#get location attributes for shorelines
Baseline.Location <-  numeric(length(Transect))
Shoreline.Location <-  numeric(length(Transect))

#get azimuth for transects
Transect.Azimuth <- numeric(length(Transect))


#get inner and outer most coordinates of shoreline's transect points (may need to be changed to get "max" for onshore baseline instead of getting "min" distance values (change for onshore at Date Distances)

Transect.Outer.Xcoord <- numeric(length(Transect))
Transect.Outer.Ycoord <- numeric(length(Transect))
Min.Dist.Position <- numeric(length(Transect))

Transect.Inner.Xcoord <- numeric(length(Transect))
Transect.Inner.Ycoord <- numeric(length(Transect))
Max.Dist.Position <- numeric(length(Transect))

# get min and max date coordinates for shoreline forecasting
Min.Date.Xcoord <- numeric(length(Transect))
Min.Date.Ycoord <- numeric(length(Transect))
Max.Date.Xcoord <- numeric(length(Transect))
Max.Date.Ycoord <- numeric(length(Transect))

#get true envelope
Transect.OEnv.Xcoord <- numeric(length(Transect))
Transect.OEnv.Ycoord <- numeric(length(Transect))
Transect.IEnv.Xcoord <- numeric(length(Transect))
Transect.IEnv.Ycoord <- numeric(length(Transect))


 #had to place an absolute reference to WorkTable1 for the EPR rate and others dealing with the DISTANCE field for R version 2.9 to fix a bug on 07-07-2009


  for (i in 1:length(Transect)) {


Baseline.Offshore[i] <- max((WorkTable1df$OFFSHORE)[WorkTable1df$TRANSECT == Transect[i]], na.rm=T)

Transect.Spacing[i]  <- max((WorkTable1df$TRANSPACE)[WorkTable1df$TRANSECT == Transect[i]])

Transect.Distance[i]  <- max((WorkTable1df$TRANDIST)[WorkTable1df$TRANSECT == Transect[i]])

Transect.StartX[i] <- max((WorkTable1df$STARTX)[WorkTable1df$TRANSECT == Transect[i]])
Transect.StartY[i] <- max((WorkTable1df$STARTY)[WorkTable1df$TRANSECT == Transect[i]])
Transect.EndX[i] <- max((WorkTable1df$ENDX)[WorkTable1df$TRANSECT == Transect[i]])
Transect.EndY[i] <- max((WorkTable1df$ENDY)[WorkTable1df$TRANSECT == Transect[i]])

Transect.Min.Date[i] <- min((WorkTable1df$DATE2)[WorkTable1df$TRANSECT == Transect[i]])

Transect.Max.Date[i] <- max((WorkTable1df$DATE2)[WorkTable1df$TRANSECT == Transect[i]])

Elapsed.Years <-  (Transect.Max.Date - Transect.Min.Date)/time.adj



Min.Date.Position[i] <- which.min((WorkTable1df$DATE2)[WorkTable1df$TRANSECT == Transect[i]])

Max.Date.Position[i] <- which.max((WorkTable1df$DATE2)[WorkTable1df$TRANSECT == Transect[i]])

Min.Date.Dist[i] <- WorkTable1df$DISTANCE[WorkTable1df$TRANSECT == Transect[i]][Min.Date.Position[i]]

Max.Date.Dist[i] <- WorkTable1df$DISTANCE[WorkTable1df$TRANSECT == Transect[i]][Max.Date.Position[i]]

Net.Change <- (Max.Date.Dist - Min.Date.Dist)

EPR <- Net.Change / Elapsed.Years

Mean.EPR.Eras[i]  <-  ifelse(length(Rates.Consec.Eras[Changes.Transects == Transect[i]]) > 0, mean(Rates.Consec.Eras[Changes.Transects == Transect[i]][-1]), 0)

StDev.EPR.Eras[i] <-  ifelse(length(Rates.Consec.Eras[Changes.Transects == Transect[i]]) > 1, sd((Rates.Consec.Eras)[Changes.Transects == Transect[i]][-1],na.rm=TRUE), 0)

Min.Date.Acc[i] <- WorkTable1df$ACCURACY[WorkTable1df$TRANSECT == Transect[i]][Min.Date.Position[i]]

Max.Date.Acc[i] <- WorkTable1df$ACCURACY[WorkTable1df$TRANSECT == Transect[i]][Max.Date.Position[i]]

Min.Date.Class1[i] <- as.character(WorkTable1df$CLASS_1[WorkTable1df$TRANSECT == Transect[i]][Min.Date.Position[i]])

Max.Date.Class1[i] <- as.character(WorkTable1df$CLASS_1[WorkTable1df$TRANSECT == Transect[i]][Max.Date.Position[i]])

Baseline.Location[i] <- as.character(WorkTable1df$BASE_LOC[WorkTable1df$TRANSECT == Transect[i]][Max.Date.Position[i]])

Shoreline.Location[i] <- as.character(WorkTable1df$SHORE_LOC[WorkTable1df$TRANSECT == Transect[i]][Max.Date.Position[i]])

Transect.Azimuth[i] <-  as.character(WorkTable1df$AZIMUTH[WorkTable1df$TRANSECT == Transect[i]][Max.Date.Position[i]])

EPR.Error <- (sqrt(((Min.Date.Acc[i])^2) + ((Max.Date.Acc[i])^2))) / Elapsed.Years


Min.Date <- as.character(format(as.POSIXct(Transect.Min.Date, origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p"))

Max.Date <- as.character(format(as.POSIXct(Transect.Max.Date, origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p"))

Number.Dates[i] <- length(WorkTable1df$DATE2[WorkTable1df$TRANSECT == Transect[i]])

Mean.EPR.Eras.L[i]  <- ifelse(Number.Dates[i] > 2, Mean.EPR.Eras[i] - StDev.EPR.Eras[i],Mean.EPR.Eras[i] - EPR.Error)

Mean.EPR.Eras.U[i]  <- ifelse(Number.Dates[i] > 2, Mean.EPR.Eras[i] + StDev.EPR.Eras[i],Mean.EPR.Eras[i] + EPR.Error)


Transect.Means[i] <- sum(Changes.Distances[Changes.Transects == Transect[i]]) / (Number.Dates[i] - 1)

Uniq.Dates[i] <- length(unique(WorkTable1df$DATE2))

Range.Distance[i] <-  (max((WorkTable1df$DISTANCE)[WorkTable1df$TRANSECT == Transect[i]])) - (min((WorkTable1df$DISTANCE)[WorkTable1df$TRANSECT == Transect[i]]))

LRR.obj <- lm((Changes.Distances)[Changes.Transects == Transect[i]] ~ (Changes.Years)[Changes.Transects == Transect[i]])

LRR.slope[i] <- (coef(LRR.obj)[2])

LRR.intercept[i] <- coef(LRR.obj)[1]

LRR.Rsquared[i] <-  ifelse(Number.Dates[i] > 2,(cor((Changes.Years)[Changes.Transects == Transect[i]], (Changes.Distances)[Changes.Transects == Transect[i]]))^2,0)

LRR.SECoef[i] <- ifelse(Number.Dates[i] > 2, (summary(LRR.obj)$coefficients)[2,2], 0)

LRR.SEResi[i] <- ifelse(Number.Dates[i] > 2, summary(LRR.obj)$sigma, 0)

LRR.Pval[i] <- ifelse(Number.Dates[i] > 2, (summary(LRR.obj)$coefficients)[2,4], 0)

LRR.CI.L[i] <- ifelse(Number.Dates[i] > 2, (confint(LRR.obj, level = ConfInter))[2,1], LRR.slope[i] - EPR.Error)

LRR.CI.U[i] <- ifelse(Number.Dates[i] > 2, (confint(LRR.obj, level = ConfInter))[2,2], LRR.slope[i] + EPR.Error)


WLR.obj <-  lm((Changes.Distances)[Changes.Transects == Transect[i]] ~ (Changes.Years)[Changes.Transects == Transect[i]], weights=(1/((WorkTable1df$ACCURACY)[Changes.Transects == Transect[i]])^2))

WLR.slope[i] <- (coef(WLR.obj)[2])

WLR.intercept[i] <- coef(WLR.obj)[1]

WLR.Rsquared[i] <-  ifelse(Number.Dates[i] > 2, summary(WLR.obj)$r.squared, 0)

WLR.SECoef[i] <- ifelse(Number.Dates[i] > 2, (summary(WLR.obj)$coefficients)[2,2], 0)

WLR.SEResi[i] <- ifelse(Number.Dates[i] > 2, summary(WLR.obj)$sigma, 0)

WLR.Pval[i] <- ifelse(Number.Dates[i] > 2, (summary(WLR.obj)$coefficients)[2,4], 0)

WLR.CI.L[i] <- ifelse(Number.Dates[i] > 2, (confint(WLR.obj, level = ConfInter))[2,1], LRR.slope[i] - EPR.Error)

WLR.CI.U[i] <- ifelse(Number.Dates[i] > 2, (confint(WLR.obj, level = ConfInter))[2,2], LRR.slope[i] + EPR.Error)

#
RLR.slope[i] <- ifelse(Number.Dates[i] > 2 & Basic.Analysis != 1,(coef(rlm((Changes.Distances)[Changes.Transects == Transect[i]] ~ (Changes.Years)[Changes.Transects == Transect[i]], weights=(1/((WorkTable1df$ACCURACY)[Changes.Transects == Transect[i]])^2),method="MM",wt.method="inv.var",psi=psi.bisquare,init="ls"))[2]),0)

#RLR.intercept[i] <- ifelse(Number.Dates[i] > 2 & Basic.Analysis != 1,coef(rlm((Changes.Distances)[Changes.Transects == Transect[i]] ~ (Changes.Years)[Changes.Transects == Transect[i]], weights= 1/((ACCURACY)[Changes.Transects == Transect[i]]) ^2,method="MM",wt.method="inv.var",psi=psi.bisquare,init="ls"))[1],0)

#RLR.Rsquared[i] <-  ifelse(Number.Dates[i] > 2 & Basic.Analysis != 1, summary(rlm((Changes.Distances)[Changes.Transects == Transect[i]] ~ (Changes.Years)[Changes.Transects == Transect[i]], weights=(1/((ACCURACY)[Changes.Transects == Transect[i]])^2),method="MM",wt.method="inv.var",psi=psi.bisquare,init="ls"))$r.squared, 0)

#RLR.SECoef[i] <- ifelse(Number.Dates[i] > 2 & Basic.Analysis != 1, (summary(rlm((Changes.Distances)[Changes.Transects == Transect[i]] ~ (Changes.Years)[Changes.Transects == Transect[i]], weights=(1/((ACCURACY)[Changes.Transects == Transect[i]])^2),method="MM",wt.method="inv.var",psi=psi.bisquare,init="ls"))$coefficients)[2,2], 0)

#RLR.SEResi[i] <- ifelse(Number.Dates[i] > 2 & Basic.Analysis != 1, summary(rlm((Changes.Distances)[Changes.Transects == Transect[i]] ~ (Changes.Years)[Changes.Transects == Transect[i]], weights=(1/((ACCURACY)[Changes.Transects == Transect[i]])^2),method="MM",wt.method="inv.var",psi=psi.bisquare,init="ls"))$sigma, 0)

#RLR.tval[i] <- ifelse(Number.Dates[i] > 2 & Basic.Analysis != 1, (summary(rlm((Changes.Distances)[Changes.Transects == Transect[i]] ~ (Changes.Years)[Changes.Transects == Transect[i]], weights=(1/((ACCURACY)[Changes.Transects == Transect[i]])^2),method="MM",wt.method="inv.var",psi=psi.bisquare,init="ls"))$coefficients)[2,3], 0)

#RLR.CI[i] <- ifelse(Number.Dates[i] > 2 & Basic.Analysis != 1, (confint(rlm((Changes.Distances)[Changes.Transects == Transect[i]] ~ (Changes.Years)[Changes.Transects == Transect[i]], weights=(1/((ACCURACY)[Changes.Transects == Transect[i]])^2),method="MM",wt.method="inv.var",psi=psi.bisquare,init="ls"), level = ConfInter))[2,2], 0)


LMS.slope[i] <- ifelse(Number.Dates[i] > 2 & Basic.Analysis != 1,(coef(lmsreg((Changes.Distances)[Changes.Transects == Transect[i]] ~ (Changes.Years)[Changes.Transects == Transect[i]]))[2]),0)

JK.prep <- lm((Changes.Distances)[Changes.Transects == Transect[i]] ~ (Changes.Years)[Changes.Transects == Transect[i]])

JK.avg[i] <- ifelse(Number.Dates[i] > 2 & Basic.Analysis != 1,jackknife.lm.avg(JK.prep),0)

JK.min[i] <- ifelse(Number.Dates[i] > 2 & Basic.Analysis != 1,jackknife.lm.min(JK.prep),0)

JK.max[i] <- ifelse(Number.Dates[i] > 2 & Basic.Analysis != 1,jackknife.lm.max(JK.prep),0)



Time.Stamp[i] <- as.character(Sys.time())



Trans.Dates[i] <- length((WorkTable1df$DATE)[WorkTable1df$TRANSECT == Transect[i]])

Transect.Flag[i] <- ifelse(Uniq.Dates[i] == Trans.Dates[i],'', 'FLAG')

Stdev.Change[i] <- sd((Changes.Distances)[Changes.Transects == Transect[i]])

Stdev.Eras.Distance[i] <- ifelse(length(Rates.Consec.Eras[Changes.Transects == Transect[i]]) > 1, sd((abs(Rates.Consec.Eras))[Changes.Transects == Transect[i]][-1],na.rm=TRUE), 0) #new added 20130916
Mean.Eras.Distance[i] <-  ifelse(length(Rates.Consec.Eras[Changes.Transects == Transect[i]]) > 1, mean((abs(Rates.Consec.Eras))[Changes.Transects == Transect[i]][-1],na.rm=TRUE), 0)#new added 20130916
CoVar.Eras.Distance[i] <- Stdev.Eras.Distance[i]/Mean.Eras.Distance[i] #new added 20130916

      


#####fixed envelope with raw data.frame(shapedata)  2/15/2013

raw.data <-  data.frame(shapedata) ##added  2/15/2013 to get the true envelope of change including shoreline switchbacks  (see next 6 lines of code)

Min.Dist.Position[i] <- ifelse(Baseline.Offshore[i] == 1, which.max((raw.data$Distance)[raw.data$Transect == Transect[i]]), which.min((raw.data$Distance)[raw.data$Transect == Transect[i]]))  #fixed 2/15/2013

Transect.Outer.Xcoord[i] <- raw.data$coords.x1[raw.data$Transect == Transect[i]][Min.Dist.Position[i]]   #fixed 2/15/2013

Transect.Outer.Ycoord[i] <- raw.data$coords.x2[raw.data$Transect == Transect[i]][Min.Dist.Position[i]]   #fixed 2/15/2013

Max.Dist.Position[i] <- ifelse(Baseline.Offshore[i] == 1, which.min((raw.data$Distance)[raw.data$Transect == Transect[i]]), which.max((raw.data$Distance)[raw.data$Transect == Transect[i]]))   #fixed 2/15/2013

Transect.Inner.Xcoord[i] <- raw.data$coords.x1[raw.data$Transect == Transect[i]][Max.Dist.Position[i]]    #fixed 2/15/2013

Transect.Inner.Ycoord[i] <- raw.data$coords.x2[raw.data$Transect == Transect[i]][Max.Dist.Position[i]]    #fixed 2/15/2013




Min.Date.Xcoord[i] <- WorkTable1df$X_COORD[WorkTable1df$TRANSECT == Transect[i]][Min.Date.Position[i]]

Min.Date.Ycoord[i] <- WorkTable1df$Y_COORD[WorkTable1df$TRANSECT == Transect[i]][Min.Date.Position[i]]

Max.Date.Xcoord[i] <- WorkTable1df$X_COORD[WorkTable1df$TRANSECT == Transect[i]][Max.Date.Position[i]]

Max.Date.Ycoord[i] <- WorkTable1df$Y_COORD[WorkTable1df$TRANSECT == Transect[i]][Max.Date.Position[i]]


#########get true envelope transects
Min.Dist.Position[i] <- ifelse(Baseline.Offshore[i] == 1, which.max((WorkTable1df$DISTANCE)[WorkTable1df$TRANSECT == Transect[i]]), which.min((WorkTable1df$DISTANCE)[WorkTable1df$TRANSECT == Transect[i]]))

Transect.OEnv.Xcoord[i] <- WorkTable1df$X_COORD[WorkTable1df$TRANSECT == Transect[i]][Min.Dist.Position[i]]

Transect.OEnv.Ycoord[i] <- WorkTable1df$Y_COORD[WorkTable1df$TRANSECT == Transect[i]][Min.Dist.Position[i]]

Max.Dist.Position[i] <- ifelse(Baseline.Offshore[i] == 1, which.min((WorkTable1df$DISTANCE)[WorkTable1df$TRANSECT == Transect[i]]), which.max((WorkTable1df$DISTANCE)[WorkTable1df$TRANSECT == Transect[i]]))

Transect.IEnv.Xcoord[i] <- WorkTable1df$X_COORD[WorkTable1df$TRANSECT == Transect[i]][Max.Dist.Position[i]]
         
Transect.IEnv.Ycoord[i] <- WorkTable1df$Y_COORD[WorkTable1df$TRANSECT == Transect[i]][Max.Dist.Position[i]]








#status update: add progress bar, estimate percent completion and map

Pcnt.Complete <-  round(((Transect[i] - Start.Transect + 1)/ length(Transect)) * 100, 0)

Pcnt.Complete2 <- paste(Pcnt.Complete," ","%",sep="")



    info <- sprintf("%d%% done", Pcnt.Complete)
    setTkProgressBar(pb, i, sprintf("AMBUR: Transect analysis (%s)", info), info)







plot(Max.Date.Xcoord[1:i], Max.Date.Ycoord[1:i], type="n", lwd= 0, col= "gray" , asp=1, xlab=paste("X ","(",Map.Units,")",sep=""),ylab=paste("Y ","(",Map.Units,")",sep=""),main=c("Transects processed:", i , "out of", length(Transect)), sub=as.character(Pcnt.Complete2),xlim=c(min(mydata$X_COORD),max(mydata$X_COORD)),ylim=c(min(mydata$Y_COORD),max(mydata$Y_COORD)))


points(mydata$X_COORD[mydata$DATE == min(mydata$DATE)], mydata$Y_COORD[mydata$DATE == min(mydata$DATE)], type="p", pch=20, lwd= 1, col= "gray")

points(Max.Date.Xcoord[1:i][Net.Change > 0], Max.Date.Ycoord[1:i][1:i][Net.Change > 0], type="p", col= "blue")

points(Max.Date.Xcoord[1:i][Net.Change <= 0], Max.Date.Ycoord[1:i][1:i][Net.Change <= 0], type="p", col= "red")


par("usr")
textalign <- par("usr")[1] - ((par("usr")[2] - par("usr")[1]) * 1.464759)

mtext("-constructing master data table...", side=3, line= 0, adj= 0, cex= 0.75, at=textalign)
mtext("AMBUR: Shoreline Change Analysis", side=3, line= 3, adj= 0, cex= 0.75, at=textalign)
mtext("-starting analyses...", side=3, line= 2, adj= 0, cex= 0.75, at=textalign)
mtext("-building data tables...", side=3, line= 1, adj= 0, cex= 0.75, at=textalign)

}

#determine morphologic change

Att.Change <- ifelse(as.character(Min.Date.Class1) != as.character(Max.Date.Class1), "CHANGE", "NO CHANGE")

#replace NA values with 0 in StDev of EPRs for Consecutive Eras

StDev.EPR.Eras <- ifelse(Mean.EPR.Eras == EPR, 0, StDev.EPR.Eras * 1)

#join results to final table
FinalTable <- cbind(Transect, Baseline.Offshore, Transect.Spacing, Transect.Distance, Transect.Flag, Transect.StartX, Transect.StartY, Transect.EndX, Transect.EndY, Transect.Inner.Xcoord, Transect.Inner.Ycoord,Transect.Outer.Xcoord, Transect.Outer.Ycoord, Min.Date.Xcoord, Min.Date.Ycoord, Max.Date.Xcoord, Max.Date.Ycoord, Number.Dates, Min.Date, Max.Date, Elapsed.Years, Transect.Means, Range.Distance, Stdev.Change, Min.Date.Position, Max.Date.Position, Min.Date.Dist, Max.Date.Dist, Min.Date.Acc, Max.Date.Acc, Net.Change, EPR, EPR.Error, Mean.EPR.Eras, StDev.EPR.Eras, Mean.EPR.Eras.L, Mean.EPR.Eras.U, LRR.slope, LRR.Rsquared, LRR.intercept, LRR.SECoef, LRR.SEResi, LRR.Pval, LRR.CI.L, LRR.CI.U, WLR.slope, WLR.Rsquared, WLR.intercept, WLR.SECoef, WLR.SEResi, WLR.Pval, WLR.CI.L, WLR.CI.U, RLR.slope, LMS.slope, JK.avg, JK.min, JK.max, Min.Date.Class1, Max.Date.Class1, Att.Change, Baseline.Location, Shoreline.Location, Transect.Azimuth, Time.Stamp,Stdev.Eras.Distance,Mean.Eras.Distance,CoVar.Eras.Distance,Transect.IEnv.Xcoord,Transect.IEnv.Ycoord,Transect.OEnv.Xcoord,Transect.OEnv.Ycoord   )


#set alternative field names for GIS compatible files

Transect <- Transect
Base_Off <- Baseline.Offshore
Tran_Spac <- Transect.Spacing
Tran_Dist <- Transect.Distance
Tran_Flag <- Transect.Flag
Start_X <- Transect.StartX
Start_Y <- Transect.StartY
End_X <- Transect.EndX
End_Y <- Transect.EndY
Inner_X <- Transect.Inner.Xcoord
Inner_Y <- Transect.Inner.Ycoord
Outer_X <- Transect.Outer.Xcoord
Outer_Y <- Transect.Outer.Ycoord
Min_DateX <- Min.Date.Xcoord
Min_DateY <- Min.Date.Ycoord
Max_DateX <- Max.Date.Xcoord
Max_DateY <- Max.Date.Ycoord
Num_Dates <- Number.Dates
Min_Date <- Min.Date
Max_Date <- Max.Date
Elp_Years <- Elapsed.Years
Tran_Mean <- Transect.Means
Range_Dst <- Range.Distance
Stdev_Chg <- Stdev.Change
MinDPos <- Min.Date.Position
MaxDPos <- Max.Date.Position
MinDDist <- Min.Date.Dist
MaxDDist <- Max.Date.Dist
MinDAcc <- Min.Date.Acc
MaxDAcc <- Max.Date.Acc
Net_Chng <- Net.Change
EPR <- EPR
EPR_Error  <- EPR.Error
EPR_MnEra  <- Mean.EPR.Eras
EPR_SDEra  <- StDev.EPR.Eras
EPR_Er_L <- Mean.EPR.Eras.L
EPR_Er_U <- Mean.EPR.Eras.U
LRR <- LRR.slope
LRR_Rsqr <- LRR.Rsquared
LRR_int <- LRR.intercept
LRR_SEcoe <- LRR.SECoef
LRR_SEres <- LRR.SEResi
LRR_Pval <- LRR.Pval
LRR_CI_L <- LRR.CI.L
LRR_CI_U <- LRR.CI.U
WLR <- WLR.slope
WLR_Rsqr <- WLR.Rsquared
WLR_int <- WLR.intercept
WLR_SEcoe <- WLR.SECoef
WLR_SEres <- WLR.SEResi
WLR_Pval  <- WLR.Pval
WLR_CI_L <- WLR.CI.L
WLR_CI_U <- WLR.CI.U
RLR <- RLR.slope
#RLR_Rsqr <- RLR.Rsquared
#RLR_int <- RLR.intercept
#RLR_SEcoe <- RLR.SECoef
#RLR_SEres <- RLR.SEResi
#RLR_tval  <- RLR.tval
#RLR_CI <- RLR.CI
LMS <- LMS.slope
JK_avg <- JK.avg
JK_min <- JK.min
JK_max <- JK.max
MinClass1 <- Min.Date.Class1
MaxClass1 <- Max.Date.Class1
Attr_Chng <- Att.Change
Base_Loc <- Baseline.Location
Shore_Loc <- Shoreline.Location
T_azimuth <- Transect.Azimuth
Time_Stmp <- Time.Stamp
SDErasDst <- Stdev.Eras.Distance
MnErasDst <- Mean.Eras.Distance
CVErasDst <- CoVar.Eras.Distance
OuterEnvX <- Transect.OEnv.Xcoord
OuterEnvY <- Transect.OEnv.Ycoord
InnerEnvX <- Transect.IEnv.Xcoord
InnerEnvY <- Transect.IEnv.Ycoord




FinalGISTable <- cbind(Transect, Base_Off, Tran_Spac, Tran_Dist, Tran_Flag, Start_X, Start_Y, End_X, End_Y, Inner_X, Inner_Y, Outer_X, Outer_Y, Min_DateX, Min_DateY, Max_DateX, Max_DateY, Num_Dates, Min_Date, Max_Date, Elp_Years, Tran_Mean, Range_Dst, Stdev_Chg, MinDPos, MaxDPos, MinDDist, MaxDDist, MinDAcc, MaxDAcc, Net_Chng, EPR, EPR_Error, EPR_MnEra, EPR_SDEra, EPR_Er_L, EPR_Er_U, LRR, LRR_Rsqr, LRR_int, LRR_SEcoe, LRR_SEres, LRR_Pval, LRR_CI_L, LRR_CI_U, WLR, WLR_Rsqr, WLR_int, WLR_SEcoe, WLR_SEres, WLR_Pval, WLR_CI_L, WLR_CI_U, RLR, LMS, JK_avg, JK_min, JK_max,MinClass1, MaxClass1, Attr_Chng, Base_Loc, Shore_Loc, T_azimuth, Time_Stmp,SDErasDst,MnErasDst,CVErasDst,OuterEnvX,OuterEnvY,InnerEnvX,InnerEnvY)


#status checkpoint
par("usr")
textalign <- par("usr")[1] - ((par("usr")[2] - par("usr")[1]) * 1.464759)

mtext("-constructing master data table...", side=3, line= 0, adj= 0, cex= 0.75, at=textalign)
mtext("-writing final data tables...", side=3, line= -1, adj= 0, cex= 0.75, at=textalign)

close(pb)


#write the final table to a csv file
write.table(FinalTable, file = "results_stats.csv", sep = ",", row.names = FALSE)

#write the final table to a R csv file for plotting
#1#write.table(FinalTable, file = "results_stats_plotting.csv", sep = ",", row.names = TRUE)

#write the final table to a csv file compatible with ArcGIS or Databases

FinalTable2 <- FinalTable

colnames(FinalTable2) <- gsub("[.]", "_", colnames(FinalTable2))

#1#write.table(FinalTable2, file = "GIS_stats_table_long.csv", sep = ",", row.names = FALSE, quote = FALSE)

write.table(FinalGISTable, file = "GIS_stats_table_short.csv", quote = FALSE, sep = ",", row.names = FALSE)

#merge FinalTable with finalchngs to create the supertable
finalchngs2 <- finalchngs[,!duplicated(colnames(finalchngs))]   #####added to remove duplicates  11/19/2012
Super.Table <- merge(FinalTable, finalchngs2, by.x = "Transect", by.y = "TRANSECT", sort= FALSE, all= TRUE)   #####added to remove duplicates  11/19/2012

write.table(Super.Table, file = "super_table.csv", quote = FALSE, sep = ",", row.names = FALSE)

#write the table for debugging of dates and distances

DebugTable <- WorkTable1a

#1# write.table(DebugTable, file = "debugging.csv", quote = FALSE, sep = ",", row.names = FALSE)

DebugTable2 <- WorkTable1

write.table(DebugTable2, file = "debugging2.csv", quote = FALSE, sep = ",", row.names = FALSE)






#status checkpoint
mtext("-plotting and saving graphics...", side=3, line= -2, adj= 0, cex= 0.75,at=textalign)
mtext("-finishing up...", side=3, line= -3, adj= 0, cex= 0.75,at=textalign)

#################################################
#Plot of histograms for EPR,LRR,WLR

pdf("PDF/histograms_1.pdf",width = 6.5, height = 6.5, bg="white", paper= "letter" )
par(mfrow=(c(3,3)))
par(pty= "m")

hist(EPR, main= "", xlab= "", ylab="", breaks= 100)
hist(LRR.slope, main= "All Transects", xlab= "", ylab="", breaks= 100)
hist(WLR.slope, main= "", xlab="", ylab="", breaks= 100)

hist(ifelse(EPR <0, EPR, 0), main= "", xlab= "", ylab="Frequency", breaks= 100)
hist(ifelse(LRR.slope <0, LRR.slope, 0), main= "Erosion Transects", xlab= "", ylab="", breaks= 100)
hist(ifelse(WLR.slope <0, WLR.slope, 0), main= "", xlab= "", ylab="", breaks= 100)


hist(ifelse(EPR >0, EPR, 0), main= "", xlab= "EPR", ylab="", breaks= 100)
hist(ifelse(LRR.slope >0, LRR.slope, 0), main= "Accretion Transects", xlab= "LRR.slope", ylab="", breaks= 100)
hist(ifelse(WLR.slope >0, WLR.slope, 0), main= "", xlab= "WLR.slope", ylab="", breaks= 100)

#turn the device off

dev.off()


#results - 6 plots net,stdev,envelope,epr,lrr,wlr

pdf("PDF/results_6_plots.pdf", width = 11, height = 8.5, bg="white")
par(mfrow=(c(3,2)))
par(pty= "m")


#Plot 1  Net Change
plot(Transect, Net.Change, type="l", lwd= 0, col= "white" , las= 1, cex.axis= 1, cex.lab= 1, xlab=expression(paste("Transect")),ylab=paste("Net Change ","(",Map.Units,")",sep=""))

points(Transect[Net.Change > 0], Net.Change[Net.Change > 0], type="h", col= "blue")

points(Transect[Net.Change <= 0], Net.Change[Net.Change <= 0], type="h", col= "red")


#Plot 2  EPR
plot(Transect, EPR, type="l", lwd= 0, col= "white" , las= 1, cex.axis= 1, cex.lab= 1, xlab=expression(paste("Transect")),ylab=paste("EPR ","(",Map.Units,"/",Time.Units,")",sep=""))

points(Transect[EPR > 0], EPR[EPR > 0], type="h", col= "blue")

points(Transect[EPR <= 0], EPR[EPR <= 0], type="h", col= "red")




#Plot 3  Standard deviation of change
plot(Transect, Stdev.Change, type="l", lwd= 0, col= "white" , las= 1, cex.axis= 1, cex.lab= 1, xlab=expression(paste("Transect")),ylab=paste("St. Dev of Change ","(",Map.Units,")",sep=""))

points(Transect[Stdev.Change > 0], Stdev.Change[Stdev.Change > 0], type="h", col= "orange")

points(Transect[Stdev.Change <= 0], Stdev.Change[Stdev.Change <= 0], type="h", col= "orange")

#Plot 4
plot(Transect, LRR.slope, type="l", lwd= 0, col= "white" , las= 1, cex.axis= 1, cex.lab= 1, xlab=expression(paste("Transect")),ylab=paste("LRR ","(",Map.Units,"/",Time.Units,")",sep=""))

points(Transect[LRR.slope > 0], LRR.slope[LRR.slope > 0], type="h", col= "blue")

points(Transect[LRR.slope <= 0], LRR.slope[LRR.slope <= 0], type="h", col= "red")



#Plot 5  Envelope of change
plot(Transect, Range.Distance, type="l", lwd= 0, col= "white" , las= 1, cex.axis= 1, cex.lab= 1, xlab=expression(paste("Transect")),ylab=paste("Envelope of Change ","(",Map.Units,")",sep=""))

points(Transect[Range.Distance > 0], Range.Distance[Range.Distance > 0], type="h", col= "orange")

points(Transect[Range.Distance <= 0], Range.Distance[Range.Distance <= 0], type="h", col= "orange")


#Plot 6 WLR
plot(Transect, WLR.slope, type="l", lwd= 0, col= "white" , las= 1, cex.axis= 1, cex.lab= 1, xlab=expression(paste("Transect")),ylab=paste("WLR ","(",Map.Units,"/",Time.Units,")",sep=""))

points(Transect[WLR.slope > 0], WLR.slope[WLR.slope > 0], type="h", col= "blue")

points(Transect[WLR.slope <= 0], WLR.slope[WLR.slope <= 0], type="h", col= "red")

dev.off()







#map detailed shoreline plot of Net Change trend
pdf("PDF/map_net_change.pdf", width = 10, height = 10, bg="white")

plot(Transect.Outer.Xcoord, Transect.Outer.Ycoord, type="l", lwd= 0, col= "gray" , asp=1, xlab=paste("X ","(",Map.Units,")",sep=""),ylab=paste("Y ","(",Map.Units,")",sep=""))


points(Transect.Outer.Xcoord[Net.Change > 0], Transect.Outer.Ycoord[Net.Change > 0], type="p", col= "blue")

points(Transect.Outer.Xcoord[Net.Change <= 0], Transect.Outer.Ycoord[Net.Change <= 0], type="p", col= "red")



text.interval1 <- pretty(min(Transect):(max(Transect)), n=6)

text.interval2 <- which(Transect %in% pretty(Transect)  > 0)

text(Transect.Outer.Xcoord[text.interval2], Transect.Outer.Ycoord[text.interval2], Transect[text.interval2], cex= 1, pos= 3, offset= 0.3, font= 2)

points(Transect.Outer.Xcoord[text.interval2], Transect.Outer.Ycoord[text.interval2], Transect[text.interval2], type="p", col= "black", pch= 13)

myLegendText <- c("Erosion", "Accretion")

myColors <- c("red", "blue")

legend("bottomright", legend= myLegendText, pch= 16, pt.cex= 1, col= myColors, bty= "n", cex= 0.7)



dev.off()





#summary report (multiplot with different size windows)
pdf("PDF/summary_report.pdf",width = 6.5, height = 9, bg="white", paper = "letter")

nf <- layout(matrix(c(4,4,2,2,4,4,3,3,1,1,5,5,1,1,6,6), 4, 4, byrow=TRUE), respect=FALSE)

#plot 1 shoreline trend
plot(Transect.Outer.Xcoord, Transect.Outer.Ycoord, type="l", lwd= 0, col= "gray" , asp=1, xlab=paste("X ","(",Map.Units,")",sep=""),ylab=paste("Y ","(",Map.Units,")",sep=""))


points(Transect.Outer.Xcoord[Net.Change > 0], Transect.Outer.Ycoord[Net.Change > 0], type="p", col= "blue", pch=20)

points(Transect.Outer.Xcoord[Net.Change <= 0], Transect.Outer.Ycoord[Net.Change <= 0], type="p", col= "red", pch=20)


text.interval1 <- pretty(min(Transect):(max(Transect)), n=6)

text.interval2 <- which(Transect %in% pretty(Transect)  > 0)

text(Transect.Outer.Xcoord[text.interval2], Transect.Outer.Ycoord[text.interval2], Transect[text.interval2], cex= 1, pos= 3, offset= 0.3, font= 2)

points(Transect.Outer.Xcoord[text.interval2], Transect.Outer.Ycoord[text.interval2], Transect[text.interval2], type="p", col= "black", pch= 13)




myLegendText <- c("Erosion", "Accretion")

myColors <- c("red", "blue")

legend("bottomright", legend= myLegendText, pch= 16, pt.cex= 1, col= myColors, bty= "n", cex= 0.7)

#Plot 2  Net Change
plot(Transect, Net.Change, type="l", lwd= 0, col= "white" , las= 1, cex.axis= 0.7, cex.lab= 0.7, xlab=expression(paste("Transect")),ylab=paste("Net Change ","(",Map.Units,")",sep=""))


points(Transect[Net.Change > 0], Net.Change[Net.Change > 0], type="h", col= "blue")

points(Transect[Net.Change <= 0], Net.Change[Net.Change <= 0], type="h", col= "red")


#Plot 3  EPR
plot(Transect, EPR, type="l", lwd= 0, col= "white" , las= 1, cex.axis= 0.7, cex.lab= 0.7,xlab=expression(paste("Transect")),ylab=paste("EPR ","(",Map.Units,"/",Time.Units,")",sep=""))

points(Transect[EPR > 0], EPR[EPR > 0], type="h", col= "blue")

points(Transect[EPR <= 0], EPR[EPR <= 0], type="h", col= "red")

#Plot 4  Summary info
dummyX <- 0
dummyY <- 0
plot(dummyX, dummyY, type="n", lwd= 0, col= "white" , xlab="",ylab="",main="", bty="n", axes=FALSE)

mtext("AMBUR SUMMARY REPORT", side=3, line= 3, adj= 0.5, cex= 0.5)

mtext("AMBUR v.1.1.6 (20130916)", side=3, line= 2, adj= 0.5, cex= 0.5)


mtext("Oldest Date:", side=3, line= 1, adj= 0, cex= 0.5)
mtext(as.character(min(WorkTable1df$DATE)), side=3, line= 1, adj= 1, cex= 0.5)


mtext("Youngest Date:", side=3, line= 0, adj= 0, cex= 0.5)
mtext(as.character(max(WorkTable1df$DATE)), side=3, line= 0, adj= 1, cex= 0.5)


mtext("Number of transects:", side=3, line= -1, adj= 0, cex= 0.5)
mtext(length(unique(Transect)), side=3, line= -1, adj= 1, cex= 0.5)

mtext("Erosion transects:", side=3, line= -2, adj= 0, cex= 0.5)
mtext(length(Transect[EPR < 0 & EPR != "NaN"]), side=3, line= -2, adj= 1, cex= 0.5)

mtext("Accretion transects:", side=3, line= -3, adj= 0, cex= 0.5)
mtext(length(Transect[EPR > 0 & EPR != "NaN"]), side=3, line= -3, adj= 1, cex= 0.5)

mtext("Flagged transects:", side=3, line= -4, adj= 0, cex= 0.5)
mtext(length(Transect[Transect.Flag == "FLAG"]), side=3, line= -4, adj= 1, cex= 0.5)

mtext(paste("Maximum Erosion ","(",Map.Units,"):",sep=""), side=3, line= -6, adj= 0, cex= 0.5)
mtext(round(min(Net.Change[Net.Change < 0 & Net.Change != "NaN"]), digits=2), side=3, line= -6, adj= 1, cex= 0.5)

mtext(paste("Mean Erosion ","(",Map.Units,"):",sep=""), side=3, line= -7, adj= 0, cex= 0.5)
mtext(round(mean(Net.Change[Net.Change < 0 & Net.Change != "NaN"]), digits=2), side=3, line= -7, adj= 1, cex= 0.5)

mtext(paste("Maximum Erosion Rate ","(",Map.Units,"/",Time.Units,"):",sep=""), side=3, line= -8, adj= 0, cex= 0.5)
mtext(round(min(EPR[EPR < 0 & EPR != "NaN"]), digits=2), side=3, line= -8, adj= 1, cex= 0.5)

mtext(paste("Mean Erosion Rate ","(",Map.Units,"/",Time.Units,"):",sep=""), side=3, line= -9, adj= 0, cex= 0.5)
mtext(round(mean(EPR[EPR < 0 & EPR != "NaN"]), digits=2), side=3, line= -9, adj= 1, cex= 0.5)

mtext(paste("Maximum Accretion ","(",Map.Units,"):",sep=""), side=3, line= -11, adj= 0, cex= 0.5)
mtext(round(max(Net.Change[Net.Change > 0 & Net.Change != "NaN"]), digits=2), side=3, line= -11, adj= 1, cex= 0.5)

mtext(paste("Mean Accretion ","(",Map.Units,"):",sep=""), side=3, line= -12, adj= 0, cex= 0.5)
mtext(round(mean(Net.Change[Net.Change > 0 & Net.Change != "NaN"]), digits=2), side=3, line= -12, adj= 1, cex= 0.5)

mtext(paste("Maximum Accretion Rate ","(",Map.Units,"/",Time.Units,"):",sep=""), side=3, line= -13, adj= 0, cex= 0.5)
mtext(round(max(EPR[EPR > 0 & EPR != "NaN"]), digits=2), side=3, line= -13, adj= 1, cex= 0.5)

mtext(paste("Mean Accretion Rate ","(",Map.Units,"/",Time.Units,"):",sep=""), side=3, line= -14, adj= 0, cex= 0.5)
mtext(round(mean(EPR[EPR > 0 & EPR != "NaN"]), digits=2), side=3, line= -14, adj= 1, cex= 0.5)

mtext(paste("Overall Mean Change ","(",Map.Units,"):",sep=""), side=3, line= -16, adj= 0, cex= 0.5)
mtext(round(mean(Net.Change[Net.Change != "NaN"]), digits=2), side=3, line= -16, adj= 1, cex= 0.5)

mtext(paste("Overall Mean Rate ","(",Map.Units,"/",Time.Units,"):",sep=""), side=3, line= -17, adj= 0, cex= 0.5)
mtext(round(mean(EPR[EPR != "NaN"]), digits=2), side=3, line= -17, adj= 1, cex= 0.5)


mtext("Analysis Date:", side=3, line= -19, adj= 0, cex= 0.5)
mtext(time.stamp1, side=3, line= -19, adj= 1, cex= 0.5)

mtext("Analyst:", side=3, line= -20, adj= 0, cex= 0.5)
mtext(Sys.info()[7], side=3, line= -20, adj= 1, cex= 0.5)

mtext("File:", side=3, line= -21, adj= 0, cex= 0.5)
mtext(path, side=3, line= -22, adj= 1, cex= 0.3)

mtext("Results:", side=3, line= -23, adj= 0, cex= 0.5)
mtext(results.path, side=3, line= -24, adj= 1, cex= 0.3)


mtext("Transect Spacing:", side=3, line= -25, adj= 0, cex= 0.5)
mtext(max(Transect.Spacing), side=3, line= -25, adj= 1, cex= 0.5)

mtext("Confidence Interval:", side=3, line= -26, adj= 0, cex= 0.5)
mtext(paste(userinput2,"%",sep=""), side=3, line= -26, adj= 1, cex= 0.5)

#Plot 5 Average Rates of Consecutive Eras
plot(Transect, Mean.EPR.Eras, type="l", lwd= 0, col= "white" , las= 1, cex.axis= 0.7, cex.lab= 0.7,xlab=expression(paste("Transect")),ylab=paste("Mean EPR Eras ","(",Map.Units,"/",Time.Units,")",sep=""))

points(Transect[Mean.EPR.Eras > 0], Mean.EPR.Eras[Mean.EPR.Eras > 0], type="h", col= "blue")

points(Transect[Mean.EPR.Eras <= 0], Mean.EPR.Eras[Mean.EPR.Eras <= 0], type="h", col= "red")


#Plot 6  stdev of EPR eras
plot(Transect, StDev.EPR.Eras, type="l", lwd= 0, col= "white" , las= 1, cex.axis= 0.7, cex.lab= 0.7,xlab=expression(paste("Transect")),ylab=paste("St.Dev. EPR Eras","(",Map.Units,")",sep=""))

points(Transect[StDev.EPR.Eras > 0], StDev.EPR.Eras[StDev.EPR.Eras > 0], type="h", col= "orange")

points(Transect[StDev.EPR.Eras <= 0], StDev.EPR.Eras[StDev.EPR.Eras <= 0], type="h", col= "orange")

#layout.show(nf)


dev.off()


#True Cumulative shoreline evolution graph to show discontinuity of shoreline data
pdf("PDF/graph_cumulative_change.pdf", width = 9, height = 6, bg="white")

x <- Changes.Transects
y <- Changes.Years
z <- Changes.Distances
my.renorm <- function (z) {
  z <- abs(z)
  z <- z*5/max(z)
  z
}
z <- my.renorm(z)

Changes.Colors <- ifelse(Changes.Distances < 0, "red", "blue")

plot(x, y, cex = z, xlab = "Transect", ylab = "Elapsed Time (years)", col= Changes.Colors, las= 1, cex.axis= 0.7, cex.lab= 0.7, bty= "l")

points(Changes.Transects[Changes.Distances == 0], Changes.Years[Changes.Distances == 0], col= "green", type= "p", pch= 39)

dev.off()




#present dates easier viewing (non time scaled)
pdf("PDF/graph_dates_present.pdf", width = 8, height = 4, bg="white")

par("usr")[1]
par(mar=c(4,7,1,1))

Changes.DatesSecs <- as.numeric(as.POSIXlt(Changes.Dates, format="%m/%d/%Y %I:%M:%S %p") - as.POSIXlt("01/01/1970 12:0:0 AM", format="%m/%d/%Y %I:%M:%S %p"),units="secs")

x <- Changes.Transects
y <- Changes.DatesSecs
z <- Changes.Distances
my.renorm <- function (z) {
  z <- abs(z)
  z <- z*5/max(z)
  z
}
z <- my.renorm(z)

ha <- sort(unique(y))
ha2 <- seq(from=1, to=length(ha), by=1)
ha3 <- cbind(ha, ha2)
ha4 <- cbind(x, y, z)
ha5 <- merge(ha3, ha4, by.x= "ha", by.y= "y")

plot(ha5$x, ha5$ha2, xlab = "Transect", ylab = "", col= "black", las= 1, cex.axis= 0.7, cex.lab= 0.7, bty= "l", yaxt= "n", pch=39)

points(ha5$x[ha5$z == 0], ha5$ha2[ha5$z == 0], col= "green", type= "p", pch= 39)

axis(2, at= ha5$ha2, labels= as.character(format(as.POSIXct(ha5$ha , origin="1970-01-01"),"%m/%d/%Y %I:%M:%S %p")), cex.axis= 0.5, las= 1)


dev.off()





#experimental shoreline plot of RANGE trend
pdf("PDF/map_range_of_change.pdf", width = 10, height = 10, bg="white")

plot(Transect.Outer.Xcoord, Transect.Outer.Ycoord, type="l", lwd= 0, col= "gray" , asp=1, xlab=paste("X ","(",Map.Units,")",sep=""),ylab=paste("Y ","(",Map.Units,")",sep=""))

NetChange.Colors <- ifelse(Net.Change < 0, "red", "blue")

segments(Transect.Outer.Xcoord, Transect.Outer.Ycoord, Transect.Inner.Xcoord, Transect.Inner.Ycoord, col= NetChange.Colors, lwd=1)

points(Transect.Inner.Xcoord[Net.Change > 0], Transect.Inner.Ycoord[Net.Change > 0], type="p", col= "blue", pch= 20, cex=0.5)

points(Transect.Inner.Xcoord[Net.Change <= 0], Transect.Inner.Ycoord[Net.Change <= 0], type="p", col= "red", pch= 20, cex=0.5)

points(Transect.Outer.Xcoord[Net.Change > 0], Transect.Outer.Ycoord[Net.Change > 0], type="p", col= "blue", pch= 20, cex=0.5)

points(Transect.Outer.Xcoord[Net.Change <= 0], Transect.Outer.Ycoord[Net.Change <= 0], type="p", col= "red", pch= 20, cex=0.5)

points(Transect.Outer.Xcoord[Transect.Flag == "FLAG"], Transect.Outer.Ycoord[Transect.Flag == "FLAG"], col= "green", pch= 20, cex= 0.5)

text.interval1 <- pretty(min(Transect):(max(Transect)), n=6)

text.interval2 <- which(Transect %in% pretty(Transect)  > 0)

text(Transect.Outer.Xcoord[text.interval2], Transect.Outer.Ycoord[text.interval2], Transect[text.interval2], cex= 1, pos= 3, offset= 0.3, font= 2)

points(Transect.Outer.Xcoord[text.interval2], Transect.Outer.Ycoord[text.interval2], Transect[text.interval2], type="p", col= "black", pch= 13)


myLegendText <- c("Erosion", "Accretion")

myColors <- c("red", "blue")

legend("bottomright", legend= myLegendText, pch= 16, pt.cex= 1, col= myColors, bty= "n", cex= 0.7)

dev.off()



# plot EPR rates & error bars

EPR.Error.L <- EPR - EPR.Error
EPR.Error.U <- EPR + EPR.Error

plotrng1 <- min(EPR.Error.L[EPR != "NaN"])
plotrng2 <- max(EPR.Error.U[EPR != "NaN"])


# plot and add the error bars:

pdf("PDF/graph_EPR.pdf", bg="white")

plot(Transect, EPR, type="p", lwd= 0, col= "black" , pch=0, cex=0, las= 1, cex.axis= 1, cex.lab= 1, xlab=expression(paste("Transect")),ylab=paste("EPR ","(",Map.Units,"/",Time.Units,")",sep=""), ylim = range(c(plotrng1,plotrng2)))

abline(h=0)

points(Transect[EPR > 0], EPR[EPR > 0], type="p", col= "blue", pch= 1, cex=0.3)

points(Transect[EPR <= 0], EPR[EPR <= 0], type="p", col= "red", pch= 1, cex=0.3)

arrows(Transect,EPR.Error.L,Transect,EPR.Error.U,length = 0.05, angle = 90, code = 3,col = "gray")

arrows(Transect,EPR.Error.L,Transect,EPR.Error.U,length = 0.05, angle = 90, code = 3,col= ifelse(EPR.Error.L < 0 & EPR.Error.U < 0, "red","transparent"))

arrows(Transect,EPR.Error.L,Transect,EPR.Error.U,length = 0.05, angle = 90, code = 3,col= ifelse(EPR.Error.L > 0 & EPR.Error.U > 0, "blue","transparent"))


dev.off()



# plot the LRR and confidence limits (CI error: lower limit & the upper limit)



plotrng1 <- min(LRR.CI.L[EPR != "NaN"])
plotrng2 <- max(LRR.CI.U[EPR != "NaN"])

if (length(DateList3) == 2) {
plotrng1 <- min(LRR.slope[EPR != "NaN"])
plotrng2 <- max(LRR.slope[EPR != "NaN"])
}

#plot and add the error bars
pdf("PDF/graph_LRR.pdf", bg="white")

plot(Transect, LRR.slope, type="p", lwd= 0, col= "black" , pch=0, cex=0, las= 1, cex.axis= 1, cex.lab= 1, xlab=expression(paste("Transect")),ylab=paste("LRR ","(",Map.Units,"/",Time.Units,")",sep=""),ylim = range(c(plotrng1,plotrng2)))

abline(h=0)


points(Transect[LRR.slope > 0], LRR.slope[LRR.slope > 0], type="p", col= "blue", pch= 1, cex=0.5)

points(Transect[LRR.slope <= 0], LRR.slope[LRR.slope <= 0], type="p", col= "red", pch= 1, cex=0.5)

if (length(DateList3) > 2) {

arrows(Transect,LRR.CI.L,Transect,LRR.CI.U,length = 0.05, angle = 90, code = 3,col = "gray")

arrows(Transect,LRR.CI.L,Transect,LRR.CI.U,length = 0.05, angle = 90, code = 3,col= ifelse(LRR.CI.L < 0 & LRR.CI.U < 0, "red","transparent"))

arrows(Transect,LRR.CI.L,Transect,LRR.CI.U,length = 0.05, angle = 90, code = 3,col= ifelse(LRR.CI.L > 0 & LRR.CI.U > 0, "blue","transparent")) }



dev.off()



# plot the WLR and confidence limits (CI error: lower limit & the upper limit)

plotrng1 <- min(WLR.CI.L[EPR != "NaN"])
plotrng2 <- max(WLR.CI.U[EPR != "NaN"])

if (length(DateList3) == 2) {
plotrng1 <- min(WLR.slope[EPR != "NaN"])
plotrng2 <- max(WLR.slope[EPR != "NaN"])
}


pdf("PDF/graph_WLR.pdf", bg="white")

plot(Transect, WLR.slope, type="p", lwd= 0, col= "black" , pch=0, cex=0, las= 1, cex.axis= 1, cex.lab= 1, xlab=expression(paste("Transect")),ylab=paste("WLR ","(",Map.Units,"/",Time.Units,")",sep=""), ylim = range(c(plotrng1,plotrng2)))

abline(h=0)


points(Transect[WLR.slope > 0], WLR.slope[WLR.slope > 0], type="p", col= "blue", pch= 1, cex=0.5)

points(Transect[WLR.slope <= 0], WLR.slope[WLR.slope <= 0], type="p", col= "red", pch= 1, cex=0.5)

if (length(DateList3) > 2) {

arrows(Transect,WLR.CI.L,Transect,WLR.CI.U,length = 0.05, angle = 90, code = 3,col = "gray")

arrows(Transect,WLR.CI.L,Transect,WLR.CI.U,length = 0.05, angle = 90, code = 3,col= ifelse(WLR.CI.L < 0 & WLR.CI.U < 0, "red","transparent"))

arrows(Transect,WLR.CI.L,Transect,WLR.CI.U,length = 0.05, angle = 90, code = 3,col= ifelse(WLR.CI.L > 0 & WLR.CI.U > 0, "blue","transparent"))
}

dev.off()



# plot the Mean EPR of Consecutive Eras and Stnd Deviations of rates as error bars (CI error: lower limit & the upper limit)

if (length(DateList3) > 2) {


plotrng1 <- min(Mean.EPR.Eras.L[EPR != "NaN"])
plotrng2 <- max(Mean.EPR.Eras.U[EPR != "NaN"])
}

if (length(DateList3) < 3) {
plotrng1 <- min(Mean.EPR.Eras[EPR != "NaN"])
plotrng2 <- max(Mean.EPR.Eras[EPR != "NaN"])
}


#plot and add the error bars
pdf("PDF/graph_EPR_Eras.pdf", bg="white")

plot(Transect, Mean.EPR.Eras, type="p", lwd= 0, col= "black" , pch=0, cex=0, las= 1, cex.axis= 1, cex.lab= 1, xlab=expression(paste("Transect")),ylab=paste("Mean EPR Eras ","(",Map.Units,"/",Time.Units,")",sep=""), ylim = range(c(plotrng1,plotrng2)))

abline(h=0)

points(Transect[Mean.EPR.Eras > 0], Mean.EPR.Eras[Mean.EPR.Eras > 0], type="p", col= "blue", pch= 1, cex=0.3)

points(Transect[Mean.EPR.Eras <= 0], Mean.EPR.Eras[Mean.EPR.Eras <= 0], type="p", col= "red", pch= 1, cex=0.3)


if (length(DateList3) > 2) {

arrows(Transect,Mean.EPR.Eras.L,Transect,Mean.EPR.Eras.U,length = 0.05, angle = 90, code = 3,col = "gray")

arrows(Transect,Mean.EPR.Eras.L,Transect,Mean.EPR.Eras.U,length = 0.05, angle = 90, code = 3,col= ifelse(Mean.EPR.Eras.L < 0 & Mean.EPR.Eras.U < 0, "red","transparent"))

arrows(Transect,Mean.EPR.Eras.L,Transect,Mean.EPR.Eras.U,length = 0.05, angle = 90, code = 3,col= ifelse(Mean.EPR.Eras.L > 0 & Mean.EPR.Eras.U > 0, "blue","transparent"))
}

dev.off()

################### Write shapefiles  #########################################################

####################set up a progress bar for the sapply function
############build progress bar for building transect lines
sapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}
########################end function






FinalGISTable2gis <- read.table("GIS_stats_table_short.csv", header=TRUE, sep=",")

time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_gisdata", showWarnings=FALSE)
setwd("AMBUR_gisdata")

#dir.create(paste(time.stamp2," ","gisdata",sep=""))
#setwd(paste(time.stamp2," ","gisdata",sep=""))

#replace characters in fields that GIS doesn't read
colnames(FinalGISTable2gis) <- gsub(".", "_", colnames(FinalGISTable2gis),fixed=TRUE)
colnames(FinalGISTable2gis) <- gsub(" ", "", colnames(FinalGISTable2gis),fixed=TRUE)


ddTable <- data.frame(FinalGISTable2gis)
coordinates(ddTable) <- data.frame((x=FinalGISTable2gis["Outer_X"]), y=(FinalGISTable2gis["Outer_Y"]))
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(ddTable) <- projectionString
writeOGR(ddTable, ".", "outer_pts", driver="ESRI Shapefile")



ddTable <- data.frame(FinalGISTable2gis) 
coordinates(ddTable) <- data.frame(x=(FinalGISTable2gis["Inner_X"]), y=(FinalGISTable2gis["Inner_Y"]))
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(ddTable) <- projectionString
writeOGR(ddTable, ".", "inner_pts", driver="ESRI Shapefile")


ddTable <- data.frame(FinalGISTable2gis) 
coordinates(ddTable) <- data.frame(x=(FinalGISTable2gis["Start_X"]), y=(FinalGISTable2gis["Start_Y"]))
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(ddTable) <- projectionString
writeOGR(ddTable, ".", "start_pts", driver="ESRI Shapefile")


ddTable <- data.frame(FinalGISTable2gis) 
coordinates(ddTable) <- data.frame(x=(FinalGISTable2gis["End_X"]), y=(FinalGISTable2gis["End_Y"]))
  projectionString <- proj4string(shapedata) # contains projection info
  proj4string(ddTable) <- projectionString
writeOGR(ddTable, ".", "end_pts", driver="ESRI Shapefile")


ddTable <- data.frame(FinalGISTable2gis) 
coordinates(ddTable) <- data.frame(x=(FinalGISTable2gis["Max_DateX"]), y=(FinalGISTable2gis["Max_DateY"]))
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(ddTable) <- projectionString
writeOGR(ddTable, ".", "max_date_pts", driver="ESRI Shapefile")


ddTable <- data.frame(FinalGISTable2gis) 
coordinates(ddTable) <- data.frame(x=(FinalGISTable2gis["Min_DateX"]), y=(FinalGISTable2gis["Min_DateY"]))
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(ddTable) <- projectionString
writeOGR(ddTable, ".", "min_date_pts", driver="ESRI Shapefile")

################################################################################################################
new_trandata <- data.frame(FinalGISTable2gis)

row.names(new_trandata) <- new_trandata$Transect

Transect.Factor <- factor(new_trandata$Transect)    #fixed 20130224 to get proper order of transects to match LineIDs with row.names of new_trandata

shape.final <- sapply_pb(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(new_trandata$Start_X[new_trandata$Transect == x], new_trandata$End_X[new_trandata$Transect == x]), y=c(new_trandata$Start_Y[new_trandata$Transect == x],new_trandata$End_Y[new_trandata$Transect == x])))), ID=(as.numeric(x))))
,simplify = TRUE)
shape.final2 <- SpatialLines(shape.final)
#edit(data.frame(getSLLinesIDSlots(shape.final2)) )
shape.final3 <- SpatialLinesDataFrame(shape.final2, new_trandata)
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(shape.final3) <- projectionString
writeOGR(shape.final3, ".", "original_transects", driver="ESRI Shapefile")


shape.final <- sapply_pb(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(new_trandata$Outer_X[new_trandata$Transect == x], new_trandata$Inner_X[new_trandata$Transect == x]), y=c(new_trandata$Outer_Y[new_trandata$Transect == x],new_trandata$Inner_Y[new_trandata$Transect == x])))), ID=(as.numeric(x))))
,simplify = TRUE)
shape.final2 <- SpatialLines(shape.final)
shape.final3 <- SpatialLinesDataFrame(shape.final2, new_trandata)
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(shape.final3) <- projectionString
writeOGR(shape.final3, ".", "envelope_transects_original", driver="ESRI Shapefile")


shape.final <- sapply_pb(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(new_trandata$OuterEnvX[new_trandata$Transect == x], new_trandata$InnerEnvX[new_trandata$Transect == x]), y=c(new_trandata$OuterEnvY[new_trandata$Transect == x],new_trandata$InnerEnvY[new_trandata$Transect == x])))), ID=(as.numeric(x))))
,simplify = TRUE)
shape.final2 <- SpatialLines(shape.final)
shape.final3 <- SpatialLinesDataFrame(shape.final2, new_trandata)
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(shape.final3) <- projectionString
writeOGR(shape.final3, ".", "envelope_transects_analysis", driver="ESRI Shapefile" )





shape.final <- sapply_pb(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(new_trandata$Min_DateX[new_trandata$Transect == x], new_trandata$Max_DateX[new_trandata$Transect == x]), y=c(new_trandata$Min_DateY[new_trandata$Transect == x],new_trandata$Max_DateY[new_trandata$Transect == x])))), ID=(as.numeric(x))))
,simplify = TRUE)
shape.final2 <- SpatialLines(shape.final)
shape.final3 <- SpatialLinesDataFrame(shape.final2, new_trandata)
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(shape.final3) <- projectionString
writeOGR(shape.final3, ".", "net_min_max__date_transects", driver="ESRI Shapefile")




################################################################################################################













#####################





#end timing the analysis and report the results
time1 <- Sys.time()
worktime <- time1 - time0


print(worktime)

#status checkpoint

mtext("...done!", side=3, line= -4, adj= -0, cex= 0.75, at=textalign)
mtext(paste("Elapsed time:"," ",format(worktime),".",sep=""), side=3, line= -5, adj= 0, cex= 0.75, at=textalign)
print("done!")

#tidy up and remove all objects
#detach(WorkTable1)
rm(list = ls())

#end the function
}

