library(googleVis)

require(tcltk)
tkmessageBox(message = "Please select the debugging3WorkTable1.csv file...")
getdata <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata, header=TRUE, sep=",")
attach(mydata)

dir.path <- dirname(getdata)
setwd(dir.path)


Setup.Date <- as.POSIXlt(mydata[ ,"DATE"], origin= "01/01/1970 12:0:0 AM", format="%m/%d/%Y %I:%M:%S %p")

Setup.Year <- as.numeric(format(Setup.Date, "%Y"))


#fix transect id to get it to plot alphabetically

fill_zeros <- paste("%","0",max(nchar(mydata[ ,"TRANSECT"])),"d",sep="")
Setup.Transect <- sprintf(fill_zeros, mydata[ ,"TRANSECT"])

#set up field for shoreline location for Jekyll Island
Setup.Location <- numeric(length(mydata[ ,"TRANSECT"]))
Setup.Location[mydata[ ,"TRANSECT"] >= 1 & mydata[ ,"TRANSECT"] <= 57] <- "south inlet"
Setup.Location[mydata[ ,"TRANSECT"] >= 58 & mydata[ ,"TRANSECT"] <= 330] <- "backbarrier"
Setup.Location[mydata[ ,"TRANSECT"] >= 331 & mydata[ ,"TRANSECT"] <= 356] <- "north inlet"
Setup.Location[mydata[ ,"TRANSECT"] >= 357 & mydata[ ,"TRANSECT"] <= 586] <- "oceanfront"

Setup.Distance <- (mydata[ , "Setup.Min.Date.Dist"]* ifelse(mydata[ ,"OFFSHORE"]==1,-1,1)) +  mydata[ ,"DISTANCE"]

test <- data.frame(mydata[ ,"X_COORD"],mydata[ ,"Y_COORD"],Setup.Transect,as.Date(Setup.Date, format="%m/%d/%Y %I:%M:%S %p"),Setup.Year, Setup.Distance,mydata[ ,"ACCURACY"], as.character(mydata[ ,"CLASS_1"]),Setup.Location)

colnames(test) <- c("X", "Y", "Transect", "Date", "Year", "Distance", "Accuracy", "Class_1", "Location")

#try to make equal axes
testdiffx <- max(mydata[ ,"X_COORD"]) - min(mydata[ ,"X_COORD"])
testdiffy <- max(mydata[ ,"Y_COORD"]) - min(mydata[ ,"Y_COORD"])
testdiff <- testdiffx - testdiffy
dum_max_x <- ifelse(testdiff > 0, max(mydata[ ,"X_COORD"]), max(mydata[ ,"X_COORD"] + abs(testdiff)))
dum_max_y <- ifelse(testdiff < 0, max(mydata[ ,"Y_COORD"]), max(mydata[ ,"Y_COORD"] + abs(testdiff)))
dum_max_tran <- "00"
dum_date <- test$Date[1]
dum_year <- min(Setup.Year)
dum_dist <- 0
dum_accuracy <- 0
dum_class1 <- "na"
dum_location <- ""

test_dum <- data.frame(dum_max_x,dum_max_y,dum_max_tran,dum_date,dum_year,dum_dist,dum_accuracy, dum_class1,dum_location)
colnames(test_dum) <- c("X", "Y", "Transect", "Date", "Year", "Distance", "Accuracy", "Class_1", "Location")

test1 <- rbind(test, test_dum)
test1$Class_1[is.na(test1$Class_1)] <- "na"

test1$X <- test1$X - min(test1$X)
test1$Y <- test1$Y - min(test1$Y)

test1$Size_1 <- 1
test1$Size_1[length(test1$Size_1)] <- 100

##test plot
ShoreMotion <- gvisMotionChart(test1, idvar="Transect", timevar="Year", options=list(height=550, width=900))
# Display chart
plot(ShoreMotion)
# Create Google Gadget


cat(createGoogleGadget(ShoreMotion), file="shoremotionchart.xml")

plot(test1$X, test1$Y)
max(test1$X) - min(test1$X)
max(test1$Y) - min(test1$Y)

#{"xZoomedDataMax":519000,"iconType":"BUBBLE","playDuration":15000,"xZoomedDataMin":517000,"iconKeySettings":[],"nonSelectedAlpha":0.4,"yZoomedDataMax":3563000,"yZoomedIn":false,"sizeOption":"_UNISIZE","xZoomedIn":false,"uniColorForNonSelected":false,"yLambda":1,"yAxisOption":"3","time":"2006","duration":{"multiplier":1,"timeUnit":"D"},"orderedByX":false,"dimensions":{"iconDimensions":["dim0"]},"orderedByY":false,"yZoomedDataMin":3561000,"xAxisOption":"2","xLambda":1,"colorOption":"4","showTrails":false}

##Example from Google
#Motion <- gvisMotionChart(Fruits, idvar="Fruit", timevar="Year", options=list(height=350, width=400))
# Display chart
#plot(Motion)
# Create Google Gadget
#cat(createGoogleGadget(Motion), file="motionchart.xml")