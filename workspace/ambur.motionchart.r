library(googleVis)

require(tcltk)
tkmessageBox(message = "Please select the ambur_motionchart_data.csv file...")
getdata <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata, header=TRUE, sep=",")
attach(mydata)

dir.path <- dirname(getdata)
setwd(dir.path)


test1 <- data.frame(mydata)

ShoreMotion <- gvisMotionChart(mydata, idvar="Transect", timevar="Year", options=list(height=550, width=900))
# Display chart
plot(ShoreMotion)
# Create Google Gadget


cat(createGoogleGadget(ShoreMotion), file="ambur_motionchart.xml")