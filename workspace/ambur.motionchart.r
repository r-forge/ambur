ambur.motionchart <-
function(id_var="Transect", time_var="Year", ht_var=550, wdth_var=900) {


#id_var="Transect"
#time_var="Year"

library(googleVis)

require(tcltk)
tkmessageBox(message = "Please select the ambur_motionchart_data.csv file...")
getdata <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata, header=TRUE, sep=",")
attach(mydata)

dir.path <- dirname(getdata)
setwd(dir.path)




ShoreMotion <- gvisMotionChart(mydata, idvar=id_var, timevar=time_var, options=list(height=550, width=900))
# Display chart
plot(ShoreMotion)

#create html file with code
ShoreMotion$html$chart
cat(ShoreMotion$html$chart, file="ambur_motionchart.html")

# Create Google Gadget
cat(createGoogleGadget(ShoreMotion), file="ambur_motionchart.xml")

}