ambur.addfields <-
function(checktype=1) {
require(tcltk)
library(foreign)

tkmessageBox(message = "Warning: Shapefile must not be empty.")

getdata <- tk_choose.files(default = "*.dbf",multi = FALSE)
shapename <- gsub(".dbf", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)
path <- getdata

mydata <- foreign::read.dbf(path)

checktype <- checktype

if (checktype == 1) shapetype <- "shorelines"
if (checktype == 2) shapetype <- "baselines"


cat("Checking",shapetype,"shapefile table...","\n",sep=" ")


if (checktype == 1) reqfields <- c("Id","DATE_","ACCURACY","SHORE_LOC","CLASS_1","CLASS_2","CLASS_3","GROUP")

if (checktype == 2) reqfields <- c("Id","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC")



#check for missing fields
fieldcheck <- toupper(colnames(mydata))

 presentfields <-  which(toupper(reqfields) %in% fieldcheck)

   missingfields <- as.character(reqfields[-presentfields])

   


cat("The following fields were missing and added:","\n")
cat(missingfields, "\n",sep = ", ")

add.fields <- data.frame(matrix(data= NA,ncol=length(missingfields),nrow=length(mydata[,1])) )



 colnames(add.fields) <- missingfields



if ("Id" %in% missingfields) add.fields$Id <- seq(1,length(mydata[,1]),1)
if ("DATE_" %in% missingfields) add.fields$DATE_ <- "mm/dd/yyyy 12:00:01 AM"
if ("ACCURACY" %in% missingfields) add.fields$ACCURACY <- 1
if ("SHORE_LOC" %in% missingfields) add.fields$SHORE_LOC <- "NA"
if ("CLASS_1" %in% missingfields) add.fields$CLASS_1 <- "NA"
if ("CLASS_2" %in% missingfields) add.fields$CLASS_2 <- "NA"
if ("CLASS_3" %in% missingfields) add.fields$CLASS_3 <- "NA"
if ("GROUP" %in% missingfields) add.fields$GROUP <- "NA"

if ("Location" %in% missingfields) add.fields$Location <- "NA"
if ("MaxBNum" %in% missingfields) add.fields$MaxBNum <- 0
if ("SHORE_LOC" %in% missingfields) add.fields$SHORE_LOC <- "NA"
if ("BaseOrder" %in% missingfields) add.fields$BaseOrder <- seq(1,length(mydata[,1]),1)
if ("OFFshore" %in% missingfields) add.fields$OFFshore <- 1
if ("CastDir" %in% missingfields) add.fields$CastDir <- -1
if ("BASE_LOC" %in% missingfields) add.fields$BASE_LOC <- "NA"



   new.table <-  cbind(mydata,add.fields)       
       foreign::write.dbf(new.table,path) 



#detach("package:foreign")
}

