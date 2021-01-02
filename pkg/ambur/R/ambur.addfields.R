ambur.addfields <-
function(checktype=1) {
#require(rgdal)
 #require(tcltk)

tkmessageBox(message = "Warning: Shapefile must not be empty.")
tkmessageBox(message = "Please select the shapefile...")

filetype <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- if(interactive()) tk_choose.files(filter = filetype)
shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)

 
   shapedata <- readOGR(getdata,layer=shapename)
   
   attrtable <- data.frame(shapedata)
   
   attrtable$Id <- seq(1,length(attrtable[,1]),1)

#checktype = 1  ####for testing


checktype <- checktype

if (checktype == 1) shapetype <- "shorelines"
if (checktype == 2) shapetype <- "baselines"


cat("Checking",shapetype,"shapefile table...","\n",sep=" ")


if (checktype == 1) reqfields <- c("Id","DATE_","ACCURACY","SHORE_LOC","CLASS_1","CLASS_2","CLASS_3","GROUP")

if (checktype == 2) reqfields <- c("Id","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC")



#check for missing fields
fieldcheck <- toupper(colnames(attrtable))

 presentfields <-  which(toupper(reqfields) %in% fieldcheck)

   missingfields <- as.character(reqfields[-presentfields])

   


cat("The following fields were missing and added:","\n")
cat(missingfields, "\n",sep = ", ")

add.fields <- data.frame(matrix(data= NA,ncol=length(missingfields),nrow=length(attrtable[,1])) )

colnames(add.fields) <- missingfields


## adjust the character field widths
fwidth <- rep(" ",30)
fwidth1 <- paste(fwidth, collapse = "")


if ("Id" %in% missingfields) add.fields$Id <- seq(1,length(attrtable[,1]),1)
if ("DATE_" %in% missingfields) add.fields$DATE_ <- "mm/dd/yyyy 12:00:01 AM"
if ("ACCURACY" %in% missingfields) add.fields$ACCURACY <- 1
if ("SHORE_LOC" %in% missingfields) add.fields$SHORE_LOC <- fwidth1
if ("CLASS_1" %in% missingfields) add.fields$CLASS_1 <- fwidth1
if ("CLASS_2" %in% missingfields) add.fields$CLASS_2 <- fwidth1
if ("CLASS_3" %in% missingfields) add.fields$CLASS_3 <- fwidth1
if ("GROUP" %in% missingfields) add.fields$GROUP <- fwidth1

if ("Location" %in% missingfields) add.fields$Location <- fwidth1
if ("MaxBNum" %in% missingfields) add.fields$MaxBNum <- 0
if ("SHORE_LOC" %in% missingfields) add.fields$SHORE_LOC <- fwidth1
if ("BaseOrder" %in% missingfields) add.fields$BaseOrder <- seq(1,length(attrtable[,1]),1)
if ("OFFshore" %in% missingfields) add.fields$OFFshore <- 1
if ("CastDir" %in% missingfields) add.fields$CastDir <- -1
if ("BASE_LOC" %in% missingfields) add.fields$BASE_LOC <- fwidth1



   new.table <-  data.frame(attrtable,add.fields)       

    
 
# finally, write a shape file (with .prj component)
  outputname <- paste("ambur_",shapename,sep="")
 

shape.final <- SpatialLinesDataFrame(shapedata, new.table)
 
   # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata) # contains projection info
  
  proj4string(shape.final) <- projectionString
  
  
   writeOGR(shape.final, ".", outputname, driver="ESRI Shapefile") 
   message("Added missing fields and created a new shapefile") 

tkmessageBox(message = paste("New shapefile --> ",outputname," <--generated with missing fields.",sep=""))


}

