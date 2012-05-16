
require(tcltk)

tkmessageBox(message = "Please select the directory...")
getdata <- tk_choose.dir(default = "", caption = "Select directory")
workingdir <- getdata
setwd(workingdir)

list.jpg <- dir(pattern = ".jpg")
repair.name <- toupper(gsub("_dd.jpg", "", list.jpg))



tkmessageBox(message = "Please select a *.csv file to compare the list against...")
getdata2 <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata2, header=TRUE, sep=",")


target.field <- select.list(colnames(mydata), multiple = FALSE, title = "Choose field:")  

field.stats <- mydata[,target.field]


#find the items not present in the directory list based on the csv file
present.files <- repair.name

all.files <- field.stats

missing.files <- field.stats[!(field.stats %in% repair.name)]

write.table(missing.files, file = "manuscripts_missing.csv", quote = FALSE, sep = ",", row.names = FALSE)
write.table(present.files, file = "manuscripts_present.csv", quote = FALSE, sep = ",", row.names = FALSE)

###########
 
