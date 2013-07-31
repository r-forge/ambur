merge.csv.files <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1



require(rgdal)
require(rgeos)
require(tcltk)


tkmessageBox(message = "Please select the directory with the CSV files...")
getdata <- tk_choose.dir(default = "", caption = "Select directory")
workingdir <- getdata
setwd(workingdir)


filenames <- list.files(path = getdata)
merge.tab <- do.call("rbind", lapply(filenames, read.csv, header = TRUE))

write.table(merge.tab, file = "merged_csv_data.csv", sep = ",", row.names = FALSE)

]