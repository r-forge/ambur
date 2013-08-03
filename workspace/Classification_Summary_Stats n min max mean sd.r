require(tcltk)
tkmessageBox(message = "Please select the results_stats.csv file...")
getdata <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata, header=TRUE, sep=",")
attach(mydata)

dir.path <- dirname(getdata)
setwd(dir.path)

filename <- gsub(".csv", "", basename(getdata))



target.field <- select.list(colnames(mydata), multiple = FALSE, title = "Choose the stats field:")  

field.stats <- mydata[,target.field]

class.field <- select.list(colnames(mydata), multiple = FALSE, title = "Choose the classification field:")  

class.stats <- mydata[,class.field]

#change from adding zero if you want to simulate a curve
#field.stats <- field.stats[field.stats > 0] + 0

gclass <- as.character(unique(class.stats) )

min.class <- numeric(length(gclass))
max.class <- numeric(length(gclass))
mean.class <- numeric(length(gclass))
sd.class <- numeric(length(gclass))
n.class <- numeric(length(gclass))

for (i in 1:length(gclass)) {

min.class[i] <- min(field.stats[mydata[class.field] == gclass[i]],na.rm=T)
max.class[i] <- max(field.stats[mydata[class.field] == gclass[i]],na.rm=T)
mean.class[i] <- mean(field.stats[mydata[class.field] == gclass[i]],na.rm=T)
sd.class[i] <- sd(field.stats[mydata[class.field] == gclass[i]],na.rm=T)
n.class[i] <- length(field.stats[mydata[class.field] == gclass[i]])

}

data.table <- data.frame(gclass,n.class,min.class,max.class,mean.class,sd.class)

data.table

outputname <- paste(class.field,target.field,filename,"summary.csv",sep="_")


write.table(data.table, file = outputname, sep = ",", row.names = FALSE)

barplot(mean.class,ylab="Mean",xlab="",axes=T,width=1,names.arg=gclass)


