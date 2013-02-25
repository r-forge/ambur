require(tcltk)
tkmessageBox(message = "Please select the results_stats.csv file...")
getdata <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata, header=TRUE, sep=",")
attach(mydata)

dir.path <- dirname(getdata)
setwd(dir.path)



target.field <- select.list(colnames(mydata), multiple = FALSE, title = "Choose the stats field:")  

field.stats <- mydata[,target.field]

class.field <- select.list(colnames(mydata), multiple = FALSE, title = "Choose the classification field:")  

class.stats <- mydata[,class.field]

#change from adding zero if you want to simulate a curve
#field.stats <- field.stats[field.stats > 0] + 0

gclass <- as.character(unique(class.stats) )


mean.class <- numeric(length(gclass))
mean.eros.class <- numeric(length(gclass))
mean.acc.class <- numeric(length(gclass))
mean.peros.class <- numeric(length(gclass))

for (i in 1:length(gclass)) {

mean.class[i] <- mean(field.stats[mydata[class.field] == gclass[i]],na.rm=T)
mean.eros.class[i] <- mean(field.stats[field.stats < 0 & mydata[class.field] == gclass[i]],na.rm=T)
mean.acc.class[i] <- mean(field.stats[field.stats  > 0 & mydata[class.field] == gclass[i]],na.rm=T)
mean.peros.class[i] <- length(field.stats[field.stats < 0 & mydata[class.field] == gclass[i]]) / length(field.stats[mydata[class.field] == gclass[i]]) * 100

}

data.table <- data.frame(gclass,mean.class,mean.peros.class,mean.eros.class,mean.acc.class)

data.table


write.table(data.table, file = "classification_summary.csv", sep = ",", row.names = FALSE)

barplot(mean.eros.class,ylab="Mean Erosion Rate",xlab="",axes=T,width=1,names.arg=gclass)


