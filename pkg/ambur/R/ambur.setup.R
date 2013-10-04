  ambur.setup <- function(setup=1) {
 packs <- installed.packages()
 exc <- names(packs[,'Package'])

 av <- names(available.packages()[,1])

 ins <- av[!av %in% exc]

 req.packs <- c("locfit", "sp", "spatial", "spatstat", "googleVis", "rgeos", "rgdal", "tcltk", "stringr")

 inst.packs <- req.packs[!req.packs %in% exc]
 
 packs.ok <- "All required packages are installed"
 packs.rem <- "All required packages have been removed"
 

 
 ifelse(setup == 1,ifelse(length(inst.packs) == 0, packs.ok,install.packages(inst.packs)),remove.packages(req.packs))
 
  
 test <- ifelse(setup == 1,packs.ok,packs.rem)
 
 test
 
 
 }
