ambur.hva.gui <-
function(){

 #require(tcltk)

       myData1 <- "path"
       myData2 <- "path"
       myData3 <- "path"
       myData4 <- "path"
       myDir1 <-  "path"
       remove.anthro <- "no"
       myLocation <- ""
  
  tt <- tktoplevel(width=800,height=800)
       tkwm.title(tt,"AMBUR-HVA")





 
       getfile <- function() {
    name <- tclvalue(tkgetOpenFile(
    filetypes = "{{Shapefile} {.shp}}"))
    if (name == "") return;
    Shapedata1 <- file.path(name)
    return(Shapedata1)
}

       getdir <- function() {
    name <- tclvalue(tkchooseDirectory(title = "Please select a directory to store AMBUR-HVA results in"))
    if (name == "") return;
    ShapeDir <- name
    return(ShapeDir)
}


       mloc <- tclVar(myLocation)
       xvar <- tclVar(myData1)
       yvar <- tclVar(myData2)
       zvar <- tclVar(myData3)
       zzvar <- tclVar(myData4)
       myDir1var <-  tclVar(myDir1)
       remove.anthrovar <- tclVar(remove.anthro)
      
       b1var <- tclVar(150)
       b2var <- tclVar(10)
       scr1avar <- tclVar(-1)
       scr1bvar <- tclVar(-0.2)
       scr1cvar <- tclVar(0.2)
       scr1dvar <- tclVar(1)
       scr2avar <- tclVar(0.5)
       scr2bvar <- tclVar(0.2)
       scr3avar <- tclVar(0.16)
       scr3bvar <- tclVar(0.074)
       soviavar <- tclVar(6.77)
       sovibvar <- tclVar(4.06)
       sovicvar <- tclVar(-1.37)
       sovidvar <- tclVar(-4.09)

       myDir1.entry <- tkentry(tt, textvariable=myDir1var)
       remove.anthro.entry <- tkentry(tt, textvariable=remove.anthrovar)
       loc.entry <- tkentry(tt, textvariable=mloc)
       x.entry <- tkentry(tt, textvariable=xvar)
       y.entry <- tkentry(tt, textvariable=yvar)
       z.entry <- tkentry(tt, textvariable=zvar)
       zz.entry <- tkentry(tt, textvariable=zzvar)
       b1.entry <- tkentry(tt, textvariable=b1var)
       b2.entry <- tkentry(tt, textvariable=b2var)
       scr1a.entry  <- tkentry(tt, textvariable=scr1avar)
       scr1b.entry <- tkentry(tt, textvariable=scr1bvar)
       scr1c.entry  <- tkentry(tt, textvariable=scr1cvar)
       scr1d.entry  <- tkentry(tt, textvariable=scr1dvar)
       scr2a.entry  <- tkentry(tt, textvariable=scr2avar)
       scr2b.entry  <- tkentry(tt, textvariable=scr2bvar)
       scr3a.entry  <- tkentry(tt, textvariable=scr3avar)
       scr3b.entry  <- tkentry(tt, textvariable=scr3bvar)
       sovia.entry  <- tkentry(tt, textvariable=soviavar)
       sovib.entry <- tkentry(tt, textvariable=sovibvar)
       sovic.entry  <- tkentry(tt, textvariable=sovicvar)
       sovid.entry  <- tkentry(tt, textvariable=sovidvar)

 
       Update1 <- function() {
         tclvalue(xvar)<- getfile()
         setwd(dirname(tclvalue(xvar)))       
        }
 
        Update2 <- function() {
         tclvalue(yvar)<- getfile()
         setwd(dirname(tclvalue(yvar)))       
        }
        
        Update3 <- function() {
         tclvalue(zvar)<- getfile()
         setwd(dirname(tclvalue(zvar)))       
        }
        
        Update4 <- function() {
         tclvalue(zzvar)<- getfile()
         setwd(dirname(tclvalue(zzvar)))       
        }
 
        Update5 <- function() {
         tclvalue(myDir1var)<- getdir()       
        }
 
 
      button.widget <- tkbutton(tt, text = "Select AMBUR Envelope Transects Shapefile", command = Update1)
      button.widget2 <- tkbutton(tt, text = "Select FEMA Q3/DFIRM Shapefile", command = Update2)
      button.widget3 <- tkbutton(tt, text = "Select NOAA/NWS/USACOE SLOSH Shapefile", command = Update3)
      button.widget4 <- tkbutton(tt, text = "Select NOAA SoVI Shapefile", command = Update4)
      button.widget5 <- tkbutton(tt, text = "Select a directory to store results", command = Update5)

 
 

        

   

 

       tkgrid(tklabel(tt,text="Enter or select your shapefiles and enter your breakpoint values (from high risk to low risk):"),columnspan=3, pady = 10)
       tkgrid(tklabel(tt,text="enter location name (optional)"), loc.entry, pady= 10, padx= 10)
       tkgrid(button.widget5, myDir1.entry, pady= 10, padx= 10)
       tkgrid(button.widget, x.entry, pady= 10, padx= 10)
       tkgrid(button.widget2, y.entry, pady= 10, padx= 10)
       tkgrid(button.widget3, z.entry, pady= 10, padx= 10)
       tkgrid(button.widget4, zz.entry, pady= 10, padx= 10)
       tkgrid(tklabel(tt,text="remove anthropogenic transects"), remove.anthro.entry, pady= 10, padx= 10)
       tkgrid(tklabel(tt,text="spatial search radius (map units)"), b1.entry, pady= 10, padx= 10)
       tkgrid(tklabel(tt,text="transect buffer (map units)"), b2.entry, pady= 10, padx= 10)
       tkgrid(tklabel(tt,text="shoreline change rate (map units/yr) breakpoints"), scr1a.entry, scr1b.entry, scr1c.entry, scr1d.entry, pady= 10, padx= 10)
       tkgrid(tklabel(tt,text="temporal shoreline change (st.dev map units/yr) breakpoints"), scr2a.entry, scr2b.entry, pady= 10, padx= 10)
       tkgrid(tklabel(tt,text="spatial shoreline change (st.dev map units/yr) breakpoints"), scr3a.entry, scr3b.entry, pady= 10, padx= 10)
       tkgrid(tklabel(tt,text="SoVI (score) breakpoints"), sovia.entry, sovib.entry, sovic.entry, sovid.entry, pady= 10, padx= 10)



 workingdir <- dirname(tclvalue(myDir1var))
setwd(workingdir)
 
 Name <- tclVar("first")
entry.Name <-tkentry(tt,width="20",textvariable=Name)  
   
OnOK <- function()
{
	NameVal <- tclvalue(Name)
	tkdestroy(tt)

workingdir <- dirname(tclvalue(myDir1var))
setwd(workingdir)

    	
#	     	msg <- paste("ambur.hva( ","loc=",tclvalue(mloc),", ","bufferdist=",as.numeric(tclvalue(b1var)),", ","bufferdist2=",as.numeric(tclvalue(b2var)),", ","removeanthro=",as.character(tclvalue(remove.anthrovar)),", ","file1=",as.character(tclvalue(xvar)),", ","file2=",as.character(tclvalue(yvar)),", ","file3=",as.character(tclvalue(zvar)),", ","file4=",as.character(tclvalue(zzvar)),", ","dir1=",as.character(tclvalue(myDir1var)),", ","scr=","c(",as.numeric(tclvalue(scr1avar)),", ",as.numeric(tclvalue(scr1bvar)),", ",as.numeric(tclvalue(scr1cvar)),", ",as.numeric(tclvalue(scr1dvar)),")", "," ,"tscr=","c(",as.numeric(tclvalue(scr2avar)),", ",as.numeric(tclvalue(scr2bvar)),")",",", "sscr=","c(",as.numeric(tclvalue(scr3avar)),",",as.numeric(tclvalue(scr3bvar)),")", ",","sovi=","c(",as.numeric(tclvalue(soviavar)),", ",as.numeric(tclvalue(sovibvar)),", ",as.numeric(tclvalue(sovicvar)),", ",as.numeric(tclvalue(sovidvar)),")", ")")
	
#  tkmessageBox(message=msg)
  
#  txt <- c("ambur.hva( ","loc=",paste('"',tclvalue(mloc),'"',sep=""),", ","bufferdist=",as.numeric(tclvalue(b1var)),", ","bufferdist2=",as.numeric(tclvalue(b2var)),", ","removeanthro=",paste('"',tclvalue(remove.anthrovar),'"',sep=""),", ","file1=",paste('"',tclvalue(xvar),'"',sep=""),", ","file2=",paste('"',tclvalue(yvar),'"',sep=""),", ","file3=",paste('"',tclvalue(zvar),'"',sep=""),", ","file4=",paste('"',tclvalue(zzvar),'"',sep=""),", ","dir1=",paste('"',tclvalue(myDir1var),'"',sep=""),", ","scr=","c(",as.numeric(tclvalue(scr1avar)),", ",as.numeric(tclvalue(scr1bvar)),", ",as.numeric(tclvalue(scr1cvar)),", ",as.numeric(tclvalue(scr1dvar)),")", "," ,"tscr=","c(",as.numeric(tclvalue(scr2avar)),", ",as.numeric(tclvalue(scr2bvar)),")",",", "sscr=","c(",as.numeric(tclvalue(scr3avar)),",",as.numeric(tclvalue(scr3bvar)),")", ",","sovi=","c(",as.numeric(tclvalue(soviavar)),", ",as.numeric(tclvalue(sovibvar)),", ",as.numeric(tclvalue(sovicvar)),", ",as.numeric(tclvalue(sovidvar)),")", ")")
  
# writeLines(txt, "parameters.txt")
    

  ambur.hva(loc=tclvalue(mloc),bufferdist=as.numeric(tclvalue(b1var)),bufferdist2=as.numeric(tclvalue(b2var)),removeanthro=as.character(tclvalue(remove.anthrovar)),file1=as.character(tclvalue(xvar)),file2=as.character(tclvalue(yvar)),file3=as.character(tclvalue(zvar)) ,file4=as.character(tclvalue(zzvar)),dir1=as.character(tclvalue(myDir1var)),scr=c(as.numeric(tclvalue(scr1avar)),as.numeric(tclvalue(scr1bvar)),as.numeric(tclvalue(scr1cvar)),as.numeric(tclvalue(scr1dvar))),tscr=c(as.numeric(tclvalue(scr2avar)),as.numeric(tclvalue(scr2bvar))),sscr=c(as.numeric(tclvalue(scr3avar)),as.numeric(tclvalue(scr3bvar))),sovi=c(as.numeric(tclvalue(soviavar)) ,as.numeric(tclvalue(sovibvar)) ,as.numeric(tclvalue(sovicvar)),as.numeric(tclvalue(sovidvar))))
  

}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.Name, "<Return>",OnOK)
tkgrid(OK.but)
tkfocus(tt)    
 
 
 
   
    }
