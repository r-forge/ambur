ambur.gui <-
function() {

# Begin AMBUR GUI: Create a "main" window with a button which activates dialog buttons

#require(tcltk)

ttMain <- tktoplevel()
tktitle(ttMain) <- "ambur"


# ambur analysis 
launchAnalysis <- function() {
tt<-tktoplevel()
tktitle(tt) <- "ambur.analysis"
Name <- tclVar("first")
entry.Name <-tkentry(tt,width="20",textvariable=Name)
tkgrid(tklabel(tt,text="first or last intersection:"))
tkgrid(entry.Name)
Conf <- tclVar(95)
entry.Conf <-tkentry(tt,width="20",textvariable=Conf)
tkgrid(tklabel(tt,text="confidence level:"))
tkgrid(entry.Conf)
Ulab <- tclVar("m")
entry.Ulab <-tkentry(tt,width="20",textvariable=Ulab)
tkgrid(tklabel(tt,text="unit label:"))
tkgrid(entry.Ulab)
Dtag <- tclVar("ambur")
entry.Dtag <-tkentry(tt,width="20",textvariable=Dtag)
tkgrid(tklabel(tt,text="folder tag (optional):"))
tkgrid(entry.Dtag)
Stran <- tclVar(1)
entry.Stran <-tkentry(tt,width="20",textvariable=Stran)
tkgrid(tklabel(tt,text="start transect:"))
tkgrid(entry.Stran)
Etran <- tclVar("all")
entry.Etran <-tkentry(tt,width="20",textvariable=Etran)
tkgrid(tklabel(tt,text="end transect:"))
tkgrid(entry.Etran)
Atype <- tclVar("basic")
entry.Atype <-tkentry(tt,width="20",textvariable=Atype)
tkgrid(tklabel(tt,text="analysis type:"))
tkgrid(entry.Atype)
Runit <- tclVar("yr")
entry.Runit <-tkentry(tt,width="20",textvariable=Runit)
tkgrid(tklabel(tt,text="Time unit for rates:"))
tkgrid(entry.Runit)

OnOK <- function()
{
	NameVal <- tclvalue(Name)
	tkdestroy(tt)
	msg <- paste("ambur.analysis( ",tclvalue(Name),", ",as.numeric(tclvalue(Conf)),", ",tclvalue(Ulab),", ",tclvalue(Dtag),", ",tclvalue(Stran),", ",tclvalue(Etran),", ",tclvalue(Atype),", ",tclvalue(Runit),")")
	tkmessageBox(message=msg)

ambur.analysis(tclvalue(Name), as.numeric(tclvalue(Conf)), tclvalue(Ulab),  tclvalue(Dtag), as.numeric(tclvalue(Stran)), ifelse(tclvalue(Etran) != "all", as.numeric(tclvalue(Etran)),tclvalue(Etran)), tclvalue(Atype), tclvalue(Runit))  
	
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.Name, "<Return>",OnOK)
tkgrid(OK.but)
tkfocus(tt)    
}

# ambur.addfields

launchAddfields <- function() {
tt<-tktoplevel()
tktitle(tt) <- "ambur.addfields"
Name <- tclVar(1)
entry.Name <-tkentry(tt,width="20",textvariable=Name)
tkgrid(tklabel(tt,text="Add fields: enter '1' for shorelines or '2' for baselines:"))
tkgrid(entry.Name)


OnOK <- function()
{
	NameVal <- tclvalue(Name)
	tkdestroy(tt)
	msg <- paste("ambur.addfields( ",tclvalue(Name),")")
	tkmessageBox(message=msg)

ambur.addfields(as.numeric(tclvalue(Name)))  
	
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.Name, "<Return>",OnOK)
tkgrid(OK.but)
tkfocus(tt)    
}


# ambur.transects

launchTransects <- function() {
tt<-tktoplevel()
tktitle(tt) <- "ambur.transects"
Space <- tclVar(50)
entry.Space <-tkentry(tt,width="20",textvariable=Space)
tkgrid(tklabel(tt,text="Transect spacing (map units):"))
tkgrid(entry.Space)
Length <- tclVar(500)
entry.Length <-tkentry(tt,width="20",textvariable=Length)
tkgrid(tklabel(tt,text="Transect length (map units):"))
tkgrid(entry.Length)



OnOK <- function()
{
	NameVal <- tclvalue(Space)
	tkdestroy(tt)
	msg <- paste("ambur.transects( ",tclvalue(Space),", ",tclvalue(Length),")")
	tkmessageBox(message=msg)

ambur.transects(as.numeric(tclvalue(Space)),as.numeric(tclvalue(Length)))  
	
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.Space, "<Return>",OnOK)
tkgrid(OK.but)
tkfocus(tt)    
}


 # ambur.filtertran

launchFiltertran <- function() {
tt<-tktoplevel()
tktitle(tt) <- "ambur.filtertran"
Winsize <- tclVar(5)
entry.Winsize <-tkentry(tt,width="20",textvariable=Winsize)
tkgrid(tklabel(tt,text="Window size (must be odd #):"))
tkgrid(entry.Winsize)
Indv <- tclVar(1)
entry.Indv <-tkentry(tt,width="20",textvariable=Indv)
tkgrid(tklabel(tt,text="Filter individual baselines: enter '1' for yes or '0' for no:"))
tkgrid(entry.Indv)



OnOK <- function()
{
	NameVal <- tclvalue(Winsize)
	tkdestroy(tt)
	msg <- paste("ambur.transects( ",tclvalue(Winsize),", ",tclvalue(Indv),")")
	tkmessageBox(message=msg)

ambur.filtertran(as.numeric(tclvalue(Winsize)),as.numeric(tclvalue(Indv)))  
	
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.Winsize, "<Return>",OnOK)
tkgrid(OK.but)
tkfocus(tt)    
}

# ambur.capture

launchCapture <- function() {
tt<-tktoplevel()
tktitle(tt) <- "ambur.capture"

tkgrid(tklabel(tt,text="Press 'OK' to capture shoreline points:"))



OnOK <- function()
{
	
	tkdestroy(tt)


ambur.capture()  
	
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)

tkgrid(OK.but)
tkfocus(tt)    
}


# ambur.statshape

launchStatshape <- function() {
tt<-tktoplevel()
tktitle(tt) <- "ambur.statshape"

tkgrid(tklabel(tt,text="Press 'OK' to build shapefiles from analysis results:"))



OnOK <- function()
{
	
	tkdestroy(tt)


ambur.statshape()  
	
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)

tkgrid(OK.but)
tkfocus(tt)    
}


# ambur.forecast

launchForecast <- function() {
tt<-tktoplevel()
tktitle(tt) <- "ambur.forecast"
Name <- tclVar(50)
entry.Name <-tkentry(tt,width="20",textvariable=Name)
tkgrid(tklabel(tt,text="Enter # of years to extrapolate shoreline position:"))
tkgrid(entry.Name)


OnOK <- function()
{
	NameVal <- tclvalue(Name)
	tkdestroy(tt)
	msg <- paste("ambur.forecast( ",tclvalue(Name),")")
	tkmessageBox(message=msg)

ambur.forecast(as.numeric(tclvalue(Name)))  
	
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.Name, "<Return>",OnOK)
tkgrid(OK.but)
tkfocus(tt)    
}




#add buttons to the main menu
launchDlg.button1 <- tkbutton(ttMain, text = "Add fields to shapefiles", command = launchAddfields)
launchDlg.button2 <- tkbutton(ttMain, text = "Construct transects", command = launchTransects)
launchDlg.button3 <- tkbutton(ttMain, text = "Filter transects", command = launchFiltertran)
launchDlg.button4 <- tkbutton(ttMain, text = "Capture shoreline & transect intersection points", command = launchCapture)
launchDlg.button5 <- tkbutton(ttMain, text = "Analyze shoreline capture points", command = launchAnalysis)
launchDlg.button6 <- tkbutton(ttMain, text = "Build shapefiles from analysis results", command = launchStatshape)
launchDlg.button7 <- tkbutton(ttMain, text = "Extrapolate future/past shoreline positions", command = launchForecast)

tkpack(launchDlg.button1,launchDlg.button2,launchDlg.button3,launchDlg.button4,launchDlg.button5,launchDlg.button6,launchDlg.button7)




###########################

}

