#MasterMix for forensim
#Hinda Haned,
#July 2009
mastermix<-function()
{
	#if(!require(tcltk)) stop("package tcltk is required")
	#tclRequire("Tktable")
	#fonts definition
	font0 <- tcltk::tkfont.create(family="times",size=35,weight="bold",slant="italic")
	font1<-tcltk::tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
	font2<-tcltk::tkfont.create(family="times",size=16,weight="bold",slant="italic")
	font3<-tcltk::tkfont.create(family="times",size=10,weight="bold")#,slant="italic")
	font4<-tcltk::tkfont.create(family="times",size=10)#,slant="italic")
	#
	tt <- tcltk::tktoplevel()
	tcltk::tkwm.title(tt,"forensic DNA mixtures resolution")	
	dudioutvar <- tcltk::tclVar()
	dudivar <- tcltk::tclVar()
	facvar <- tcltk::tclVar()
	nfvar <- tcltk::tclVar()

	scannfvar <- tcltk::tclVar(1)
#
# Title
#
	TFrame <- tcltk::tkframe(tt, relief="groove")
	labh <- tcltk::tklabel(TFrame)
	font0 <- tcltk::tkfont.create(family="times",size=35,weight="bold",slant="italic")
	font1<-tcltk::tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
	font2<-tcltk::tkfont.create(family="times",size=16,weight="bold",slant="italic")
	
	tcltk::tkgrid(tcltk::tklabel(TFrame,text="MasterMix", font=font0, foreground="red"), labh)
	tcltk::tkgrid(tcltk::tklabel(TFrame,text="Two-person DNA mixtures resolution using allele peak height/ or area information", font=font2, foreground="red"), labh)
	#tcltk::tkbind(labh, "<Button-1>", function() print(help("between")))
	tcltk::tkgrid(TFrame)

	



#____________________________________________
	
RCSFrame <- tcltk::tkframe(tt, relief="groove")
A2.but <- tcltk::tkbutton(RCSFrame, text="Two-allele model", font=font1, command=A2.simu)
tcltk::tkbind(A2.but, "<Button-3>", function() print(help("A2.simu")))
A3.but <- tcltk::tkbutton(RCSFrame, text="Three-allele model",font=font1, command=A3.simu)
tcltk::tkbind(A3.but, "<Button-3>", function() print(help("A3.simu")))
A4.but <- tcltk::tkbutton(RCSFrame, text="Four-allele model", font=font1,command=A4.simu)
tcltk::tkbind(A4.but, "<Button-3>", function() print(help("A4.simu")))
tcltk::tkgrid(A2.but, A3.but,A4.but,ipadx=20)	
tcltk::tkgrid(RCSFrame)

}




	# tcltk::tkbind(getdata.but, "<Button-3>", function() print(help("data")))








