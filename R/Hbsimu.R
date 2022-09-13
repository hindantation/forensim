# Hinda Haned, May 2010- Lyon
"Hbsimu" <- function()
{
	#if(!require(tcltk)) stop("package tcltk is required")
	#if(!require(tkrplot)) stop("package tkrplot is required")
	tcltk::tclRequire("Tktable")
	font0 <- tcltk::tkfont.create(family="courrier",size=35,weight="bold",slant="italic")
	font1<-tcltk::tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
	font2<-tcltk::tkfont.create(family="times",size=16,weight="bold",slant="italic")
	font3<-tcltk::tkfont.create(family="times",size=12)#,slant="italic")
	font4<-tcltk::tkfont.create(family="courrier",size=14)#,slant="italic")
	font5<-tcltk::tkfont.create(family="courrier",size=13,weight="bold")#,slant="italic")
	font6<-tcltk::tkfont.create(family="times",size=13)#tcltk::tkframe entries labels
	tf <- tcltk::tktoplevel()
	tcltk::tkwm.title(tf,"Hbdemo: heterozygote balance simulation")
	
	
	
	
	#function Hb

	done <- tcltk::tclVar(0)
		
	frame1 <- tcltk::tkframe(tf, relief="groove", borderwidth=4)
	labh <- tcltk::tklabel(tf)
	#labh <- tcltk::tklabel(frame1)
	#tkbind(labh, "<Button-1>", function() 'hh')
	tcltk::tkgrid(tcltk::tklabel(frame1,text="   Heterozygote balance simulator  ", font="courrier 22", 
	foreground="darkblue"),labh)
	#tkbind(frame1)
	tcltk::tkpack(frame1)#,padx=10)#,pady=10)

	frame2 <- tcltk::tkframe(tf, relief="flat", borderwidth=2)#,bg="white")	
	frame3 <- tcltk::tkframe(tf, relief="groove", borderwidth=2)#,bg="white")	
	# tab.entry <- tcltk::tkentry(frame2, textvariable=tabnamevar)
	# file.entry <- tcltk::tkentry(frame2, textvariable=filenamevar)
	# choosefile.but <- tcltk::tkbutton(frame2, text="Set", command=function() print('hh'))
	#tcltk::tkgrid(tcltk::tklabel(frame1, text="- Input parameters -",font=font3,foreground="blue"), columnspan=6)
	# tcltk::tkgrid(tcltk::tklabel(frame2,text=" "), file.entry)
	# tcltk::tkgrid(tcltk::tklabel(frame2,text="Dataframe to receive the data : "), tab.entry)
	# varnames.cbut <- tkcheckbutton(frame2,text="Variables names on the first row of data file", variable=varnames)
	# tcltk::tkgrid(varnames.cbut, columnspan=2, sticky="w")
	
	inFrame <- tcltk::tkframe(frame2, relief="groove", borderwidth=2)
	inFrame2 <- tcltk::tkframe(frame2, relief="groove", borderwidth=2)
	#inFrame3 <- tcltk::tkframe(frame2, relief="groove", borderwidth=2)

	#define new entries for input parameters
	#
	# Variables for text fields
	# Input parameters
	#PCR FRAME
	ncells<-tcltk::tclVar(5)
	varncells<-tcltk::tclVar(1)
	Pextrac <- tcltk::tclVar(0.6)
	Paliquot<-tcltk::tclVar(0.30)
	PPCR<-tcltk::tclVar(0.80)
	Tdrop<-tcltk::tclVar(2e+07)#threshold of detection for each allele signal
	Tcyc<-tcltk::tclVar(28)#of PCR cycles
	#ncells.entry <- tcltk::tkentry(inFrame, textvariable=ncells, width=4,relief="solid",justify="center")
	ncells.entry<-tcltk::tkentry(inFrame,textvariabl=ncells,width=4,relief="solid",justify="center")
	varncells.entry<-tcltk::tkentry(inFrame,textvariabl=varncells,width=4,relief="solid",justify="center")

	Pextrac.entry <-tcltk::tkentry(inFrame, textvariable=Pextrac, width=4,relief="solid",justify="center")
	Paliquot.entry<-tcltk::tkentry(inFrame, textvariable=Paliquot, width=4,relief="solid",justify="center")
	#repl.entry<-tcltk::tkentry(frame2, textvariable=repl, width=5)
	PPCR.entry<-tcltk::tkentry(inFrame2, textvariable=PPCR, width=4,relief="solid",justify="center")
	Tcyc.entry<-tcltk::tkentry(inFrame2, textvariable=Tcyc, width=4,relief="solid",justify="center")
	Tdrop.entry<-tcltk::tkentry(inFrame2, textvariable=Tdrop, width=8,relief="solid",justify="center")
	
	#tcltk::tkgrid(tcltk::tklabel(inFrame, text="________________",font=font3,foreground="blue"), columnspan=9)
	tcltk::tkgrid(tcltk::tklabel(inFrame, text="pre-PCR parameters",font=font4,foreground="blue"), columnspan=9)
	#tcltk::tkgrid(tcltk::tklabel(inFrame, text="Number of cells",font=font3,foreground="blue"), columnspan=9)
	#tcltk::tkgrid(tcltk::tklabel(inFrame,text="#Cells",font=font6), ncells.entry,sticky="we")#title varibale ncells.entry
	tcltk::tkgrid(tcltk::tklabel(inFrame,text="Mean #cells",font=font6), ncells.entry,sticky="we")#title varibale ncells.entry
	tcltk::tkgrid(tcltk::tklabel(inFrame,text="Variance #cells",font=font6), varncells.entry,sticky="we")#title varibale ncells.entry
	tcltk::tkgrid(tcltk::tklabel(inFrame,text="Prob. extraction",font=font6), Pextrac.entry,sticky="we")
	tcltk::tkgrid(tcltk::tklabel(inFrame,text="Prob. surviving aliquot",font=font6), Paliquot.entry,sticky="we")
	#tcltk::tkgrid(tcltk::tklabel(inFrame, text="________________",font=font3,foreground="blue"), columnspan=9)
	
	tcltk::tkgrid(tcltk::tklabel(inFrame2, text="PCR parameters",font=font4,foreground="blue"), columnspan=9)
	tcltk::tkgrid(tcltk::tklabel(inFrame2,text="Prob. PCR efficiency",font=font6), PPCR.entry,sticky="n")
	tcltk::tkgrid(tcltk::tklabel(inFrame2,text="#PCR cycles",font=font6), Tcyc.entry,sticky="n")#title varibale ncells.entry
	tcltk::tkgrid(tcltk::tklabel(inFrame2,text="Allele detection threshold",font=font6), Tdrop.entry,sticky="n")
	#tcltk::tkgrid(tcltk::tklabel(inFrame2, text="________________",font=font3,foreground="blue"), columnspan=9)
	
	# tcltk::tkgrid(tcltk::tklabel(inFrame3, text="Cells ploidy",font=font4,foreground="blue"), columnspan=9)
	# dip<-tcltk::tclVar(TRUE)#ploidy variable, TRUE cells are diploid
	# hap<-tcltk::tclVar(0.50)#probability of encountering allele of type A in haploid cells
	# dip.entry<-tcltk::tkentry(inFrame3, textvariable=hap, width=4,state="disabled",relief="solid",justify="center")
	
	# tcltk::tkgrid(tkcheckbutton(inFrame3, text="Diploid", variable=dip,font=font6, 
	# command=function() if (!as.logical(tclObj(dip))) tkconfigure(dip.entry, state="normal") else tkconfigure(dip.entry, state="disabled") ))
	# tcltk::tkgrid(tcltk::tklabel(inFrame3,text="prob. allele A in haploid cells",font=font6), dip.entry,sticky="n")

    repFrame<-tcltk::tkframe(frame3, relief="flat", borderwidth=2)
    repl<-tcltk::tclVar(100)
	repl.entry<-tcltk::tkentry(repFrame,textvariable=repl,width=10,highlightthickness=2,relief="solid",justify="center")
	tcltk::tkgrid(tcltk::tklabel(repFrame,text="# Replicate simulations",font=font6),repl.entry,sticky="we")
    tcltk::tkgrid(inFrame, inFrame2, padx=12,pady=8 )
	tcltk::tkgrid(repFrame,padx=8)
	tcltk::tkpack(frame2,frame3, pady=18, padx=28)#,side="top")
	
	#now, the following functions check that the user entered the correct values
	
	
	ncells.verif<-function()
	{
		
		q<-tcltk::tclvalue(ncells)
		if(q=='NULL'){return(NULL)}
		else{
			if(!is.null(q) & q<0){tcltk::tkmessageBox(message="Invalid value for the number of cells",icon="error",type="ok")}
			else{return(as.numeric(q))}
		}
	}
	
	varncells.verif<-function()
	{
		
		q<-tcltk::tclvalue(varncells)
		if(q=='NULL'){return(NULL)}
		else{
			if(!is.null(q) & q<0){tcltk::tkmessageBox(message="Invalid value for the variance of the number of cells",icon="error",type="ok")}
			else{return(as.numeric(q))}
		}
	}
	Pextrac.verif<-function()
	{
		p2<-as.numeric(tcltk::tclvalue(Pextrac))
		if(p2<0 || p2>1){tcltk::tkmessageBox(message="Invalid value for the probability of extraction",icon="error",type="ok")
		}
		else{return(p2)}
	}
	
	Paliquot.verif<-function()
	{
		p3<-as.numeric(tcltk::tclvalue(Paliquot))
		if(p3<0 || p3>1){tcltk::tkmessageBox(message="Invalid value, please enter a probability ([0,1])",icon="error",type="ok")
		}
		else{return(p3)}
	}
	
		
	PPCR.verif<-function()
	{
		p4<-as.numeric(tcltk::tclvalue(PPCR))
		if(p4<0 || p4>1){tcltk::tkmessageBox(message="Invalid value for the probability of PCR efficiency",icon="error",type="ok")
		}
		else{return(p4)}
	}
	Tcyc.verif<-function()
	{
		p5<-as.numeric(tcltk::tclvalue(Tcyc))
		if(p5<=0) {tcltk::tkmessageBox(message="Invalid value for the number of PCR cycles",icon="error",type="ok")
		}
		else{return(p5)}
	}
	repl.verif<-function()
	{
		p6<-as.numeric(tcltk::tclvalue(repl))
		if(p6<1) {tcltk::tkmessageBox(message="At least one replicate is needed",icon="error",type="ok")
		}
		else{return(p6)}
	}
	Tdrop.verif<-function()
	{
		p6<-as.numeric(tcltk::tclvalue(Tdrop))
		if(p6<10^7) {tcltk::tkmessageBox(message="Detection threshold, in number of molecules, must be at last 10^7",icon="error",type="ok")
		}
		else{return(p6)}
	}
	########### function Hb.local in tcltk 
	# main frame
	Hb.local<-function()
	{
		#first, get the parameters, check their validity
		
		ncells<-ncells.verif()
		varncells<-varncells.verif()
		#n<-ncells.verif()
		p1<-Pextrac.verif()
		p2<-Paliquot.verif()
		p3<-PPCR.verif()
		cyc<-Tcyc.verif()
		repl<-repl.verif()
		Tdrop<-Tdrop.verif()#allele detection threshold
		# an then get the results from the pcrsim function
		# all replicates are concatenated in a singla table
		
		
		#diploid case
		nsim<-as.integer(pmax(rnorm(repl,mean=ncells,sd=sqrt(varncells)),1))
		tab1<-t(sapply(nsim,function(i)simPCR2(i,probEx=p1,probAlq=p2,probPCR=p3,cyc=cyc,dip=TRUE) ))
		Hbsim<-na.omit(unlist(tab1[,5]))
		#hist(p1,col="gray",prob=T,xlab="Hb",ylab="dF(x)",cex.lab=1.3,main="Heterozygous balance, mu=1,sd=1\nProbability density function")
		
		dd <- tcltk::tktoplevel()
		frameC <- tcltk::tkframe(dd, relief="flat", borderwidth=4)
		tcltk::tkwm.title(dd,"Hb distribution")
		Myhscale <- 1    # Horizontal scaling
		Myvscale <- 1

		# Vertical scaling
		plotHb<-function()
		{
			params <- par(bg="white")
			hist(Hbsim,col="gray",prob=TRUE,xlab="Hb",ylab="f(Hb)",cex.lab=1.3,main="Heterozygote balance\n Probability density function",xlim=c(0,1))
			par(params)
		}	

		img <- tcltk::tkgrid(frameC,fun= plotHb,hscale=Myhscale,vscale=Myvscale)
		CopyToClip <- function()
		{
			tkrplot::tkrreplot(img)
		}
		copy.but <- tcltk::tkbutton(frameC,text="Copy to Clipboard",font="courrier 10",fg="darkblue",command=CopyToClip)
		tcltk::tkgrid(img)
		tcltk::tkgrid(copy.but)
			
		tcltk::tkpack(frameC)
	}
		
		
		
	

	#Bottom buttons #
	ok.but <- tcltk::tkbutton(tf, text="Simulate!", font=font5,fg="darkblue",command=function() Hb.local())
	cancel.but <- tcltk::tkbutton(tf, text="Dismiss", font=font5,fg="darkblue", command=function() tcltk::tkdestroy(tf))
	tcltk::tkpack(cancel.but, ok.but, side="left", fill="x", expand=1)

}
