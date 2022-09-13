
#________________________________	
#Two-allele model simulation
#________________________________
"A2.simu" <- function()
{	
	#if(!require(tcltk)) stop("package tcltk is required")
	tcltk::tclRequire("Tktable")
	#fonts definition
	font0 <- tcltk::tkfont.create(family="times",size=35,weight="bold",slant="italic")
	font1<-tcltk::tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
	font2<-tcltk::tkfont.create(family="times",size=16,weight="bold",slant="italic")
	font3<-tcltk::tkfont.create(family="times",size=10,weight="bold")#,slant="italic")
	font4<-tcltk::tkfont.create(family="times",size=12)#,slant="italic")
	#Two-allele model implementation
	twoAmod<- function(Mx,peak1,peak2)
	{
		expA <- c(0.5+ Mx/2, 0.5,Mx,1-(Mx/2),1-Mx,Mx/2,(1-Mx)/2)
		expB <- c((1-Mx)/2, 0.5,1-Mx,Mx/2,Mx,1-Mx/2,0.5+Mx/2)
		obsA <- peak1/(peak1+peak2)
		obsB <- peak2/(peak1+peak2)
		resid <- (expA-obsA)^2+(expB-obsB)^2
		genotypes<-c('AA,AB','AB,AB','AA,BB','AB,AA','BB,AA','AB,BB','BB,AB')
		Mx.Conditioned<-abs(c((peak1-peak2)/(peak1+peak2),NA,peak1/(peak1+peak2),
		2*peak2/(peak1+peak2),peak2/(peak1+peak2),2*peak1/(peak1+peak2),(peak2-peak1)/(peak1+peak2)))
		#to be completed
		#Hb1<-
		#Hb2<-
		result<-data.frame(genotypes,expA,expB,obsA,obsB,resid,Mx.Conditioned)
		return(result)
	}
	
	#Two-allele model: graphical simulation
	twoAmodG<-function(obsMx,peak1,peak2)
	{
		Mx<-seq(0.1,0.9,by=0.05)
		res<-sapply(Mx, function(i) twoAmod(i,peak1,peak2)$resid)
		par(mar = c(4,4,2,3)+0.1,oma = c(1,2,0.3,5))
		#log="y"
		plot(Mx,(res[1,]),type='n',lab=c(6,5,1),xlab="Mx (mixture proportion)",ylab="Residuals",cex.lab=1.2,las=1,ylim=c(min(res),max(res)))
		title("Two-allele model simulations")
		pchref<-c(18,15,17,4,8,16,3)
		col<-c("lightgreen","magenta","purple","cyan","red","darkblue","gold")
		for(i in 1:7)#7 genotype combinations
		{
			points(Mx,(res[i,]),pch=pchref[i],col=col[i])
			lines(Mx,(res[i,]),col=col[i])
		}
		#print(obsMx)
		abline(v=obsMx,col="gray",lty=2)
		#text(obsMx,median(res),"Observed Mx",   adj = c(0, -.1))
		legend(0.95,max(res),c('AA,AB','AB,AB','AA,BB','AB,AA','BB,AA','AB,BB','BB,AB','Obs. Mx'),
		,pch =c(pchref,32),lty=2,col=c(col,'gray'),bg='white',xpd=NA,cex=1.2,
		 box.lwd=2)
	}
		
	#Two-allele model slikelihood results
	twoAmodT<-function(peak1,peak2)
	{
		Mx<-seq(0.1,0.9,by=0.05)
		geno<-c('AA,AB','AB,AB','AA,BB','AB,AA','BB,AA','AB,BB','BB,AB')
		res<-sapply(Mx, function(i) signif(twoAmod(i,peak1,peak2)$resid,digits=3))
		row.names(res)<-c(geno)
		colnames(res)<-paste('Mx=',Mx,sep='')
		minX<-which(res==min(res),arr.ind=TRUE)
		minC<-signif(Mx[minX[,2]],digits=2)
		minY<-row.names(minX)	
		return(list(minY,minC,res))
		#return(res)
	}
	# Two-allele model main frame
	tt <- tcltk::tktoplevel()
	tcltk::tkwm.title(tt,"Two-allele model simulations")
	frame1 <- tcltk::tkframe(tt, relief="groove", borderwidth=2)
	frame2 <- tcltk::tkframe(tt, relief="groove", borderwidth=2)
 
	xyframe <- tcltk::tkframe(frame1, relief="groove", borderwidth=2)
	labframe <- tcltk::tkframe(frame1, relief="groove", borderwidth=2)
	limframe <- tcltk::tkframe(frame2, relief="groove", borderwidth=2)
	posframe <- tcltk::tkframe(frame2, relief="groove", borderwidth=2)
	legframe <- tcltk::tkframe(frame2, relief="groove", borderwidth=2)

	#
	# Variables for text fields
	#
	xy1var <- tcltk::tclVar(899)
	xy2var <- tcltk::tclVar(2183)
	mixvar <- tcltk::tclVar(0.70)
	
		

	# Title
	#
	TFrame <- tcltk::tkframe(tt, relief="groove")
	labh <- tcltk::tklabel(TFrame)
	tcltk::tkgrid(tcltk::tklabel(TFrame,text="Two-allele model", font=font2, foreground="red"), labh)
	
	tcltk::tkpack(TFrame)
	
	xy1.entry <- tcltk::tkentry(xyframe, textvariable=xy1var, width=8)
	xy2.entry <- tcltk::tkentry(xyframe, textvariable=xy2var, width=8)
		
	peak1<-function()
	{
		p1<-as.numeric(tcltk::tclvalue(xy1var))
		if(p1<0){tcltk::tkmessageBox(message="Invalid value for the peak height of allele #1",icon="error",type="ok")}
		else{return(p1)}
	}
		
		
		
	peak2<-function()
	{
		p2<-as.numeric(tcltk::tclvalue(xy2var))
		if(p2<0){tcltk::tkmessageBox(message="Invalid value for the peak height of allele #2",icon="error",type="ok")
		}
		else{return(p2)}
	}
	
		
	# choosexy1.but <- tcltk::tkbutton(xyframe, text="Enter", command=function()peak1())
	# choosexy2.but <- tcltk::tkbutton(xyframe, text="Enter", command=function() peak2())
	tcltk::tkgrid(tcltk::tklabel(xyframe, text="- Peak heights (rfu) -",font=font3, foreground="blue"), columnspan=5)
	tcltk::tkgrid(tcltk::tklabel(xyframe,text="Allele #1: "), xy1.entry)#, choosexy1.but)
	tcltk::tkgrid(tcltk::tklabel(xyframe,text="Allele #2: "), xy2.entry)#, choosexy2.but)
	#tcltk::tkgrid(tcltk::tklabel(xyframe,text="X axis col. #: "), nx.entry)
	#tcltk::tkgrid(tcltk::tklabel(xyframe,text="Y axis col. #: "), ny.entry)
	#
	# Labels frame
	#
	lab.entry <- tcltk::tkentry(labframe, textvariable=mixvar, width=8)
	prop<-function() 
	{
		Mx<-as.numeric(tcltk::tclvalue(mixvar))
		if(Mx>1 || Mx<0)
		{
			tcltk::tkmessageBox(message="Mx is the mixture proportion, it must be comprised in the interval [0,1]",
				icon="error",type="ok")

		}
		
		else{return(Mx)}
			#print(Mx)
			
	}
	#chooselab.but <- tcltk::tkbutton(labframe, text="Set", command= function() prop())#function() print(tcltk::tclvalue(mixvar)))
	tcltk::tkgrid(tcltk::tklabel(labframe, text="- Mixture proportion  -",font=font3, foreground="blue"), columnspan=3)
	#lab.entry pour dire que c'est l'entree lab.entry, et le label est celui de "set"
	tcltk::tkgrid(tcltk::tklabel(labframe,text="Mx : "), lab.entry)#, chooselab.but)
	#tcltk::tkgrid(tcltk::tklabel(labframe,text="1-Mx : "), lab2.entry)
	tcltk::tkpack(xyframe, labframe, side="left")
	tcltk::tkpack(frame1)

	# 
	RCSFrame <- tcltk::tkframe(tt, relief="groove")
	plotFunction<-function()
	{
		p1<-peak1()
		p2<-peak2()
		obsMx<-prop()
		twoAmodG(obsMx,p1,p2)
	}

	propFunction<-function()
	{
		
		p1<-peak1()
		p2<-peak2()
		res<-twoAmodT(p1,p2)
		resdata<-res[[3]]
		resgeno<-res[[1]]
		resmix<-res[[2]]
		
		myRarray<- c("Genotype","Mx=0.1","Mx=0.15","Mx=0.20","Mx=0.25","Mx=0.30",
		"Mx=0.35","Mx=0.40","Mx=0.45","Mx=0.50","Mx=0.55","Mx=0.60","Mx=0.65",
		"Mx=0.70","Mx=0.75","Mx=0.80","Mx=0.85","Mx=0.90",
		"AA,AB",resdata[1,],"AB,AB",resdata[2,],"AA,BB",resdata[3,],"AB,AA",resdata[4,],"BB,AA",resdata[5,],"AB,BB",
		resdata[6,],"BB,AB",resdata[7,])

		dim(myRarray) <- c(18,8)
		tclarray <- tcltk::tclArray()
		tclarray2 <- tcltk::tclArray()
		tclarray3 <- tcltk::tclArray()

		for(i in 0:17){
		for (j in (0:7)){
		  tclarray[[i,j]] <- myRarray[i+1,j+1]}}

		myRarray2<-resgeno
		dim(myRarray2)<-c(1,length(resgeno))
		myRarray3<-resmix
		dim(myRarray3)<-c(1,length(resmix))
		#myRarray3

		tclarray2 <- tcltk::tclArray()
		for(i in 0:(length(resgeno)-1))
		{
			tclarray2[[0,i]] <- myRarray2[1,i+1]
		}
				
		tclarray3 <- tcltk::tclArray()
		
		for(i in 0:(length(resmix)-1))
		{
			tclarray3[[0,i]] <- myRarray3[1,i+1]
		}
				
		
		
		save1<-function(filename1="simulation2.txt",filename2="likelihood2.txt")
		{
			#myRarray
			write.table(myRarray,file=filename1,row.names=FALSE,col.names=FALSE)
			loctab<-data.frame(resmix,resgeno)
			colnames(loctab)<-c('Mixture proportion','Genotype combination')
			#writeLines('\n',filename)
			write.table(loctab,file=filename2,row.names=FALSE)
			#write.table(myRarray3,"filter.txt",row.names=FALSE,col.names=FALSE)
		}
		saveFunction<-function()
		{
			
			ss<-tcltk::tktoplevel()
			SSframe <- tcltk::tkframe(ss, relief="groove",width=35)
			tcltk::tkwm.title(ss,"")
			filevar1 <- tcltk::tclVar("simulation2.txt")
			filevar1.entry <- tcltk::tkentry(SSframe, textvariable=filevar1, width=12)
			filevar2 <- tcltk::tclVar("likelihood2.txt")
			filevar2.entry <- tcltk::tkentry(SSframe, textvariable=filevar2, width=12)
			#filevar2.entry <- tcltk::tkentry(SSframe, textvariable=filevar, width=12)
			#tcltk::tkgrid(tcltk::tklabel(SSframe, text="- Enter filenames -",font=font1, foreground="blue"), columnspan=15)
			save1.butt<-tcltk::tkbutton(ss, text="Enter", font=font3,command=function() save1(tcltk::tclvalue(filevar1),tcltk::tclvalue(filevar2)))
			tcltk::tkgrid(tcltk::tklabel(SSframe,text="Simulations results",font=font4), filevar1.entry)
			tcltk::tkgrid(tcltk::tklabel(SSframe,text="Maximum likelihood",font=font4), filevar2.entry)
			tcltk::tkgrid(filevar2.entry, save1.butt)		
			tcltk::tkpack(SSframe)
		}
		
		
		tt<-tcltk::tktoplevel()
		tcltk::tkwm.title(tt,"Most likely genotypes combination")
		table1 <- tcltk::tkwidget(tt,"table",variable=tclarray,rows=18,colwidth=8,cols=8,titlerows=1,background="white")
		table2 <- tcltk::tkwidget(tt,"table",variable=tclarray2,
		cols=length(resgeno),selectmode="extended",colwidth=10,rows=1,background="lightblue")
		table3 <- tcltk::tkwidget(tt,"table",variable=tclarray3,
		cols=length(resmix),selectmode="extended",colwidth=10,rows=1,background="lightblue")
		tit1<-tcltk::tkwidget(tt,"label",text="Matrix of the residuals",font=font1,foreground="blue")	
		tit2<-tcltk::tkwidget(tt,"label",text="Maximum likelihood estimation results",font=font1,foreground="blue")
		tit3<-tcltk::tkwidget(tt,"label",text="Most likely genotype combinations",font=font1,foreground="blue")
		tit4<-tcltk::tkwidget(tt,"label",text="Corresponding mixture proportions",font=font1,foreground="blue")
		#tcltk::tkgrid(tcltk::tklabel(tt,text="File name"), filevar.entry, save1.butt)
		filelab<-tcltk::tkwidget(tt,"label",text="-Save the results-",font=font3,foreground="blue")
		save.butt<-tcltk::tkbutton(tt, text="Save", font=font3,command=saveFunction)
		tcltk::tkpack(tit1,table1,tit3,table2,tit4,table3,save.butt)
		
	#tcltk::tkpack(filevar.entry,save1.butt,side="left")
	}


	filterFunction<-function()
	{
		
		p1<-peak1()
		p2<-peak2()
		obsMx<-prop()
		tmp1<-twoAmod(obsMx,p1,p2)
		Mx.C<-signif(tmp1$Mx.Conditioned,digits=2)
		genotypes<-as.character(tmp1$genotypes)
		tab1<-c("Genotype",genotypes,"Mx conditioned",Mx.C)
		dim(tab1)<-c(8,2)
		tab1array <- tcltk::tclArray()

		for(m in 0:7){
		for (n in (0:1)){
		  tab1array[[m,n]] <- tab1[m+1,n+1]}}
		
		
		saveHH<-function(name1="filter2.txt")
		{
			#myRarray
			write.table(tab1,name1,row.names=FALSE,col.names=FALSE)
				
		}		
			
		saveFunction2<-function()
		{
			
			hh<-tcltk::tktoplevel()
			HHframe<- tcltk::tkframe(hh, relief="groove")
			tcltk::tkwm.title(hh,"Filenames")
			filtervar<- tcltk::tclVar("filter2.txt")
			filtervar.entry <- tcltk::tkentry(HHframe, textvariable=filtervar, width=12)
			#filevar2.entry <- tcltk::tkentry(SSframe, textvariable=filevar, width=12)
			#tcltk::tkgrid(tcltk::tklabel(HHframe, text="- Enter filename -",font=font1, foreground="blue"), columnspan=15)
			saveHH.butt<-tcltk::tkbutton(hh, text="Enter", font=font3,command=function() saveHH(tcltk::tclvalue(filtervar)))
			tcltk::tkgrid(tcltk::tklabel(HHframe,text="Genotypes filter",font=font4), filtervar.entry)
			#tcltk::tkgrid(tcltk::tklabel(SSframe,text="Maximum likelihood",font=font4), filevar2.entry)
			tcltk::tkgrid(filtervar.entry, saveHH.butt)		
			tcltk::tkpack(HHframe)
		}
		
		
		
		tt2<-tcltk::tktoplevel()
		tcltk::tkwm.title(tt2,"Genotypes filter")
		#tab1.tit<-tcltk::tkwidget(tt2,"label",text="Genotypes filter",font=font1,foreground="blue")
		save2.butt<-tcltk::tkbutton(tt2, text="Save", font=font3,command=saveFunction2)
		tab1.tcl<-tcltk::tkwidget(tt2,"table",variable=tab1array,rows=8,colwidth=18,cols=2,titlerows=1,background="white")
		tcltk::tkpack(tab1.tcl,save2.butt)
		#tcltk::tkpack(tab1.tit,tab1.tcl,save2.butt)
		
		
		
			
	}

	A1.but <- tcltk::tkbutton(RCSFrame, text="Plot simulations",font=font3, command=plotFunction)#twoAmod())
	A2.but <- tcltk::tkbutton(RCSFrame, text="Simulation details", font=font3,command=propFunction)
	A3.but <- tcltk::tkbutton(RCSFrame, text="Genotype filter", font=font3,command=filterFunction)
	tcltk::tkgrid(A1.but,A2.but,A3.but,ipadx=20)	
	tcltk::tkpack(RCSFrame)
}


