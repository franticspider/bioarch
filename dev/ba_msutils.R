
#Globals


#
generatefilenames <- function(directory,sampleindexes){

	append = F
	
	if(exists("filenames"))
		rm(filenames)
	
	for(i in 1:length(sampleindexes)){
	
		name <- sprintf("%s%s.txt",directory,sampleindexes[i])
	
		if(!file.exists(name)){
			message(sprintf("File %s does not exist",name))
			return (F)
		}
			
		if(exists("filenames"))
			filenames<- c(filenames,name)
		else
			filenames<-name	
	}

	return (filenames)
}



#
validmaldidata <- function(filenames){

	if(exists("maldidata"))
		rm(maldidata)
		

	for(i in 1:length(filenames)){
		goodread = T
		e1 <- read.table(filenames[i])
		
		#It must be a data frame
		if(!is.data.frame(e1))
			goodread = F
		
		#It must have two columns
		if(ncol(e1) != 2)
			goodread = F
		
		if(!goodread){
			message(sprintf("Data in file %s is not in maldi format",filenames[i]))
			return (F)
		}
	}

	return (maldidata)
}






#convert a character to an index value for googlesheets 
ba_chartoint<-function(cc,subA=T){

	if(subA)
		mf = 64
	else
		mf = 0
  	return (as.integer(charToRaw(cc))-mf)
}



#get the column index from the column name in google sheets (e.g. ANU = 1061)
ba_colidx_from_name <- function (name,verbose=F){

	if(verbose){
		message(sprintf("Length of name %s is %d",name,length(name)))
	}
	
	i=nchar(name)
	
	idx=0;
	mypow =0;
	while(i>0){
		cc = substr(name,i,i)
		cv = ba_chartoint(cc)
		val =  ((26 ** mypow) * (cv))
		idx = idx + val
		
		if(verbose){
			message(sprintf("char %d is %s, contributing (26 ^ %d) *  %d = %d to idx, which now = %d",i,cc,mypow,cv,val,idx))
			
		}
		i = i-1
		mypow = mypow+1
	}
	return (idx)
	
}








ba_msload <- function(){



}






ms_subrange <- function(ms,lbl,ubl){
	subms <- ms[  
		ms[,1] <= ubl & ms[,1] > lbl
		,]		
}





ba_ms_align <- function(ts,data,txlim,doplot=F){

	if(doplot){
		#save the old par
    	op <- par(no.readonly = TRUE)
		#set up a new plot window to show the alignment
		dev.new()
    	#set new par
		par(mar=c(0.9,2.3,2.9,.3), mfrow = c(3,1), oma=c(5,0,2,0))
	}



	#Now we need to generate the correlations data by normalising the sampling frequency

	#STEPSIZE & MAX LAG
	myby <- 0.005 #125
	#the lagmax parameter is expressed as the number of samples (so it's scale-free)
	mylagmax <- 1/myby



	#create an interpolation (isodists is accurate to 2 decimal places)
	
	#Generate the x values
	xout = seq(from = txlim[1], to = txlim[2], by = myby)
	#Resample against xout. 
	yii <- approx(x=data[,1], y=data[,2], xout=xout, method="linear", rule = 2)
	#renormalise this segment
	yii$y = yii$y/max(yii$y)
	
	#TODO: Pass this data out so we can plot it elsewhere
	if(doplot){
		plot(yii,type="l",col="red")
	}

	#Now resample the theoretical data: 
	yri <- approx(x=ts$mass,y=ts$prob,xout=xout, method="linear", rule = 2)
	#set yvals to zero
	yri$y[] <-0
	#go through each peak
	for(i in 1:length(ts$prob)){
		idx <- which.min(abs(yri$x-ts$mass[i]))
		yri$y[idx] <- ts$prob[i]
	}
	
	#TODO: Pass this data out so we can plot it elsewhere
	if(doplot){
		lines(yri)
		title(sprintf("Data resampled to resolution %0.2f Da",myby), line = "-2")
	}

	#Now we can do the cross-correlation: 
	ccd <- ccf(yri$y,yii$y,ylim=c(-0.5,1.0),plot=doplot,axes=F, lag.max = mylagmax)
		#message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
		#plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )
		
		
	cor = ccd$acf[,,1] 
	lag = ccd$lag[,,1] 
	res = data.frame(cor,lag) 
	res_max = res[which.max(res$cor),] 
	
	out <- data.frame(
		cor = res_max$cor,
		lag = res_max$lag * myby
		)

	td <- sprintf("Cross-Correlation\n max c=%.2f, at lag=%0.3f",res_max$cor,res_max$lag*myby)
	message(td)
	
	if(doplot){
		title(td, line = "-2")
		axis(1, at = c(-50,0,50), labels = c(-0.5,0,0.5))
		axis(2)
		#text(0,-0.4,td)
	}


	#Plot the calculated peaks as probabilities 
	#until this is in the library, we have to load it:
	#source("plotseqpeaks.R")
	if(doplot){
	
		#plot the peaks from the sequence
		ba_plotseqpeaks(ts,txlim)
		#####plotseqpeaks(ocow,myxlim)
	
		message(sprintf("Lag is %0.3f",res_max$lag*myby))
	
	
	
		mycol="red"
		if(res_max$cor > 0.1)
			if(abs(res_max$lag*myby) < 0.4)
				mycol="green"
			else{
				message(sprintf("Lag too great: %0.3f",res_max$lag*myby))
			}
		else{
			message("Weak correlation - ignore")
		}
		
		lines(x=data[,1],    y = data[,2]/max(data[,2]),col="grey50")
		lines(x=data[,1]+res_max$lag*myby,y = data[,2]/max(data[,2]),col=mycol)
		
	
		
		readline("hit <return> to close the plot window and carry on")
		dev.off()
		
		#dev.set(dev.prev()) # go back to first
		#reset the par
		op
	}

	return(out)

}




ms_peaklineplot <- function(sms,ms,mycol){
	lines(sms[,1],sms[,2]/max(ms[,2]),col=mycol)
}


ms_offset_peaklineplot <- function(ms,offset,mycol){
	lines(x=(ms[,1]+offset),ms[,2]/max(ms[,2]),col=mycol)
}




ts_index <- function(sheet,spp){

	#TODO: There is surely a more efficient way of doing this...
	 found <- F
	 spidx <- 0
	 for(i in 1:nrow(sheet)){
	 	if(grepl(spp,sheet[i,1],ignore.case=TRUE)){
	 		if(!found){ 
		 		message(sprintf("FOUND: index is %d, species name is %s",i,sheet[i,1]))
		 		found <- T
		 		spidx <- i
		 	}
	 		else{
		 		message(sprintf("Further match found at  index %d, species name is %s - this will be ignored",i,sheet[i,1]))
	 		}
	 	}
	 }
	 if(!found){
	 	message("Match not found, exiting")
	 	return (-1)
	 }
	 return (spidx)

}



###########################################################################################
# The Next function, 161109, Generates the data for three replicates. Plan is to do this for 
# All five of the norfolk samples, and analyse BY HAND - see where we get and what we need
# So we build a function that helps the manual analysis first - see where that takes us...


seqvmaldi <- function(){


	#TODO: generate the ouput file names!

	#TODO: really, we want to open the file and close it - the following is a bit crude
	system("rm analysis1.dat")

	#ranked_alignment_of_mass_spectrum

	#par(mar=c(0.9,2.3,2.9,.3), mfrow = c(6,1), oma=c(5,0,2,0))
	par(mar=c(3.8,3.8,3,1), mfrow = c(1,1), oma=c(2,2,2,2))
	pdf(file = "C1.pdf",w=6,h=4)
		al <- rams(sheet,"human",e1,e2,e3,dopause=F)
	dev.off()

	readline("analysis done. Hit <return> to see results")

	#load the data from the analysis just done:
	x<-read.table("analysis1.dat",sep=",")
	colnames(x) <- c("sequence","position","length","fragno","plotno","cor","lag","pm1","pm5","energy","efrac")
	x$plotno <- as.numeric(x$plotno)

	#here's the list of plots we think are 'true'¬
	#originally
	#gt<-c(3,5,8,9,10,12,14,23,31,32,38,39,42,44,45,51,52,59,64,68,69,72,79,85,86,87)
	#removing any 'suspect' plots:
	gt<-c(3,5,8,9,10,12,14,23,31,32,38,39,42,44,45,51,52,59,64,68,69,72,79,85,86,87)

	t <- x[x$plotno %in% gt,]
	f <- x[!(x$plotno %in% gt),]

	myyint = 1

	pdf(file = "C1lagvcor.pdf",w=8,h=6)
	par(mar=c(4,3.8,3,1), mfrow = c(1,1), oma=c(1,1,1,1))
	plot(f$lag,f$cor,col="red",pch=19,main="Lag vs Correlation\nC1 sample 1 (G7)",xlab="Lag (Da)",ylab="Correlation")
	legend(x=0.25,y=0.2,c("good match","poor match"),pch = c(19,19), col = c("green","red"),y.intersp=myyint)
	points(t$lag,t$cor,col="green",pch=19)
	dev.off()

	readline("hit <return> for next plot (mass v lag)")

	pdf(file = "C1massvlag.pdf",w=8,h=6)
	plot(x=f$pm1,y=f$lag,col="red",pch=19,
	 xlab="Mass of first peak (Da)",
	 ylab="Lag (Da)",
	 main="Lag vs Correlation\nC1 sample 1 (G7)")
	 
	legend(x=2500,y=0.5,c("good match","poor match"),pch = c(19,19), col = c("green","red"),y.intersp=myyint)
	points(x=t$pm1,y=t$lag,col="green",pch=19)
	dev.off()



	readline("hit <return> to look at same for sheep")

	par(mar=c(3.8,3.8,3,1), mfrow = c(1,1), oma=c(2,2,2,2))
	pdf(file = "C1sheep.pdf",w=6,h=4)
		alsheep <- rams(sheet,"sheep",e1,e2,e3,dopause=T)
	dev.off()


	#1: get the TS
	#2: for each fragment:
	#	figure out nglut and nhyd
	#	do an alignment
	#	get a score

	##########################################################
}


load.analysis <- function(fn,gt,tf=T){
	x<-read.table(fn,sep=",")
	colnames(x) <- c("sequence","position","length","fragno",
	                 "plotno","ndeam","nhyd","cor","lag","pm1","pm5","energy","efrac")
	x$plotno <- as.numeric(x$plotno)

	if(tf){
		d <- x[x$plotno %in% gt,]
	}
	else{
		d <- x[!(x$plotno %in% gt),]
	}
	
	return (d)
}



# scode = the sample code (e.g. C1)
# gt = ground truth - list of plots that are good matches
gtanalysis <- function(scode,gt){

	fn <- sprintf("%sanalysis.dat",scode)
	message(sprintf("Opening %s",fn))


	#load the data from the analysis just done:
	x<-read.table(fn,sep=",")
	colnames(x) <- c("sequence","position","length","fragno",
	                 "plotno","ndeam","nhyd","cor","lag","pm1","pm5","energy","efrac")
	x$plotno <- as.numeric(x$plotno)

	t <- x[x$plotno %in% gt,]
	f <- x[!(x$plotno %in% gt),]


	mytlabs <- sprintf("%s/%s/%s",t$position,t$ndeam,t$nhyd)
	myflabs <- sprintf("%s/%s/%s",f$position,f$ndeam,f$nhyd)


	myyint = 1
	
	message("Created t and f from gt")

	fn <- sprintf("%slagvcor.pdf",scode)
	#str(f)
	message(sprintf("lag: %d",length(f$lag) ))
	message(sprintf("cor: %d",length(f$cor) ))
	
	
	
	pdf(file = fn ,w=8,h=6)
	par(mar=c(4,3.8,3,1), mfrow = c(1,1), oma=c(1,1,1,1))
	plot(f$lag,f$cor,col="red",pch=19,main=sprintf("Lag vs Correlation\n%s sample 1",
											scode),xlab="Lag (Da)",ylab="Correlation")
	legend(x=0.25,y=0.2,c("good match","poor match"),pch = c(19,19), 
			col = c("green","red"),y.intersp=myyint)
			
	text(x=f$lag,y=f$cor,labels=myflabs, cex=0.25, pos = 3, srt = 0, col = "grey")
	points(t$lag,t$cor,col="green",pch=19)
	text(x=t$lag,y=t$cor,labels=mytlabs, cex=0.25, pos = 3, srt = 0)
	dev.off()

#####readline("hit <return> for next plot (mass v lag)")

	#create a polynomial
	set.seed(20)
	message("Making model")
	ts <- t[ order(t$pm1), ]
	model <- lm(ts$lag ~ poly(ts$pm1,3))
	pi <- predict(model,data.frame(x=ts$pm1,interval='confidence',level=0.99))
#> plot(ts$pm1,ts$lag,col="red",pch=19)
#> lines(ts$pm1,pi)

	

	message("massvlag")

	fn <- sprintf("%smassvlag.pdf",scode)
	pdf(file = fn,w=8,h=6)
	plot(x=f$pm1,y=f$lag,col="red",pch=19,
	 xlab="Mass of first peak (Da)",
	 ylab="Lag (Da)",
	 main=sprintf("Mass vs Lag\n%s sample 1",scode)
	 ) 
	text(x=f$pm1,y=f$lag,labels=myflabs, cex=0.25, pos = 3, srt = 0, col = "grey")
	
	legend(x=2500,y=0.5,c("good match","poor match"),pch = c(19,19), col = c("green","red"),y.intersp=myyint)
	
	points(x=t$pm1,y=t$lag,col="green",pch=19)
	#Let's hav text for the 'true' points
	text(x=t$pm1,y=t$lag,labels=mytlabs, cex=0.25, pos = 3, srt = 0)
	lines(ts$pm1,pi,col='orange')
	dev.off()





	fn <- sprintf("%smassvlag_zoom.pdf",scode)
	pdf(file = fn,w=8,h=6)
	plot(x=f$pm1,y=f$lag,col="red",pch=19,
	 xlab="Mass of first peak (Da)",
	 ylab="Lag (Da)",
	 #this is the only extra line on the zoom plot!
	 ylim=c(-0.25,0.1),
	 main=sprintf("Mass vs Lag\n%s sample 1",scode)
	 ) 
	points(x=t$pm1,y=t$lag,col="green",pch=19)
	lines(ts$pm1,pi,col='orange')
	text(x=f$pm1,y=f$lag,labels=myflabs, cex=0.25, pos = 3, srt = 0, col = "grey")
	
	legend(x=2500,y=0.5,c("good match","poor match"),pch = c(19,19), col = c("green","red"),y.intersp=myyint)
	
	#Let's hav text for the 'true' points
	text(x=t$pm1,y=t$lag,labels=mytlabs, cex=0.25, pos = 3, srt = 45)
	dev.off()





	fn <- sprintf("%smassvlag_zoom_zoom.pdf",scode)
	pdf(file = fn,w=8,h=6)
	plot(x=f$pm1,y=f$lag,col="red",pch=19,
	 xlab="Mass of first peak (Da)",
	 ylab="Lag (Da)",
	 #this is the only extra line on the zoom plot!
	 ylim=c(-0.1,0.1),
	 xlim=c(1425,1900),
	 main=sprintf("Mass vs Lag\n%s sample 1",scode)
	 ) 
	
	points(x=t$pm1,y=t$lag,col="green",pch=19)
	lines(ts$pm1,pi,col='orange')
	text(x=f$pm1,y=f$lag,labels=myflabs, cex=0.25, pos = 3, srt = 45, col = "grey")
	
	legend(x=2500,y=0.5,c("good match","poor match"),pch = c(19,19), col = c("green","red"),y.intersp=myyint)
	
	points(x=t$pm1,y=t$lag,col="green",pch=19)
	lines(ts$pm1,pi,col='orange')
	#Let's hav text for the 'true' points
	text(x=t$pm1,y=t$lag,labels=mytlabs, cex=0.25, pos = 3, srt = 45)
	dev.off()



#####readline("hit <return> to see 

	fn <- sprintf("%sseqvionct.pdf",scode)
	pdf(file = fn,w=8,h=6)
	plot(x=t$position,y=t$energy,col="red",pch=19,
	 xlab="Position on sequence",
	 ylab="Ion Count",
	 main=sprintf("Seqpos vs Ion Count\n%s sample 1",scode)
	 )
	text(x=t$position,y=t$energy,labels=mytlabs, cex=0.7, pos = 3)
	dev.off()

	#par(mar=c(3.8,3.8,3,1), mfrow = c(1,1), oma=c(2,2,2,2))
	#pdf(file = "C1sheep.pdf",w=6,h=4)
	#	alsheep <- rams(sheet,"sheep",e1,e2,e3,dopause=T)
	#dev.off()



}






