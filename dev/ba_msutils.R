
#Globals



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

