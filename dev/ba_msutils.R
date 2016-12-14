
#Globals




ba_peptides_from_sequence <- function (sheet,spp="human",verbose=T,moff = 1){


	 spidx<-ts_index(sheet,spp)


	 endcol<-ncol(sheet)	
	 	
	 message("Calculating sequences now")
	 start<-4
	 count<-1
	 
	 #seqset <- data.frame(seq=character(),mass=vector(),prob=vector())
	 
	 gotdata = F
	 
	 for(j in start:endcol){
 		#TODO: reject non a-a entries
 		if(grepl("K|R",sheet[spidx,j])){
 		
 			count<-count+1
 			message(sprintf("%s\n%d:\t",sheet[spidx,j],count),appendLF=F)
 			end=j
 			#to make the string, do something like:
 			#cc <- paste0(shit[85,4:12],collapse="")

			sequence <- paste0(sheet[spidx,start:end],collapse="")
			
			nhyd <- str_count(sequence,"P")
			
			nglut <- str_count(sequence,"Q")
			
			message(sprintf("Generating spectrum for sequence %s ...",sequence,appendLF=F))

			cd1 <- q2e_tpeaks(sequence)
			lbl <- min(cd1$mass) - moff
			ubl <- max(cd1$mass) + moff

			if(!gotdata){
				seqset <- data.frame(seq=as.character(cd1$infile),start=start, mass = t(cd1$mass), prob = t(cd1$prob))
				gotdata = T
			}
			else{
				nr <- data.frame(seq=as.character(cd1$infile),start=start, mass = t(cd1$mass), prob = t(cd1$prob))
				seqset <- rbind(
				  seqset,nr)
			}
			
			#readline("hit <return> to continue...\n")
			start = j+1
			
			
		}
		else{
			message(sprintf("%s",sheet[spidx,j]),appendLF=F)
		}
	}
	 	
	#tell R that the 'seq' column is strings:
	seqset$seq <- as.character(seqset$seq)

	return(seqset)
}











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

#Get the column name from the column index in googlesheets
ba_name_from_colidx <- function(idx,verbose=F){

	


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
	#mylagmax gives the 'reverse scaling' of the stepsize - useful when comparing etc. 
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

	######################################################################################
	#Now we can do the cross-correlation:                                 4*mylagmax
	ccd <- ccf(yri$y,yii$y,ylim=c(-0.1,0.5),plot=doplot,axes=F, lag.max = mylagmax)
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
	
	#let's set a limit on the permissable lag: 
	laglim = 0.3
	
	lag_range = c(-laglim,laglim)
	lagset = res[ (res$lag > -(laglim/myby)) & (res$lag < (laglim/myby)),]
	
	res_lrmax = lagset[which.max(lagset$cor),]
	
	message("\nComparing max correlation with within-range correlation:")
	message(sprintf("  cor = %0.2f, lag = %0.2f\nlrcor = %0.2f, lrlag = %0.2f\n",out$cor,out$lag,res_lrmax$cor,res_lrmax$lag * myby))
	
	
	
	
	
	
	
	#Here's where we can do a more detailed analysis
	message(sprintf("max cor = %0.2f at lag %0.2f",out$cor, out$lag))
	if(doplot){
		points(x=res_max$lag,y=res_max$cor,pch=19,col="red")
		points(x=res_lrmax$lag,y=res_lrmax$cor,pch=10,col="green",cex = 3)
	}
	
	message(sprintf("There are %d points in the correlation from %0.2f to %0.2f (scaled to %0.2f to %0.2f)",nrow(res),res$lag[1],res$lag[nrow(res)],res$lag[1]*myby,res$lag[nrow(res)]*myby))
	
	
	
	
	
	

	td <- sprintf("Cross-Correlation\n max c=%.2f, at lag=%0.3f\n max inrange c = %.2f at lag %.3f",res_max$cor,res_max$lag*myby,res_lrmax$cor,res_lrmax$lag * myby)
	
	message(td)
	#readline(sprintf("Hit <return> for %s",td))
	
	
	labelvals = c(-1,-laglim,0,laglim,1)
	
	if(doplot){
		title(td, line = "-2")
		#axis(1, at = c(-mylagmax/2,0,mylagmax/2), labels = c(-0.5,0,0.5))
		axis(1, at = mylagmax * labelvals, labels = labelvals)
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
		plot.new()
		#dev.set(dev.prev()) # go back to first
		#reset the par
		#op
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
		 		message(sprintf("FOUND: index is %d, search term is \"%s\"",i,sheet[i,1]))
		 		found <- T
		 		spidx <- i
		 	}
	 		else{
		 		message(sprintf("Further match found at  index %d, search term is \"%s\" - this entry will be ignored",i,sheet[i,1]))
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

	#TODO: This is the wrong structure! use load.analysis()
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



load.analysis <- function(fn){

	x<-read.table(fn,sep=",",, stringsAsFactors=FALSE)
	colnames(x) <- c("sequence","position","length","fragno",
	                 "plotno","ndeam","nhyd","cor","lag","pm1","pm5","energy","efrac")
	x$plotno <- as.numeric(x$plotno)
	
	return(x)
}


#Used to be just 'load.analysis
load.analysis.tf <- function(fn,gt,tf=T){

	x <- load.analysis(fn)

	if(tf){
		d <- x[x$plotno %in% gt,]
	}
	else{
		d <- x[!(x$plotno %in% gt),]
	}
	
	return (d)
}


#Load the hydroxylation probabilities
load.phyd <- function(fn){

 	#h <- read.table("N1_hprobs.dat",fill=T,sep=",")
 	h <- read.table(fn,fill=T,sep=",")
 	
	#h[,2:ncol(h)] <- as.numeric(h[,2:ncol(h)])
	h[,1] <- as.character(h[,1])
	
	#get rid of the nas
	h[is.na(h)]<-0
	
	
	#Name the columns
	nc <- ncol(h)-5
	colnames(h) <- c("sequence","start","len","nh",paste("nh",0:nc,sep=""))


	
	return(h)

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















ms_fit <- function(sheet,species,froot,scode,spots,doyn=T){

	#Pick the 'ground truth' as small lag and high correlation
	cor.threshold = 0.1
	lag.threshold = 0.2	
													# b  l  t  r	
	par(mar=c(3.8,3.8,2.4,3.8), mfrow = c(1,1), oma=c(.1,.1,.2,.2))
	layout(matrix(c(1,1,
	                2,2,
	                3,3,
	                4,5
	                ), 4, 2, byrow = TRUE))
	#                5,5), 4, 2, byrow = TRUE))


	x <- load.analysis(sprintf("%sanalysis.dat",scode))
	
	t <- x[which(x$cor > cor.threshold & abs(x$lag) < lag.threshold),]
	f <- x[which(x$cor < cor.threshold | abs(x$lag) > lag.threshold),]
	
	xseqs  <- x[!duplicated(x$sequence),]

	#TODO: We need to get this data from somewhere more reliable!
	#hp <- load.phyd(sprintf("hprobs.dat",scode))
	
	
	message("Loading master plot")
	eraw <- read.table(sprintf("%s%s.txt",froot,spots[1]))
	#Make a copy that we can manipulate freely
	e1<-eraw
	
	message("Loading hydroxylation data")	
	#hprobs <- load.phyd(sprintf("%s%s_hprobs.dat",froot,spots[1]))
	hp <- load.phyd("hprobs.dat")
	
	#readline(sprintf("hp$sequence[1] is %s",hp$sequence[1]))
	
	
	
	#create the matched peak set array
	matchedpeaks <- e1
	matchedpeaks[,2]=0
	
	gxlim = c(800,3200)		

	finished = F
	
	
	
	#create a polynomial
	set.seed(20)
	message("Making model")
	ts <- t[ order(t$pm1), ]
	if(nrow(ts)<5){
		message(sprintf("ts has %d points - less than polynomial degree!",nrow(ts)))
		}
	else{
		model <- lm(ts$lag ~ poly(ts$pm1,3))
		pi <- predict(model,data.frame(x=ts$pm1,interval='confidence',level=0.99))
		}
	
	#We need this scale because we are going to remove peaks from e1
	ms_yscale <-max(e1[,2])
	
	
	
	#ok now lets go through each peak in the ms
	while(!finished){

		if(nrow(t)==0){
			finished = T
		}
		else{
			message(sprintf("%d theoretical peaksets remaining",nrow(t)))

			cca <- t[which(t$energy == max(t$energy)),]
		
			message(sprintf("Found %d peaks with highest ion count",nrow(cca)))
		
			#Go through each result with this score...
			for(ii in 1:nrow(cca)){
		
				cc <- cca[ii,]
			
				message(sprintf("Current highest peak has range %0.2f..%0.2f",cc$pm1,cc$pm5))
		
				#ff =  ALL peaks for this sequence
				ff <- t[which(t$position == cc$position),]
		
				#remove cc from t
				t <- t[-which(t$energy == max(t$energy)),]
		
				#To check we've removed only one peak: 
				#message(sprintf("           %d  left....",nrow(t)))
		
				#Get the 'raw' peaks from the sequence
				# - we'll add deamidations and hydroxylations later: 
			
				message(sprintf("cc seqeunce is %s",cc$sequence))
				#readline("Press<return> to get the peak data from the sequence")
				cd1 <- q2e_tpeaks(cc$sequence)
				#message("Counting the number of hydroxylations:")
				nhyd <- str_count(cc$sequence,"P")
				nglut <- str_count(cc$sequence,"Q") + str_count(cc$sequence,"N")

				#TODO: fixed - but why was I using seq anyway? 
				#Get the hydroxylation probabilities
				#shp <- hp[which(hp$sequence == seq),]
				shp <- hp[which(hp$sequence == cc$sequence),]


				# Calculating ccxlim: originally, we only used the peaksets that were 'correct'
				# but it may be better to see the *whole* range of potential changes to mass
				#ccxlim = c(min(ff$pm1)-5,max(ff$pm1)+10)
		
				ccxlim = c(cd1$mass[1]-5, cd1$mass[1] + (nglut*0.984015)+(nhyd*16) + 10)




				#TODO: We are assuming an analyis has been done and saved in C1analysis.dat!
				#PLOT 1: THE RAW MASS SPECTRUM
				plot(e1,type="l",xlab="Mass (Da)",ylab="Ion Count",xlim=gxlim,mgp=c(.6,0.6,.6),
					ylim = c(0,ms_yscale),
					main = sprintf("Global view. Target is for peptide %s, nDeam = %d, nHyd = %d",
						cc$sequence,cc$ndeam,cc$nhyd
						)
				)	
				lines(matchedpeaks,col="red",xlim=gxlim)
				#Mass V lag - true only - overlaid on the mass plot
				par(new=T)
				plot(t$pm1+t$lag,t$lag,col="green",pch=20,
						xlim=gxlim,
						xaxt="n",yaxt="n",xlab="",ylab="cor",ylim=c(-0.25,0.25))	
						
				if(nrow(ts)>5){	
					lines(ts$pm1,pi,col='orange')
					}
				mtext("Cor",side=4,line=3)			
				points(ff$pm1,ff$lag,col = "blue",pch=10)
				points(cc$pm1,cc$lag,col = "red",pch=10,cex=3)
		









			
			
				#PLOT 2: Show all candidate peaks..
				#calculate the range from the hydroxylations
				mass_min <- cc$pm1
				mass_max <- cc$pm1+10
				nseqs=0;
				for(xrow in 1:nrow(xseqs)){
					c_nhyd  <- str_count(xseqs$sequence[xrow],"P")
					c_nglut <- str_count(xseqs$sequence[xrow],"Q") + str_count(xseqs$sequence[xrow],"N")
					message(sprintf("Sequence %s has %d hydroxylations and %d deamitation candidates",
							xseqs$sequence[xrow],c_nhyd,c_nglut))
					thiscc <-  q2e_tpeaks(cc$sequence)
					#This next doesn't work because 
					#message(sprintf("thiscc has %d rows",nrow(thiscc)))
				
					tmin <- thiscc$mass[1] 
					tmax <- thiscc$mass[1] +(c_nglut*0.984015)+(c_nhyd*16) + 10
					#if(exists("cc$pm1") & exists("tmin") & exists("tmax")){
						if(tmin < cc$pm1 & 
						   tmax > cc$pm1){
						   	mass_min <- tmin  
						   	mass_max <- tmax
						   	nseqs <- nseqs + 1
						   	message(sprintf("DEBUG: mass_min is %f, mass_max is %f",mass_min,mass_max))
						}
					#}
					#else{
					#	message(sprintf("one of tmax,tmin,or cc$pm1 doesn't exist"))
					#	#finished = T
					#	break
					#}
				}
				#	
				plot(NA,xlab="Mass (Da)",ylab="Ion Count",xlim=c(mass_min,mass_max),mgp=c(.3,0.3,.3),ylim=c(0,1),
					main = sprintf("ALL Theoretical spectra around the target peak (random colours)")
					)
			
				#readline(sprintf("nseqs is %d",nseqs))
				#colours <- rainbow(nseqs)
				cccolor = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
				colours=sample(cccolor, nseqs)
			
				cseqs=0 
				for(xrow in 1:nrow(xseqs)){
					message(sprintf("DEBUG - iteration %d ",xrow))
					thiscc <-  q2e_tpeaks(xseqs$sequence[xrow])
					c_nhyd  <- str_count(xseqs$sequence[xrow],"P")
					c_nglut <- str_count(xseqs$sequence[xrow],"Q") + str_count(xseqs$sequence[xrow],"N")
				
					tmin <- thiscc$mass[1] 
					tmax <- thiscc$mass[1] +(c_nglut*0.984015)+(c_nhyd*16) + 10
					if(tmin < cc$pm1 & 
					   tmax > cc$pm1){
				
						message(sprintf("DEBUG: seq = %s ,x$pm1 = %0.3f",#, thiscc$mass[1] = %0.3f",
										xseqs$sequence[xrow]
										, xseqs$pm1[xrow] #-(x$nglut[xx]*0.984015) -(x$nhyd[xx]*16),
										#thiscc$mass[1])
										))
								
						cseqs <- cseqs+1
						ccol <- colours[cseqs]
						#readline(sprintf("Cseqs is %d",cseqs))
						for(e in 0:c_nglut){
							for( p in 0:c_nhyd){
							
								xx <- thiscc$mass +(e*0.984015)+(p*16)
								yy <- thiscc$prob * as.numeric(shp[5+p])
								message(sprintf("xrow[1] is %0.2f",xrow[1]))
								segments(x0=xx, y0=yy, y1=0, col=ccol)
								points(xx, yy, pch=21, col=ccol, bg=ccol)
							}
						}
					}
				}
				lines(x=e1[,1],y=e1[,2]/ms_yscale)
				lines(x=matchedpeaks[,1],y=matchedpeaks[,2]/ms_yscale,col="red")
	
				points(ff$pm1,ff$lag+0.5,col = "blue",pch=10,cex=2)	
				points(cc$pm1,cc$lag+0.5,col = "red",pch=10,cex=3)
		
				#DEBUG:readline("2nd plot should now be visible")
		
				#PLOT 3: Now zoom in to the target
				#calculate the range from the hydroxylations
				#	
				plot(NA,type="l",xlab="Mass (Da)",ylab="Ion Count",xlim=ccxlim,mgp=c(.3,0.3,.3),ylim=c(0,1),
					main = sprintf("Aligned view. Target is for peptide %s, nDeam = %d, nHyd = %d, peak at %0.1f\n significant peaksets are pink (others grey)",
						cc$sequence,cc$ndeam,cc$nhyd,cc$pm1
						))
				#First draw all possible theoretical peaks...
				for(e in 0:nglut){
					for( p in 0:nhyd){
						xx <- cd1$mass +(e*0.984015)+(p*16)
						yy <- cd1$prob  * (as.numeric(shp[5+p])/max(shp[5:ncol(shp)]))
						segments(x0=xx, y0=yy, y1=0, col="grey80")
						points(xx, yy, pch=21, col="grey80", bg="grey80")
					}
				}
				#NOW Draw the GOOD theoretical peak matches:
				message(sprintf("nrow(ff) = %d",nrow(ff)))
				if(nrow(ff)>0){
					for( i in 1:nrow(ff)){
						xx <- cd1$mass +(ff$ndeam[i]*0.984015)+(ff$nhyd[i]*16)
						yy <- cd1$prob #* (as.numeric(shp[5+ff$nhyd])/max(shp[5:ncol(shp)]))
						segments(x0=xx, y0=yy, y1=0, col="pink")
						points(xx, yy, pch=21, col="pink", bg="pink")
					}
				}
				lines(x=e1[,1],y=e1[,2]/ms_yscale)
				lines(x=matchedpeaks[,1],y=matchedpeaks[,2]/ms_yscale,col="red",xlim=gxlim)
				par(new=T)
				plot(t$pm1,t$lag,col="green",pch=20,
						xlim=ccxlim,
						xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(-0.25,0.25))	
				if(nrow(ts)>5){	
					lines(ts$pm1,pi,col='orange')
					}
				points(ff$pm1,ff$lag,col = "blue",pch=10,cex=2)	
				points(cc$pm1,cc$lag,col = "red",pch=10,cex=3)
		
	










		
		
	
				#PLOT Lag v cor - true AND false
				plot(f$lag,f$cor,xlab="Lag",ylab="cor",col="red",pch=20,
						main="Lag vs Correlation\nall seqences",
						xlim=c(-1,1),ylim=c(0,0.3),
						#xlim=c(min(x$lag),max(x$lag)),
						#ylim=c(min(x$cor),max(x$cor))
						)		
				
				points(t$lag,t$cor,col="green",pch=20)
		
				points(ff$lag,ff$cor,col = "blue",pch=10,cex=2)
				points(cc$lag,cc$cor,col = "red",pch=10,cex=3)





				#PLOT Lag V cor - true only
				plot(t$lag,t$cor,xlab="Lag",ylab="cor",col="green",pch=20,
						xlim=c(-0.3,0.3),ylim=c(0.1,0.3),
						main="Lag vs Correlation\ngood candidates",
						)	
				points(ff$lag,ff$cor,col = "blue",pch=10,cex=2)
				points(cc$lag,cc$cor,col = "red",pch=10,cex=3)
		
		
				answer = "q"
				if(doyn){
					answer <- ba_ynq("Does this theoretical peak match the data?")
				}
	
		
				if(answer == "Y" || answer == "y" || doyn ==F){
		
					#Move the MS intensities to the 
					lbl <- min( cd1$mass) + (cc$ndeam*0.984015)+(cc$nhyd*16) - 0.5
					ubl <- max( cd1$mass) + (cc$ndeam*0.984015)+(cc$nhyd*16) + 0.5
		
		
					#message(sprintf("cc$ndeam is %f, cc
					message(sprintf("min(cd1$mass) is %f, max is %f",min( cd1$mass),max( cd1$mass)))
					message(sprintf("lbl is %f, ubl is %f",lbl,ubl))
		
					#In case of overlap, we have to *ADD* e1 to matchedpeaks...not *REPLACE*
					matchedpeaks[     which(e1[,1] > lbl & e1[,1] < ubl),2] <-
						#matchedpeaks[ which(e1[,1] > lbl & e1[,1] < ubl),2] 
						#+
						eraw[           which(e1[,1] > lbl & e1[,1] < ubl),2]
		
					#Set the remaining peakset to zero
		 			e1[ which(e1[,1] > lbl & e1[,1] < ubl),2] <- 0
			
			
				}
				else{
					answer <- ba_ynq("Do you want to quit?")
					if(answer == "Y" || answer == "y"){
						finished = T
					}
				}
			
				if(doyn)
					readline("Hit <return> to do the next iteration\n\n")	
		#		readline("Press <return> to do the next most intense peak")		
			}	
		}
	}
}



#we know start and end, so we can calculate n_hyds in this range
#phydp <- get_hydroxylations(sheet,start,end)
get_hydroxylations <- function(sheet,start,end,dopause=T){


	hidx<-ts_index(sheet,"hydroxylation")
	if(hidx<0){
	readline("Couldn't find the hydroxylation row - check the sheet!\nhit <return> to continue")
	}
	
	message("pos\tprob")
	#TODO: this is redundant if the method below works...
	for(i in start:end){
		message(sprintf("%d\t%f",i,sheet[hidx,i]))
	}

	#TODO: Check this is a better method:
	d <- as.numeric(sheet[hidx,start:end])
	i <- (start-4):(end-4)
	
	hoffset <- -16
	output<-data.frame(helixpos=i[which(d>0)]+hoffset,pos=i[which(d>0)],prob=d[which(d>0)])
	
	message(sprintf("start=%d; end=%d\n",start,end))
	if(nrow(output)>0)
		print(output)
	else
		message("No ",appendLF=F)
	
	message("Hydroxylation probs found")
	#\nhit <return> to continue")

	return(output)

}











#ranked_alignment_of_mass_spectrum
rams <- function(sheet, spp, ms1, ms2, ms3, dopause=F, scode, doccplot=F){

	aa <- list()


	message(sprintf("Ranking alignement with theoretical spectrum for %s",spp))
	message(sprintf("Searching for %s...",spp))

	spidx<-ts_index(sheet,spp)
	if(spidx<0){
		return
	}
	
	

	

	###################################################
	message(sprintf("Match for %s found at index %d, now calculating spectrum", spp, spidx))

	endcol<-ncol(sheet)	

	message("Calculating sequences now")
	start<-4
	count<-1

	moff = 1.5
	
	plotno=0;
	
	#TODO: Use this variable to check whether the expression rate is single or double
	#There's a marker in the mammalian collagen file.
	times_expressed = 1;
	
	#TODO: figure out the best way to look this up in the Mammalin Collagen Sequences sheet
	ANU = 1061 
	 
	 
	#outfn1 = ("analysis1.dat")
	#outfn2 = ("analysis2.dat")
	#outfn3 = ("analysis3.dat")
	
	
	outfn1 <- sprintf("%sanalysis.dat",scode[1])
	outfnhyd <- sprintf("%s_hprobs.dat",scode[1])
	
	#TODO: split this out of this function 
	#Create a new table for the hydroxylations..
	
	if(file.exists(outfnhyd))
		file.remove(outfnhyd)
	
	#write_header(outfn1)
	##ad1 <- NA #new_align_df()
	#ad2 <- new_align_df()
	#ad3 <- new_align_df()
	 	 
 	message(sprintf("%d:\t",count),appendLF=F)
	for(j in start:endcol){
		#TODO: reject non a-a entries
		if(grepl("K|R",sheet[spidx,j])){

			mybg = 1 #Set 'standard' pinhead colour to red

			count<-count+1
			message(sprintf("%s\n%d:\t",sheet[spidx,j],count),appendLF=F)
			end=j
			#to make the string, do something like:
			#cc <- paste0(shit[85,4:12],collapse="")
			
			sequence <- paste0(sheet[spidx,start:end],collapse="")

			nhyd <- str_count(sequence,"P")
			nglut <- str_count(sequence,"Q") + str_count(sequence,"N")
			
			
			###########################################################################
			#TODO: We only need to calculate this ONCE per available sequence - so we should do this elsewhere...
			#we know start and end, so we can calculate n_hyds in this range
			phydp <- get_hydroxylations(sheet,start,end)
			str(phydp)
			
			
			hdat <- data.frame(sequence,start-4,end-start,nhyd)
			
			if(nrow(phydp) >0 ){
				pnh <- ba_nhyd(phydp$prob)
				#print(pnh)
				message(sprintf("nhyd\tprob"))
				message(sprintf("%d\t%0.4f\n",pnh$nhyd,pnh$prob))
				
				#TODO: find a better way do create out data for the nhyds file

				hdat <- data.frame(hdat,t(pnh$prob))
		
				#readline("Calculated Nhyd prob\nhit<return> to continue")
			}
			else{
				message("No hydroxylation probabilities for this sequence")
				#TODO: later on - check if there any P's with no probs available
				
				hdat <- data.frame(hdat,1)
			}
			
			


			#################################################			
			#WRITE THE HYDROXY PROBS TO FILE
			
			write.table(
				hdat,
				file = outfnhyd,append=T,sep = ",", row.names = F, col.names = F)
						







			message(sprintf("Generating spectrum for sequence %s ...",sequence,appendLF=F))

			cd1 <- q2e_tpeaks(sequence)

			###########################################################
			#RESTRICTION: Only do this for the range we have data for #
			if(max(cd1$mass) > 800 && min(cd1$mass) < 3500){
			
				message(sprintf("Sequence is %s\nThere are %d prolines in this segment",sequence,nhyd))
				message(sprintf("There are %d glutamines  in this segment",nglut))
	
				message(sprintf("Mass range is %f to %f",min(cd1$mass),max(cd1$mass)))
			
				cc=0;
				for(e in 0:nglut){
					for( p in 0:nhyd){
					
					
#grab the data for the range we are intested in 
#xm <-       testdata[[spot]]@mass[myxlim[1] <= testdata[[spot]]@mass & testdata[[spot]]@mass < myxlim[2] ]
#yi <-  testdata[[spot]]@intensity[myxlim[1] <= testdata[[spot]]@mass & testdata[[spot]]@mass < myxlim[2] ]
			
						lbl <- min(cd1$mass) - moff + (e*0.984015)+(p*16)
						ubl <- max(cd1$mass) + moff + (e*0.984015)+(p*16)
								
						subms1 <- ms_subrange(ms1,lbl,ubl)
						subms2 <- ms_subrange(ms2,lbl,ubl)
						subms3 <- ms_subrange(ms3,lbl,ubl)
					
						#Only bother if the energy is there	
						#if(max(subms1[,2]) > (0.1*max(ms1[,2])) ){
						
						mir <- 0.05
						
						#TODO: Need to use ALL THREE SAMPLES 	
						if(max(subms1[,2]) > (mir*max(ms1[,2])) ){
							message(sprintf("Max intensity  ratio sufficient in this segment (%f > %f) ", max(subms1[,2]), mir*max(ms1[,2]) ))
					
							plotno = plotno+1
			
							message(sprintf("\nPlot number %d\nSegment at position %d of %d",
									plotno,j-4,endcol-4))
							#0.984015 - if Q changes to E - add this much....
			
							myxlim = c(lbl,ubl)
							#TODO: start-4 is used a few times - so it needs setting as a variable
							if(nrow(phydp)>0){		
								mymain <- sprintf(
								"plot %d, seqpos %d\n%s\nnglut = %d/%d, nhyd = %d/%d, hp=%0.4f (max = %0.4f)",
								plotno,start-4,sequence,e,nglut,p,nhyd,pnh$prob[which(pnh$nhyd == p)],max(pnh$prob))
							}
							else{	
								mymain <- sprintf(
								"plot %d, seqpos %d\n%s\nnglut = %d/%d, nhyd = %d/%d, hp=UNKNOWN",
								plotno,start-4,sequence,e,nglut,p,nhyd)
							}							
							
							
							plot(1, type="n", xlab="Mass", ylab = "Probability", 
									xlim=myxlim, ylim=c(0, 1), main=mymain)
							
							legend('topright',c(scode[2],scode[3],scode[4],'pre-align'), lty = c(1,1,1,1),
								   col=c('red','green','blue','grey'),ncol=1,bty ="n")
			
							#ba_plotseqpeaks(cd1,myxlim)
				
							cc = cc+1;
							x <- cd1$mass +(e*0.984015)+(p*16)
							y <- cd1$prob
							
							cdshift <-cd1
							cdshift$mass <- cd1$mass +(e*0.984015)+(p*16)
							
							segments(x0=x, y0=y, y1=0, col=8)
							points(x, y, pch=21, col=1+e, bg=mybg+p)
						
							message(sprintf("Calculating alignment for sequence %s, #%d of %d, 1stMass: %0.2f nglut: %d, nhyd %d\n",
									sequence, cc,(nglut+1)*(nhyd+1),x[1],e,p))
							
							#plot the peaks
							#ms_peaklineplot <- function(sms,ms,col)
							
							#/Calculate the shift (if any) */
						
							
							
							ms_offset_peaklineplot(ms1,0,"grey")
							ms_offset_peaklineplot(ms2,0,"grey")
							ms_offset_peaklineplot(ms3,0,"grey")
							
										
							#readline(" Press <return> to do the alignment")
						
							#align1 <- ba_ms_align(cdshift,subms1,myxlim,doplot=T)
							align1 <- ba_ms_align(cdshift,subms1,myxlim)
							message(sprintf("red:   cor = %0.3f, lag = %0.3f",align1$cor,align1$lag))
							#sequence,position,length,fragno,plotno,cor,lag,pm1,energy,efrac
							
							subms1[,1] = subms1[,1] + align1$lag
							#ms_peaklineplot(subms1,ms1,"red")
							ms_offset_peaklineplot(ms1,align1$lag,"red")
						
							write.table(
								t(c(
									#sequence,position,length,fragno,plotno,
									sequence,start-4,end-start,count,plotno,
									#ndean,nhyd
									e,p,
									#cor,      lag,        
									align1$cor,align1$lag,
									#pm1,pm5,energy,efrac
									x[1],x[5],sum(subms1[,2]), sum(subms1[,2]) / sum(ms1[,2])
									)),
								file = outfn1,append=T,sep = ",", row.names = F, col.names = F)
						
							
							align2 <- ba_ms_align(cdshift,subms2,myxlim)
							message(sprintf("green: cor = %0.3f, lag = %0.3f",align2$cor,align2$lag))
							subms2[,1] = subms2[,1] + align2$lag
							ms_offset_peaklineplot(ms2,align2$lag,"green")
							
							
							align3 <- ba_ms_align(cdshift,subms3,myxlim)
							message(sprintf("blue:  cor = %0.3f, lag = %0.3f",align3$cor,align3$lag))
							subms3[,1] = subms3[,1] + align3$lag
							ms_offset_peaklineplot(ms3,align3$lag,"blue")
							
							
							
							
							if(dopause){		
								readline("hit <return> to look at the next segment\n\n")
							}
						}
						else{
							message(sprintf("Max intensity  ratio too small in this segment (%f/%f)", max(subms1[,2]), max(ms1[,2]) ))
						}
					}
				}
				
				
				
				#plot(d1[,1],d1[,2]/max(d1[,2]),type="l",xlim=myxlim, ylim=c(0,1), bty='l', main = "Method K: 1,Norfolk (black); 6,Italy (red)")
				#lines(i1[,1],i1[,2]/max(i1[,2]), col = "red")
				
				#plot(d2[,1],d2[,2]/max(d2[,2]),type="l",xlim=myxlim, ylim=c(0,1), bty='l', main = "Method I: 1,Norfolk (black); 6,Italy (red)")
				#lines(i2[,1],i2[,2]/max(i2[,2]), col = "red")
				
				#plot(d3[,1],d3[,2]/max(d3[,2]),type="l",xlim=myxlim, ylim=c(0,1), bty='l', main = "Method C: 1,Norfolk (black), 6,Italy (red)")
				#lines(i3[,1],i3[,2]/max(i3[,2]), col = "red")
				
				#readline("hit <return> to look at the next segment")

			}
			message("done")
			


			#readline("hit <return> to continue...\n")
			start = j+1
 		}
 		else{
 			message(sprintf("%s",sheet[spidx,j]),appendLF=F)
 		}	 	
	 }
	 
	 #return(ad1)
}



generate_alignment_pdf <- function(froot,scode,spots,species = "human"){

	
	anfn <- sprintf("%sanalysis.dat",scode)
	#Delete the analyisis file
	if(file.exists(anfn))
		file.remove(anfn)
	
	e1 <- read.table(sprintf("%s%s.txt",froot,spots[1]))
	e2 <- read.table(sprintf("%s%s.txt",froot,spots[2]))
	e3 <- read.table(sprintf("%s%s.txt",froot,spots[3]))
	
	par(mar=c(3.8,3.8,3,1), mfrow = c(1,1), oma=c(2,2,2,2))
	pdf(file = sprintf("%s.pdf",scode),w=6,h=4)
		al <- rams(sheet,species,e1,e2,e3,dopause=F,c(scode,spots),doccplot=T)
	dev.off()

	#TODO: We'll need to do an absoulute fit, but this is NOT practical here...
	#readline("Press <return> to do absolute fit now")
	#ms_fit(sheet,species,froot,scode,spots)

}







