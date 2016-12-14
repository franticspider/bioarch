

#load the utils file
source("~/git/bioarch/dev/ba_msutils.R")

#load the hyd probs data
hp <- load.phyd("C5_hprobs.dat")


xlim=c(500,3250)
myylim=c(0,1.2)
plot(1, type="n", xlab="Mass", ylab = "Probability", 
			xlim=xlim, ylim=myylim, main="Weighted theoretical spectrum for human collagen\nnumber of hydroxylations is indicated via 'rainbow' colouring: red = 0, violet = 11 hydroxylations")


hlim = (ncol(hp)-5)

mycol = rainbow(hlim+1)


for(s in 1: nrow(hp)){

	##set a sequence
	#seq <- "GLPGPPGAPGPQGFQGPPGEPGEPGASGPMGPR"
	#
	##get the hyd data for this sequence
	#shp <- hp[which(hp$sequence == seq),]

	shp <- hp[s,]
	seq <- shp$sequence


	#get the mass, nhyd and ndeam values:
	require(stringr)
	nhyd_P <- str_count(seq,"P")
	ndeam <- str_count(seq,"Q") + str_count(seq,"N")
			
	#get the theoretical peaks
	require(q2e)
	cd1 <- q2e_tpeaks(seq)



	#We need to count the number of hydroxylations *that have probabilities* - otherwise we could
	#try to reference outside the calculated values.  
	#swrk <- shp[5:ncol(shp)]
	#swrk <- swrk[which(swrk>1e-12)]
	#nhyd_withprobs <- length(swrk)#-1 #we have to -1 because otherwise we count the h=0 twice

	message(sprintf("nhdyd_withprobs is %d",nhyd_withprobs))	


	mymain <- sprintf("Seq %d: %s\npos %d, #P = %d, #p2 = %d",s,seq,shp$start,nhyd_P,shp$nh)

	#TODO: make this a switch
	##ok, now we have everything we need to plot: 
	#xlim <- c(cd1$mass[1]-3,cd1$mass[1]+(ndeam*0.984015)+(nhyd_P*16)+4+3)	
	#
	#if(shp$nh == 0)
	#	myylim <- c(0,1)
	#else
	#	myylim <- c(0, as.numeric(max(shp[5:ncol(shp)],0.05)))
	#
	#plot(1, type="n", xlab="Mass", ylab = "Probability", 
	#		xlim=xlim, ylim=myylim, main=mymain)



	message(sprintf("myylim: %0.3f-%0.3f",myylim[1],myylim[2]))


	for(d in 0:ndeam){
		for( h in 0:hlim){#nhyd_withprobs){
	
	
	
	
	
	
			x <- cd1$mass +(d*0.984015)+(h*16)
			#NB the as.numeric function stops the cd1$prob collapsing to a single value!
			if(nhyd_P==0){
				y <- cd1$prob
			}
			else{
				y <- cd1$prob * as.numeric(shp[5+h])
			}
				
			text(cd1$mass[1],1.2,sprintf("%0.1f: %s",cd1$mass[1],seq),srt=90,cex=0.5,adj=c(1,0),font=3)
							
			segments(x0=x, y0=y, y1=0, col="black")
			points(x, y, pch=21, col=mycol[h+1], bg=mycol[h+1], cex = 0.25)
			
			
			#segments(x0=x, y0=y, y1=0, col="black")
			#points(x, y, pch=21, col=d+3, bg=d+3, cex = 0.5)
				
				
			#For sanity, put a point at the first calculated mass:	
			#points(cd1$mass[1],-0.005,col = "red",pch=10,cex=3)
	
		}
	}
	
	#readline("hit <return> for the next plot")
}
