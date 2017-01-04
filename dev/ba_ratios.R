

# Building on 'ratios.R' This is a generic approach to comparing two samples - for whatever reason
# Impractical to do ratio(3hyd): ratio(2hyd) - you'd have to compare all possible hydroxylations
# Better perhaps to aim at 'fit' of the observed peaks to the expected distribution. 
# Let's look at that now..

require("stringr")
library("googlesheets")
require(q2e)
require(bioarch)

source("~/git/bioarch/dev/ba_msutils.R")
######################################################################

compare_two_spectra <- function(tss, froot1, sample1, spots1, froot2, sample2, spots2){


	generate_alignment_pdf(froot1, sample1, spots1)
	pdf(file = sprintf("%s_hits.pdf",sample1),w=8,h=8)
	ms_fit(sheet,"human",froot1,sample1,spots1,doyn=F)
	dev.off()


	generate_alignment_pdf(froot2, sample2, spots2)
	pdf(file = sprintf("%s_hits.pdf",sample2),w=8,h=8)
	ms_fit(sheet,"human",froot2,sample2,spots2,doyn=F)
	dev.off()


	#generate_alignment_pdf(froot, sample1, spots1)
	#pdf(file = "C1_hits.pdf",w=8,h=8)
	#ms_fit(sheet,"human",froot,sample1,spots1,doyn=F)
	#dev.off()


	#generate_alignment_pdf(froot, sample2, spots2)
	#pdf(file = "C1.SE_hits.pdf",w=8,h=8)
	#ms_fit(sheet,"human",froot,sample2,spots2,doyn=F)
	#dev.off()


	#readline("Finished initial pass")
	###################################################################
	#Load data sets: 
	#TODO: These have been loaded already!


	#Load data set 1
	n1.1 <- read.table(sprintf("%s%s.txt",froot1,spots1[1]))
	n1.2 <- read.table(sprintf("%s%s.txt",froot1,spots1[2]))
	n1.3 <- read.table(sprintf("%s%s.txt",froot1,spots1[3]))
	#TODO: analysis must be available - build the generation of this into the system?
	n1.a <- load.analysis(sprintf("%sanalysis.dat",sample1))


	#Load data set 2
	i1.1 <- read.table(sprintf("%s%s.txt",froot2,spots2[1]))
	i1.2 <- read.table(sprintf("%s%s.txt",froot2,spots2[2]))
	i1.3 <- read.table(sprintf("%s%s.txt",froot2,spots2[3]))
	i1.a <- load.analysis(sprintf("%sanalysis.dat",sample2))


	###################################################################
	#Load the probability of hydroxylation: 
	#TODO: these do not change! Have a default file (or table) in the package!
	phyd<- load.phyd("hprobs.dat")



	#####################################################################
	#Global parameters

	moff = 0.5 #offset value in mass - daltons either side of a peak

	#Set the plot parameters
	pdf(file=sprintf("%s_vs_%s_ratios.pdf",sample1,sample2),w=6,h=8)
	par(mar=c(0.9,2.3,2.9,.3), mfrow = c(5,1), oma=c(5,0,2,0))

	tss <- tss[order(tss$mass.1),]

	#GO THROUGH EVERY PEPTIDE...
	for(i in 1:nrow(tss)){

		#################################################################
		#GET THE MATCHES
		hits.n <- n1.a[which(n1.a$sequence == tss$seq[i] & abs(n1.a$lag)  < 0.3 & n1.a$cor >0.1),]

		hits.i <- i1.a[which(i1.a$sequence == tss$seq[i] & abs(i1.a$lag)  < 0.3 & i1.a$cor >0.1),]

		#################################################################
		#GATHER THE DATA
		#Get the deam and hydr levels
		nhyd <- str_count(tss$seq[i],"P")
		nglut <- str_count(tss$seq[i],"Q") + str_count(tss$seq[i],"N")
		message(sprintf("nglut = %d, nhyd = %d",nglut,nhyd))

		#Theoretical mass data: 
		mass <- t(tss[i,3:7])
		prob <- t(tss[i,8:12])	
		hpdata <- phyd[which(phyd$sequence == tss$seq[i]),]



		#Get the h-level probs

		################################################################
		#PLOT THE RAW SPECTRUM FOR SAMPLE 1
		mymain <- sprintf("%s raw data",sample1)
		plot(x=n1.1[,1],y=n1.1[,2]/max(n1.1[,2]), xlim = c(tss$mass.1[i]-moff,tss$mass.5[i]+ (nglut*0.984015)+(nhyd*16)+moff),type="l",main=mymain,ylim = c(0,1),col="red")
		lines(x=n1.2[,1],y=n1.2[,2]/max(n1.2[,2]),col="green")
		lines(x=n1.3[,1],y=n1.3[,2]/max(n1.3[,2]),col="blue")



		################################################################
		#TODO: atm we are basing all this on the first spectrum - we need to align and 
		#      then merge the hits - not sure how best to do this..... 
		################################################################
		#PLOT THE MATCHED PEAKS FOR SAMPLE1
		if(nrow(hits.n)>0){
			mymain <- sprintf("%s aligned matches",sample1)
			plot(NA, xlim = c(tss$mass.1[i]-moff,
					tss$mass.5[i]+ (nglut*0.984015)+(nhyd*16)+moff),
					type="l",main=mymain,ylim = c(0,1))
				
			#Draw the mass for each hydroxylation (to line things up) - no deamidation
			for( p in 0:nhyd){
				xx <- mass + (p*16)
				segments(x0=xx, y0=1, y1=0, col="grey")
			}	
				
			for(j in 1:nrow(hits.n)){
				hms <- n1.1[which(n1.1[,1] > hits.n$pm1[j]-0.5 & n1.1[,1] < hits.n$pm5[j]+0.5),]
				lines(x=hms[,1],y=hms[,2]/max(hms[,2]), col = "red")
			
				hms <- n1.2[which(n1.2[,1] > hits.n$pm1[j]-0.5 & n1.2[,1] < hits.n$pm5[j]+0.5),]
				lines(x=hms[,1],y=hms[,2]/max(hms[,2]), col = "green")
			
				hms <- n1.3[which(n1.3[,1] > hits.n$pm1[j]-0.5 & n1.3[,1] < hits.n$pm5[j]+0.5),]
				lines(x=hms[,1],y=hms[,2]/max(hms[,2]), col = "blue")
		
			}		
				
			#lines(x=i2[,1],y=i2[,2]/max(i2[,2]),col="green")
			#lines(x=i3[,1],y=i3[,2]/max(i3[,2]),col="blue")
		}
		else{
			#Empty plot....
			plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), xaxt='n',yaxt='n',bty='n',pch='')
			text(x=5,y=5,labels=sprintf("No matches for sequence %s",tss$seq[i]))
		}

		#TODO: OVERLAY THE RATIO DATA


		################################################################
		#PLOT THE THEORETICAL SPECTRUM
	
		#message(mass)
		message(hpdata)
	
		ccol <- "red"
		plot(NA,xlim = c(tss$mass.1[i]-moff,tss$mass.5[i]+ (nglut*0.984015)+(nhyd*16)+moff),
			ylim = c(0, 1  ),
			main = "Theoretical spectrum" )
		for(e in 0:nglut){
			message(sprintf("glut %d:",e))
			for( p in 0:nhyd){
				xx <- mass +(e*0.984015)+(p*16)
				yy <- prob * as.numeric(hpdata[5+p]) / max(as.numeric(hpdata[5:(5+nhyd)]))
				#message(sprintf("mass[1] is %f, xx[1] is %f",mass[1],xx[1]))
				message(sprintf("hpdata[5+%d] is %f, prob[%d] is %f",p,hpdata[5+p],p,prob[p]))
				#message(sprintf("xrow[1] is %0.2f",xrow[1]))
				segments(x0=xx, y0=yy, y1=0, col="red")
				points(xx, yy, pch=21, col=ccol, bg=ccol)
			}
		}

	
		################################################################
		#PLOT THE MATCHED PEAKS FOR SAMPLE2
		if(nrow(hits.i)>0){
			mymain <- sprintf("%s aligned matches",sample2)
			plot(NA, xlim = c(tss$mass.1[i]-moff,
					tss$mass.5[i]+ (nglut*0.984015)+(nhyd*16)+moff),
					type="l",main=mymain,ylim = c(0,1))
				
			#Draw the mass for each hydroxylation (to line things up) - no deamidation
			for( p in 0:nhyd){
				xx <- mass + (p*16)
				segments(x0=xx, y0=1, y1=0, col="grey")
			}	
			
			for(j in 1:nrow(hits.i)){
				hms <- i1.1[which(i1.1[,1] > hits.i$pm1[j]-0.5 & i1.1[,1] < hits.i$pm5[j]+0.5),]
				lines(x=hms[,1],y=hms[,2]/max(hms[,2]), col = "red")
			
				hms <- i1.2[which(i1.2[,1] > hits.i$pm1[j]-0.5 & i1.2[,1] < hits.i$pm5[j]+0.5),]
				lines(x=hms[,1],y=hms[,2]/max(hms[,2]), col = "green")
			
				hms <- i1.3[which(i1.3[,1] > hits.i$pm1[j]-0.5 & i1.3[,1] < hits.i$pm5[j]+0.5),]
				lines(x=hms[,1],y=hms[,2]/max(hms[,2]), col = "blue")
			}		
				
			#lines(x=i2[,1],y=i2[,2]/max(i2[,2]),col="green")
			#lines(x=i3[,1],y=i3[,2]/max(i3[,2]),col="blue")
		}
		else{
			#Empty plot....
			plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), xaxt='n',yaxt='n',bty='n',pch='')
			text(x=5,y=5,labels=sprintf("No matches for sequence %s",tss$seq[i]))
		}
		################################################################
		#PLOT THE RAW SPECTRUM FOR SAMPLE2
		mymain <- sprintf("%s raw data",sample2)
		plot(x=i1.1[,1],y=i1.1[,2]/max(i1.1[,2]), xlim = c(tss$mass.1[i]-moff,tss$mass.5[i]+ (nglut*0.984015)+(nhyd*16)+moff),type="l",main=mymain,ylim = c(0,1),col="red")
		lines(x=i1.2[,1],y=i1.2[,2]/max(i1.2[,2]),col="green")
		lines(x=i1.3[,1],y=i1.3[,2]/max(i1.3[,2]),col="blue")




	
		#readline("Hit <return> for the next plot")
		mtext(sprintf("%s vs %s for seq %s",sample1,sample2,tss$seq[i]), outer = TRUE, cex = 1.2)
	}
	
	
	dev.off()
}



site_ratios <- function(){

	nsamples <- 10	
	Csample<-list()
	Nsample<-list()
	
	maxhyd <- 0
	maxeratio <- 0
	
	for(i in 1:nsamples){
		af <- sprintf("C%danalysis.dat",i)
		Csample[[i]] <- load.analysis(af)
		maxhyd <- max(maxhyd,Csample[[i]]$nhyd)
		maxeratio <- max(maxeratio,Csample[[i]]$efrac)
		
		
		af <- sprintf("N%danalysis.dat",i)
		Nsample[[i]] <- load.analysis(af)
		maxhyd <- max(maxhyd,Csample[[i]]$nhyd)
		maxeratio <- max(maxeratio,Csample[[i]]$efrac)
	}
	
	sd <- load.phyd("C1_hprobs.dat")
	
	pdf(file="NIratios.pdf",w=8,h=8)
	
	par(mar=c(3.8,3.8,3,1), mfrow = c(2,1), oma=c(2,2,2,2))
	
	plot(NA,ylim=c(0.001,maxeratio),xlim=c(0,maxhyd+0.5),ylab="Fraction of total ion count",xlab="Number of hydroxylations",
	#log="y",
	main="method C")
	
	#Now let's try and plot some lines...
	CNh <- rep(list(NA),(maxhyd+1))#length=maxhyd+1)
	NNh <- rep(list(NA),(maxhyd+1))#length=maxhyd+1)
	CIh <- rep(list(NA),(maxhyd+1))#length=maxhyd+1)
	NIh <- rep(list(NA),(maxhyd+1))#length=maxhyd+1)
	
	CNhmean = vector(length=length(CNh))
	NNhmean = vector(length=length(CNh))
	CIhmean = vector(length=length(CNh))
	NIhmean = vector(length=length(CNh))
	
	for(i in 1:5){
		points(x=Csample[[i]]$nhyd-0.1,y=Csample[[i]]$efrac,pch=1,col="red")
		for(h in 0:maxhyd){
			CNh[[(h+1)]] <- c(CNh[[(h+1)]],Csample[[i]]$efrac[which(Csample[[i]]$nhyd == h)])
		}	
	}
	
	for(i in 6:10){
		points(x=Csample[[i]]$nhyd+0.1,y=Csample[[i]]$efrac,pch=1,col="green")
		for(h in 0:maxhyd){
			CIh[[h+1]] <- c(CIh[[h+1]],Csample[[i]]$efrac[which(Csample[[i]]$nhyd == h)])
		}
	}
	
	for(gg in 0:maxhyd+1){
	                  #mean(CNh[[1]][2:length(CNh[[1]])])
		CNhmean[gg] <- mean(CNh[[gg]][2:length(CNh[[gg]])])
		CIhmean[gg] <- mean(CIh[[gg]][2:length(CIh[[gg]])])
		
	}
	
	
	lines(x=c(0:6)-0.1,y=CNhmean)
	points(x=c(0:6)-0.1,cex=1.5,y=CNhmean)

	lines(x=c(0:6)+0.1,y=CIhmean,col="blue")
	points(x=c(0:6)+0.1,cex=1.5,y=CIhmean,col="blue")
	
	legend('topright', col= c("red","green","black","blue"), lty = c(0,0,1,1),pch=1,
		legend=c("Norfolk","Italy","mean Norfolk","mean Italy"),cex=0.75)
	
	plot(NA,ylim=c(0.001,0.12),xlim=c(0,maxhyd+0.5),ylab="Fraction of total ion count",xlab="Number of hydroxylations",
	#log="y",
	main="method N")
	
	#Now let's try and plot some lines...
	CNh <- rep(list(NA),(maxhyd+1))#length=maxhyd+1)
	NNh <- rep(list(NA),(maxhyd+1))#length=maxhyd+1)
	CIh <- rep(list(NA),(maxhyd+1))#length=maxhyd+1)
	NIh <- rep(list(NA),(maxhyd+1))#length=maxhyd+1)
	
	
	CNhmean = vector(length=length(CNh))
	NNhmean = vector(length=length(CNh))
	CIhmean = vector(length=length(CNh))
	NIhmean = vector(length=length(CNh))
	
	for(i in 1:5){
	
		points(x=Nsample[[i]]$nhyd-0.1,y=Nsample[[i]]$efrac,pch=2,col="red")
		
		for(h in 0:maxhyd){
			NNh[[(h+1)]] <- c(NNh[[(h+1)]],Nsample[[i]]$efrac[which(Nsample[[i]]$nhyd == h)])
		}
		
	}
	
	
	
	#italy
	for(i in 6:10){
		points(x=Nsample[[i]]$nhyd+0.1,y=Nsample[[i]]$efrac,pch=2,col="green")
		
		for(h in 0:maxhyd){
			NIh[[h+1]] <- c(NIh[[h+1]],Nsample[[i]]$efrac[which(Nsample[[i]]$nhyd == h)])
		}
	}
	
	
	for(gg in 0:maxhyd+1){
	                  #mean(CNh[[1]][2:length(CNh[[1]])])
		NNhmean[gg] <- mean(NNh[[gg]][2:length(NNh[[gg]])])
		NIhmean[gg] <- mean(NIh[[gg]][2:length(NIh[[gg]])])		
	}
			
	lines(x=c(0:6)-0.1,y=NNhmean)
	points(x=c(0:6)-0.1,cex=1.5,y=NNhmean,pch=2)

	lines(x=c(0:6)+0.1,y=NIhmean,col="blue")
	points(x=c(0:6)+0.1,cex=1.5,y=NIhmean,col="blue",pch=2)

	legend('topright', col= c("red","green","black","blue"), lty = c(0,0,1,1),pch=2,
		legend=c("Norfolk","Italy","mean Norfolk","mean Italy"),cex=0.75)
			
	mtext("Hydroxylation levels, Norfolk vs Italy for methods C and N", outer = TRUE, cex = 1.2)
	
	dev.off()
}






