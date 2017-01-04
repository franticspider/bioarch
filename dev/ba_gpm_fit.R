###############################################################################
###############################################################################
###############################################################################
###############################################################################
# OBSOLETE FILE! - SEE PACKAGE 'bacollite' FOR WORKING VERSION
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################



library(stringr)
source("~/git/bioarch/dev/ba_msutils.R")

###############################################################################
#
# read the gpm file
#####################

#TODO: get the info from google sheets....
load.gpm <- function(fn="~/tmp/bioarch_keri/161208meeting/gmp_human_collagen.dat"){


	gpm <- read.table(fn,
	fill=T,header=T,sep="\t",comment.char="")

	gpm2 <- data.frame(seq=as.character(gpm$sequence), 
		nhyd = str_count(gpm$modifications,"P"), 
		nglut = (str_count(gpm$modifications,"Q")+str_count(gpm$modifications,"N")), 
		mass1=0, prob1=0 )

	gpm2 <- unique(gpm2)

	for(i in 1:nrow(gpm2)){

		#get the mass data...
	
		cd1 <- q2e_tpeaks(gpm2$seq[i])
		gpm2$mass1[i] = cd1$mass[1] +  (gpm2$nglut[i]*0.984015)+(gpm2$nhyd[i]*16)
		gpm2$prob1[i] = cd1$prob[1]
	}

	#remove outliers
	gpm2 <- gpm2[which(gpm2$mass1>795 & gpm2$mass1 <3500),]

	#Sort by mass
	gpm2 <- gpm2[order(gpm2$mass1),]
	
	return(gpm2)
	
}



########################################################################
#
# load the spectrum data
##########################

count = 0;

#Load the data

froot  <- "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_"
spots  <- c("G7","G10","G13")
sample <- "C1"
	
ms1 <- read.table(sprintf("%s%s.txt",froot,spots[1]))
ms2 <- read.table(sprintf("%s%s.txt",froot,spots[2]))
ms3 <- read.table(sprintf("%s%s.txt",froot,spots[3]))




########################################################################
#
# carry out the analysis
##########################

#set up the pdf:
par(mar=c(3.8,3.8,3,1), mfrow = c(1,1), oma=c(2,2,2,2))
pdf(file="gpm.pdf",w=6,h=4)


#initialise:
plotno <- 0
moff <- 0.5


for(i in 1:nrow(gpm2)){

	#get the mass data...
	#TODO: has this been done? 
	cd1 <- q2e_tpeaks(gpm2$seq[i])
	
	#readline(sprintf("Mass of sequence %s is %0.2f, nhyd = %d, nglut = %d, hit <return>",gpm2$seq[i], cd1$mass[1], gpm2$nhyd[i], gpm2$nglut[i]))
	
	#TODO: consolidate with the _rams_ function...
	###########################################################
	#RESTRICTION: Only do this for the range we have data for #
	if(max(cd1$mass) > 800 && min(cd1$mass) < 3500){
	
		message(sprintf("Sequence is %s\nThere are %d prolines in this segment",
						gpm2$seq[i],gpm2$nhyd[i]))
		#message(sprintf("There are %d glutamines  in this segment",nglut))
		message(sprintf("Mass range is %f to %f",min(cd1$mass),max(cd1$mass) )) 
		count = count +1	
		
		lbl <- min(cd1$mass) - moff + (gpm2$nglut[i]*0.984015)+(gpm2$nhyd[i]*16)
		ubl <- max(cd1$mass) + moff + (gpm2$nglut[i]*0.984015)+(gpm2$nhyd[i]*16)
				
		subms1 <- ms_subrange(ms1,lbl,ubl)
		subms2 <- ms_subrange(ms2,lbl,ubl)
		subms3 <- ms_subrange(ms3,lbl,ubl)
		
		mir <- 0.05
		
		if(max(subms1[,2]) > (mir*max(ms1[,2])) ){
			message(sprintf("Max intensity  ratio sufficient in this segment (%f > %f) ", max(subms1[,2]), mir*max(ms1[,2]) ))
			
			cdshift <-cd1
			cdshift$mass <- cd1$mass + (gpm2$nglut[i]*0.984015)+(gpm2$nhyd[i]*16)
			
			myxlim = c(lbl,ubl)			

			align1 <- ba_ms_align(cdshift,subms1,myxlim)
			align2 <- ba_ms_align(cdshift,subms2,myxlim)
			align3 <- ba_ms_align(cdshift,subms3,myxlim)
			
			message(sprintf("align1$lag is %0.3f, cor is %0.3f",align1$lag,align1$cor))
			message(sprintf("align2$lag is %0.3f, cor is %0.3f",align2$lag,align2$cor))
			message(sprintf("align3$lag is %0.3f, cor is %0.3f",align3$lag,align3$cor))
			
			doplot=F
			if(abs(align1$lag) < 0.3 & align1$cor > 0.15)
			{ doplot=T}
			if(abs(align2$lag) < 0.3 & align2$cor > 0.15)
			{ doplot=T}
			if(abs(align3$lag) < 0.3 & align3$cor > 0.15)
			{ doplot=T}
			
			if(doplot){
			
				plotno = plotno+1

				message(sprintf("\nPlot number %d\nSegment at row %d of %d",
						plotno,i,nrow(gpm2)))
				#0.984015 - if Q changes to E - add this much....

				#TODO: start-4 is used a few times - so it needs setting as a variable
				#if(nrow(phydp)>0){		
				#	mymain <- sprintf(
				#	"plot %d, seqpos %d\n%s\nnglut = %d/%d, nhyd = %d/%d, hp=%0.4f (max = %0.4f)",
				#	plotno,start-4,sequence,e,nglut,p,nhyd,pnh$prob[which(pnh$nhyd == p)],max(pnh$prob))
				#}
				#else{	
				#	mymain <- sprintf(
				#	"plot %d, seqpos %d\n%s\nnglut = %d/%d, nhyd = %d/%d, hp=UNKNOWN",
				#	plotno,start-4,sequence,e,nglut,p,nhyd)
				#}	
			
				mymain <- sprintf(
				"plot %d, entry %d, mass %0.1f\n%s\nndean = %d, nhyd = %d lag: %0.2f,%0.2f,%0.2f cor: %0.2f,%0.2f,%0.2f",
				plotno,i,gpm2$mass1[i],
				gpm2$seq[i],gpm2$nglut[i],gpm2$nhyd[i],
				align1$lag,align2$lag,align3$lag,
				align1$cor,align2$cor,align3$cor)						
			
				myxxlim <- c(lbl-1,ubl+1)
				plot(1, type="n", xlab="Mass", ylab = "Probability", 
						xlim=myxxlim, ylim=c(0, 1), main=mymain)
			
				#legend('topright',spots,'pre-align'), lty = c(1,1,1,1),
				#	   col=c('red','green','blue','grey'),ncol=1,bty ="n")

				#ba_plotseqpeaks(cd1,myxlim)

				#cc = cc+1;
				x <- cd1$mass + (gpm2$nglut[i]*0.984015)+(gpm2$nhyd[i]*16)
				y <- cd1$prob
			
				segments(x0=x, y0=y, y1=0, col="black")
				points(x, y, pch=21, col="black", bg="red")
			
				#####################
			
							
				ms_offset_peaklineplot(ms1,0,"grey")
				ms_offset_peaklineplot(ms2,0,"grey")
				ms_offset_peaklineplot(ms3,0,"grey")
			
			
				ms_offset_peaklineplot(ms1,align1$lag,"red")
				ms_offset_peaklineplot(ms2,align2$lag,"green")
				ms_offset_peaklineplot(ms3,align3$lag,"blue")
							
				#readline("hit <return> to continue")
			}
		}
		
	}
	else{
		message(sprintf("Mass %0.2f out of range",cd1$mass[1]))
	}


}

message(sprintf("%d theoretical peptides in range",count))	

dev.off()


#Now try to do the equivalent of ms_fit...






















