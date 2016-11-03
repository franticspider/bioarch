



# useful links:
# Basic tutorial
# http://blog.revolutionanalytics.com/2015/09/using-the-googlesheets-package-to-work-with-google-sheets.html


##Generate the theroretical masses
#- Use googlesheets to access the data
#- Find the row corresponding to the species
#- Generate the sequence fragments by split at K
#- For each subsequence
#	- Record the sequence
#    - then check for Xs (unknowns) 
#    - if there's not too many Xs, check for prob of hydroxylation and deamidation
#- Plot the theoretical masses that work. Adjust the relative masses
library("googlesheets")
library("bioarch")
library("q2e")

#we need this to search for letters (P) etc. 
require("stringr")

#TODO: functions from this file to be included in the bioarch package
source("~/git/bioarch/dev/ba_msutils.R")



######################################################################################

generate_diff_spectra <- function(sheet,spp1,spp2){

	sp1idx = ts_index(sheet,spp1)


}


generate_theroretical_spectrum <- function(sheet, spp){

	
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

	moff = 1
	
	plot(1, type="n", xlab="Mass", ylab = "Probability", xlim=c(800, 3500), ylim=c(0, 1), main=sheet[spidx,1],)
	 
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
			
			message(sprintf("Generating spectrum for sequence %s ...",sequence,appendLF=F))

			cd1 <- q2e_isodists(sequence)
			lbl <- min(cd1$mass) - moff
			ubl <- max(cd1$mass) + moff
			
			if(max(cd1$mass) > 800){
			
			
				myxlim = c(lbl,ubl)
			
				#ba_plotseqpeaks(cd1,myxlim)
				x <- cd1$mass
				y <- cd1$prob
				segments(x0=x, y0=y, y1=0, col=8)
				points(x, y, pch=21, col=1, bg=mybg+nhyd)

			}
			message("done")
			


			#readline("hit <return> to continue...\n")
			start = j+1
 		}
 		else{
 			message(sprintf("%s",sheet[spidx,j]),appendLF=F)
 		}
	 	
	 }
}
###################################################################






generate_theroretical_spectrum_and_ms <- function(sheet, spp, masses, intensities){

	
	message(sprintf("Searching for %s...",spp))
	
	 spidx<-ts_index(sheet,spp)
	 if(spidx<0){
	 	return
	 }

	 ###################################################
	 readline(sprintf("Match for %s found at index %d, now calculating spectrum", spp, spidx))
	 	
	 endcol<-ncol(sheet)	
	 	
	 message("Calculating sequences now")
	 start<-4
	 count<-1

	moff = 1
	
	plotno=0;
	 
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
			
			nglut <- str_count(sequence,"Q")
			
			message(sprintf("Generating spectrum for sequence %s ...",sequence,appendLF=F))

			cd1 <- q2e_isodists(sequence)
			lbl <- min(cd1$mass) - moff
			ubl <- max(cd1$mass) + moff
			
			if(max(cd1$mass) > 800){
			
				plotno = plotno+1
			
				message(sprintf("\n\nPlot number %d\nSegment at position %d of %d",plotno,j-4,endcol-4))
				message(sprintf("Sequence is %s\nThere are %d prolines in this segment",sequence,nhyd))
				message(sprintf("There are %d glutamines  in this segment",nglut))
			
			
			
				message(sprintf("Mass range is %f to %f",min(cd1$mass),max(cd1$mass)))
			
			    #0.984015 - if Q changes to E - add this much....
			
				myxlim = c(lbl,ubl+(nhyd*16))
				
				
				plot(1, type="n", xlab="Mass", ylab = "Probability", xlim=myxlim, ylim=c(0, 1), main=sequence,)
			
				#ba_plotseqpeaks(cd1,myxlim)
				
				for(e in 0:nglut){
					for( p in 0:nhyd){
				
						x <- cd1$mass +(e*0.984015)+(p*16)
						y <- cd1$prob
						segments(x0=x, y0=y, y1=0, col=8)
						points(x, y, pch=21, col=1+e, bg=mybg+p)
					}
				}
				
				
				lines(masses,intensities/max(intensities))
				
				readline("hit <return> to look at the next segment")

			}
			message("done")
			


			#readline("hit <return> to continue...\n")
			start = j+1
 		}
 		else{
 			message(sprintf("%s",sheet[spidx,j]),appendLF=F)
 		}
	 	
	 }
}
























##################################################################################################
#load the test data
if(!exists("testdata")){
	testdata <- bioarch_loadBrukerXML("/home/sjh/Desktop/sjh/bioarch/SF/SFData/20140307_SF_UPenn164-181/")
}
# Let's pick an example spot
spot <- 6
message(sprintf("Working on spot %d",spot))
#now we have to dismantle the function plotclass_v2 to get what we are after. 
#set up the plot area: (1 column, 4 rows
par(mar=c(0.9,2.3,0.9,.3), mfrow = c(2,1), oma=c(5,0,2,0))
#Let's do some setting up:
myxlim = c(800,4000)#c(1192,1198)

#grab the data for the range we are intested in 
xm <-       testdata[[spot]]@mass[myxlim[1] <= testdata[[spot]]@mass & testdata[[spot]]@mass < myxlim[2] ]
yi <-  testdata[[spot]]@intensity[myxlim[1] <= testdata[[spot]]@mass & testdata[[spot]]@mass < myxlim[2] ]
#Plot the calculated peaks as probabilities 
#until this is in the library, we have to load it:
#source("plotseqpeaks.R")
#plotseqpeaks(cd1,myxlim)
plot(xm,yi/max(yi),type="l",col= "red",xlim=myxlim)

##Compare with sample data


#Fetch the sheet if it isn't available...
if(!exists("sheet")){

	message("Fetching the sheet of mammalian collagen sequences")
	mcs <- gs_title("Mammalian Collagen Sequences v0.0.1")
	print(mcs)

	sheet <- gs_read(mcs)
}
generate_theroretical_spectrum(sheet,"cow")



#generate_diff_spectra(sheet,"cow","sheep")


readline("Hit return to see this spectrum compared with a cow sample")

par(mar=c(0.9,2.3,0.9,.3), mfrow = c(1,1), oma=c(5,0,2,0))
generate_theroretical_spectrum_and_ms(sheet,"cow",xm,yi)



