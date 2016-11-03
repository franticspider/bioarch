
source("~/git/bioarch/dev/ba_msutils.R")


#step 1: Get the human spectra



###################################
#So, folder one (20160803_Keri12) contains the spectra for methods B and N (both non-destructive). Folder two (20160909_Keri13) contains the spectra for methods K, I and C (K and C are destructive methods involving acid at 4 degrees, and method I is destructive, involving acid at -20 degrees). Everything is spotted in triplicate as usual.

#The files I've shared with you have already been converted so that they can be opened in mMass - if you need the raw files as they are when they come off the Ultraflex that's no problem, I can go and get those first thing tomorrow. (Sorry, I should be able to remember whether you need converted or unconverted files but my brain's been like a sieve lately!)

#Samples 1-5 are from a poor population in England (and therefore might have been suffering from scurvy). Samples 6-11 are from an Italian population, and samples 12-20 are from a Portuguese population, so we would expect that none of these individuals would have had scurvy. Samples 21-24 are animals and so can be ignored, as they'll have different collagen sequences to the humans and definitely won't have had scurvy. All of the samples are archaeological.



#	REVERSE PLATEMAP FOR KERI'S DATA
#   |-----ENGLAND-----|------ITALY-------------|------------PORTUGAL--------------|	
#	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20
#-------------------------------FOLDER ONE: 20160803------------------------------|
#			
#B	A1	A2	A3  B1	B3	C1	C2	C3	A10	A11	A12	B10	B12	C10 C11	C12	A19	A20	A21	B19		
#	A4	A5	A6	B4	B6	C4	C5	C6	A13	A14	A15	B13	B15	C13	C14	C15	A22	A23	A24	B22	
#	A7	A8	A9	B7	B9	C7	C8	C9	A16	A17	A18	B16	B18	C16	C17	C18	D1	D2	D3	E20	
#																						
#N	D4	D5	D6	E4	E5	F4	F5	F6	D13 D14 D15	E13	E15	F13	F14	F15	D22	D23 D24	E22					
#	D7	D8	D9	E7	E9	F7	F8	F9	D16	D17 D18	E16	E18	F16	F17	F18	G1	G2	G3	H1			
#	D10	D11	D12	E10	E12	F10	F11	F12	D19	D20	D21	E19	E21	F19	F20	F21	G4	G5	G6	H2		
#																						
#-------------------------------FOLDER TWO: 20160909------------------------------|
#	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20
#	
#K 	A1	A2	A3	B1	B3	C1	C2	C3	A10	A11	A12	B10	B12	C10	C11	C12	A19	A20	A21 B19
#	A4	A5	A6	B4	B6	C4	C5	C6	A13	A14	A15	B13	B15	C13	C14	C15	A22	A23	A24	B22
#	A7	A8	A9	B7	B9	C7	C8	C9	A16	A17	A18	B16	B18	C16	C17	C18	D1	D2	D3	E1	
#
#I	D4	D5	D6	E4	E6	F4	F5	F6	D13	D14	D15	E13	E15	F13	F14	F15	D22	D23	D24	E22	
#	D7	D8	D9	E7	E9	F7	F8	F9	D16	D17	D18	E16	E18	F16	F17	F18 G1	G2	G3	H1		
#	D10	D11	D12	E10	E12	F10	F11	F12	D19	D20	D21	E19	E21	F19	F20	F21	G4	G5	G6	H4	
#
#C	G7	G8	G9	H7	H9	I7	I8	I9	G16	G17	G18	H16	H18	I16	I17	I18	J1	J2	J3	K1
#	G10	G11	G12	H10	H12	I10	I11	I12	G19	G20	G21	H19	H21	I19 I20 I21	J4	J5	J6	K4				
#	G13	G14	G15	H13	H15	I13	I14	I15	G22	G23	G24	H22	H24	I22	I23	I24	J7	J8	J9	K7				



gths <- function(sheet, spp, d1, d2, d3, d4, d5, i1, i2, i3, i4, i5){

	
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
				
				mymain <- sprintf("%s, nglut = %d, nhyd = %d",sequence,nglut,nhyd)
				
				plot(1, type="n", xlab="Mass", ylab = "Probability", xlim=myxlim, ylim=c(0, 1), main=mymain)
			
				#ba_plotseqpeaks(cd1,myxlim)
				
				for(e in 0:nglut){
					for( p in 0:nhyd){
				
						x <- cd1$mass +(e*0.984015)+(p*16)
						y <- cd1$prob
						segments(x0=x, y0=y, y1=0, col=8)
						points(x, y, pch=21, col=1+e, bg=mybg+p)
					}
				}
				
				
				#plot(masses,intensities/max(intensities))
				plot(d1[,1],d1[,2]/max(d1[,2]),type="l",xlim=myxlim, ylim=c(0,1), bty='l', main = "Method K: 1,Norfolk (black); 6,Italy (red)")
				lines(i1[,1],i1[,2]/max(i1[,2]), col = "red")
				
				plot(d2[,1],d2[,2]/max(d2[,2]),type="l",xlim=myxlim, ylim=c(0,1), bty='l', main = "Method I: 1,Norfolk (black); 6,Italy (red)")
				lines(i2[,1],i2[,2]/max(i2[,2]), col = "red")
				
				plot(d3[,1],d3[,2]/max(d3[,2]),type="l",xlim=myxlim, ylim=c(0,1), bty='l', main = "Method C: 1,Norfolk (black), 6,Italy (red)")
				lines(i3[,1],i3[,2]/max(i3[,2]), col = "red")
				
				plot(d4[,1],d4[,2]/max(d4[,2]),type="l",xlim=myxlim, ylim=c(0,1), bty='l', main = "Method B: 1,Norfolk (black), 6,Italy (red)")
				lines(i4[,1],i4[,2]/max(i4[,2]), col = "red")
				
				plot(d5[,1],d5[,2]/max(d5[,2]),type="l",xlim=myxlim, ylim=c(0,1), bty='l', main = "Method N: 1,Norfolk (black), 6,Italy (red)")
				lines(i4[,1],i4[,2]/max(i4[,2]), col = "red")
				
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
}






#ranked_alignment_of_mass_spectrum
rams <- function(sheet, spp, ms1, ms2, ms3){

	aa <- list()


	message(sprintf("Ranking alignement with theoretical spectrum for %s",spp))
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

	moff = 1.5
	
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

			###########################################################
			#RESTRICTION: Only do this for the range we have data for #
			if(max(cd1$mass) > 800 && max(cd1$mass) < 3500){
			
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
					
							
						if(max(subms1[,2]) > (0.1*max(ms1[,2])) ){
							message(sprintf("Max intensity  ratio sufficient in this segment (%f > %f) ", max(subms1[,2]), 0.1*max(ms1[,2]) ))
					
							plotno = plotno+1
			
							message(sprintf("Plot number %d\nSegment at position %d of %d",
									plotno,j-4,endcol-4))
							#0.984015 - if Q changes to E - add this much....
			
							myxlim = c(lbl,ubl)
									
							mymain <- sprintf("%s, nglut = %d, nhyd = %d",sequence,nglut,nhyd)
							plot(1, type="n", xlab="Mass", ylab = "Probability", 
									xlim=myxlim, ylim=c(0, 1), main=mymain)
			
							#ba_plotseqpeaks(cd1,myxlim)
				
							cc = cc+1;
							x <- cd1$mass +(e*0.984015)+(p*16)
							y <- cd1$prob
							segments(x0=x, y0=y, y1=0, col=8)
							points(x, y, pch=21, col=1+e, bg=mybg+p)
						
						
							message(sprintf("Calculating alignment for sequence %s, #%d/%d, nglut: %d, nhyd %d",sequence,cc,nglut*nhyd,e,p))
							
							#plot the peaks
							#ms_peaklineplot <- function(sms,ms,col)
							
							#/Calculate the shift (if any) */
							
							
							ms_peaklineplot(subms1,ms1,"red")
							ms_peaklineplot(subms2,ms2,"green")
							ms_peaklineplot(subms3,ms3,"blue")
						
						
						
						
							readline("hit <return> to look at the next segment\n\n")
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
}















#NB: you need to have run theoretical_spectra.R before calling this!


fn <- "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_A1.txt"
e1 <- read.table(fn)
fn <- "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_D4.txt"
e2 <- read.table(fn)
fn <- "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_G7.txt"
e3 <- read.table(fn)
fn <- "~/tmp/bioarch_keri/20160803_Keri12/20160803_Keri12_0_A1.txt"
e4 <- read.table(fn)
fn <- "~/tmp/bioarch_keri/20160803_Keri12/20160803_Keri12_0_D4.txt"
e5 <- read.table(fn)

fn <- "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_C1.txt"
i1 <- read.table(fn)
fn <- "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_F4.txt"
i2 <- read.table(fn)
fn <- "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_I7.txt"
i3 <- read.table(fn)
fn <- "~/tmp/bioarch_keri/20160803_Keri12/20160803_Keri12_0_C1.txt"
i4 <- read.table(fn)
fn <- "~/tmp/bioarch_keri/20160803_Keri12/20160803_Keri12_0_F4.txt"
i5 <- read.table(fn)


par(mar=c(0.9,2.3,0.9,.3), mfrow = c(1,1), oma=c(5,0,2,0))

#This makes the pdf that was used in the meeting on 161101. 
	pdf(file="keri2.pdf",w=8,h=12)
	#generate_theroretical_spectrum_and_ms
	par(mar=c(0.9,2.3,2.9,.3), mfrow = c(6,1), oma=c(5,0,2,0))

	readline("hit <return> to generate global spectra")

	generate_theroretical_spectrum(sheet,"human")


	col1 = rgb(1,0,0,0.4)
	col2 = rgb(0,0,1,0.4)

	plot(e1[,1],e1[,2],type="l", xlim=c(800, 3500), main = "Method K",col=col2)
	lines(i1[,1],i1[,2],col=col1)
	legend("topright",inset=c(0,0),legend=c("Sample 1, Norfolk","Sample 6, Italy"),col=c(col2,col1),lty=c(1),lwd=c(2))
	#legend("topright",inset=c(0,0.5),legend=lnames,col=col,lty=lty,lwd=lwd,y.intersp=myyint)
	
	plot(e2[,1],e2[,2],type="l", xlim=c(800, 3500), main = "Method I",col=col2)
	lines(i2[,1],i2[,2],col=col1)
	plot(e3[,1],e3[,2],type="l", xlim=c(800, 3500), main = "Method C",col=col2)
	lines(i3[,1],i3[,2],col=col1)
	plot(e4[,1],e4[,2],type="l", xlim=c(800, 3500), main = "Method B",col=col2)
	lines(i4[,1],i4[,2],col=col1)
	plot(e5[,1],e5[,2],type="l", xlim=c(800, 3500), main = "Method N",col=col2)
	lines(i5[,1],i5[,2],col=col1)




	#readline("hit <return> to generate fragment spectra")
	

	#generate_theroretical_spectrum_and_ms
	#gths(sheet,"human",e1,e2,e3,e4,e5,i1,i2,i3,i4,i5)
	dev.off()


#ranked_alignment_of_mass_spectrum
#al <- rams(sheet,"human",e1,e2,e3)

#1: get the TS
#2: for each fragment:
#	figure out nglut and nhyd
#	do an alignment
#	get a score






