#/* Copyright (C) 2014 Simon Hickinbotham, Matthew Collins, Keri Roswell  */
#/* When you use this, send an email to: simon.hickinbotham@york.ac.uk    */
#/* with an appropriate reference to your work.                           */

#/* This file is part of exportBruker version 1.0.2  		  	  */

#/* exportBruker is free software: you can redistribute it and/or modify  */
#/* it under the terms of the GNU General Public License as published by  */
#/* the Free Software Foundation, either version 3 of the License, or     */
#/* (at your option) any later version.                                   */

#/* This program is distributed in the hope that it will be useful,       */
#/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
#/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
#/* GNU General Public License for more details.                          */

#/* You should have received a copy of the GNU General Public License     */
#/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */


library("MALDIquantForeign")
library("bioarch")
#source("exportBruker.R")
#source("flamePlot_3.R")







#' Plot a classification of a Bruker data set
#' 
#' @param fn the pdf file name
#' @param onepdf whether to print the analysis into a single file or not
#' @param refdata the reference data that the classification is based on
#' @param testdata the data that is being classified
#' @param threshold classification threshold
#' @param uselag whether to calculate the lag
#' @param thislb lower bound 1
#' @param thisub upper bound 1
#' @param lbl lower bound limit
#' @param ubl upper bound limit
#' @param classnames the names of the classes 
#' @keywords question
#' @export
#' @examples
#' plot_zooms_v1       (fn,onepdf=T,refdata,t1l,     0,        uselag=TRUE,LB,    UB,    lbl,ubl,classnames){
#PROOF OF CONCEPT - CCF AGAINST A SET OF SPOTS (WE'LL NEED REFERENCE DATA TO MAKE THIS WORK)
plot_zooms_v1<-function(fn,onepdf=F,refdata,testdata,threshold,uselag=TRUE,thislb,thisub,lbl,ubl,classnames){


	#create data structures to hold the correlation and lag values:
	cordata <- lbl
	lagdata <- lbl




}




















getclass <- function(idx,fn){
	name = "unknown"
	t <- read.table(fn,sep='\t')
	#'b' %in% v
	if(idx %in% t[,1]){
		name<-t[grep(idx,t[,1]),3]
		message(sprintf("found a %s entry at %s",name,idx))
	}
	return(name)

}




moff = -1
mrange = 4

message("Loading known cow data")
#fallowG9<- import('../FallowDeer/Fallow deer_0_G9.txt')
#fallowH9<- import('../FallowDeer/Fallow deer_0_H9.txt')
#fallowI9<- import('../FallowDeer/Fallow deer_0_I9.txt')

if(!exists("cow"))
cow <- import("/home/sjh/Desktop/sjh/bioarch/SF/SFData/averaged/Averaged Cow Spectrum.txt")



message("Loading ref vals for sheep, goat and cow")
massTable<-read.table('/home/sjh/Desktop/sjh/bioarch/SF/SFData/refmasses.tsv')

cowTable<-massTable[grep("Cattle",massTable[,1]),]
cowVals<-as.integer(cowTable[,3:ncol(cowTable)])
cowLB<-cowVals+moff
cowUB<-cowVals+mrange


message("Creating combined data values")

LB = sort(unique(cowLB))
UB = sort(unique(cowUB))
lbl = list(cowLB)
ubl = list(cowUB)

 
if(!exists("testdata"))
testdata <- bioarch_loadBrukerXML("/home/sjh/Desktop/sjh/bioarch/SF/SFData/20140307_SF_UPenn164-181/")

t1 <- testdata[[1]]
t1l <- list(t1)

#refdata <- list(cow)

#Need to sort the mapping from aa sequence to chemical formula - 
#Is the charge always 1?
#When does deamidation occur?
#When does oxidation occur?
#Do we generate formulae for ALL these, or is there an obvious subset? 



sequences=c(
#1191 
"IGQPGAVGPAGIR",
#1207 
"IGQPGAVGPAGIR", #need to add an oxygen...
#1426 
"GIPGEFGLPGPAGAR",
#1579 
"GEPGPAGAVGPAGAVGPR",
#2852 
"", #this is missing!
#3016 
   "GPSGEPGTAGPPGTPGPQGLLGAPGFLGLPGSR",
#or GPPGASGAPGPQGFQGPPGEPGEPGQTGPAGAR
#3032
   "GPSGEPGTAGPPGTPGPQGLLGAPGFLGLPGSR"
#or GPPGASGAPGPQGFQGPPGEPGEPGQTGPAGAR
)





fn = sprintf("/home/sjh/plotcow_theory.pdf",t1l[[1]]@metaData$name)

pdf(file=fn,width = 14,height=8)

for(ii in 1:10){#length(testdata)){
	

	t1<-testdata[[ii]]
	t1l<-list(t1)


	#get the platemap location index
	iname <- substr(t1l[[1]]@metaData$fullName,nchar(t1l[[1]]@metaData$sampleName)+2,nchar(testdata[[1]]@metaData$sampleName)+4)
	inum <-strtoi(substr(iname,2,3), base = 0L)
	in2 <- sprintf("%s%02d",substr(iname,1,1),inum)
	message(sprintf("iname = %s, inum = %d, in2 = %s",iname,inum,in2))

	class <- getclass(in2,"/home/sjh/Desktop/sjh/bioarch/SF/SFData/platemap_classes.txt")
	message(sprintf("Outside of getclass, name for %s is %s",in2,class))


	ffn = sprintf("%s",t1l[[1]]@metaData$name)
	message(sprintf("%02d: Processing %s...",ii,ffn))

	fn = sprintf("Classification of %s (%s)",t1l[[1]]@metaData$name,class)


	classnames=c("Sheep","Goat","Cattle")

	plot_zooms_v1(fn,onepdf=T,refdata,t1l,0,uselag=TRUE,LB,UB,lbl,ubl,classnames)

}
dev.off()





#spots <- c("C4","C9","C4","C8")
#plotRuminants("20131024_TedG1_A2_C9_C4_C8.pdf",data,spots)
#plotccfs("20131024_TedG1_A2_C9_C4_C8_ccf.pdf",data,spots,100)
#plotflame("20131024_TedG1_A2_C9_C4_C8_flame.pdf",data,spots,100)


#spots <- c("C4","C8")
#plotcomp("20131024_TedG1_C4vsC8.pdf",fallowG9,list(testdata[[1]]),spots,150)
























