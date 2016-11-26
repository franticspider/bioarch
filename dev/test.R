

require("stringr")
library("googlesheets")
library("bioarch")
library("q2e")

source("~/git/bioarch/dev/ba_msutils.R")



#copied from theoretical_spectra.R... 
#Fetch the sheet if it isn't available...
if(!exists("sheet")){

	message("Fetching the sheet of mammalian collagen sequences")
	mcs <- gs_title("Mammalian Collagen Sequences v0.0.1")
	print(mcs)

	sheet <- gs_read(mcs)
}

####################

par(mar=c(3.8,3.8,2.4,3.8), mfrow = c(3,1), oma=c(.1,.1,.2,.2))


#testing 

cor.threshold = 0.15
lag.threshold = 0.225

x <- load.analysis(sprintf("C1analysis.dat"))

t <- x[which(x$cor > cor.threshold & abs(x$lag) < lag.threshold),]

cca <- t[which(t$energy == max(t$energy)),]

cc <- cca[1,]


#############################
# Now need to test how peaksets are moved to the 'true' position


eraw <- read.table(sprintf("%s%s.txt","~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_","G7"))

e1<-eraw

#create the matched peak set array
matchedpeaks <- e1
matchedpeaks[,2]=0

#lbl1 = 1563
#lbl2 = 1562
#ubl1 = 1568
#ubl2 = 1567
lbl1 = 1562.3
lbl2 = 1561.3
ubl1 = 1567.3
ubl2 = 1566.3

plot(e1,type = "l", xlim = c(1560,1570))

########################################################################
#In case of overlap, we have to *ADD* e1 to matchedpeaks...not *REPLACE*
matchedpeaks[  which(e1[,1] > lbl1 & e1[,1] < ubl1),2] <-
#matchedpeaks[  which(e1[,1] > lbl1 & e1[,1] < ubl1),2] 
#+
eraw[          which(e1[,1] > lbl1 & e1[,1] < ubl1),2]

#Set the remaining peakset to zero
e1[            which(e1[,1] > lbl1 & e1[,1] < ubl1),2] <- 0
			

plot(e1,type = "l", xlim = c(1560,1570))
lines(matchedpeaks,col="red")

########################################################################
#In case of overlap, we have to *ADD* e1 to matchedpeaks...not *REPLACE*
matchedpeaks[  which(e1[,1] > lbl2 & e1[,1] < ubl2),2] <-
#matchedpeaks[  which(e1[,1] > lbl2 & e1[,1] < ubl2),2] 
#+
eraw[          which(e1[,1] > lbl2 & e1[,1] < ubl2),2]

#Set the remaining peakset to zero
e1[            which(e1[,1] > lbl2 & e1[,1] < ubl2),2] <- 0


plot(e1,type = "l", xlim = c(1560,1570))
lines(matchedpeaks,col="red")



##########################################################################
#testing unigue - to get unique sequences out...
