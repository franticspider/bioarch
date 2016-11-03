library("q2e")
library("bioarch")

#need to have the following loaded before loading this: ~/git/bioarch/tests/plotseqpeaks.R
#source("plotseqpeaks.R")

#Specify the `margins' of the sample 
moff = 1


#load the cow sample
cow <- isodists("../q2e/q2e/testdata/cowPeptides")
lbl <- min(cow$mass) - moff
ubl <- max(cow$mass) + moff

#load the q23 library
require(q2e)


#Now load the line data: 
message("Generating peaks from 


3 sequences...")
cd1 <- isodists("/home/sjh/Desktop/Dropbox/bioarch/peptideData/cowpeptides")
cd2 <- isodists("/home/sjh/Desktop/Dropbox/bioarch/peptideData/cowpeptides2")
cd3 <- isodists("/home/sjh/Desktop/Dropbox/bioarch/peptideData/cowpeptides3")

#load the test data
if(!exists("testdata"))
testdata <- loadBrukerXML("/home/sjh/Desktop/sjh/bioarch/SF/SFData/20140307_SF_UPenn164-181/")


# Let's pick an example spot
spot <- 6
message(sprintf("Working on spot %d",spot))

#now we have to dismantle the function plotclass_v2 to get what we are after. 

#set up the plot area: (1 column, 7 rows
par(mar=c(0.9,2.3,0.9,.3), mfrow = c(4,1), oma=c(5,0,2,0))

#Let's do some setting up:
myxlim = c(1192,1198)



#grab the data for the range we are intested in 
		xm <-       testdata[[spot]]@mass[myxlim[1] <= testdata[[spot]]@mass & testdata[[spot]]@mass < myxlim[2] ]
		yi <-  testdata[[spot]]@intensity[myxlim[1] <= testdata[[spot]]@mass & testdata[[spot]]@mass < myxlim[2] ]




#Plot the calculated peaks as probabilities 
#until this is in the library, we have to load it:
#source("plotseqpeaks.R")
plotseqpeaks(cd1,myxlim)
		lines(xm,yi/max(yi),type="l",col= "red")
		
title("Raw data (red), isodists (grey with red dots)", line = "-2")
#plotseqpeaks(cd2,myxlim)
#		lines(xm,yi/max(yi),type="l",col= "red")
#plotseqpeaks(cd3,myxlim)
#		lines(xm,yi/max(yi),type="l",col= "red")

#Plot the data we are testing (i.e. the middle row plot in plotclass_v2)
		#yi[yi<threshold]<-0
		#plot(xm,yi,type="l",col= "red")

#Now we need to generate the correlations data by normalising the sampling frequency


myby = 0.01

#create an interpolation (isodists is accurate to 2 decimal places)
xout = seq(from = myxlim[1], to = myxlim[2], by = myby)


yii <- approx(x=xm,y=yi,xout=xout, method="linear", rule = 2)
yii$y = yii$y/max(yii$y)
plot(yii,type="l",col="red")


#This doesn't work - the interpolation is all wrong - gives convex hull of peaks, but useful for getting structure:
yri <- approx(x=cd1$mass,y=cd1$prob,xout=xout, method="linear", rule = 2)
#set yvals to zero
yri$y[] <-0
#go through each peak
for(i in 1:length(cd1$prob)){
	idx <- which.min(abs(yri$x-cd1$mass[i]))
	yri$y[idx] <- cd1$prob[i]
}
lines(yri)

title(sprintf("Data resampled to resolution %0.2f Da",myby), line = "-2")

mylagmax <- 90

ccd <- ccf(yri$y,yii$y,ylim=c(-0.5,1.0),plot=T,axes=F, lag.max = mylagmax)
	#message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
	#plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )
cor = ccd$acf[,,1] 
lag = ccd$lag[,,1] 
res = data.frame(cor,lag) 
res_max = res[which.max(res$cor),] 

td <- sprintf("Cross-Correlation\n max c=%.2f, at lag=%0.3f",res_max$cor,res_max$lag*myby)
message(td)
title(td, line = "-2")
axis(1, at = c(-50,0,50), labels = c(-0.5,0,0.5))
axis(2)
#text(0,-0.4,td)



#Plot the calculated peaks as probabilities 
#until this is in the library, we have to load it:
#source("plotseqpeaks.R")
plotseqpeaks(cd1,myxlim)
		lines(xm,yi/max(yi),type="l",col= "gray80")
		lines(xm+(res_max$lag*myby),yi/max(yi),type="l",col= "green")
title("Shifted spectrum (in green), original (grey)",line = "-2")








dummy <- function (){

	############################################


	#below we pick a spot with good data:
	tdi <- 24

	#not sure which of these we'll need:
	t1<-testdata[[tdi]]
	t1l<-list(t1)

	#get the subset of the data
	testi <-  t1l[[1]]@intensity[lbl <= t1l[[1]]@mass & t1l[[1]]@mass < ubl ]                   
	testm <-  t1l[[1]]@mass[lbl <= t1l[[1]]@mass & t1l[[1]]@mass < ubl ]

	#Plot the raw data 
	plot(x=testm,y=testi)

	#create an interpolation
	xout = seq(from = lbl, to = ubl, by = 0.02)

	#message("bang!")
	#testint <- approx(x=xi,y=yi,xout=xout, method="linear", rule = 2)
	yri <- approx(x=testm,y=testi,xout=xout, method="linear", rule = 2)

	#plot the interpolated data
	lines(x=yri[[1]],y=yri[[2]])

}
