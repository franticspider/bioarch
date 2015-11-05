#/* Copyright (C) 2014 Simon Hickinbotham	                          */
#/* When you use this, send an email to: simon.hickinbotham@york.ac.uk    */
#/* with an appropriate reference to your work.                           */

#/* This file is part of flamePlot version 1.0.0  		  	  */

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





#load the library
library("MALDIquant")

#source("exportBruker.R")



mydata<-NULL

#TODO: write this function
printUsage_flame <- function(){
message("=============================================================")
message("This file is part of flamePlot version 1.0.0")
message("Typical usage is:\n")

message("  > plotRuminants(PATH,SPOT)")
message("\n..where PATH is the file path to the Bruker source directory\n")
message("This program has created a subdirectory called \'csv\',\nwhere the exported text files will be placed")
message("=============================================================")
}


lb <-c(1105,1180,1192,1196,1208,1427,1580,1648,2131,2792,2853)
ub <-c(1108,1184,1198,1202,1212,1431,1584,1652,2135,2796,2857)


plotRuminants<-function(fn,data,spot){

	pdf(fn,width=14,height=10)

		np <- length(lb)
		ns <- length(spot)
		#par(mfrow=c(ns,np))
		par(mar=c(1.7,1.7,0.3,.3), mfrow = c(ns,np), oma=c(5,0,2,0))

		for(s in 1:ns){

			idx <- indexFromPlatePos(data,spot[s])

			for(i in 1:np){
		
				xm <-       data[[idx]]@mass[lb[i] <= data[[idx]]@mass & data[[idx]]@mass < ub[i] ]
				yi <-  data[[idx]]@intensity[lb[i] <= data[[idx]]@mass & data[[idx]]@mass < ub[i] ]
		
				plot(xm,yi,type="l",ylim=c(0,max(1000,max(yi)))  )

	
			}
		}

	dev.off()

}



Find_Max_CCF<- function(a,b) 
{ 
 d <- ccf(a, b, plot = FALSE) 
 cor = d$acf[,,1] 
 lag = d$lag[,,1] 
 res = data.frame(cor,lag) 
 res_max = res[which.max(res$cor),] 
 return(res_max) 
}



#PROOF OF CONCEPT - CCF AGAINST A SET OF SPOTS (WE'LL NEED REFERENCE DATA TO MAKE THIS WORK)
plotccfs<-function(fn,data,spot,threshold){

	pdf(fn,width=14,height=10)

		np <- length(lb)
		ns <- length(spot)
		#par(mfrow=c(ns,np))
		par(mar=c(0.3,0.3,0.3,.3), mfrow = c(ns,np), oma=c(5,0,2,0))

		#ASSUME SPOT 1 IS THE REFERENCE
		idr <- indexFromPlatePos(data,spot[1])		

		for(s in 1:ns){

			idx <- indexFromPlatePos(data,spot[s])

			for(i in 1:np){
		
				#inefficient, but simple:
				yr <-  data[[idr]]@intensity[lb[i] <= data[[idr]]@mass & data[[idr]]@mass < ub[i] ]

				yi <-  data[[idx]]@intensity[lb[i] <= data[[idx]]@mass & data[[idx]]@mass < ub[i] ]
		

				yr[yr<threshold]<-0
				yi[yi<threshold]<-0


				#plot(xm,yi,type="l",ylim=c(0,max(1000,max(yi)))  )

				if(max(yr) > 0 && max(yi) >0){
					ccd <- ccf(yr,yi,ylim=c(-0.5,1.0))
					cor = ccd$acf[,,1] 
					lag = ccd$lag[,,1] 
					res = data.frame(cor,lag) 
					res_max = res[which.max(res$cor),] 

					td <- sprintf("c=%0.3f, lag=%0.0f",res_max$cor,res_max$lag)
					message(td)
					text(0,-0.4,td)

					message(sprintf("Between mass %d and %d, correlation for %s and %s is %f, lag is %f",lb[i],ub[i],spot[1],spot[s],res_max$cor,res_max$lag))
				}
				else{
					plot.new()
					text(0,-0.4,"empty array")
				}
			}
		}

	dev.off()

}



#PROOF OF CONCEPT - CCF AGAINST A SET OF SPOTS (WE'LL NEED REFERENCE DATA TO MAKE THIS WORK)
plotflame<-function(fn,data,spot,threshold,uselag=TRUE){

	pdf(fn,width=14,height=10)

		np <- length(lb)
		ns <- length(spot)
		#par(mfrow=c(ns,np))
		par(mar=c(0.6,0.3,0.6,.3), mfrow = c(ns,np), oma=c(5,0,2,0))

		#ASSUME SPOT 1 IS THE REFERENCE
		idr <- indexFromPlatePos(data,spot[1])		

		for(s in 1:ns){

			idx <- indexFromPlatePos(data,spot[s])

			for(i in 1:np){
		
				#inefficient, but simple:
				yr <-  data[[idr]]@intensity[lb[i] <= data[[idr]]@mass & data[[idr]]@mass < ub[i] ]

				yi <-  data[[idx]]@intensity[lb[i] <= data[[idx]]@mass & data[[idx]]@mass < ub[i] ]
		

				yr[yr<threshold]<-0
				yi[yi<threshold]<-0


				

				if(max(yr) > 0 && max(yi) >0 && length(yr) == length(yi) ){
					ccd <- ccf(yr,yi,ylim=c(-0.5,1.0),plot=FALSE)
					message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
					plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )
					cor = ccd$acf[,,1] 
					lag = ccd$lag[,,1] 
					res = data.frame(cor,lag) 
					res_max = res[which.max(res$cor),] 

					td <- sprintf("c=%0.3f, lag=%0.0f",res_max$cor,res_max$lag)
					message(td)
					text(0,-0.4,td)

					message(sprintf("Between mass %d and %d, correlation for %s and %s is %f, lag is %f",lb[i],ub[i],spot[1],spot[s],res_max$cor,res_max$lag))
				}
				else{
					plot.new()
					text(0,-0.4,"empty array")
				}
			}
		}

	dev.off()

}



#PROOF OF CONCEPT - CCF AGAINST A SET OF SPOTS (WE'LL NEED REFERENCE DATA TO MAKE THIS WORK)
plotcomp<-function(fn,refdata,testdata,threshold,uselag=TRUE,thislb,thisub){

	pdf(fn,width=14,height=10)

		np <- length(thislb)
		ns <- 5# 5 rows: ref, data, ccf, flame, text
		#par(mfrow=c(ns,np))
		par(mar=c(0.9,0.3,0.9,.3), mfrow = c(ns,np), oma=c(5,0,2,0))

		#ASSUME SPOT 1 IS THE REFERENCE
		#idr <- indexFromPlatePos(data,spot[1])		

		#PLOT THE REFERENCE DATA
		#idx <- indexFromPlatePos(data,spot[1])
		message("Plotting the reference peaks...")
		for(i in 1:np){
	
			xm <-  refdata[[1]]@mass[thislb[i] <= refdata[[1]]@mass & refdata[[1]]@mass < thisub[i] ]
			yi <-  refdata[[1]]@intensity[thislb[i] <= refdata[[1]]@mass & refdata[[1]]@mass < thisub[i] ]
	
			plot(xm,yi,type="l",ylim=c(0,max(1000,max(yi)))  )
		}
		message("...done")

		#PLOT THE SAMPLE DATA
		#idx <- indexFromPlatePos(data,spot[2])
		for(i in 1:np){
	
			xm <-       testdata[[1]]@mass[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]
			yi <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]
	
			yi[yi<threshold]<-0
			plot(xm,yi,type="l",ylim=c(0,max(1000,max(yi)))  )

		}


		#PLOT THE CORRELATION
		for(i in 1:np){
	
			#inefficient, but simple:
			yr <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]

			yi <-  refdata[[1]]@intensity[thislb[i] <= refdata[[1]]@mass & refdata[[1]]@mass < thisub[i] ]
	
			yr[yr<threshold]<-0
			yi[yi<threshold]<-0

			if(max(yr) > 0 && max(yi) >0){
				ccd <- ccf(yr,yi,ylim=c(-0.5,1.0),plot=TRUE)
				message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
				#plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )
				cor = ccd$acf[,,1] 
				lag = ccd$lag[,,1] 
				res = data.frame(cor,lag) 
				res_max = res[which.max(res$cor),] 

				td <- sprintf("c=%0.3f, lag=%0.0f",res_max$cor,res_max$lag)
				message(td)
				text(0,-0.4,td)

				message(sprintf("Between mass %d and %d, correlation for ref and test is %f, lag is %f",thislb[i],thisub[i],res_max$cor,res_max$lag))
			}
			else{
				plot.new()
				text(0,-0.4,"empty array")
			}
		}


		for(i in 1:np){
	
			#inefficient, but simple:
			yr <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]

			yi <-  refdata[[1]]@intensity[thislb[i] <= refdata[[1]]@mass & refdata[[1]]@mass < thisub[i] ]
	

			yr[yr<threshold]<-0
			yi[yi<threshold]<-0

			yr = yr/max(yr)
			yi = yi/max(yi)

			minl = min(length(yr),length(yi))
			yr = yr[1:minl]
			yi = yi[1:minl]


			#ccd <- ccf(yr,yi,ylim=c(-0.5,1.0),plot=FALSE)
			message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
			plot(yr,yi,type="l",xlim=c(0,1),ylim=c(0,1)  )
		}

		#SHIFT THE FLAMEPLOT ACCORDING TO THE LAG
		for(i in 1:np){
	
			#inefficient, but simple:
			yr <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]

			yi <-  refdata[[1]]@intensity[thislb[i] <= refdata[[1]]@mass & refdata[[1]]@mass < thisub[i] ]
	

			yr[yr<threshold]<-0
			yi[yi<threshold]<-0

			if(max(yr) > 0 && max(yi) >0){
				yr = yr/max(yr)
				yi = yi/max(yi)

				minl = min(length(yr),length(yi))
				yr = yr[1:minl]
				yi = yi[1:minl]


				ccd <- ccf(yr,yi,ylim=c(-0.5,1.0),plot=FALSE)
				cor = ccd$acf[,,1] 
				lag = ccd$lag[,,1] 
				res = data.frame(cor,lag) 
				res_max = res[which.max(res$cor),] 

				lag = res_max$lag
				message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
				if(lag<0){
					message(sprintf("Negative lag: %f",lag))
					lag = abs(lag)
					lyi = length(yi)

					yi=yi[(lag+1):lyi]
					bbb = length(yr)-(lag)

					message(sprintf("bbb = %f",bbb))

					yr=yr[1:bbb]
					#message(sprintf("yi from %d to %d,\nyr from %d to %d",lag+1,lyi,1,bbb))
				
				}
				else if(lag>0){
					message("Positive lag")
					lag = abs(lag)
					lyr = length(yr)

					yr=yr[(lag+1):lyr]

					bbb = length(yi)-(lag)

					message(sprintf("bbb = %f",bbb))

					yi=yi[1:bbb]
					#message(sprintf("yi from %d to %d,\nyr from %d to %d",lag+1,lyi,1,bbb))
			
				}


				message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
				plot(yr,yi,type="l",xlim=c(0,1),ylim=c(0,1)  )
			}
			else{
				plot.new()
				text(0,-0.4,"empty array")
			}
		}


		

	dev.off()
	message("Finished!")

}








#PROOF OF CONCEPT - CCF AGAINST A SET OF SPOTS (WE'LL NEED REFERENCE DATA TO MAKE THIS WORK)
plotclass<-function(fn,refdata,testdata,threshold,uselag=TRUE,thislb,thisub,lbl,ubl){

	pdf(fn,width=14,height=10)

		np <- length(thislb)
		ns <- 5# 5 rows: ref, data, ccf, flame, text
		#par(mfrow=c(ns,np))
		par(mar=c(0.9,0.3,0.9,.3), mfrow = c(ns,np), oma=c(5,0,2,0))

		#ASSUME SPOT 1 IS THE REFERENCE
		#idr <- indexFromPlatePos(data,spot[1])		

		#PLOT THE REFERENCE DATAhttps://www.facebook.com/
		#idx <- indexFromPlatePos(data,spot[1])
		message("Plotting the reference peaks...")
		for(b in 1:length(lbl)){
			for(i in 1:np){
				if(thislb[i] %in% lbl[[b]]){
					xm <-  refdata[[1]]@mass[thislb[i] <= refdata[[1]]@mass & refdata[[1]]@mass < thisub[i] ]
					yi <-  refdata[[1]]@intensity[thislb[i] <= refdata[[1]]@mass & refdata[[1]]@mass < thisub[i] ]
	
					plot(xm,yi,type="l",ylim=c(0,max(1000,max(yi)))  )
				}
				else{
					plot.new()
				}
			}
		}
		message("...done")

		#PLOT THE SAMPLE DATA
		#idx <- indexFromPlatePos(data,spot[2])
		for(i in 1:np){
	
			xm <-       testdata[[1]]@mass[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]
			yi <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]
	
			yi[yi<threshold]<-0
			plot(xm,yi,type="l",ylim=c(0,max(1000,max(yi)))  )

		}


		#PLOT THE CORRELATION
		for(i in 1:np){
	
			#inefficient, but simple:
			yr <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]

			yi <-  refdata[[1]]@intensity[thislb[i] <= refdata[[1]]@mass & refdata[[1]]@mass < thisub[i] ]
	
			yr[yr<threshold]<-0
			yi[yi<threshold]<-0

			if(max(yr) > 0 && max(yi) >0){
				ccd <- ccf(yr,yi,ylim=c(-0.5,1.0),plot=TRUE)
				message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
				#plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )
				cor = ccd$acf[,,1] 
				lag = ccd$lag[,,1] 
				res = data.frame(cor,lag) 
				res_max = res[which.max(res$cor),] 

				td <- sprintf("c=%0.3f, lag=%0.0f",res_max$cor,res_max$lag)
				message(td)
				text(0,-0.4,td)

				message(sprintf("Between mass %d and %d, correlation for ref and test is %f, lag is %f",thislb[i],thisub[i],res_max$cor,res_max$lag))
			}
			else{
				plot.new()
				text(0,-0.4,"empty array")
			}
		#}


	
	}
	dev.off()
	message("Finished!")
}













#USAGE: 
#data <- loadBruker("20131024_TedG1/20131024_TedG1")
#spots <- c("C4","C9","C4","C8")
#plotRuminants("20131024_TedG1_A2_C9_C4_C8.pdf",data,spots)
#plotccfs("20131024_TedG1_A2_C9_C4_C8_ccf.pdf",data,spots,100)
#plotflame("20131024_TedG1_A2_C9_C4_C8_flame.pdf",data,spots,100)


#spots <- c("C4","C8")
#plotcomp("20131024_TedG1_C4vsC8.pdf",data,spots,150)



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
#' plotclass_v2() [TODO: fill this example in!]
#PROOF OF CONCEPT - CCF AGAINST A SET OF SPOTS (WE'LL NEED REFERENCE DATA TO MAKE THIS WORK)
plotclass_v2<-function(fn,onepdf=F,refdata,testdata,threshold,uselag=TRUE,thislb,thisub,lbl,ubl,classnames){

	mylagmax = 50


	if(onepdf==F){
		pdf(fn,width=14,height=10)
	}

	#create data structures to hold the correlation and lag values:
	cordata <- lbl
	lagdata <- lbl

	#get the number of masses we are testing against:
	np <- length(thislb) 

	#readline(sprintf("np = %d, hit <return> to continue",np))


	#GET THE CORRELATION DATA
	for(b in 1:length(lbl)){
		#create an index into the data list:
		j<-1

		for(i in 1:np){
			##inefficient, but simple:
			#message(sprintf("Getting correlation %d for %d",i,j))

#			yr <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]

#			yi <-  refdata[[b]][[1]]@intensity[thislb[i] <= refdata[[b]][[1]]@mass & refdata[[b]][[1]]@mass < thisub[i] ]

			yr <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]                   
			xr <-  testdata[[1]]@mass[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]

			yi <-  refdata[[b]][[1]]@intensity[thislb[i] <= refdata[[b]][[1]]@mass & refdata[[b]][[1]]@mass < thisub[i] ]
			xi <-  refdata[[b]][[1]]@mass[thislb[i] <= refdata[[b]][[1]]@mass & refdata[[b]][[1]]@mass < thisub[i] ]

			yr[yr<threshold]<-0
			yi[yi<threshold]<-0

			#set a noise floor
			yr[yr<threshold]<-0
			yi[yi<threshold]<-0

			
			#message("biff!")
			#create an interpolation
			xout = seq(from = thislb[i], to = thisub[i], by = 0.02)

			#message("bang!")
			yii <- approx(x=xi,y=yi,xout=xout, method="linear", rule = 2)
			yri <- approx(x=xr,y=yr,xout=xout, method="linear", rule = 2)


			if(thislb[i] %in% lbl[[b]]){
				#message("doing ccd now...")
				#ccd <- ccf(yr,yi,ylim=c(-0.5,1.0),plot=F, lag.max = mylagmax)
				ccd <- ccf(yri$y,yii$y,ylim=c(-0.5,1.0),plot=F, lag.max = mylagmax)

				cor = ccd$acf[,,1] 
				lag = ccd$lag[,,1] 
				res = data.frame(cor,lag) 
				res_max = res[which.max(res$cor),] 
				#message(res_max)
				td <- sprintf("c=%.2f, lag=%0.0f",res_max$cor,res_max$lag)
				#message(td)

				#message(sprintf("Between mass %d and %d, correlation for ref and test is %0.2f, lag is %0.f",thislb[i],thisub[i],res_max$cor,res_max$lag))
				#flush.console()
				#TODO: test if res_max has data!
				cordata[[b]][j] <- res_max$cor
				lagdata[[b]][j] <- res_max$lag

				j<-j+1
			}
			else{
				#message("skipping this range")
			}
		}
	}



	#readline(sprintf("finished calculating correlations, hit <return> to continue"))



	#get the number of masses we are testing against:
	np <- length(thislb)



	ns <- 7# 5 rows: ref, data, ccf, flame, text
	#par(mfrow=c(ns,np))
	par(mar=c(0.9,0.3,0.9,.3), mfrow = c(ns,np+1), oma=c(5,0,2,0))

	#ASSUME SPOT 1 IS THE REFERENCE
	#idr <- indexFromPlatePos(data,spot[1])		



	#PLOT THE REFERENCE DATA
	#idx <- indexFromPlatePos(data,spot[1])
	#message("Plotting the reference peaks...")
	for(b in 1:length(lbl)){
		for(i in 1:np){
			#message(sprintf("plotting %d freq %d",b,i))
			if(thislb[i] %in% lbl[[b]]){
				xm <-  refdata[[b]][[1]]@mass[     thislb[i] <= refdata[[b]][[1]]@mass & refdata[[b]][[1]]@mass < thisub[i] ]
				yi <-  refdata[[b]][[1]]@intensity[thislb[i] <= refdata[[b]][[1]]@mass & refdata[[b]][[1]]@mass < thisub[i] ]

				#plot(xm,yi,type="l",ylim=c(0,max(1000,max(yi))) , axes=F )
				plot(xm,yi,type="l", axes=F )
			}
			else{
				plot.new()
			}
			if(i==1){
				message(sprintf("Printing %s",classnames[b]))
				#text(xm+1,max(1000,max(yi))/2, classnames[b], col="red")
				mtext(classnames[b], side = 2,  col="red", padj = 1)
			}
		}
		plot(x=lbl[[b]],y=cordata[[b]],type='b',ylim=c(0,1))
	}
	#message("...done")


	mtext(fn, outer = TRUE, cex = 1.5)


	#PLOT THE SAMPLE DATA
	#idx <- indexFromPlatePos(data,spot[2])
	for(i in 1:np){

		xm <-       testdata[[1]]@mass[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]
		yi <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]

		yi[yi<threshold]<-0
		plot(xm,yi,type="l",col= "red")#, ylim=c(0,max(1000,max(yi)))  )

	}


	suppressWarnings(boxplot(cordata[[1]],cordata[[2]],cordata[[3]],notch=T,ylim=c(0,1),names=classnames,las=2))


	#PLOT THE CORRELATION
	for(b in 1:length(lbl)){
		for(i in 1:np){

			#inefficient, but simple:
#			yr <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]

#			yi <-  refdata[[b]][[1]]@intensity[thislb[i] <= refdata[[b]][[1]]@mass & refdata[[b]][[1]]@mass < thisub[i] ]

			yr <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]                   
			xr <-  testdata[[1]]@mass[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]

			yi <-  refdata[[b]][[1]]@intensity[thislb[i] <= refdata[[b]][[1]]@mass & refdata[[b]][[1]]@mass < thisub[i] ]
			xi <-  refdata[[b]][[1]]@mass[thislb[i] <= refdata[[b]][[1]]@mass & refdata[[b]][[1]]@mass < thisub[i] ]

			yr[yr<threshold]<-0
			yi[yi<threshold]<-0

			#message("biff!")
			#create an interpolation
			xout = seq(from = thislb[i], to = thisub[i], by = 0.02)

			#message("bang!")
			yii <- approx(x=xi,y=yi,xout=xout, method="linear", rule = 2)
			yri <- approx(x=xr,y=yr,xout=xout, method="linear", rule = 2)



			if(max(yr) > 0 && max(yi) >0 && thislb[i] %in% lbl[[b]]){
				#ccd <- ccf(yr,yi,ylim=c(-0.5,1.0),plot=TRUE,axes=F, lag.max = mylagmax)#, lag.max = 10)
				ccd <- ccf(yri$y,yii$y,ylim=c(-0.5,1.0),plot=T,axes=F, lag.max = mylagmax)
				#message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
				#plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )
				cor = ccd$acf[,,1] 
				lag = ccd$lag[,,1] 
				res = data.frame(cor,lag) 
				res_max = res[which.max(res$cor),] 

				td <- sprintf("c=%.2f, lag=%0.0f",res_max$cor,res_max$lag)
				#message(td)
				text(0,-0.4,td)

				#message(sprintf("Between mass %d and %d, correlation for ref and test is %f, lag is %f",thislb[i],thisub[i],res_max$cor,res_max$lag))
			}
			else{
				plot.new()
				text(0,-0.4,"empty array")
			}
		}
		plot(x=lbl[[b]],y=lagdata[[b]],type='l')
	}

	if(onepdf==F){
		dev.off()
		message("Finished!")
	}

	return(cordata)
}









#PROOF OF CONCEPT - CCF AGAINST A SET OF SPOTS (WE'LL NEED REFERENCE DATA TO MAKE THIS WORK)
plotclass_v2bits<-function(fn,onepdf=F,refdata,testdata,threshold,uselag=TRUE,thislb,thisub,lbl,ubl,classnames){


	if(onepdf==F){
		pdf(fn,width=14,height=10)
	}

	#create data structures to hold the correlation and lag values:
	cordata <- lbl
	lagdata <- lbl

	#get the number of masses we are testing against:
	np <- length(thislb) 

	#readline(sprintf("np = %d, hit <return> to continue",np))



	#readline(sprintf("finished calculating correlations, hit <return> to continue"))



	#get the number of masses we are testing against:
	np <- length(thislb)



	ns <- 7# 5 rows: ref, data, ccf, flame, text
	#par(mfrow=c(ns,np))
	par(mar=c(0.9,0.3,0.9,.3), mfrow = c(ns,np+1), oma=c(5,0,2,0))

	#ASSUME SPOT 1 IS THE REFERENCE
	#idr <- indexFromPlatePos(data,spot[1])		



	#PLOT THE REFERENCE DATA
	#idx <- indexFromPlatePos(data,spot[1])
	#message("Plotting the reference peaks...")
	for(b in 1:length(lbl)){
		for(i in 1:np){
			#message(sprintf("plotting %d freq %d",b,i))
			if(thislb[i] %in% lbl[[b]]){
				xm <-  refdata[[b]][[1]]@mass[     thislb[i] <= refdata[[b]][[1]]@mass & refdata[[b]][[1]]@mass < thisub[i] ]
				yi <-  refdata[[b]][[1]]@intensity[thislb[i] <= refdata[[b]][[1]]@mass & refdata[[b]][[1]]@mass < thisub[i] ]

				plot(xm,yi/max(yi),type="l", axes=F)#ylim=c(0,max(1000,max(yi))))
			}
			else{
				plot.new()
			}
			if(i==1){
				message(sprintf("Printing %s",classnames[b]))
				#text(xm+1,max(1000,max(yi))/2, classnames[b], col="red")
				mtext(classnames[b], side = 2,  col="red", padj = 1)
			}
		}
		plot(x=lbl[[b]],y=cordata[[b]],type='b',ylim=c(0,1))
	}
	#message("...done")


	mtext(fn, outer = TRUE, cex = 1.5)


	#PLOT THE SAMPLE DATA
	#idx <- indexFromPlatePos(data,spot[2])
	for(i in 1:np){

		xm <-       testdata[[1]]@mass[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]
		yi <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]

		yi[yi<threshold]<-0
		plot(xm,yi,type="l",col= "red", ylim=c(0,max(1000,max(yi)))  )

	}


	suppressWarnings(boxplot(cordata[[1]],cordata[[2]],cordata[[3]],notch=T,ylim=c(0,1),names=classnames,las=2))


	#PLOT THE CORRELATION
	for(b in 1:length(lbl)){
		for(i in 1:np){

			#inefficient, but simple:
#			yr <-  testdata[[1]]@intensity[thislb[i] <= testdata[[1]]@mass & testdata[[1]]@mass < thisub[i] ]

#			yi <-  refdata[[b]][[1]]@intensity[thislb[i] <= refdata[[b]][[1]]@mass & refdata[[b]][[1]]@mass < thisub[i] ]

#			yr[yr<threshold]<-0
#			yi[yi<threshold]<-0

			#if(max(yr) > 0 && max(yi) >0 && thislb[i] %in% lbl[[b]]){
			#	ccd <- ccf(yr,yi,ylim=c(-0.5,1.0),plot=TRUE,axes=F)#, lag.max = 10)
			#	#message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
			#	#plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )
			#	cor = ccd$acf[,,1] 
			#	lag = ccd$lag[,,1] 
			#	res = data.frame(cor,lag) 
			#	res_max = res[which.max(res$cor),] 
#
#				#td <- sprintf("c=%.2f, lag=%0.0f",res_max$cor,res_max$lag)
#				#message(td)
#				text(0,-0.4,td)
#
#				#message(sprintf("Between mass %d and %d, correlation for ref and test is %f, lag is %f",thislb[i],thisub[i],res_max$cor,res_max$lag))
#			}
#			else{
				plot.new()
#				text(0,-0.4,"empty array")
#			}
		}
		plot.new()#plot(x=lbl[[b]],y=lagdata[[b]],type='l')
	}

	if(onepdf==F){
		dev.off()
		message("Finished!")
	}

	return(cordata)
}


















