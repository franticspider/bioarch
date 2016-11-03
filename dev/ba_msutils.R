


ba_msload <- function(){



}




ms_subrange <- function(ms,lbl,ubl){
	subms <- ms[  
		ms[,1] <= ubl & ms[,1] > lbl
		,]		
}




ms_peaklineplot <- function(sms,ms,mycol){
	lines(sms[,1],sms[,2]/max(ms[,2]),col=mycol)
}




ts_index <- function(sheet,spp){

	#TODO: There is surely a more efficient way of doing this...
	 found <- F
	 spidx <- 0
	 for(i in 1:nrow(sheet)){
	 	if(grepl(spp,sheet[i,1],ignore.case=TRUE)){
	 		if(!found){ 
		 		message(sprintf("FOUND: index is %d, species name is %s",i,sheet[i,1]))
		 		found <- T
		 		spidx <- i
		 	}
	 		else{
		 		message(sprintf("Further match found at  index %d, species name is %s - this will be ignored",i,sheet[i,1]))
	 		}
	 	}
	 }
	 if(!found){
	 	message("Match not found, exiting")
	 	return (-1)
	 }
	 return (spidx)

}

