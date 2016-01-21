
# Function to get the index of a platemap location in a MassSpectrum Object
# Assume the platemap location is after the '.' in data[[i]]@metaData$fullName
# e.g. L9 in '20131016_Oxford_Samples.L9'

#' Get the index of a platemap location in a MALDIquant data object
#' 
#' @param data the MALDIquant object holding the bruker data
#' @param posString the spot postition, e.g. 'C9'
#' @param verbose verbose output if True
#' @keywords MALDIquant
#' @export
#' @examples
#' ba_indexFromPlatePos("/home/anon/brukerdata/20151202",'C9')
ba_indexFromPlatePos <- function(data,posString,verbose=F){
	idx<--1
	#idx <- sapply(1:length(data), function(i){
	for ( i in 1:length(data)){
		#graphics:::lines(x=times(x), y=(e[,i])^2, col=cols[i])
		sname <- data[[i]]@metaData$sampleName
		name  <- data[[i]]@metaData$fullName
		lenpos <- nchar(posString)
		subname <- substring(name,nchar(sname)+2,nchar(sname)+4)
		if(verbose)message(sprintf("Entry %d, posString is %s sname is %s, name is %s,lenpos is %d, subname is %s = %d char",i,posString,sname,name,lenpos,subname,nchar(subname)))
		if(nchar(subname)==nchar(posString)){
			if(posString==subname){
				idx<-i
				if(verbose)
					message(sprintf("Match found at %d: %s %s",i,name,posString))
				break	
			}
		}
	}
	
	if(verbose)
		message(sprintf("Returning %d",idx))

	return(idx)
}

