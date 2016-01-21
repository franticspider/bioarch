
# Function to get the index of a platemap location in a MassSpectrum Object
# Assume the platemap location is after the '.'
# e.g. L9 in '20131016_Oxford_Samples.L9'
indexFromPlatePos <- function(data,posString){
	idx<--1
	#idx <- sapply(1:length(data), function(i){
	for ( i in 1:length(data)){
		#graphics:::lines(x=times(x), y=(e[,i])^2, col=cols[i])
		sname <- data[[i]]@metaData$sampleName
		name  <- data[[i]]@metaData$fullName
		lenpos <- nchar(posString)
		subname <- substring(name,nchar(sname)+2,nchar(sname)+4)
		#message(sprintf("sname is %s, name is %s, lenpos is %d, subname is %s = %d char",sname,name,lenpos,subname,nchar(subname)))
		if(nchar(subname)==nchar(posString)){
			if(posString==subname){
				idx<-i
				message(sprintf("Match found at %d: %s %s",i,name,posString))
				break	
			}
		}
	}

	message(sprintf("Returning %d",idx))

	return(idx)
}

