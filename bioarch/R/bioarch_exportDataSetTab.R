
#' Export Bruker data to files, with lot Number as part of the filename
#' this function will create files either with the file extension .txt, or using the fullname in the Bruker format.
#'
#' @param path the path to the Bruker data set
#' @keywords Bruker 
#' @export
#' @examples
#'
#' 1: to export directly from a Bruker directory:
#' 
#'   > bioarch_exportDataSetTab(bioarch_loadBruker("path/to/bruker/data"))
#'   
#' 2: to export from a Bruker data set that has been loaded into an 
#' object called 'bdata':
#' 
#'   > bioarch_exportDataSetTab(bdata)
#'
#' 3: to export from a Bruker object called 'bdata' to files named
#' after the Bruker 'fullname':
#' 
#'   > bioarch_exportDataSetTab(bdata,txt=T)
bioarch_exportDataSetTab <-function(data,txt=TRUE){

	for(i in 1:length(data)){
		outfn <- sprintf("csv/%s",data[[i]]@metaData$fullName)
		if(txt){
			#replace the '.' character with '_'
			outfn <- gsub("[.]","_",outfn)

			#append the '.txt'
			outfn <- sprintf("%s.txt",outfn)
		}
		message(sprintf("Exporting %s to %s",data[[i]]@metaData$fullName,outfn))
		MALDIquantForeign::exportTab(data[[i]], file = outfn, force=TRUE, col.names = FALSE, sep = "\t")
	}
}


