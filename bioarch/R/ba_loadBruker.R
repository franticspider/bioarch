



#' load a bruker directory structure as a dataset
#' 
#' @param xml the path to the xml file on the computer you are working at
#' @keywords bruker
#' @export
#' @examples
#' ba_loadBruker("/home/anon/brukerdata/20151202")


ba_loadBruker <- function(xml){
	mydata = MALDIquantForeign::importBrukerFlex(xml)
	return (mydata)
}
