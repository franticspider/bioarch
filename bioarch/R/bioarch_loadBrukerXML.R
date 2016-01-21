#TODO: Not sure if there's a better way to do this..
#load the library
require("MALDIquant")


#' load Bruker dataset via xml
#' 
#' @keywords MALDIquant platemap
#' @export
#' @examples
#' data <- bioarch_loadBrukerXML("Bruker.xml")
bioarch_loadBrukerXML <- function(xml){
	mydata = MALDIquantForeign::importBrukerFlex(xml)
	return (mydata)
}
