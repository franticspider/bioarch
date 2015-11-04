#/* Copyright (C) 2015 Simon Hickinbotham	                          */
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


#' Get an arbitrary answer to a question
#' 
#' @param data the path to the Bruker data set
#' @keywords question
#' @export
#' @examples
#' ba_checkyn("how are you feeling?")
exportLotDataSet <-function(path,lotno=0,txt=TRUE){


    message("Loading Bruker Plate Data from file: ",path)
    data<-loadBrukerXML(path)

   
    a<-ba_ynq(sprintf("Is the lot number %06d correct?",lotno))
    if(a=="y"){
    }


    message("\nNow we need to specify which spots belong to lot number ",lotno)



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



































