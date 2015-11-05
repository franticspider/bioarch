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


#' Export Bruker data to files, with lot Number as part of the filename
#' 
#' @param path the path to the Bruker data set
#' @param lotno the lot number of the data set (checked in the function)
#' @keywords Bruker 
#' @export
#' @examples
#' exportLotData("/home/anon/Bruker/20151105sample",123456)
exportLotData <-function(path,lotno=0,txt=TRUE){


    message("Loading Bruker Plate Data from file: ",path)
    data<-loadBrukerXML(path)

   
    a<-ba_ynq(sprintf("Is the lot number %06d correct?",lotno))
    if(a=="n"){
    	lotno <- ba_checkyn("Enter the lot number (eg \"123456\")")
    	#TODO: Check that this is valid
    }


    message("\nNow we need to specify which spots belong to lot number ",lotno)
    
    
    message("Please enter at least one spot position")
    finspot = F
    nspots = 0
    while(!finspot){
    
    	if(nspots==0){#INITIALISE
    		sp <- ba_checkyn("Enter the first spot position (eg \"M9\")")
    		spots = c(sp)
    	}
    	else{
    		sp <- ba_checkyn("Enter the next spot position (eg \"M10\")")
    		spots = c(spots,sp)	
    	}
    	nspots = nspots+1
    
    	message("so far we have ")
    	ss <- print(spots)
    	a<-ba_ynq("Is that all the spots")
    	if(a=="y")
    		finspot=T
    }

	#make sure we have a csv directory...
	dir.create(file.path(getwd(), "csv"), showWarnings = FALSE)


	#TODO: make a function out of this, pass spots in...
	for(spot in 1:length(spots)){
	
		i <- indexFromPlatePos(data,spots[spot])
	
		if(i==-1){
			message(sprintf("Plate position %s was not found in this data set",spots[spot]))
		}
		else{
	
			message(sprintf("Index of plate position %s is %d",spots[spot],i))
			outfn <- sprintf("csv/%s_lot%s",data[[i]]@metaData$fullName,lotno)
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

}



































