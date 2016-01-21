#/* Copyright (C) 2014 Simon Hickinbotham	                          */
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





#load the library
require("MALDIquant")

#mydata<-NULL

#' DEPRECIATED print the usage of exportBruker
#' 
#' @keywords MALDIquant platemap Bruker
#' @examples
#' printUsage()
printUsage <- function(){
    message("=============================================================")
    message("This file is part of exportBruker version 1.0.2")
    message("Typical usage is:\n")

    message("  > exportDataSetTab(loadBruker(PATH))")
    message("\n..where PATH is the file path to the Bruker source directory\n")
    message("This program has created a subdirectory called \'csv\',\nwhere the exported text files will be placed")
    message("=============================================================")
}





###EXAMPLE USAGE###
#	mydata<-loadBruker('20131024_TedG1')
#	exportPlateMap(mydata, 'B2','plateB2.txt')


#TODO: want to be able to export a WHOLE bruker dataset to multiple csv files 
#TODO: This needs to be called from within the functions that need it!
#message('Creating directory \'csv\'')
#system('mkdir -p csv')

listData <- function(data){
	for(i in 1:length(data)){
		message(sprintf("entry %d has sample name %s and full name %s",i,data[[i]]@metaData$sampleName,data[[i]]@metaData$fullName))

	}
}

#' Export Bruker data to files, with lot Number as part of the filename
#' 
#' @param path the path to the Bruker data set
#' @keywords Bruker 
#' @export
#' @examples
#' exportDataSetTab("/home/anon/Bruker/20151105sample")
# exportDataSetTab USAGE:
#	this function will create files either with the file extension .txt, or using the fullname in the Bruker format.
#
#	to export to .txt files, do:
# 	> exportDataSetTab(data)
#
#	to export to files using the Bruker fullname, do:
# 	> exportDataSetTab(data,txt=FALSE)
exportDataSetTab <-function(data,txt=TRUE){

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



exportDataSetCsv <-function(data){

	for(i in 1:length(data)){
		outfn <- sprintf("csv/%s",data[[i]]@metaData$fullName)
		message("Exporting %s",data[[i]]@metaData$fullName)
		MALDIquantForeign::exportCsv(data[[i]], file = outfn, force=TRUE, col.names = FALSE)
	}
}





































