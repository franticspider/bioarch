#/* Copyright (C) 2014 Simon Hickinbotham, Matthew Collins, Keri Roswell  */
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


require("MALDIquantForeign")
#source("exportBruker.R")
#source("flamePlot_3.R")


#' Create a structure that holds the q2e parameters
#' 
#' @param spot the spot address (e.g. "A12")
#' @keywords platemap
#' @export
#' @examples
#' ba_logspot("A12")
getclass <- function(idx,fn){
	name = "unknown"
	t <- read.table(fn,sep='\t')
	#'b' %in% v
	if(idx %in% t[,1]){
		name<-t[grep(idx,t[,1]),3]
		message(sprintf("found a %s entry at %s",name,idx))
	}
	return(name)

}
