#' DEPRECIATED: Log a new spot position
#' 
#' @param spot the spot address (e.g. "A12")
#' @keywords platemap
#' @examples
#' ba_logspot("A12")
ba_logspot <- function(spot){

# 	https://docs.google.com/a/york.ac.uk/spreadsheet/pub?key=0Ag83-Eguk23ldDYzSzRocWZpMkFtQWxHNW9ySGhzMWc&output=html


# My test spreadsheet is here: 
#https://docs.google.com/spreadsheets/d/1N8gnNY4TO8qoYWoM-25BlOjXpXB21Ol5nnpCZYWqZG8/edit#gid=1474463154


	answer <- ba_ynq("Do you want to create a new platemap?")

	mygs <- NA

	if(answer == "y"){
		tag<-ba_checkyn("What is the tag of the new platemap file?\n")
		#message("checking to see if this exists...\n")
		#sterm <- sprintf("^%s$%s$",baname,tag)
		
		message("ok, we'll create that sheet for you now...\n")

		sheet <- ba_newplate(tag)
	
	}
	else{
		tag<-ba_checkyn("ok, please enter part or all of the platemap name and we'll find it for you")
		sterm <- sprintf("^%s$%s$",b384_master,tag)
		spshts <- gs_ls(sterm)

		message("we were looking for:")
		message(sterm)

		#TODO: Access the correct sheet
	}


	#Now we have the sheet, let's see if there's any data to add.
	answer <- ba_ynq("Do you have data for this platemap?")
	if(answer=="y"){

		row<-ba_checkyn("please enter the row (A-P)")
		col<-ba_checkyn("please enter the column (1-24)")
		message("<TODO: error checking on these entries>")
		date<- ba_checkyn("please enter the date in format YYYYMMDD")
		message("<TODO: check that the date is the same for all spots>")
		lotno<-ba_checkyn("please enter the lot number")
		mnemonic<-ba_checkyn("please provide a keyphrase for this spot (less than 12 characters)")
		author<-ba_checkyn("please enter the researcher name")

		rownum<- as.integer(charToRaw(row))
		colnum<- strtoi(col)+1

		maprow <- rownum-63
		mapcol <- colnum

		logrow <- ((rownum-65)*24)+colnum

		

		cell<-(sprintf("%s%d",intToUtf8(mapcol+64),maprow))
		message(sprintf("Adding data to row %d and column %d of the platemap sheet (cell %s)",maprow,mapcol,cell))

		sheet <- sheet %>% gs_edit_cells(ws="MAP",anchor=cell,input=author)

		message(sprintf("Adding data to row %d and column %d of the platemap sheet (cell %s)",maprow,mapcol,cell))


		data <- c(date,lotno,mnemonic,author)
		logdata <- t(as.data.frame(data))
		cell<-(sprintf("C%d",logrow))
		message(sprintf("Adding data to row C and column %d of the platemap sheet (cell %s)",logrow,cell))

		sheet <- sheet %>% gs_edit_cells(ws="LOG",anchor=cell,input=logdata)
		

	}

}
