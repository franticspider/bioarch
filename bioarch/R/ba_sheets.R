

# INSTALL NOTES
# devtools::install_github("ProjectMOSAIC/mosaic")
# Tried to use 'fetchGoogle.R' which I downloaded from the github mosaic repository listed above. Got the following message:
# fetchGoogle() is no longer actively supported.
# See the `googlesheets' package, which is dedicated to the R/Google interface.
#
# tried install.packages("googlesheets"), and got:
#
# Warning message:
# package ‘googlesheets’ is not available (for R version 3.0.2) 
#
# tried devtools::install_github("jennybc/googlesheets"). and got:
# ERROR: this R is version 3.0.2, package 'googlesheets' requires R >=  3.1.1
# Error: Command failed (1)
# In addition: Warning messages:
# 1: packages ‘tidyr’, ‘xml2’ are not available (for R version 3.0.2) 
# 2: In utils::install.packages(deps, dependencies = NA, Ncpus = threads) :
#  installation of package ‘dplyr’ had non-zero exit status
# 
# follwed this link to upgrade R to version 3.2: http://stackoverflow.com/questions/10476713/how-to-upgrade-r-in-ubuntu

# Then was able to proceed, as shown in the next section

#########################################
# HOW TO INSTALL:

# install devtools: 	install.packages("devtools"). Be patient - this takes a while.
# then  do: 		devtools::install_github("jennybc/googlesheets")

# some of the functions below use the '%>%' syntax (I know. don't worry). To follow this, you also need to install dplyr:
# then do: 		install.packages("dplyr")

#########################################

#Following the tutorial at http://htmlpreview.github.io/?https://raw.githubusercontent.com/jennybc/googlesheets/master/vignettes/basic-usage.html

require("googlesheets")
require("dplyr")

#Register as a user: (you only have to do this once)
#gs_auth(new_user = TRUE)

#########################################

#Following the tutorial at http://htmlpreview.github.io/?https://raw.githubusercontent.com/jennybc/googlesheets/master/vignettes/basic-usage.html

tutorial1 <-function(){

	# Copy the gapminder sheet to follow this tutorial:
	message("Copying gapminder into the Gapminder spreadsheet in your account")
	gs_gap() %>% gs_copy(to = "Gapminder")


	# list your sheets: 
	readline("hit <return>, and we'll list the sheets that are available to you")
	my_sheets <- gs_ls()
	my_sheets


	# register a sheet, so we can work on it (Beware making multiple copies - this function will create a new sheet with the same name):
	gap <- gs_title("Gapminder")

	#find all sheets matching a certain name:
	just_gap <- gs_ls("^Gapminder$")
	just_gap$sheet_key

	#then use the key to register a sheet (see how sheet_key[1] is passed into gs_key()
	ss2 <- just_gap$sheet_key[1] %>% gs_key()
	#alternative (same functionality)
	ss1 <- gs_key(just_gap$sheet_key[2])

}


########################################



########################################
# new_spots creates a new worksheet with a specified name, based on "Master"
new_spots <- function(sheet,ws_name,m_name="Master"){




}
##################################
# remove test data and others 
cleanup <- function(namefrag){



}

checkyn <- function(question){

	answer="n"

	while(answer != "y"){

		tag <- readline(question)

		check <- sprintf("Your answer is \"%s\", is that correct? (y/n)\n",tag)
		answer <- readline (check)
	}
	
	return (tag)

}

ynq <- function(question){

	q <- paste (question, "(y/n)\n", sep=" ")
	a <- "bibble"
	found = F

	while((a != "y" ) &&  (a != "n")){
		a <- readline(q)
		
		if(a == "y" ) found = T
		if(a == "n" ) found = T
		if(!found){
			message("try again, answer y or n\n\t")
		}
		else{
			message("ok then!")
		}
	}
	return (a)

}

platemap_new <- function(tag){

	baname = "BIOARCH-PLATEMAP-DATA-DO-NOT-DELETE"
	b384_master = "BIOARCH_BRUKER384_PLATEMAP_MASTER"
	sh_new_name <- sprintf("%s%s",baname,tag)

	#get the master sheet

	master <- gs_title(b384_master)	

	#copy it to the new sheet
	gs_copy(master,to=sh_new_name)	
	sheet <- gs_title(sh_new_name)

	return(sheet)

}

rownumfromchar <- function(c){
 return(strtoi(charToRaw(c),16L)-63)
}



log_bioarch <- function(){
	answer <- ynq("Do you want to create a new platemap?")

	mygs <- NA

	if(answer == "y"){
		tag<-checkyn("What is the tag of the new platemap file?\n")
		#message("checking to see if this exists...\n")
		#sterm <- sprintf("^%s$%s$",baname,tag)
		
		message("ok, we'll create that sheet for you now...\n")

		sheet <- platemap_new(tag)
	
		
		
		
	}
	else{
		tag<-checkyn("ok, please enter part or all of the platemap name and we'll find it for you")
		sterm <- sprintf("^baname$tag$")
		spshts <- gs_ls(sterm)
		#TODO: Access the correct sheet
	}


	#Now we have the sheet, let's see if there's any data to add.
	answer <- ynq("Do you have data for this platemap?")
	if(answer=="y"){

		lotno<-checkyn("please enter the lot number")
		author<-checkyn("please enter the author")
		row<-checkyn("please enter the row (A-P)")
		col<-checkyn("please enter the column (1-24)")




 		as.integer(charToRaw(c))


		colnum<- strtoi(col)+1

		message("Adding data to row %d and column %d of the platemap sheet")

	}


}


# get a list of all sheets with the phrase "Spot Positions" in the title
#spots <- gs_ls("^Spot Positions$")

# now *register* the first sheet in that list:
#spotssheets <- gs_title(spots$sheet_title[1])


testname = "THIS_IS_A_TEST_PLEASE_DELETE"








