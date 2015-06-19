#' Get an arbitrary answer to a yes/no question
#' 
#' @param question the question for which you want a yes/no answer
#' @keywords question
#' @export
#' @examples
#' ba_ynq("are you well?")

ba_ynq <- function(question){

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
			#message("ok then!")
		}
	}
	return (a)

}

