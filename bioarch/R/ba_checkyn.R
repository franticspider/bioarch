#' Get an arbitrary answer to a question
#' 
#' @param question the question for which you want an answer
#' @keywords question
#' @export
#' @examples
#' ba_checkyn("how are you feeling?")

ba_checkyn <- function(question){

	answer="n"
	question = sprintf("%s\n",question)

	while(answer != "y"){

		tag <- readline(question)

		check <- sprintf("Your answer is \"%s\", is that correct? (y/n)\n",tag)
		answer <- readline (check)
	}
	
	return (tag)

}
