

ba_newplate <- function(tag){

	sh_new_name <- sprintf("%s%s",b384_root,tag)

	#get the master sheet

	master <- gs_title(b384_master)

	#copy it to the new sheet
	gs_copy(master,to=sh_new_name)	
	sheet <- gs_title(sh_new_name)

	return(sheet)

}


