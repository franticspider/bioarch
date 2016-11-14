#TDD on automated analysis of MALDI files, using Keri Rowsell's data


source("~/git/bioarch/dev/ba_msutils.R")






test_human <- function(directory,sampleindexes){

	#Check the inputs
	if(!generatefilenames(directory,sampleindexes)){
		message("ERROR: sample filenames are incorrect")
		return
		}
	
	
		
	if(!validmaldidata(filenames)){
		message("ERROR: sample filenames are in the wrong format")
		return
		}

}







############################################################

test_human("~/tmp/bioarch_keri/20160909_Keri13_0_",c("G7","G10","G13"))
