

#TODO: Need to have sheet loaded (in human_spectra.R)


library("googlesheets")
library("bioarch")
library("q2e")

source("~/git/bioarch/dev/ba_msutils.R")



#copied from theoretical_spectra.R... 
#Fetch the sheet if it isn't available...
if(!exists("sheet")){

	message("Fetching the sheet of mammalian collagen sequences")
	mcs <- gs_title("Mammalian Collagen Sequences v0.0.1")
	print(mcs)

	sheet <- gs_read(mcs)
}




#copied from theoretical_spectra.R... 
#Fetch the sheet if it isn't available...
if(!exists("sheet")){

	message("Fetching the sheet of mammalian collagen sequences")
	mcs <- gs_title("Mammalian Collagen Sequences v0.0.1")
	print(mcs)

	sheet <- gs_read(mcs)
}



#readline("Press <return> to do interactive analysis of sample C1")

pdf(file = "N1byhand_v2.pdf",w=8,h=8)
ms_fit(sheet,"human","~/tmp/bioarch_keri/20160803_Keri12/20160803_Keri12_0_","N1",c("D4","D7","D10"),doyn=F)
dev.off()


#ms_fit(sheet,"human","~/tmp/bioarch_keri/20160803_Keri12/20160803_Keri12_0_","N2",c("D5","D8","D11"))
