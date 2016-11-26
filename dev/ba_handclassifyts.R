
#	REVERSE PLATEMAP FOR KERI'S DATA
#   |-----ENGLAND-----|------ITALY-------------|------------PORTUGAL--------------|	
#	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20
#-------------------------------FOLDER ONE: 20160803------------------------------|
#			
#B	A1	A2	A3  B1	B3	C1	C2	C3	A10	A11	A12	B10	B12	C10 C11	C12	A19	A20	A21	B19		
#	A4	A5	A6	B4	B6	C4	C5	C6	A13	A14	A15	B13	B15	C13	C14	C15	A22	A23	A24	B22	
#	A7	A8	A9	B7	B9	C7	C8	C9	A16	A17	A18	B16	B18	C16	C17	C18	D1	D2	D3	E20	
#																						
#N	D4	D5	D6	E4	E5	F4	F5	F6	D13 D14 D15	E13	E15	F13	F14	F15	D22	D23 D24	E22					
#	D7	D8	D9	E7	E9	F7	F8	F9	D16	D17 D18	E16	E18	F16	F17	F18	G1	G2	G3	H1			
#	D10	D11	D12	E10	E12	F10	F11	F12	D19	D20	D21	E19	E21	F19	F20	F21	G4	G5	G6	H2		
#																						
#-------------------------------FOLDER TWO: 20160909------------------------------|
#	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20
#	
#K 	A1	A2	A3	B1	B3	C1	C2	C3	A10	A11	A12	B10	B12	C10	C11	C12	A19	A20	A21 B19
#	A4	A5	A6	B4	B6	C4	C5	C6	A13	A14	A15	B13	B15	C13	C14	C15	A22	A23	A24	B22
#	A7	A8	A9	B7	B9	C7	C8	C9	A16	A17	A18	B16	B18	C16	C17	C18	D1	D2	D3	E1	
#
#I	D4	D5	D6	E4	E6	F4	F5	F6	D13	D14	D15	E13	E15	F13	F14	F15	D22	D23	D24	E22	
#	D7	D8	D9	E7	E9	F7	F8	F9	D16	D17	D18	E16	E18	F16	F17	F18 G1	G2	G3	H1		
#	D10	D11	D12	E10	E12	F10	F11	F12	D19	D20	D21	E19	E21	F19	F20	F21	G4	G5	G6	H4	
#
#C	G7	G8	G9	H7	H9	I7	I8	I9	G16	G17	G18	H16	H18	I16	I17	I18	J1	J2	J3	K1
#	G10	G11	G12	H10	H12	I10	I11	I12	G19	G20	G21	H19	H21	I19 I20 I21	J4	J5	J6	K4				
#	G13	G14	G15	H13	H15	I13	I14	I15	G22	G23	G24	H22	H24	I22	I23	I24	J7	J8	J9	K7				


#TODO: Need to have sheet loaded (in human_spectra.R)



require("stringr")
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



#readline("Press <return> to do interactive analysis of sample C1")

###################################################################################
#pdf(file = "C1byhand_v2.pdf",w=8,h=8)
#ms_fit(sheet,"human","~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_","C1",c("G7","G10","G13"),doyn=F)
#dev.off()

#froot = "~/tmp/bioarch_keri/20160803_Keri12/20160803_Keri12_0_"
#spots = c("D4","D7","D10")
#sample = "N1"
#generate_alignment_pdf(froot, sample, spots)
#pdf(file = "N1byhand_v2.pdf",w=8,h=8)
#ms_fit(sheet,"human",froot,sample,spots,doyn=F)
#dev.off()

###################################################################################
#froot = "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_"
#spots = c("G7","G10","G13")
#sample = "C1"
#generate_alignment_pdf(froot, sample, spots, species = "sheep")
#pdf(file = "C1byhand_v2_sheep.pdf",w=8,h=8)
#ms_fit(sheet,"sheep",froot,sample,spots,doyn=F)
#dev.off()


#froot = "~/tmp/bioarch_keri/20160803_Keri12/20160803_Keri12_0_"
#spots = c("D4","D7","D10")
#sample = "N1"
#generate_alignment_pdf(froot, sample, spots, species = "sheep")
#pdf(file = "N1byhand_v2_sheep.pdf",w=8,h=8)
#ms_fit(sheet,"sheep",froot,sample,spots,doyn=F)
#dev.off()

###################################################################################
#froot = "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_"
#spots = c("I7","I10","I13")
#sample = "C6"
#generate_alignment_pdf(froot, sample, spots)
#pdf(file = "C6byhand_v2.pdf",w=8,h=8)
#ms_fit(sheet,"human",froot,sample,spots,doyn=F)
#dev.off()


froot = "~/tmp/bioarch_keri/20160803_Keri12/20160803_Keri12_0_"
spots = c("F4","F7","F10")
sample = "N1"
generate_alignment_pdf(froot, sample, spots, species = "sheep")
pdf(file = "N6byhand_v2.pdf",w=8,h=8)
ms_fit(sheet,"human",froot,sample,spots,doyn=F)
dev.off()









#ms_fit(sheet,"human","~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_","C2",c("G8","G11","G14"))
