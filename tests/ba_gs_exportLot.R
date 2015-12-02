

require(bioarch)
require(googlesheets)

#String constants for the googlesheets columns
SPOTCOL = 'Spot.no.'
LOTCOL = 'Lot.no.'

message("Testing the gsExportLot function")

bpath = "/home/sjh/Desktop/sjh/bioarch/20151106_AshleySophySarah"

message(sprintf("Bruker data is at:\n\t%s",bpath))

tag<-ba_checkyn("ok, please enter part or all of the platemap name and we'll find it for you")
sterm <- sprintf("*%s*",tag)
message("searching....")
spshts <- gs_ls(sterm,ignore.case=T)
message(sprintf("Found %d sheets matching that query:",length(spshts$sheet_title)))
message("\tINDEX\tTITLE")
for(i in 1:length(spshts$sheet_title)){
	message(sprintf("\t%d\t%s",i,spshts$sheet_title[i]))
}

idx <- strtoi(ba_checkyn("Please enter the index of the platemap"))

#ok, now we can fetch the sheet...

message(sprintf("OK, fetching data from the sheet titled \"%s\"",spshts$sheet_title[idx]))

#Get a list of lot numbers:
sheet <- gs_title(spshts$sheet_title[idx])
lot <- gs_read(sheet,"List of samples")
lnos <-unique(lot[,LOTCOL])
lnos <- lnos[!is.na(lnos),1]
message("Found the following lot numbers:")
for(i in 1:nrow(lnos)){
	message(sprintf("\t%d\t%s",i,lnos[i,1]))
}
	
lidx <- strtoi(ba_checkyn("Please enter the index of the lot number"))


message(sprintf("Lot number %s has spots at the following positions:",lnos[lidx,1]))

found1 = F
for(i in 1:nrow(lot)){
	if(identical(lot[i,LOTCOL],lnos[lidx,1])){
		message(sprintf("\t%d\t%s",i,lot[i,SPOTCOL]))
		if(!found1){
			found1=T
			spotidx=i
		}
		else{
			spotidx=c(spotidx,i)
		}
	}
}

ans <- ba_ynq("Do you want to export these spots?")

if(identical(ans,"y")){
	mydata <- ba_loadBruker(bpath)
	message("Exporting these spots now")
	for(i in 1:length(spotidx)){
		message(sprintf("Exporting spot %d: position %s...",i,lot[spotidx[i],SPOTCOL]))
		j <- ba_indexFromPlatePos(mydata,lot[spotidx[i],SPOTCOL])
		#message(sprintf("entry %d has sample name %s and full name %s",j,mydata[[j]]@metaData$sampleName,mydata[[j]]@metaData$fullName))
		
		#Make the filename
		lotname = gsub(" ", "_", lnos[lidx,1])
		
		fn <- sprintf("%s_%s.txt",lotname,mydata[[j]]@metaData$fullName)
		message(sprintf("\t ...exporting to file named '%s'",fn))
		MALDIquantForeign::exportCsv(mydata[[j]], file = fn, force=TRUE, col.names = FALSE)
	}
}



#bioarch_GetSpotList()


#bioarch_ExportLotSpots(bpath,123456,spots)


#exportLot(path,lotno=0,txt=TRUE)
