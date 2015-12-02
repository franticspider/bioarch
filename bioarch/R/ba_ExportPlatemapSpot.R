

ba_ExportPlateMapSpot <- function(mydata, cell, outfn){
	ixz1 <- indexFromPlatePos(mydata,cell)
	MALDIquantForeign::exportCsv(mydata[[ixz1]], file = outfn, force=TRUE, col.names = FALSE)

}

