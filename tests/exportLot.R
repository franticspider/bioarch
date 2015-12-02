

require(bioarch)

message("Testing the exportLot function")

bpath = "/home/sjh/Desktop/sjh/bioarch/SF/SFData/20140307_SF_UPenn164-181"

bpath = "/home/sjh/Desktop/sjh/bioarch/20151106_AshleySophySarah"



bioarch_GetSpotList()


bioarch_ExportLotData(bpath,123456)


#exportLot(path,lotno=0,txt=TRUE)
