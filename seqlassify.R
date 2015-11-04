library("q2e")
library("bioarch")


#Specify the `margins' of the sample 
moff = 1


#load the cow sample
cow <- isodists("../q2e/q2e/testdata/cowPeptides")
lbl <- min(cow$mass) - moff
ubl <- max(cow$mass) + moff

#load the test data
if(!exists("testdata"))
testdata <- loadBrukerXML("/home/sjh/Desktop/sjh/bioarch/SF/SFData/20140307_SF_UPenn164-181/")

#below we pick a spot with good data:
tdi <- 24

#not sure which of these we'll need:
t1<-testdata[[tdi]]
t1l<-list(t1)

#get the subset of the data
testi <-  t1l[[1]]@intensity[lbl <= t1l[[1]]@mass & t1l[[1]]@mass < ubl ]                   
testm <-  t1l[[1]]@mass[lbl <= t1l[[1]]@mass & t1l[[1]]@mass < ubl ]

#Plot the raw data 
plot(x=testm,y=testi)

#create an interpolation
xout = seq(from = lbl, to = ubl, by = 0.02)

#message("bang!")
#testint <- approx(x=xi,y=yi,xout=xout, method="linear", rule = 2)
yri <- approx(x=testm,y=testi,xout=xout, method="linear", rule = 2)

#plot the interpolated data
lines(x=yri[[1]],y=yri[[2]])


