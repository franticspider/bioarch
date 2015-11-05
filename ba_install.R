

#The idea is that we run this whenever we update the q2e code - not sure if it'll work!
#http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/

#Tried doing this, but it doesn't ececute the rest of the file
#Assuming you have this in your .Rprofile
#makeActiveBinding("refresh", function() { system("R"); q("no") }, .GlobalEnv)
#refresh



#These are the libraries you need to build a package:
require("devtools")
require(roxygen2)

#Try to build the help files:
setwd("bioarch")


document()
roxygenize(".")

setwd("..")


install("bioarch")





############################
# Want to Upload to github?
# DO THIS OUTSIDE R!
# git add .
# git commit
# git push
############################

#########################################
# HOW TO INSTALL:
#
# install devtools: 	install.packages("devtools"). Be patient - this takes a while.
# then  do: 		devtools::install_github("jennybc/googlesheets")
#
# some of the functions below use the '%>%' syntax (I know. don't worry). To follow this, you also need to install dplyr:
# then do: 		install.packages("dplyr")
#
#########################################

############################
# want to build a windows binary? 
#
# instructions from http://win-builder.r-project.org/
# Switch to the directory above the package root (the same directory this file resides in)
#
# 
#
# First, build the tarball that'll be checked
# $ R CMD build q2e
#
# Then check *THE TARBALL* (Not the directory - that'll give errors)
# $ R CMD check q2e_VERSION.tgz
# (Where VERSION is the version number used in the DESCRIPTION FILE
# There may be errors, which need to be fixed before the next step
# *MAKE SURE YOU REBUILD THE HELP BY RUNNING THIS SCRIPT WHEN YOU FIX DOCUMENTATION ERRORS!*
#
############################


#Then to download:
#install_github('franticspider/q2e')


#TODO: Somehow check the outputs with the outputs of the original C code...


#now run the error tests (see 'q2e_tests.R')

#message("Testing q2eloadfiles:\n\t testing non-existent file:")
#try(q2eloadfiles("bobbins"))

