install.packages("devtools")
install.packages("roxygen2")
install.packages("usethis")

setwd('C:/Users/deepa/Dropbox (GaTech)/scripts/GitHub/FeatureCorr')

'''
Prepare a text file called DESCRIPTION with folwoing details
Package: CorrClust
Type: Package
Title: Finding associations through correlation of varied features
Version: 0.0.1.0
'''

#Setting up enviroment to debug, document and create a package
library(usethis)
library(devtools)
library(roxygen2)
library(usethis)
#devtools::create("FeatureCorr")
#use_gpl_license()
#use_git()


#saving usable data in the package funcitons
GTEXv7<-read.csv('C:/Users/deepa/Dropbox (GaTech)/scripts/GitHub/GTEX v7 Tissues - Copy.csv', header = TRUE, sep = ",", quote = "")
GTEXv7<-GTEXv7[-1]
GTEX<-GTEXv7[sample(rownames(GTEXv7), 8000), ]
save(GTEX, file = "GTEX_8000.RData")



TCGA<-read.csv('C:/Users/deepa/Dropbox (GaTech)/scripts/GitHub/TCGA 40genes.csv', header = TRUE, sep = ",", quote = "")
TCGA40<-t(TCGA[,3:ncol(TCGA)])
TCGA40<-TCGA40[,sample(seq(1:5000),2000)]
logTCGA40<-logTCGA40[,sample(seq(1:5000), 2000)]
save(GTEX, file = "GTEX_10000.RData")
save(GTEXv7, file = "GTEX.RData")
save(TCGA40, file = "TCGA40.RData")
save(logTCGA40, file = "logTCGA40.RData")
save(pairwisedf, file="TCGA.Rdata")
tools::checkRdaFiles('./data/GTEX.RData')
tools::resaveRdaFiles('./data/GTEX.RData', compress = "auto")
tools::checkRdaFiles('./data/logTCGA40.RData')
tools::resaveRdaFiles('./data/logTCGA40.RData', compress = "auto")
tools::checkRdaFiles('./data/TCGA40.RData')
tools::resaveRdaFiles('./data/TCGA40.RData', compress = "auto")
tools::checkRdaFiles('./data/GTEX_8000.RData')
tools::resaveRdaFiles('./data/GTEX_8000.RData', compress = "auto")

#Everytime you update your package, use the following command
devtools::document()
roxygenise()
devtools::load_all()
devtools::check()

#Create R markdown
use_readme_rmd(open = rlang::is_interactive()) #edit this file
https://github.com/SoftwareImpacts/SIMPAC-2021-53/blob/master/README.Rmd

#Update Readme
devtools::build_readme()


#Submit to CRAN
devtools::spell_check()
check_rhub(".", platforms = c("windows-x86_64-release","ubuntu-gcc-release","macos-highsierra-release")) #Check package on R-hub for CRAN submissions
devtools::release()

#Push this to check on github
usethis::use_github_action("check-standard")
devtools::release()

devtools::build_readme()

rhub::platforms()
check_rhub(".", platforms = c("windows-x86_64-release","ubuntu-gcc-release","macos-highsierra-release")) #Check package on R-hub for CRAN submissions

library(FeatureCorr)
