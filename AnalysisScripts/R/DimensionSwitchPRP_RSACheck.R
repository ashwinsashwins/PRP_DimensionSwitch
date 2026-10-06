# DimensionSwitchPRP
# Small checking script for RSA db
# NOTE:
# =============================================================

graphics.off()
rm(list = ls(all = TRUE))

#Setting up directory
#Setting up directory
fsep<-.Platform$file.sep;
userName <- "Ashwin Srinivasan"; # Atsushi Kikumoto/Ashwin Srinivasan
if (userName == "Atsushi Kikumoto") {Dir_R<-path.expand("~/Dropbox/w_ONGOINGRFILES/w_OTHERS")
}else if(userName == "Ashwin Srinivasan"){Dir_R<-(path.expand("/users/asrini17/data/asrini17/R_Helper_Scripts"))
}else{print("uh oh! check for typos")}
Dir_HERE<-dirname(rstudioapi::getSourceEditorContext()$path)
Dir_BDATA<-file.path(dirname(Dir_HERE),"Matlab","BEH","z_ALLGRAND(EEG,6sessions)")
Dir_EDATA<-file.path(dirname(Dir_HERE),"Matlab","EEG")
Dir_GRAND<-paste0(Dir_EDATA,"/w_ALLGRAND")
Dir_RSAMODEL<-paste0(Dir_EDATA,"/w_RSAMODELS")

# Load libraries
library(data.table)
library(dplyr)
library(RSQLite)


setwd(Dir_GRAND)
fileN<-"DimensionSwitchPRP"
dbName<-paste0("FIT_RSA_TEST_CTRR")
dbCon=dbConnect(dbDriver("SQLite"), paste0(fileN,'_',dbName))# Open database connection
dscheck=as.data.table(dbGetQuery(dbCon,paste0('SELECT * FROM ',dbName, ' LIMIT 15')))


unique_subids <- dbGetQuery(dbCon, paste0("
  SELECT DISTINCT SUBID_S
  FROM ", dbName, "
"))

unique_subids