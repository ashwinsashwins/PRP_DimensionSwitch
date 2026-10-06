# DimensionSwitchPRP(6 session, Study2) Cross Condition Behavior merging
# To do 
# 1. think of more analyses for switching of features!
# 2. Re-code information in meaningful manner
# 3. Exclude flipped order responses...(how to detect this??)!
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))
gc()

# Load libraries
library(data.table)
library(tidyverse)
library(broom)
library(effects)
library(afex)
library(psych)
library(stringr)
library(ggplot2)
library(here)


#Setting up directory
fsep<-.Platform$file.sep;
userName <- "Ashwin Srinivasan"; # Atsushi Kikumoto/Ashwin Srinivasan
if (userName == "Atsushi Kikumoto") {Dir_R<-path.expand("~/Dropbox/w_ONGOINGRFILES/w_OTHERS")
}else if(userName == "Ashwin Srinivasan"){Dir_R<-(path.expand("/users/asrini17/data/asrini17/R_Helper_Scripts"))
}else{print("uh oh! check for typos")}
Dir_HERE<-dirname(rstudioapi::getSourceEditorContext()$path)
Dir_BDATA<-file.path(dirname(Dir_HERE),"Matlab","BEH","z_ALLGRAND(EEG,6sessions)")
Dir_EDATA<-file.path(dirname(Dir_HERE),"Matlab","EEG")


# Source Files
setwd(Dir_R)
source('basic_theme.R')
source('basic_lib.R')
# source('anovakun_480.txt')
theme_set(theme_bwL(base_size = 10))#32/28

#======================================================================================================
# Load data
#======================================================================================================
# READING FILES
setwd(Dir_BDATA)#Change directory
ds_b<-fread("DimensionSwitchPRP_BehPP.txt") 
ds_b[is.nanM(ds_b)]<-NA
setDT(ds_b)


#keep cutRT1 and 2 from original template (1 value per subject per session) to plug back in later
#RT_cuts <- unique(ds_b[,c("SUBID", "SUBID_S", "cutRT1", "cutRT2")]) #"BLOCK", "TRIAL", 

idVs <- names(ds_b)[!grepl("_2|_1", names(ds_b))];print(idVs)

respCols <- names(ds_b)[grepl("_2|_1", names(ds_b))];print(respCols) #grabs all columns ending in _2 or _1
mesVs <- unique( sub("_1|_2", "", respCols));print(mesVs)

mes_list <- lapply(mesVs,
                   function(x) paste0(x, "_", 1:2)) # list of _1/_2 columns for each measure
names(mes_list) <- mesVs
mes_list[["TASKSET_CORSIDEc"]] <- c("TASKSET_CORSIDE_1c","TASKSET_CORSIDE_2c")

#form <- paste0("CLASS~",paste0(idVs, collapse="+"))

ds_long <- melt(ds_b,
             measure.vars = mes_list,
             variable.name = "STIMNUM", 
             value.name     = names(mes_list),
             variable.factor = FALSE)
ds_long[,RT:=as.numeric(RT)]
ds_long[,TASKTS:=as.numeric(TASKTS)]

ds_long[, cutRT := ifelse(STIMNUM == 1, cutRT1, cutRT2)]

#TASKSET_CORSIDE variable (for RSA_)
#ds_long[, TASKSET_CORSIDEc := paste(ds$TASKSET, ds$CORSIDE, sep = "_")]
#ds_long[, TASKSET_CORSIDE := as.numeric(as.factor(ds$TASKSET_CORSIDEc))]# create numeric equivalents for the 8 different combos

ds_long[,cutRT1:=NULL] 
ds_long[,cutRT2:=NULL]

#write it out!
write.table(ds_long,paste0(Dir_BDATA,'/DimensionSwitchPRP_BehPP_LONG.txt'), sep="\t",row.names=FALSE)
