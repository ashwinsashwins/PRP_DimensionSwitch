# DimensionSwitchPRP_AggRegress
# load pred file (contains single-trial accuracy & confidence)
# then produce aggregated data (by TRIAL!!) for specified conditions (with few attributes)
# SOA Realignment
# 1. add timeLK_SOA column, using correct timeL 
# 2. Use following settings
# timeL<-seq(-400,350,4);
# timeS<-seq(-400,350,20);
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))
gc()

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
Dir_REPO <- path.expand("/oscar/home/asrini17/data/asrini17/DimensionSwitch_PRP/PRP_DimensionSwitch/figures")

# Load libraries
library(data.table)
#library(tidyverse)
library(tidyr)
library(dplyr)
library(broom)
library(rhdf5)
library(foreach)
library(zoo)
library(doMC)
library(binhf)
library(RColorBrewer)
library(RSQLite)
library(lemon)
library(DBI)
library(lme4)
library(sjPlot)
library(cluster)
#library(dbplyr)

# Source Files
setwd(Dir_R)
source('basic_lib.R')
source('basic_theme.R')

#Define basic info
saveOn<-T;# save results or not
fileN<-"test"
#timeL<-seq(-200,700,4); # labels of time for time indexes
#timeS<-seq(-200,700,12); # labels of time you want

#if sessN == "nb_S2":
timeL<-seq(-500,700,4); # labels of time for time indexes
timeS<-seq(-500,700,12);

refRLV <- NA;#NA, SOA, Realignment: re-epoch data referencing refRLV timing
srate<-4
if (!is.na(refRLV)){adj<-which(timeL==0)}else{adj<-0}
cutF <- function(t){as.numeric(cut(t,breaks=timeS,include.lowest=T))}
idxF<-function(i,adj) which.min(abs(timeL-i))-adj # normal

#load in behavioral data in order to include task relevance - test
setwd(Dir_BDATA)
ds_beh<-fread("DimensionSwitchPRP_BehPP.txt"); #loading in in order to get TASKSET_2/TASKSET_1
ds_beh[is.nanM(ds_beh)]<-NA


# Database= PP: REG / RSAFits: FIT_RSA_TSCONJ,FIT_RSA_L1TSCONJ,FIT_RSA_TSCONJ_RL,FIT_RSA_TCTSCONJ
setwd(Dir_GRAND)
expN<-"DimensionSwitchPRP";
modelV <- "TASKSET_CORSIDE_2"
sessN <- "nb_S2"
regs <- "BIS" #include RT/ACC nuisance regressors? leave blank if no, or for balanced integration, "BIS"
tableN<-"FIT_RSA_TEST_CTRR";# 
rsaValueN <- "betas"

# Dependent variables
# depV<-c("ACCRESP_acc","RESP_ED_acc");# REG
#depV<-sprintf(c("LH_RSA%s","OE_RSA%s","BF_RSA%s","RB_RSA%s"),"");#
depV<-sprintf(c("TASKSET_RSA%s","CORSIDE_RSA%s", "CONJ_RSA%s"),"")

# Grouping variables (keep first two columns fixed!)
grpV<-c("SUBID","SESSION", "BLOCK","TRIAL" ,"NUMTASK", "SOA", "time") #, "SOA"
taskV <- c("TASKSET_1","TASKSET_2") #"TASKSET_1", 

# Open connections
dbCon=dbConnect(dbDriver("SQLite"),paste0(expN,'_',modelV ,'_',regs,'_',tableN))# Open database connection
dscheck=as.data.table(dbGetQuery(dbCon,paste0('SELECT * FROM ',tableN,'  LIMIT 10000')));
dbCL<-dbListFields(dbCon,dbListTables(dbCon))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Realignment: add timeLK_[event] variable, idx in same unit as time
# timeL will be arbitrary length, and 0 = onset of event
# it willadd timeLK_[ev] column in the dataset (when it does not exist)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# if (!is.na(refRLV)){
#   # Define RL eqn and update few variables for later process
#   timeRLV<-paste0("timeLK_",refRLV)
#   grpV<-gsub("time",timeRLV,grpV);
# 
#     # update original database (only when it does not exist before)
#   if (!any(dbCL %in% timeRLV)){
#     zeroIDX <- max(which(timeL<=0))
#     tleqn<-sprintf("-1 *(floor(%s/%d) + %d - time)",refRLV,srate,zeroIDX);
#     #eqn<-function(RT,time){-1 *(round(RT/srate) - time + adjIDX)}
#     print(paste0('Realignment is requested:',timeRLV," = ",tleqn))
#     dbExecute(dbCon, sprintf("ALTER TABLE %s ADD COLUMN %s NUMBER",tableN,timeRLV))
#     dbExecute(dbCon,sprintf("UPDATE %s SET %s = %s",tableN,timeRLV,tleqn))
#     }
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AGGREGATION: grouping SUBID x gvar in database with flexible query sentences:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Aggregation with SQ language
#tIDX<-sapply(range(timeS),idxF,adj);whereQ="";print(tIDX)
if (sessN == "nb_S2"){tIDX<- c(1,301)} else if (sessN == "nb_S1"){tIDX <- c(1,226)}; print(tIDX) #cut to time windows of interest
selectQ<-paste0('SELECT ',paste(unique(c(grpV)),collapse = ", "))
compQ<-paste(sprintf(", AVG(%s) AS %s",depV,depV),collapse = "")
if (!any(is.na(timeL))){whereQ <- sprintf(" WHERE %s >= %d AND %s <= %d",grpV[7],tIDX[1],grpV[7],tIDX[2])}
groupQ<-paste0(' GROUP BY ',paste(grpV,collapse = " , "))
querySTW<-paste0(selectQ,compQ,sprintf(" FROM %s",tableN),whereQ,groupQ) %>% print()
tic();dsA=as.data.table(dbGetQuery(dbCon,querySTW));toc()
#Aggregation with dbplyr package
#empty_name <- tbl(dbCon, tableN) #creates local reference for db
#querySTWb <- empty_name %>%
#  filter(time >=1& time<=250)%>%
#  group_by(SUBID, time, SESSION, NUMTASK) %>%
#  summarise(LH_RSA = mean(LH_RSA),
#            OE_RSA = mean(OE_RSA),
#            BF_RSA = mean(BF_RSA),
#            RB_RSA = mean(RB_RSA))
#querySTWb %>% show_query() #check query
#tic();dsB<- as.data.table(dbGetQuery(dbCon,querySTWb))();toc()


# Keep consistant time labels (this should be fixed in decoding step...)!
timeC<-colnames(dsA)[grepl("time",colnames(dsA))]
if (!any(is.na(timeL)) & diff(unique(dsA[,get(timeC)]))[1]!=srate){
  #tIDX<-sapply(range(timeS),idxF,0)# no adjustment needed here
  timeLL<- timeL[tIDX[1]:tIDX[2]] %>% print();
  if (is.na(refRLV)){dsA[,(timeC):=lapply(.SD,function(t)timeL[t]),.SDcols=timeC]}
  if (!is.na(refRLV)){dsA[,(timeC):=timeLL[as.numeric(as.factor(dsA[,get(timeC)]))]]}
}

#Ashwin addition -- grouping RSA by relevant features
# #tIDX<-sapply(range(timeS),idxF,adj);whereQ="";print(tIDX)
# if (modelV == "S2ID"){tIDX<- c(125,225)} else if (modelV == "S1ID"){tIDX <- c(50,150)}; print(tIDX)
# selectQT<-paste0('SELECT ',paste(unique(c(grpV, taskV)), collapse = ", "))
# compQT<-paste(sprintf(", AVG(%s) AS %s",depV,depV),collapse = "")
# if (!any(is.na(timeL))){whereQ <- sprintf(" WHERE %s >= %d AND %s <= %d",grpV[6],tIDX[1],grpV[6],tIDX[2])}
# groupQT<-paste0(' GROUP BY ',paste(c(grpV[1:6],taskV),collapse = ", "))
# querySTWb<-paste0(selectQT,compQT,sprintf(" FROM %s",tableN),whereQ,groupQT) %>% print()
# tic();dsB=as.data.table(dbGetQuery(dbCon,querySTWb));toc();
# 
# # Now for dsB
# timeC<-colnames(dsB)[grepl("time",colnames(dsB))]
# if (!any(is.na(timeL)) & diff(unique(dsB[,get(timeC)]))[1]!=srate){
#   tIDX<-sapply(range(timeS),idxF,0)# no adjustment needed here
#   timeLL<- timeL[tIDX[1]:tIDX[2]] %>% print();
#   if (is.na(refRLV)){dsB[,(timeC):=lapply(.SD,function(t)timeL[t]),.SDcols=timeC]}
#   if (!is.na(refRLV)){dsB[,(timeC):=timeLL[as.numeric(as.factor(dsB[,get(timeC)]))]]}
# }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Segmentation: averaging data over time window (aggregation with time intervals)
# Need to be adjusted fro DYN, CROSS, JMB, etc...(those already labeled vs. not)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (diff(timeS[c(1,2)])>srate & !any(is.na(timeL))){
  # Define averaging window then aggregate
  dsA[,(timeC):=(lapply(.SD,cutF)),.SDcols=timeC]
  dsA<-dsA[,lapply(.SD, mean,na.rm=T),by=grpV]# time labels will be mid-point
  dsA[,(timeC):=lapply(.SD,function(t)rollmean(timeS,2)[t]),.SDcols=timeC]
}

#now for dsB
# if (diff(timeS[c(1,2)])>srate & !any(is.na(timeL)) & !any(is.na(taskV))){
#   # Define averaging window then aggregate
#   dsB[,(timeC):=(lapply(.SD,cutF)),.SDcols=timeC]
#   grpV2<- paste(c(grpV,taskV))
#   dsB<-dsB[,lapply(.SD, mean,na.rm=T),by=grpV2]# time labels will be mid-point
#   dsB[,(timeC):=lapply(.SD,function(t)rollmean(timeS,2)[t]),.SDcols=timeC]
# }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Averaging over time to get trial-wise estimates (late & early)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
timeValues <- unique(dsA$time); print(timeValues)
if(sessN == "nb_S2"){
  timeEarly <- timeValues[which(timeValues %in% 0:350)]; timeLate <-timeValues[which(timeValues %in% 351:700)]} else if (sessN == "nb_S1"){
    timeEarly <- timeValues[which(timeValues %in% 0:350)]; timeLate <- timeValues[which(timeValues %in% 351:700)]}

dsEarly <- dsA %>% #there is definitely a better way to this that doesn't require creating 2 whole new dataframes...
  filter(time %in% timeEarly)%>%
  group_by(SUBID, SESSION, BLOCK, TRIAL, NUMTASK, SOA) %>%
  summarise(earlyTASKSET_RSA = mean(TASKSET_RSA),
            earlyCORSIDE_RSA = mean(CORSIDE_RSA),
            earlyCONJ_RSA = mean(CONJ_RSA))
dsLate<- dsA %>% 
  filter(time %in% timeLate)%>%
  group_by(SUBID, SESSION, BLOCK, TRIAL, NUMTASK, SOA) %>%
  summarise(lateTASKSET_RSA = mean(TASKSET_RSA),
            lateCORSIDE_RSA = mean(CORSIDE_RSA),
            lateCONJ_RSA = mean(CONJ_RSA))
dsA <- dsEarly %>%
  left_join(dsLate, by = c(grpV[1:6]))
setDT(dsA)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summary
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save session information as attributes
setwd(Dir_GRAND)
#attributes(dsA)$timeL<-unique(dsA[,get(timeC)])
#attributes(dsA)$timeS<-timeS
attributes(dsA)$grpV<-grpV
attributes(dsA)$depV<-depV
if(grepl("time",grpV[2])){code<-"timeSERIES"}else if (grepl("time",grpV[7])){code<-"byTrial"}
if (saveOn){saveRDS(dsA,paste0(expN,"_",modelV,"_",regs,"_",rsaValueN,"_",code,".rds"))}
#if (saveOn){write.table(dsA,paste0('DimensionSwitchPRP_S2ID_SOA.txt'), sep="\t",row.names=FALSE)} #be careful with naming 

# checking
print(dsA)

