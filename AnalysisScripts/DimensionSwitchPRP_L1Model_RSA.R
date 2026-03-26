# DimensionSwitchPRP_L1Model_RSA - amended by Ashwin 10.21.25 for merged data
# load pred file (contains single-trial accuracy & confidence)
# fit RSA models to decoding resulf of the conjunction variable
# _pred could be from: CROSS,JMB,SEG, and DYN!
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))
gc()

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
library(broom)
library(rhdf5)
library(caret)
library(foreach)
library(doMC)
library(tidyr)
library(binhf)
#library(GGally)
library(RColorBrewer)
library(RSQLite)
library(lazyeval)
library(stringr)
library(psych)
library(lme4)
library(RcppEigen)
library(rio)

# Source Files
setwd(Dir_R)
source('basic_lib.R')
registerDoMC(6)

# Setting for analysis session
fileN<-"DimensionSwitchPRP"
modelVs<-"TASKSET_CORSIDE_1"
sessN <- "nb_S1"
f2Load<-paste0("merged_",modelVs,'_',sessN, '_pred.rds'); #for merged data

# Setting for output database
# dbName<-paste0("FIT_RSA_",modelVs,'_S2_','_CTRR')
dbName<-paste0("FIT_RSA_TEST_CTRR")
rsaext = switch(modelVs, "TASKSET_CORSIDE_1"= "")  #confused about what this line does but keeping it....
saveOn<-T;

# # Load beh data
setwd(Dir_BDATA)
ds_beh<-fread("DimensionSwitchPRP_BehPP.txt");
ds_beh[is.nanM(ds_beh)]<-NA
#ds_beh[,(c("TASKSET_1","TASKSET_2")):=NULL] why add this line Atsushi?
subs<-list.files(path=Dir_EDATA,pattern = "^A6\\d{2}$") #grabs only "A6xx", since merged files are there

# # Load Individuals' Conjunction RT (TSCONJ,DSCONJ,SSCONJ)
# ds_cj<-fread(sprintf("DimensionSwitchPRP_%s_CONTROL.txt", modelVs))

#===============================================================================
#Load RSA Models (adjusted for different decoding variables)
#===============================================================================
# # Adjust RSA model sources!
setwd(paste0(Dir_RSAMODEL,"/",modelVs[1]))
mlist <- list.files(pattern = ".txt")
mL<-lapply(as.list(mlist),function(f){read.table(f,header = F)})
names(mL)<-gsub("*.txt","",mlist) # # List of RSA models
for (m in 1:length(mlist)){assign(gsub(".txt","",mlist)[m],mL[[m]],envir=.GlobalEnv)}
#===============================================================================
#Open database & do single-trial RSA
#===============================================================================

# Open Database  (erase old database of the same datatype, if it exists!!)
setwd(Dir_GRAND);
if (file.exists(paste0(fileN,'_',modelVs,'_',dbName))){file.remove(list.files(pattern=paste0(dbName,"$")))}
dbCon=dbConnect(dbDriver("SQLite"), paste0(fileN,'_',modelVs,'_',dbName))# Open database connection


s<-1#length(subs)
cc<-1

for (s in 1:c(length(subs))){#c(1:length(subs))
  # Load data & set properties
  Dir_Data_i<-paste0(Dir_EDATA,fsep,subs[s],fsep,"DATASETS");
  setwd(Dir_Data_i);f2LoadL<-paste0(subs[s],"_",f2Load);
  subN <- gsub('[[:alpha:]]+', '', subs[s])
  s2_subN <- sub("[0-9]", "2", subN)
  s3_subN <- sub("[0-9]", "3", subN)
  
  # Transform trial-by-trial confusion pattern
  pred<-data.table(readRDS(f2LoadL));#pred<-rio::import(f2LoadL)
  mesVs<-colnames(pred)[colnames(pred) %in% LETTERS_ex(50)]
  idVs<-colnames(pred)[!colnames(pred) %in% c(mesVs,"acc")]
  nmbVs<-idVs[!idVs %in% c("FBAND","ELEC","CV")]# avoid character columns in idVs
  timeVs<-nmbVs[grepl('time',nmbVs)]# time coding columns
  #castF1<-paste0("CLASS~",paste0(idVs[!idVs %in% c("SUBID")],collapse="+"))
  castF1<-paste0("CLASS~",paste0(idVs, collapse="+"))
  castF2<-paste0(paste0(idVs,collapse="+"),"~vars")
  
  # Step1: Melt data from wide (A:L) to long format (BLOCK,TRIAL,time) #confusion profile for all Classes for all timepoints for all trials (all 3 sessions)
  pred<-data.table::melt(pred,#this is way faster than gather/spread 
                         id.vars = idVs,
                         measure.vars=mesVs,
                         variable.name="CLASS",value.name="PP")
  
  # Step1':Log/Logit transform
  pred[PP==0,PP:=1e-10];pred[PP==1,PP:=1-1e-10]
  #pred$PP<-psych::logit(pred$PP);
  pred$PP<-log(pred$PP);
  
  # Check if all classes have complete cases 
  if (any(c(is.infinite(pred$PP) | is.nan(pred$PP)))){
    pred<-pred[!c(is.infinite(pred$PP) | is.nan(pred$PP)),]
    pred[,cmp:=uniqueN(CLASS),by=c("BLOCK","TRIAL",timeVs)]
    pred<-pred[cmp==length(mesVs),];pred[,cmp:=NULL]
  }
  
  # Step2: Add condition of Modeled variable from behavior (if necessary!)
  if (!any(grepl("obs",names(pred)))){
    ds_behi<-ds_beh[SUBID_S==subN|SUBID_S==s2_subN|SUBID_S==s3_subN,c(c("SUBID_S","BLOCK","TRIAL"),modelVs),with=F]
    pred<-merge(pred,ds_behi,by=c("SUBID_S","BLOCK","TRIAL"));setnames(pred,modelVs,"obs")
    } #generates a 'key' -- correct modelV class for all trials 
  
  # Step3-1:Long format (BLOCK,TRIAL,time) to another wide format (A:L to all combination of BLOCK,TRIAL,time)
  # Then, add model terms (e.g.,TASK,STIM,RESP,CONJ)
  if(is.factor(pred$obs)){pred[,obs:=as.numeric(obs)]}
  clist<-which(mesVs %in% LETTERS_ex(50))
  
  # Step (4): Prepare control vector of RT
  #CJCTR_RT<-scale(ds_cj[ds_cj$SUBID==subN,]$RT)[,1];
  #CJCTR_ACC<-scale(ds_cj[ds_cj$SUBID==subN,]$ACC)[,1];
  
  # Fit regression models condition-wise all at once
  results<-
    foreach(cc=clist) %dopar% {
      
      # Get condition specific data
      d<-pred[obs==cc,];
      predF<-data.table::dcast(d,castF1, value.var="PP")
      
      # Add regressors (inflexible but fast!)
      #predF$COLOR_RSA<-F_COLOR_2_M[,cc] # COLOR
      #predF$NUMBER_RSA<-F_NUMBER_2_M[,cc] # NUMBER
      predF$TASKSET_RSA <- TASKSET_1_M[,cc] #TASKSET
      predF$CORSIDE_RSA <- CORSIDE_1_M[,cc] #CORSIDE
      # predF$LH_RSA<-FV_LH_2_M[,cc] # LOW_HIGH
      # predF$OE_RSA<-FV_OE_2_M[,cc] # ODD_EVEN
      # predF$RB_RSA<-FV_RB_2_M[,cc] # RED_BLUE
      # predF$BF_RSA<-FV_BF_2_M[,cc] # BOLD_FAINT
      #predF$CTRRT_RSA<-CJCTR_RT;#CONTROL RT MODEL
      #predF$CTRACC_RSA<-CJCTR_ACC;#CONTROL MODEL
      

      # Step3-2:Regression all at once!
      dvIdx<-grepl("CLASS|*._RSA$|*._RSA_CFR",names(predF))
      Y<-as.matrix(predF[,which(!dvIdx),with=FALSE]);#DVs(this indexing takes time...)
      X<-as.matrix(predF[,which(dvIdx)[-1],with=FALSE]);#IVs(first var is CLASS!)
      # r<-coeff(.lm.fit(cbind(1,X),Y));#fastest, but only gives unstandardized coefficients
      r<-ls.print(lsfit(X,Y),print.it=F);#includes intercept!

      # Step3-3:Summarize results!
      names(r$coef.table)<-colnames(predF)[!dvIdx]
      estList<-rownames(r$coef.table[[1]]);#list of vars for estimates
      r<-rbindlist(lapply(r$coef.table,as.data.frame),idcol=TRUE)
      setnames(r,old=c("t-value","Pr(>|t|)","Std.Err"),new=c("tvalue","pvalue","SE"))
      
      #Rename variables and add basic identifier variables
      r[,c(modelVs):=list(cc)]
      r[,vars:=rep(estList,dim(r)[1]/length(estList))]
      
      #Restore labels from list names, then convert to wide format
      bltrtime<-str_split_fixed(r$.id,"_",n=Inf);#no assumption for .id tokens
      r[,c(idVs):=narray::split(bltrtime,along=2)]# assume SUBID is the first column
      r[,(nmbVs):=lapply(.SD,as.numeric),.SDcols=nmbVs]
      r<-r[vars!="Intercept",c(idVs,"vars","tvalue","Estimate"),with=F]
      
      # Need wide format(dcast) // long-wide(dcast) & wide-long(melt)
      r_wide<-dcast(r,castF2, value.var=c("tvalue")) #result is t-value of each predictor vector (stim feature) per timepoint per trial of modelV value
      return(r_wide)
    }
  
  # Summarize all results!
  results <- rbindlist(results)
  setorderv(results,nmbVs)
  depV<-colnames(results)[!colnames(results) %in% c(idVs,"CTRRT_RSA","CTRACC_RSA")]
  
  # Detect outliers within each time point and replace with NA
  for (cc in depV){
    results[,std:=lapply(.SD,sd),by=c(timeVs),.SDcols=cc]
    results[,bad:=abs(results[[cc]]) > std*4] # outlier based on std
    results[,(cc):=ifelse(bad,NA,results[[cc]])]
    print(sprintf("For %s, %f2 percent of rows were removed",cc,(length(which(results$bad))/dim(results)[1])*100))
    results<-results[bad!=T,] # remove bad trials altogether??
    #results[,bad:=outlier(results[[cc]],logical=T)] # outlier based on package(only one value??)
  }
  
  # Merge to behavioral template and put into DB
  results[,c("std","bad"):=NULL]# reduce data!
  predB<-merge(results,ds_beh,by=c("SUBID","SUBID_S","BLOCK","TRIAL"));#use ds_beh(containing all vars)
  dbWriteTable(dbCon,name=dbName,value=predB,row.names=FALSE,append=TRUE);#Try to append!
}

 # Checking
# dscheck=as.data.table(dbGetQuery(dbCon,paste0('SELECT * FROM ',dbName,'  LIMIT 5')));
# beepr::beep()
