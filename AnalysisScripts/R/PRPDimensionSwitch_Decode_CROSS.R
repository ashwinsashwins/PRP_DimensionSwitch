# HAL2017_ERROR02_Decode_CROSS
# peforms basic time-series decoding with mcca denoising and cross-validation
# NOTE:
# -
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))
gc()

#Setting up directory
fsep<-.Platform$file.sep;
Dir_R<-path.expand("/users/asrini17/data/asrini17/R_Helper_Scripts")
Dir_BDATA<-path.expand("/users/asrini17/data/asrini17/DimensionSwitch_PRP/DataAnalysis/Matlab/BEH/z_ALLGRAND(EEG,6sessions)")
Dir_EDATA<-path.expand("/users/asrini17/data/asrini17/DimensionSwitch_PRP/DataAnalysis/Matlab/EEG")
Dir_GRAND<-paste0(Dir_EDATA,"/w_ALLGRAND")
# fsep<-.Platform$file.sep;
# userName <- "Ashwin Srinivasan"; # Atsushi Kikumoto/Ashwin Srinivasan
# if (userName == "Atsushi Kikumoto") {Dir_R<-path.expand("~/Dropbox/w_ONGOINGRFILES/w_OTHERS")
# }else if(userName == "Ashwin Srinivasan"){Dir_R<-(path.expand("/users/asrini17/data/asrini17/R_Helper_Scripts"))
# }else{print("uh oh! check for typos")}
# Dir_HERE<-dirname(rstudioapi::getSourceEditorContext()$path)
# Dir_BDATA<-file.path(dirname(Dir_HERE),"Matlab","BEH","z_ALLGRAND(EEG,6sessions)")
# Dir_EDATA<-file.path(dirname(Dir_HERE),"Matlab","EEG")
# Dir_GRAND<-paste0(Dir_EDATA,"/w_ALLGRAND")

# Load libraries
library(data.table)
#library(tidyverse)
# library(narray) # don't load
library(dplyr)
library(tidyr)
library(collapse)
library(broom)
library(rhdf5)
library(caret)
library(foreach)
library(doMC)
library(binhf)
library(pROC)
library(RColorBrewer)
library(scales)
library(forcats)
library(stringr)
library(abind)
library(sjmisc)

# Source Files
setwd(Dir_R)
source('basic_theme.R')
source('basic_lib.R')
#try(source('noiseTools_lib.R'))
registerDoMC(4);

#Load behavioral data
setwd(Dir_BDATA)
ds_b<-fread("DimensionSwitchPRP_BehPP_LONG.txt")
#ds_b <-fread("DimensionSwitchPRP_BehPP.txt")
ds_b[is.nanM(ds_b)]<-NA;


#Analysis Setting
setwd(Dir_EDATA)
subs<-paste0("A",unique(ds_b$SUBID));
sessN<-  c("nb_CROSS") #"nbB_CROSS";#nb,nbRLLONG,nbRA,nbRLONGRA = regular decoding, nppr = shuffled label
eegV<-"eegpower" #eegraw,eegpower
modelV<-"F_COLOR"#TSCONJ,TRCONJ,ERCONJ,WRCONJ,CFRCONJ,RESP_ED
crossV<-"STIMNUM"

ds_b[, balance := paste0(TASKSET_CORSIDE, NUMTASK)]
balanceV<-"balance" #RESP_ED, ACCRESP
cutV<-NA;# 1288 for HAL2017_BI
evLock <- "RESP" # "STIM"

saveACC<-T;
savePRED<-T;
saveCM<-F;
saveIMP<-F;
saveROC<-F;
saveCOEF<-T;

# [X]RLED, [X]RESPLONG, RL


#Analysis Setting(time)
#timeL<-seq(-200,700,4);timeS<-seq(-200,700,12) #AS changed 3.4. to get 226 time points Stim_LK/RespLK
if(evLock == "RESP"){timeL<-seq(-700,200,4);timeS<-seq(-700,200,12)
} else if(evLock == "STIM"){timeL<-seq(-200,700,4);timeS<-seq(-200,700,12)}
#if (grepl("._RL$ | ._BRL$",sessN)){timeL<-seq(-800,200,4);timeS<-seq(-800,200,12);} #Resp_LK
#if (grepl("._RLED$ | ._BRLED$",sessN)){timeL<-seq(-800,200,4);timeS<-seq(-800,200,12);} #EDResp_LK
timeSS<-cbind(timeS[1:length(timeS)-1],timeS[2:length(timeS)])
timeSSL<-apply(timeSS,1,mean)
timeAV<-as.numeric(cut(timeL,timeS,include.lowest=T))
timeAV[is.na(timeAV)] <- 0

# Classification Setting
#"lda" = linear discriminant analysis
#"knn" = k nearest means
#"svmRadial" = support vector machine
#"nb" = naive bayes
#"rf" = random forest
#"svmLinear3" = L2 Regularized Support Vector Machine (dual) with Linear Kernel
#"pda" = penalized linear discriminant analysis
method<-'pda';# dda, lda, pda, naive_bayes
formula<-as.formula(paste(modelV,' ~ .'))
metric<-"Accuracy";
control<-trainControl(method="repeatedcv",
                      number=5,repeats=50,
                      selectionFunction = "oneSE",
                      sampling = "down",# rose,smote,down,up: how to balance # of observations for labels
                      classProbs=TRUE,allowParallel=TRUE,savePredictions=TRUE)

#Feature labels
freqL=c("Delta","Theta","Alpha","Beta","Gamma");
#elecL<-fread(paste0(Dir_R,fsep,"chanlocs_32MR_BROWN.txt"))
elecL<-fread(paste0(Dir_R,fsep,"chanlocs_32E_RIKEN.txt"))
elecL<-elecL[!Elec %in% c("A1","GND","EOG"),]# A2 was included 
varL = expand.grid(elec = as.vector(elecL$Elec),freq = freqL)
varL = str_c(varL$freq,"_",varL$elec)

# MERGE two MRI session!
#s<-"A201"
s <- paste0("A",Sys.getenv("SLURM_ARRAY_TASK_ID"))
tt<-1
trn_i<-"1"
which(subs %in% "A543")


# IN THE LOOP FOR SUBJECT!! 
for (s in s){ # completed upto A819
  # STEP 1: Merging data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load EEG data via HDF5 file
  # Dir_Data_i<-paste0(Dir_EDATA,fsep,s,fsep,"DATASETS")
  # setwd(Dir_Data_i);H5close();
  # dataF<-'*FFB_STIMLK.h5'# by deault use STIMLK data
  # if (grepl(".RL|RLRA",sessN)){dataF<-"*FFB_RESPLK.h5"}
  # if (grepl(".RLLONG|RLLONGRA",sessN)){dataF<-"*FFB_RESPLONGLK.h5"}
  # if (grepl(".RLED",sessN)){dataF<-"*FFB_RESPEDLK.h5"}
  # if (grepl(".GO",sessN)){dataF<-"*FFB_GOLK.h5"}
  
  Dir_Data_i<-file.path(Dir_EDATA,s,"DATASETS")
  #getwd()
  #print(Dir_Data_i)
 
  ses2 <- sub("[0-9]", "3", s)
  ses3 <- sub("[0-9]", "6", s)
  sessions <- c(s, ses2,ses3)
  bl_list <- vector("list",3)
  eeg_list <- vector("list",3)
  
  #Trying by loading in all EEG data and then running prep
  for (ses in sessions){ # just loading in EEG data in this loop
    Dir_Data_i<-paste0(Dir_EDATA,fsep,ses,fsep,"DATASETS") #sets directory for specified SUBID. copy this line for linked subids?
    #print(Dir_Data_i)
    setwd(Dir_Data_i);H5close();
    # if (evLock == "STIM") {
    #   if (grepl(".S2$",sessN)){dataF<-"*FFB_STIM2LK.h5"}else{dataF<-'*FFB_STIM1LK.h5'}} else if (evLock == "RESP"){
    #     if (grepl(".S2$",sessN)){dataF<-"*FFB_RESP2LK.h5"}else{dataF<-'*FFB_RESP1LK.h5'}}
    
    
     if (evLock == "STIM") {dataF1 <- "*FFB_STIM1LK.h5"; dataF2 <- "*FFB_STIM2LK.h5"}
     if (evLock == "RESP"){dataF1 <- "*FFB_RESP1LK.h5"; dataF2 <- "*FFB_RESP2LK.h5"}
     f1load<-list.files(path=Dir_Data_i,pattern=dataF1,full.names=TRUE)
     f2load<-list.files(path=Dir_Data_i,pattern=dataF2,full.names=TRUE)
     ds_eeg<-lapply(list(f1load,f2load), function(f){h5read(f,eegV)});#ds_info<-h5ls(f2load);
     ds_eeg <- abind::abind(ds_eeg, along = 1) # merge list into one big array with both stim periods
    
     btIDX<-lapply(list(f1load,f2load), function(f){data.table(h5read(f,"IDX")) %>%dplyr::rename(BLOCK=V1,TRIAL=V2)});
     #btIDX <- bind_rows(`names<-`(btIDX, unique(ds_b$EXP)), .id = "EXP") #EXP, in this case, is trials(?) in S1 vs S2 period #what is this line doing? 
     btIDX <- bind_rows(`names<-`(btIDX, unique(ds_b$STIMNUM)), .id = "STIMNUM")
     btIDX$STIMNUM <- as.integer(btIDX$STIMNUM)
     
     # Get Individual data to match to EEG data (this will reflect EEG artifact rejection)
     subN<-as.numeric(gsub("A","",ses));# Filter data and keep indexes
     dimL <- str_split(h5readAttributes(f2load,"/")$dim_label,"_") %>% unlist()
     timeIDX <- 1:dim(ds_eeg)[str_detect(dimL,"time")];
     elecIDX <- 1:dim(ds_eeg)[str_detect(dimL,"chan")]
     if (any(dimL=="freq")){freqIDX <- 1:dim(ds_eeg)[str_detect(dimL,"freq")]} 
     ds_bl<-merge(btIDX,ds_b[SUBID_S==subN],by=c('STIMNUM','BLOCK','TRIAL'))
   
    
    id <- match(ses,sessions)
    #print(dim(ds_eeg))
    bl_list[[id]] <- ds_bl
    eeg_list[[id]] <- ds_eeg
  }    
  #bind lists together
  ds_bl  <- data.table::rbindlist(bl_list, use.names = TRUE, fill = TRUE)
  ds_eeg <- abind::abind(eeg_list, along = 1) # merge along trials dimension
  rm(eeg_list); rm(bl_list); gc() #clear up memory!
  
  #Old method for single-session decoding ======================================
  # # 
  #   setwd(Dir_Data_i);H5close();
  # #if (evLock == "STIM") {dataF1 <- "*FFB_STIM1LK.h5"; dataF2 <- "*FFB_STIM2LK.h5"}
  #  if (evLock == "RESP"){dataF1 <- "*FFB_RESP1LK.h5"; dataF2 <- "*FFB_RESP2LK.h5"}
  #  f1load<-list.files(path=Dir_Data_i,pattern=dataF1,full.names=TRUE)
  #  f2load<-list.files(path=Dir_Data_i,pattern=dataF2,full.names=TRUE)
  #  ds_eeg<-lapply(list(f1load,f2load), function(f){h5read(f,eegV)});#ds_info<-h5ls(f2load);
  #  #ds_eeg[[2]] <- ds_eeg[[2]] [,1:226,,] # for testing - reduce timepoints in S2 period - ok this works but we lose a lot of time so....
  #  ds_eeg <- abind::abind(ds_eeg, along=1) 
  # # 
  # #have not yet merged sessions.... do that here#
  # # 
  #  btIDX<-lapply(list(f1load,f2load), function(f){data.table(h5read(f,"IDX")) %>%dplyr::rename(BLOCK=V1,TRIAL=V2)});
  #  #btIDX <- bind_rows(`names<-`(btIDX, unique(ds_b$EXP)), .id = "EXP") #EXP, in this case, is trials(?) in S1 vs S2 period #what is this line doing? 
  #  btIDX <- bind_rows(`names<-`(btIDX, unique(ds_b$STIMNUM)), .id = "STIMNUM")
  #  btIDX$STIMNUM <- as.integer(btIDX$STIMNUM)
  # 
  # 
  # # Get Individual data to match to EEG data (this will reflect EEG artifact rejection)
  # # ds_bl is typically smaller than original ds_b because btIDX reflects AR...
  #  subN<-as.numeric(gsub("A","",s));# Filter data and keep indexes
  #  dimL <- str_split(h5readAttributes(f2load,"/")$dim_label,"_") %>% unlist()
  #  timeIDX <- 1:dim(ds_eeg)[str_detect(dimL,"time")];
  #  elecIDX <- 1:dim(ds_eeg)[str_detect(dimL,"chan")]
  #  if (any(dimL=="freq")){freqIDX <- 1:dim(ds_eeg)[str_detect(dimL,"freq")]} 
  #  ds_bl<-merge(btIDX,ds_b[SUBID_S==subN],by=c('STIMNUM','BLOCK','TRIAL'))
  #=============================================================================
  
  
  # Applies 1) filtering, 2) cutting, 3) balancing, and 4) factorization
  prepSet <- list(
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$ACC_EDVALID==T & ds_bl$BLOCK>1 & ds_bl$OMIT_all==0 & (ds_bl$PRELATE==0 | is.na(ds_bl$PRELATE)),# exclude premature responses & complete omission!
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$ACC_EDVALID==T & ds_bl$BLOCK>1 & ds_bl$OMIT_all==0 & (ds_bl$PRELATE==0 | is.na(ds_bl$PRELATE)) & ds_bl$RESPALLOWED==1,
    filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$RESPORDER_K==1 & ds_bl$RESPOMIT==0 & ds_bl$RT < ds_bl$cutRT &  ds_bl$ACC==1,
    modelV = modelV, balanceV = balanceV, cutV = cutV, crossV=crossV, controlL=control)
  
  g(ds_bl, ds_eeg, filtIDX, foldsL) %=% preprop_decode(ds_bl,ds_eeg,prepSet)
  #ds_bl[,rowIndex:=.I] # add rowindex for later merging
  
  # Adjust dimensions of data without reducing time dimensions
  # d will be shifted from (obs, time, elec, fb) -> (time, obs, feature)
  if (eegV=="eegraw"){ds_eeg<-ds_eeg[,,elecIDX]} # index by dim?
  if (eegV=="eegpower"){ds_eeg<-ds_eeg[,,freqIDX,elecIDX];}
  # str_detect(dimL,"rpt|time"), str_detect(dimL,"rpt|time")
  dim(ds_eeg)<-c(dim(ds_eeg)[c(1:2)],prod(dim(ds_eeg)[3:4],na.rm=T))# obs, time, features 
  

  # STEP (2):Window Average~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Average over requested time windows
  ds_eeg<-narray::split(ds_eeg,which(str_detect(c("time"),dimL)),subsets=timeAV)
  ds_eeg[names(ds_eeg) %in% "0"]<-NULL;# remove unnecessary time
  ds_eeg<-lapply(ds_eeg, function(d){means.along(d,2)})
  ds_eeg<-aperm(do.call("abind",c(ds_eeg,list(along=3))),c(1,3,2))
  
  
  # STEP 3: Cross-condition Classification~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  registerDoMC(detectCores()-2);
  dimsTRN<-dim(ds_eeg);# always trial x time x feature!
  
  # Cross-condition settings (trn_tst)
  # ds_BL and ds_EL are prepared for test sets...?
  crossV_L <- as.character(unique(ds_bl[[crossV]]))
  crossV_P <- levels(interaction(crossV_L,crossV_L,sep="_"))
  
  # Focus on training on ACCRESP==0 (error commission) trials
  if (crossV=="ACCRESP" & modelV=="RESP_ED"){crossV_L <- "0";}
  
  results <-
    foreach(tt = 1:dimsTRN[2]) %dopar% {
      foreach(trn_i = unique(crossV_L)) %do% {
        # #Assign label (specific to cross-condition decoding)
        time_now <- mean(timeL[timeAV == tt])
        idx_trn <- foldsL[[trn_i]]
        trnIDX <- ds_bl[[crossV]] == as.numeric(trn_i) #separates row indices by crossV class
        crossIDX <- ds_bl[[crossV]] != as.numeric(trn_i)
        
        # Get target labels and data
        C_w <- ds_bl[[modelV]][trnIDX]
        C_c <- ds_bl[[modelV]][crossIDX]
        #if (grepl("nppr", sessN)) {C_w <- sample(C_w);C_c <- sample(C_c)}
        
        d_w <- data.table(ds_eeg[trnIDX, tt, ])[, (modelV) := C_w]
        d_c <- data.table(ds_eeg[crossIDX, tt, ])[, (modelV) := C_c]
        
        # Meta level function used later
        add_meta <- function(dt,tst,f)
        dt[,`:=`(SUBID=subN,time=time_now,TRN=as.numeric(trn_i),
                 TST=tst,WITHIN=f,CROSS=paste0(trn_i,"_",tst))]
        
        # Custom repeat loop to retain "cross" decoding results
        r_rep <-foreach(k = 1:control$repeats) %do% {
          # #Run classification
          rep_trn <- grepl(sprintf("Rep%d$",k),names(idx_trn))
          control_i <- control;if (!is_empty(foldsL)) {control_i$index <- idx_trn[rep_trn]}
          
          # # Train a model
          m <- train(formula, data = d_w, method = method, metric = metric, trControl = control_i)
          #print("trained!")
          
          # # Summarize results for within conditions
          r <- summary_decode(m, varL, list(PRED = T, CM = saveCM, IMP = saveIMP, ROC = saveROC, COEF = saveCOEF))
          r <- r[sapply(r, function(x) dim(x)[1]) > 0]  # remove empty slots
          r <- lapply(r, add_meta, tst = as.numeric(trn_i), f = 1L) 
          r$r_prob$rowIndex <- ds_bl$rowIndex[trnIDX][r$r_prob$rowIndex]
          
          # # Summarize results for cross conditions
          r_cross <- summary_cross_decode(m, d_c, ds_bl[crossIDX, .(rowIndex, obs = get(modelV))])
          tst_cross <- ds_bl[[crossV]][crossIDX] # label of all crossed conditions
          r_cross$r_prob <- add_meta(r_cross$r_prob, tst_cross, 0L)
          
          # # Bind all together
          list(r_acc = r$r_acc[,parameter:=NULL], r_prob = rbind(r$r_prob, r_cross$r_prob, use.names=TRUE))
        } 
        
        # # Average over repeats and merge
        grp_cols <- c("rowIndex","TRN","TST","WITHIN","CROSS")
        avg_rep <- function(lst,name,grp)rbindlist(lapply(lst,`[[`,name),idcol="REP")[,lapply(.SD,mean,na.rm=TRUE),by=grp,.SDcols=-"REP"]
        r_acc  <- avg_rep(r_rep,"r_acc",str_removeM(grp_cols,"rowIndex"))
        r_prob <- avg_rep(r_rep,"r_prob",grp_cols)
        r_f <- list(r_acc = r_acc, r_prob = merge(r_prob,ds_bl[, .(rowIndex, STIMNUM, SESSION, BLOCK, TRIAL, NUMTASK, SOA)], by="rowIndex",sort=FALSE))
        return(r_f)
      }
    }
  
  # Summarize all results
  rALL <- unlist(results, recursive = FALSE, use.names = TRUE)
  r_accG<-rbindlist(rALL$r_acc)
  
  # Result 1:Overall accuracy
  r_accG <- rbindlist(lapply(rALL, `[[`, "r_acc"))
  if (saveACC){saveRDS(r_accG,paste0(s,"_",modelV,"_",sessN,"_",crossV,"_accG.rds"))}
  
  # Result 2:Single-trial accuracy and confidence
  r_prob <- rbindlist(lapply(rALL, `[[`, "r_prob"))
  rAGG<-r_prob[,list(acc=mean(acc,na.rm=T),accSE=sd(acc, na.rm=TRUE) / sqrt(sum(!is.na(acc)))),by=c("time","CROSS","WITHIN")] #%>% print()
  if (savePRED){saveRDS(r_prob[,rowIndex:=NULL],paste0(s,"_",modelV,"_bal_",balanceV,"_CROSS_",crossV,"_",evLock,"_",sessN,"_pred.rds"))}

  # Checking
  r_prob[time==-190,.N,by=c("CROSS","WITHIN")]
  
  # Result 4:Importance map
  r_imp <- rbindlist(lapply(rALL, `[[`, "r_imp"))
  if (saveIMP){saveRDS(r_imp,paste0(s,"_",modelV,"_",sessN,"_",crossV,"_imp.rds"))}
  
  # Result 5:Coefficient(weight)
  #r_coef <- rbindlist(lapply(rALL, `[[`, "r_coef"))
  #if (saveCOEF){saveRDS(r_coef,paste0(s,"_",modelV,"_",sessN,"_",crossV,"_coef.rds"))}
  
  # Qucik check
  theme_set(theme_bw(base_size = 20))#32/28
  CSCALE_YlOrRd = rev(brewer.pal(9,"YlOrRd"));
  chance<-1/length(unique(ds_bl[[modelV]]))
  
  #quartz(width=7.5,height=4.5)
p <- ggplot(data=rAGG,aes(x=time,y=acc,ymin=acc-accSE,ymax=acc+accSE,group=CROSS,color=CROSS,fill=CROSS)) +
          geom_vline(xintercept=51,linetype=1,size=1)+annotate("text",x=150, y=chance-0.05,label="Stimulus + Cue")+
          geom_hline(yintercept=chance,linetype=1,size=1)+
          geom_ribbon(alpha=0.2,linetype=0)+#
          ggtitle(s)+
          geom_line(size=1.5)

#  quartz.save(paste0(s,"_",modelV,"_",sessN,".pdf"), type="pdf");
  
  
  outfile <- file.path(Dir_EDATA, paste0(s, "_merged_",modelV, "_", sessN, "_",crossV, ".png"))
  png(outfile, width = 6.5, height = 5, units = "in", res = 300)
  print(p)
  dev.off()
  
  message("Saved plot to: ", outfile)
  
}

#beepr::beep()
