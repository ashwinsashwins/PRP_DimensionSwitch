# DimensionSwitchPRP
# peforms basic time-series decoding
# NOTE:
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

# Load libraries
library(data.table)
#library(tidyverse)
library(tidyr)
library(dplyr)
library(broom)
library(rhdf5)
library(purrr)
library(caret)
library(foreach)
library(doMC)
library(binhf)
# library(pROC)
library(RColorBrewer)
library(scales)
# library(forcats)
library(stringr)


# Source Files
setwd(Dir_R)
source('basic_theme.R')
#source('basic_libS.R')
source('basic_libOLD.R')

#Load behavioral data
setwd(Dir_BDATA)
ds_b<-fread("DimensionSwitchPRP_BehPP.txt") 
ds_b[is.nanM(ds_b)]<-NA;

#Analysis Setting
setwd(Dir_EDATA)
#subs<-list.files( pattern="^A[[:digit:]]{3}"); don't use with slurm
eegV<-"eegpower" #eegraw
#modelV<-"CORSIDE_1"
modelV<-"TASKSET_CORSIDE_2"
#sessN<-"nb_S1";#nb = regular decoding, nppr = shuffled label
sessN<-"nb_S2";#nb = regular decoding, nppr = shuffled label # for stimLK 2


#Ashwin edits 9.22.26-----
ds_b[, balance := paste0(TASKSET_CORSIDE_2, NUMTASK)] #need to make this more robust... build into behavioral file i think
balanceV<- "balance"
evLock <- "RESP" # "STIM"
#----------
cutV<-NA;# 1288 for HAL2017_BI
saveACC<-F;
savePRED<-T;
saveCM<-F;
saveIMP<-F;
saveROC<-F;

# Classification Setting
#"lda" = linear discriminant analysis
#"knn" = k nearest means
#"svmRadial" = support vector machine
#"nb" = naive bayes
#"rf" = random forest
#"svmLinear3" = L2 Regularized Support Vector Machine (dual) with Linear Kernel
#"pda" = penalized linear discriminant analysis
method<-'pda';
formula<-as.formula(paste(modelV,' ~ .'))
metric<-"Accuracy";
control<-trainControl(method="repeatedcv",
                      number=5,repeats=10,
                      selectionFunction = "oneSE",
                      sampling = "down",# rose,smote,down,up: how to balance # of observations for labels
                      classProbs=TRUE,allowParallel=TRUE,savePredictions=TRUE)

#Feature labels
freqL=c("Delta","Theta","Alpha","Beta","Gamma");
elecL<-fread(paste0(Dir_R,fsep,"chanlocs_32E_RIKEN.txt"))
elecL<-elecL[!Elec %in% c("A1","GND","EOG"),]# A2 was included 
varL = expand.grid(elec = as.vector(elecL$Elec),freq = freqL)
varL = str_c(varL$freq,"_",varL$elec)

#s<-subs[1]
#subid is defined by array in SLURM
set.seed(612)
s <- paste0('A',Sys.getenv('SLURM_ARRAY_TASK_ID'))
t<-1

# IN THE LOOP FOR SUBJECT!!
#for (s in subs[c(1)]){ don't use with SLURM
for (s in s){
  # STEP 1: Merging data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load EEG data via HDF5 file
  Dir_Data_i<-paste0(Dir_EDATA,fsep,s,fsep,"DATASETS")
  setwd(Dir_Data_i);H5close();
  #if (grepl(".S2$",sessN)){dataF<-"*FFB_STIM2LK.h5"}else{dataF<-'*FFB_STIM1LK.h5'}
  if (grepl(".S2$",sessN)){dataF<-"*FFB_RESP2LK.h5"}else{dataF<-'*FFB_RESP1LK.h5'}
  
  f2load<-list.files(path=Dir_Data_i,pattern=dataF,full.names=TRUE)
  ds_eeg<-h5read(f2load,eegV);#ds_info<-h5ls(f2load);
  btIDX<-data.table(h5read(f2load ,"IDX"))%>%dplyr::rename(BLOCK=V1,TRIAL=V2);
  
  # Get Individual data to match to EEG data (this will reflect EEG artifact rejection)
  # ds_bl is typically smaller than original ds_b because btIDX reflects AR...
  subN<-as.numeric(gsub("A","",s));# Filter data and keep indexes
  ds_bl<-merge(btIDX,ds_b[SUBID_S==subN],by=c('BLOCK','TRIAL')) # use unique ID 
  
  # Applies 1) filtering, 2) cutting, 3) balancing, and 4) factorization
  prepSet <- list(
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$RESPORDER_K==1 & ds_bl$RESPOMIT==0 & ds_bl$RT_1 < ds_bl$cutRT1 & ds_bl$RT_2 < ds_bl$cutRT2,
    filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$RESPORDER_K==1 & ds_bl$RESPOMIT==0 & ds_bl$RT_1 < ds_bl$cutRT1 & ds_bl$RT_2 < ds_bl$cutRT2 & ds_bl$ACC_1==1 & ds_bl$ACC_2==1,
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$RESPORDER_K==1 & ds_bl$RESPOMIT==0 & ds_bl$RT_1 < ds_bl$cutRT1 & ds_bl$RT_2 < ds_bl$cutRT2 & ds_bl$ACC_1==1 & ds_bl$ACC_2==1 & ds_bl$NUMTASK=="double",
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$ACC_1==1 & ds_bl$ACC_2==1 ,
    modelV = modelV, balanceV = balanceV, cutV = cutV, controlL=control)
  
  g(ds_bl, ds_eeg, filtIDX, foldsL) %=% preprop_decode(ds_bl,ds_eeg,prepSet)
  print(sprintf("SUBID: %s DATA size: %d",s, dim(ds_bl)[1]))
  
  # STEP (2):Aggregation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # STEP 2: Classification~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  registerDoMC(8);
  ptm<- proc.time();
  dimS <- dim(ds_eeg)
  dimL <- str_split(h5readAttributes(f2load,"/")$dim_label,"_") %>% unlist()
  timeIDX <- 1:dimS[str_detect(dimL,"time")]
  elecIDX <- 1:dimS[str_detect(dimL,"chan")]
  if (any(dimL=="freq")){freqIDX <- 1:dimS[str_detect(dimL,"freq")]}
  #sprintf("Classification starts!! Retained data are:")
  print(sprintf(paste0(dimL,"= %d"),dimS))
  
 # shortTimeIDX <- timeIDX[1:]
  
  results<-
    foreach(t=timeIDX) %dopar% {
      #Get data of one point!
      if (eegV=="eegraw"){d<-ds_eeg[,t,elecIDX]} # index by dim?
      if (eegV=="eegpower"){d<-ds_eeg[,t,freqIDX,elecIDX];}
      dim(d)<-c(dim(d)[1],prod(dim(d)[2:3],na.rm=T))#adjust dimensions !
      
      #Assigne label
      C <- ds_bl[[modelV]];if(grepl("nppr",sessN)){C<-sample(C)} 
      d<-data.table(C=C,data.table(d));#make sure to data.table(d) before appending!
      setnames(d,"C",modelV)
      
      #Run classification
      if (is_empty(foldsL)){control[['index']] <- NULL}
      if (!is_empty(foldsL)){control[["index"]]<-foldsL}
      #if (!is.null(foldsL[1])){control[["index"]]<-foldsL}# removing this line allowed train() to run
      m<-train(formula, data=d, method=method, metric=metric, trControl=control);
      print("trained")
      
      # #Get all results
      r <- summary_decode(m,varL,list(PRED=T,CM=saveCM,IMP=saveIMP,ROC=saveROC))
      r <- r[sapply(r, function(x) dim(x)[1]) > 0];# remove empty slots
      r_f <- lapply(r, cbind,SUBID=unique(ds_bl$SUBID),SUBID_S=unique(ds_bl$SUBID_S),time=timeIDX[t])# don't lapply r_f repeatedly
      #print(class(ds_bl))
      #print(length(ds_bl$BLOCK))
      r_f$r_prob <- r_f$r_prob %>% mutate(BLOCK=ds_bl$BLOCK,TRIAL=ds_bl$TRIAL) # should be safe!
      return(r_f)
    }
  proc.time() - ptm
  
  # Summarize all results
  rALL<-as.list(as.data.frame(do.call(rbind, results)))
  
  # Result 1:Overall accuracy
  r_accG<-rbindlist(rALL$r_acc);
  if (saveACC){saveRDS(r_accG,paste0(s,"_unmerged_",modelV,"_",sessN, "_",evLock,"_",balanceV,"_accG.rds"))}
  
  # Result 2:Single-trial accuracy and confidence
  # Both acc and prob would be 0~1 since its cross-validated!
  r_pred<-rbindlist(rALL$r_prob) %>% dplyr::select(-rowIndex);
  if (savePRED){saveRDS(r_pred,paste0(s,"_unmerged_",modelV,"_",sessN, "_",evLock,"_",balanceV,"_pred.rds"))}
  
  # Result 3:Confusion matrix
  r_cm<-rbindlist(rALL$r_cm)
  if (saveCM){saveRDS(r_cm,paste0(s,"_",modelV,"_",sessN,"_cm.rds"))}
  
  # Result 4:Importance map
  r_imp<-rbindlist(rALL$r_imp)
  if (saveIMP){saveRDS(r_imp,paste0(s,"_",modelV,"_",sessN,"_cm.rds"))}
  
  # STEP 3: Quickly plotting?~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot related setting (for quick check!)
  theme_set(theme_bw(base_size = 20))#32/28
  CSCALE_YlOrRd = rev(brewer.pal(9,"YlOrRd"));
  chance<-1/length(unique(ds_bl[[modelV]]))
  
  # Overall Accuracy
  #quartz(width=7.5,height=4.5);theme_set(theme_bw(base_size = 20))#32/28
  p<-ggplot(data=r_accG,aes(x=time,y=Accuracy,ymin=Accuracy-AccuracySD,ymax=Accuracy+AccuracySD)) +
          geom_vline(xintercept=c(126),linetype=1,size=1)+annotate("text",x=150, y=chance-0.05,label="Stimulus + Cue")+
          geom_hline(yintercept=chance,linetype=1,size=1)+
          geom_ribbon(alpha=0.2,linetype=0,fill="red")+#
          ggtitle(s)+
          geom_line(size=1.5,color="red")
  
  # Save figure
  

  outfile <- file.path(Dir_EDATA, paste0(s, "_", modelV, "_", sessN, ".png"))
  
  png(outfile, width = 6.5, height = 5, units = "in", res = 300)
  print(p)
  dev.off()
  
  message("Saved plot to: ", outfile)
  #ggsave(filename = file.path(Dir_EDATA, paste0(s,"_",modelV,"_",sessN,".png")), plot = p, width = 6.5, height = 5)
}


# Send email to me 
# notifyM(paste0("HAL2017_",date()))
