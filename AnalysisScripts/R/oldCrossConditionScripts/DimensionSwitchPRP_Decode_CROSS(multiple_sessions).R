# NEW METHOD DimensionSwitchPRP_Decode_CROSS(OSCAR)
# peforms basic time-series decoding with mcca denoising and cross-validation
# NOTE:
# -
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

# Load libraries
library(data.table)
library(tidyverse)
library(tidyr)
library(dplyr)
# library(narray) # don't load
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

# Source Files
setwd(Dir_R)
source('basic_theme.R')
source('basic_lib.R')
# source('noiseTools_lib.R')


#Load behavioral data
setwd(Dir_BDATA)
# ds_bs<-fread("DimensionSwitchPRP_BehPP.txt") 
ds_b<-fread("DimensionSwitchPRP_BehPP_LONG.txt")
ds_b[is.nanM(ds_b)]<-NA;
# registerDoMC(4);
# ds_b1<-fread(file.path(Dir_BDATA1,"HAL2017_ERROR_BehPP.txt"))
# ds_b1[,OMIT_ED:=NA] # not defined in ERROR1 study
# ds_b2<-fread(file.path(Dir_BDATA2,"HAL2017_ERROR2_BehPP.txt"))
# ds_b<-rbind(ds_b1,ds_b2)
# ds_b[is.nanM(ds_b)]<-NA;
# ds_b[,SUBID_K:=ifelse(EXP=="ERROR_EEG1",SUBID,SUBID-200)]

# # Identify subjects who have both sessions
# subDIFF <- setdiff(unique(ds_b1$SUBID),unique(ds_b2$SUBID)-200)
# print(paste0("Subjects missing both sessions:",paste(subDIFF,collapse = "_")))
 dsSUBLIST<-ds_b[,.(SUBID = unique(SUBID),SUBID_S=unique(SUBID_S)), by = SESSION]
# dsSUBLIST<-dsSUBLIST[!SUBID_K %in% subDIFF,]
 subs<-str_c("A",as.character(unique(dsSUBLIST$SUBID_S)));

#Analysis Setting
eegV<-"eegpower" #eegraw
modelV<- "CORSIDE"  #"RESP_ED"
sessN<-   "nb_S2" #"nbBRL_CROSS";#nb = regular decoding, nppr = shuffled label
crossV<- "STIMNUM" # figure out what crossV is used for... #"ACCRESP"
cv_iter<-?1
balanceV<- "NUMTASK" #"ACC_EDc" # RESP_ED, ACC_EDc
cutV<-NA
saveOn<-TRUE;
saveACC<-F;
savePRED<-T;
saveMCCA<-T;
saveCM<-F;
saveIMP<-F;
saveROC<-F;
saveCOEF<-F

#Analysis Setting(time)
timeL<-seq(-800,200,4);timeS<-seq(-800,200,30) #Stim_LK
if (grepl("._RL$",sessN)){timeL<-seq(-800,200,4);timeS<-seq(-800,200,4);} #Resp_LK
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
                      number=5,repeats=20,
                      selectionFunction = "oneSE",
                      sampling = "down",# rose,smote,down,up: how to balance # of observations for labels
                      classProbs=TRUE,allowParallel=TRUE,savePredictions=TRUE)

#Feature labels
freqL=c("Delta","Theta","Alpha","Beta","Gamma");
elecL<-fread(paste0(Dir_R,fsep,"chanlocs_32E_RIKEN.txt"))
#elecL<-fread(paste0(Dir_R,fsep,"chanlocs_32MR_BROWN.txt"))
elecL<-elecL[!Elec %in% c("A1","GND","EOG"),]# A2 was included 
varL = expand.grid(elec = as.vector(elecL$Elec),freq = freqL)
varL = str_c(varL$freq,"_",varL$elec)

# MERGE two MRI session!
#s<-"A201"
set.seed(612);
s <- paste0("A",Sys.getenv("SLURM_ARRAY_TASK_ID"))
tt<-1 #time index?
trn_i<-"1" #what the hell is this?

# IN THE LOOP FOR SUBJECT!! 
for (s in s){ # completed upto A819
  # STEP 1: Merging data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load EEG data via HDF5 file
  
  # dataF<-'*FFB_STIMLK.h5'# by deault use STIMLK data
  # if (grepl(".RL|RLRA",sessN)){dataF<-"*FFB_RESPLK.h5"}
  # if (grepl(".RLLONG|RLLONGRA",sessN)){dataF<-"*FFB_RESPLONGLK.h5"}
  # if (grepl(".GO",sessN)){dataF<-"*FFB_GOLK.h5"}
  Dir_Data_i<-file.path(Dir_EDATA,s,"DATASETS")
  setwd(Dir_Data_i);H5close();
  dataF1 <- "*FFB_STIM1LK.h5"; dataF2 <- "*FFB_STIM2LK.h5"
  f1load<-list.files(path=Dir_Data_i,pattern=dataF1,full.names=TRUE)
  f2load<-list.files(path=Dir_Data_i,pattern=dataF2,full.names=TRUE)
  
  subNK <- dsSUBLIST[SUBID_S==as.numeric(gsub("A","",s))]
  subNK_L <- str_c("A",subNK$SUBID)
  # f2load1<-list.files(path=file.path(Dir_EEG1,subNK_L[1],"DATASETS"),pattern=dataF,full.names=TRUE) 
  # f2load2<-list.files(path=file.path(Dir_EEG2,subNK_L[2],"DATASETS"),pattern=dataF,full.names=TRUE) 
  ds_eeg<-lapply(list(f1load,f2load), function(f){h5read(f,eegV)});#ds_info<-h5ls(f2load);
  ds_eeg[[2]] <- ds_eeg[[2]] [,1:226,,] # for testing - reduce timepoints in S2 period - ok this works but we loose a lot of time so....
  ds_eeg <- abind::abind(ds_eeg, along=1)  #issue here is dimensions do not match -- # of trials is different & # of timepoints is different (226 v 301)
  
  # Prepare labels of dimensions and merging keys
  # Bind multiple files (e.g., SESSIONs) together
  dimL <- str_split(h5readAttributes(f1load,"/")$dim_label,"_") %>% unlist()
  timeIDX <- 1:dim(ds_eeg)[str_detect(dimL,"time")];
  elecIDX <- 1:dim(ds_eeg)[str_detect(dimL,"chan")]
  if (any(dimL=="freq")){freqIDX <- 1:dim(ds_eeg)[str_detect(dimL,"freq")]} 
  btIDX<-lapply(list(f1load,f2load), function(f){data.table(h5read(f,"IDX")) %>%dplyr::rename(BLOCK=V1,TRIAL=V2)});
  #btIDX <- bind_rows(`names<-`(btIDX, unique(ds_b$EXP)), .id = "EXP") #EXP, in this case, is trials(?) in S1 vs S2 period #what is this line doing? 
  btIDX <- bind_rows(`names<-`(btIDX, unique(ds_b$STIMNUM)), .id = "STIMNUM")
  btIDX$STIMNUM <- as.integer(btIDX$STIMNUM)
  
  # Get Individual data to match to EEG data (this will reflect EEG artifact rejection)
  # ds_bl is typically smaller than original ds_b because btIDX reflects AR...
  ds_bl<-bind_rows(base::split(ds_b[SUBID_S==subNK$SUBID_S[1]],by="SUBID"))
  ds_bl<-merge(btIDX,ds_bl,by=c('BLOCK','TRIAL', 'STIMNUM')) # use unique ID #added 'stimnum' because otherwise it creates a duplicate column
  if (all(na.omit(unique(ds_bl[[crossV]])) %in% c(0,1))){ds_bl[[crossV]]<-ds_bl[[crossV]]+1}
  
  # Applies 1) filtering, 2) cutting, 3) balancing, and 4) factorization
  # filtIDX_C = ! is.na(ds_bl[[modelV]]) & ds_bl$ACC_EDVALID==T & ds_bl$BLOCK>1 & ds_bl$OMIT_all==0 & (ds_bl$PRELATE==0 | is.na(ds_bl$PRELATE)) & ds_bl$RESPALLOWED==1
  # if (modelV=="PTSRCONJ"){filtIDX_C=!is.na(ds_bl[[modelV]]) & ds_bl$ACC==1  & ds_bl$PTSRGROUP==as.numeric(str_extract(sessN,"\\d"))}
  # prepSet <- list(filtIDX = filtIDX_C, modelV = modelV, balanceV = balanceV, cutV = cutV, crossV = crossV, controlL = control)
  prepSet <- list(
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$RESPORDER_K==1 & ds_bl$RESPOMIT==0 & ds_bl$RT_1 < ds_bl$cutRT1 & ds_bl$RT_2 < ds_bl$cutRT2,
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$RESPORDER_K==1 & ds_bl$RESPOMIT==0 & ds_bl$RT_1 < ds_bl$cutRT1 & ds_bl$RT_2 < ds_bl$cutRT2 & ds_bl$ACC_1==1 & ds_bl$ACC_2==1,
    filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$RESPORDER_K==1 & ds_bl$RESPOMIT==0 & ds_bl$RT < ds_bl$cutRT &  ds_bl$ACC==1,
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$RESPORDER_K==1 & ds_bl$RESPOMIT==0 & ds_bl$RT_1 < ds_bl$cutRT1 & ds_bl$RT_2 < ds_bl$cutRT2 & ds_bl$ACC_1==1 & ds_bl$ACC_2==1 & ds_bl$NUMTASK=="double",
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$ACC_1==1 & ds_bl$ACC_2==1 ,
    modelV = modelV, balanceV = balanceV, cutV = cutV, controlL=control) # need to add crossV in here somewhere
  
  # Prepare for decoding (ds_eeg: trials x time x freq x elec)
  g(ds_bl, ds_eeg, filtIDX, foldsL) %=% preprop_decode(ds_bl,ds_eeg,prepSet) #preprop_decode function does not have crossV at all -- does this need to change?
  ds_bl[,rowIndex:=.I] # add rowindex for later merging
  
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
  
  results <-
    foreach(tt = 1:dimsTRN[2]) %dopar% {
      foreach(trn_i = unique(crossV_L)) %do% {
        # #Assign label
        C <- ds_bl[[modelV]]; if (grepl("nppr", sessN)) {C <- sample(C)}
        d <- data.table(ds_eeg[, tt, ])[, (modelV) := C]
        
        # #Run classification
        control_i <- control
        if (!is_empty(foldsL)) {control_i$index <- lapply(foldsL, `[[`, trn_i) }
        #if (!is_empty(is_empty(foldsL))){control_i$indexOut<-lapply(foldsL,`[[`, tst_i)} # only needed to constrain tst within cv
        m <- train(formula, data = d, method = method, metric = metric, trControl = control_i)
        
        # # Summarize results
        r <- summary_decode(m, varL, list(PRED = T, CM = saveCM, IMP = saveIMP, ROC = saveROC, COEF = saveCOEF))
        r <- r[sapply(r, function(x) dim(x)[1]) > 0]  # remove empty slots
        r_f <- lapply(r, cbind, SUBID_K = subNK$SUBID_K[[1]], time = mean(timeL[timeAV == tt]))
        r_f$r_prob <- merge(ds_bl[, .(rowIndex, EXP, BLOCK, TRIAL)], r_f$r_prob, by = "rowIndex", sort = FALSE, all.x = TRUE)
        
        # # Add information about crossing conditions
        testV <- ds_bl[[crossV]]
        r_f$r_prob[, `:=`(TRN = (trn_i), TST = testV), ]
        r_f$r_prob[, `:=`(WITHIN = as.integer((trn_i) == testV), CROSS = paste0((trn_i), "_", testV)), ]
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
  rAGG<-r_prob[,list(acc=mean(acc,na.rm=T),accSE=sd(acc, na.rm=TRUE) / sqrt(sum(!is.na(acc)))),by=c("time","CROSS","WITHIN")] %>% print()
  if (savePRED){saveRDS(r_prob[,rowIndex:=NULL],paste0(s,"_",modelV,"_",sessN,"_",crossV,"_pred.rds"))}

  # Result 4:Importance map
  r_imp <- rbindlist(lapply(rALL, `[[`, "r_imp"))
  if (saveIMP){saveRDS(r_imp,paste0(s,"_",modelV,"_",sessN,"_",crossV,"_imp.rds"))}
  
  # Result 5:Coefficient(weight)
  r_coef <- rbindlist(lapply(rALL, `[[`, "r_coef"))
  if (saveCOEF){saveRDS(r_coef,paste0(s,"_",modelV,"_",sessN,"_",crossV,"_coef.rds"))}
  
  # Qucik check
  theme_set(theme_bw(base_size = 20))#32/28
  CSCALE_YlOrRd = rev(brewer.pal(9,"YlOrRd"));
  chance<-1/length(unique(ds_bl[[modelV]]))
  
  quartz(width=7.5,height=4.5)
  print(ggplot(data=rAGG,aes(x=time,y=acc,ymin=acc-accSE,ymax=acc+accSE,group=CROSS,color=CROSS,fill=CROSS)) +
          geom_vline(xintercept=51,linetype=1,size=1)+annotate("text",x=150, y=chance-0.05,label="Stimulus + Cue")+
          geom_hline(yintercept=chance,linetype=1,size=1)+
          geom_ribbon(alpha=0.2,linetype=0)+#
          ggtitle(s)+
          geom_line(size=1.5))
  
  quartz.save(paste0(s,"_",modelV,"_",sessN,".pdf"), type="pdf");
  
  
}

beepr::beep()
