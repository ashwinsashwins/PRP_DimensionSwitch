# DimensionSwitchPRP_Decode_CROSS(OSCAR)
# cross-temporal decoding (train all models at once method)
# same as DYN but saves _pred file as .rdf so use large time segments!
# 2-steps cross-validation
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))
gc()

# Libraries
library(data.table)
library(tidyverse)
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
library(viridis)
library(RSQLite)
library(here)
library(abind)
library(rlist)

# # #Setting up directory
# fsep <- .Platform$file.sep;
# Dir_S <- dirname(here())
# Dir_R <- paste("","gpfs","data","dbadre","akikumot","w_Support",sep=fsep)
# Dir_BDATA <- paste(Dir_S,"Data",sep=fsep)
# Dir_EDATA <- paste(Dir_S,"DataAnalysis",sep=fsep)

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


# Source Files
setwd(Dir_R)
source('basic_theme.R')
source('basic_lib.R')

#Load behavioral data
setwd(Dir_BDATA)
# ds_b<-fread("DimensionSwitchPRP_BehPP.txt") 
ds_b<-fread("DimensionSwitchPRP_BehPP_LONG.txt")
ds_b[is.nanM(ds_b)]<-NA;

#Analysis Setting
setwd(Dir_EDATA)
subs<-list.files( pattern="^A[[:digit:]]{3}");
eegV<-"eegpower" #eegraw
modelV<-"CORSIDE"
#stimN <- 2
#sessN<-"nb_S1";#nb = regular decoding, nppr = shuffled label
sessN<-"nb_S2";#nb = regular decoding, nppr = shuffled label
crossV<-""
cv_iter<-50;# number of repeated cv
cv_prop<-0.80;# proportion of data for cv
balanceV<-"NUMTASK" #NA
cutV<-NA;# 1288 for HAL2017_BI

# Saving results
saveACC<-F;
savePRED<-T;
saveCM<-F;
saveIMP<-F;
saveROC<-F;

#Analysis Setting(time)
# timeL<-seq(-200,800,4);timeS<-seq(-200,800,8) #Stim_LK
# if (grepl("._RL",sessN)){timeL<-seq(-800,200,4);timeS<-seq(-800,200,8);} #Resp_LK
# timeSS<-cbind(timeS[1:length(timeS)-1],timeS[2:length(timeS)])
# timeSSL<-apply(timeSS,1,mean)
# timeAV<-as.numeric(cut(timeL,timeS,include.lowest=T))
# timeAV[is.na(timeAV)] <- 0

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
                      sampling = "down",
                      classProbs=TRUE,allowParallel=FALSE,savePredictions=TRUE)

#Feature labels
freqL=c("Delta","Theta","Alpha","Beta","Gamma");
elecL<-fread(paste0(Dir_R,fsep,"chanlocs_32E_RIKEN.txt"))
elecL<-elecL[!Elec %in% c("A1","GND","EOG"),]# A2 was included 
varL = expand.grid(elec = as.vector(elecL$Elec),freq = freqL)
varL = str_c(varL$freq,"_",varL$elec)
# elecL<-fread(paste0(Dir_R,fsep,"chanlocs_20E_Oregon.txt")) #switch to "chanlocs_32E_RIKEN.txt"
# varL = expand.grid(elec = as.vector(elecL$Elec),freq = freqL)
# varL = str_c(varL$freq,"_",varL$elec)

# Functions for temporal generalization
genFF<-function(mm,dd,C_tst){
  # class, prob, obs, acc
  p<-predict(mm,newdata=data.table(dd),type="prob")
  p<-cbind(p,C_tst[,c("rowIndex","BLOCK","TRIAL")]) %>% data.table()
  p[,`:=`(prob=do.call(pmax,.SD),class=names(.SD)[max.col(.SD)]),.SDcols=names(p)%in%LETTERS]
  p[,`:=`(obs=C_tst[[modelV]],acc=(class==C_tst[[modelV]])*1)]
}

# SUBJECT is defined by ArrayJob in run_BI_DIM_rh.sh
set.seed(612);
s <- paste0("A",Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(s);t<-1;i<-1;#s<-"A201"

for (s in s){
  # STEP 1: Merging data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load EEG data via HDF5 file
  Dir_Data_i<-file.path(Dir_EDATA,s,"DATASETS")
  setwd(Dir_Data_i);H5close();
  #if (grepl(".S2$",sessN)){dataF<-"*FFB_STIM2LK.h5"}else{dataF<-'*FFB_STIM1LK.h5'}
  #if (grepl(".S2$",sessN)){dataF<-"*FFB_RESP2LK.h5"}else{dataF<-'*FFB_RESP1LK.h5'}
  dataF1 <- "*FFB_STIM1LK.h5"; dataF2 <- "*FFB_STIM2LK.h5" #from stim 1 and stim 2 period
  f1load<-list.files(path=Dir_Data_i,pattern=dataF1,full.names=TRUE)
  f2load<-list.files(path=Dir_Data_i,pattern=dataF2,full.names=TRUE)

  #ds_eeg<-h5read(f2load,"eegpower");ds_info<-h5ls(f2load);
  #eeg_1<-h5read(f1load,eegV);eeg_2<-h5read(f2load,eegV)
  ds_eeg<-lapply(list(f1load,f2load), function(f){h5read(f,eegV)});
  ds_eeg <- abind::abind(ds_eeg, along = 1); gc() #issue here is dimensions do not match -- # of trials is different & # of timepoints is different (226 v 301)
  
  btIDX<-data.table(h5read(f2load ,"IDX"))%>%dplyr::rename(BLOCK=V1,TRIAL=V2);
  #dimL <- str_split(h5readAttributes(f2load,"/")$dim_label,"_") %>% unlist()
  
  # Get Individual data to match to EEG data (this will reflect EEG artifact rejection)
  # ds_bl is typically smaller than original ds_b because btIDX reflects AR...
  subN<-as.numeric(gsub("A","",s));# Filter data and keep indexes
  ds_bl<-merge(btIDX,ds_b[SUBID_S==subN],by=c('BLOCK','TRIAL'))
  
  # Block/Trial order differs between ds_eeg and ds_bl. Fix!
  ds_bl <- ds_bl %>% arrange(STIMNUM)
  
  # Applies 1) filtering, 2) cutting, 3) balancing, and 4) factorization
  prepSet <- list(
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$RESPORDER_K==1 & ds_bl$RESPOMIT==0 & ds_bl$RT_1 < ds_bl$cutRT1 & ds_bl$RT_2 < ds_bl$cutRT2,
    filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$RESPORDER_K==1 & ds_bl$RESPOMIT==0 & ds_bl$RT_1 < ds_bl$cutRT1 & ds_bl$RT_2 < ds_bl$cutRT2 & ds_bl$ACC_1==1 & ds_bl$ACC_2==1, #issue. none exist
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$RESPORDER_K==1 & ds_bl$RESPOMIT==0 & ds_bl$RT_1 < ds_bl$cutRT1 & ds_bl$RT_2 < ds_bl$cutRT2 & ds_bl$ACC_1==1 & ds_bl$ACC_2==1 & ds_bl$NUMTASK=="double",
    #filtIDX =! is.na(ds_bl[[modelV]]) & ds_bl$ACC_1==1 & ds_bl$ACC_2==1 ,
    modelV = modelV, balanceV = balanceV, cutV = cutV)
  
  g(ds_bl, ds_eeg) %=% preprop_decode(ds_bl,ds_eeg,prepSet);
  ds_eeg <- dimAdj(ds_eeg, which(str_detect(c("freq|chan"),dimL)));
  
  # STEP (2):Window Average~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Average over requested time windows
  ds_eeg<-narray::split(ds_eeg,which(str_detect(c("time"),dimL)),subsets=timeAV)
  ds_eeg[names(ds_eeg) %in% "0"]<-NULL;# remove unnecessary time
  ds_eeg<-lapply(ds_eeg, function(d){means.along(d,2)})
  ds_eeg<-aperm(do.call("abind",c(ds_eeg,list(along=3))),c(1,3,2))
  
  # STEP 2: Classification~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  registerDoMC(6);
  rRESALL <- list()
  dimsTRN<-dim(ds_eeg);# always trial x time x feature!
  #paste(dimL,str_split(ds_info$dim[2],"x")[[1]],sep=" =") %>% paste(collapse = "")
  sprintf("Dyn-Classification starts!! %d percent of %d trials with %d iterations",cv_prop*100,dimsTRN[1],cv_iter)
  
  #  Loop over each time point
  for (t in 1:dimsTRN[2]){
    rRES<-
      foreach(i = 1:cv_iter) %dopar% {
        # Partition training and test set (1st = train, 2nd = test!!)
        # Reflecting crossing variables....!
        ds_bl[,gIDX:=2-(.I %in% sample(1:dimsTRN[1],floor(dimsTRN[1]*cv_prop)))]
        ds_bl[,gIDX:=interaction(get(crossV),gIDX,sep="_")] # include crossV: ds_bl[,list(.N),by=test]
        ds_BL<-split(ds_bl, list(factor(ds_bl$gIDX)))
        ds_EL<-narray::split(ds_eeg,1,subsets=ds_bl$gIDX)
        ds_BL<-ds_BL[order(names(ds_BL))]#necessary to align list index...
        ds_EL<-ds_EL[order(names(ds_EL))]#necessary to align list index...
        crossP<-names(ds_EL);crossL<-list();# all pattern of crossing
        
        # Loop over different levels of crossV
        for (trn_i in which(str_detect(crossP,"_1"))){
          # STEP1:Accecss data and do fitting attime point "t" (loop over)!
          # Do this for one specific level of crossV!
          C_trn<-ds_BL[[trn_i]][[modelV]]
          if(grepl("nppr",sessN)){C_trn<-sample(C_trn)} # Add correct label
          d<-data.table(C=C_trn,data.table(ds_EL[[trn_i]][,t,]));setnames(d,"C",modelV)
          m<-train(formula,data=d,method=method,trControl=control)
          
          for (tst_i in which(str_detect(crossP,"_2"))){
            # STEP2:Make predictions for each time points,then summarize!
            C_tst<-ds_BL[[tst_i]]# give entire data set!
            d_tst<-narray::split(ds_EL[[tst_i]],2,subsets=1:dimsTRN[2])
            rlist<-lapply(d_tst,function(d){genFF(m,d,C_tst)})# can this be faster??
            names(rlist) <- timeSSL;# label of time for tst samples
            
            # STEP3:Summarize results so far
            r <- bind_rows(rlist,.id="time_g") %>% data.table()
            r[,`:=`(time_g=as.numeric(time_g), CROSS=paste0(gsub("_1","",crossP[trn_i]),"_",gsub("_2","",crossP[tst_i])))]
            crossL<-rlist::list.append(crossL,r)
          }
        }
        return(rbindlist(crossL))
      } %>% bind_rows(.,.id="CV") %>% data.table()
    
    # Summarize results across cv iterations
    rRES[,`:=`(time_f=timeSSL[t],SUBID = subN)]
    depV <- c("acc","prob",names(rRES)[names(rRES)%in%LETTERS])
    gvar <- c("SUBID","BLOCK","TRIAL","obs","time_f","time_g","CROSS")
    rRESALL[t] <-list(rRES[,lapply(.SD,mean,na.rm=T),.SD=depV,by=gvar])
  } 
  
  # Summarize results over time samples
  r_ALL <- bind_rows(rRESALL) %>% data.table()
  r_AGG <- r_ALL[,list(acc=mean(acc)),by=c("time_f","time_g","CROSS")]
  
  # Save all results (r_AGG as .rds, r_PRED as h5 file)
  if (saveOn){
    # Aggregated results
    #saveRDS(r_AGG,paste0(s,"_",paste(modelV,collapse="_"),"_",sessN,"_DYN.rds"))
    
    # Save as h5 ?? or .RDS? 
    if (savePRED){
      #saveRDS(r_PRED,paste0(s,"_",modelV,"_",sessN,"_pred_DYN.rds"))
      saveRDS(r_ALL,paste0(s,"_",modelV,"_",crossV,"_",sessN,"_pred_CROSS.rds"))#for RSA and DB 
    }
    
    # # # Save as h5
    # h5fname<-paste0(s,"_",modelV,"_",sessN,"_DYN_pred.h5")
    # if (file.exists(h5fname)){file.remove(list.files(pattern=h5fname))}
    # h5createFile(h5fname)
    # timeDim <- n_distinct(r_PRED$time_f)
    # trialDim <- dim(r_PRED)[1]/(timeDim^2)
    # dimS<-c(trialDim,timeDim,timeDim)
    # for (cc in c("prob","acc")){h5write(array(r_PRED[[cc]],dim=dimS),h5fname,cc)}
    # #for (cc in colnames(r_pred)){h5write(array(r_pred[[cc]],dim=dimS),h5fname,cc)}
  }
  
  # Plot result
  theme_set(theme_bw(base_size = 16))#32/28
  actDF<-r_AGG;
  chanceV<-1/n_distinct(ds_bl[[modelV]])
  actDF$DV<-actDF$acc-chanceV
  
  #quartz(width=6,height=5.5);
  print(ggplot() +
          facet_wrap(~CROSS)+
          geom_raster(data=actDF,aes(time_f, time_g,fill=DV), interpolate = F)+
          #geom_tile(aes(fill = DV))+
          #geom_contour(data=actDF,aes(time_f,time_g,z=v),colour='white',alpha=1,bins=1) +
          #Axis + scales-------------------
          coord_cartesian(xlim=range(timeL),ylim=range(timeL)) +
          scale_x_continuous(expand=c(0,0),breaks=c(-300,0,800,1600))+scale_y_continuous(expand=c(0,0),breaks=c(-300,0,800,1600))+
          scale_fill_gradientn(limits=c(0,0.03),colours=inferno(30),oob=squish,breaks=c(0.33,0.35,0.37,0.39))+
          #Labels--------------------------
          ggtitle((s))+
          ylab("Testing time (ms)")+xlab("Training time (ms)")+
          geom_abline(intercept=0, slope = 1,size=0.8,color="white",linetype="dashed",alpha=1)+
          annotate("text",x=900,y=650,label="Distractor",size=6,fontface="bold",color="white")+
          annotate("text",x=130,y=-150,label="Stimulus",size=6,fontface="bold",color="white")+
          annotate("rect", xmin=0, xmax=Inf, ymin=0, ymax=Inf,linetype="dashed",color="white",size=0.8,alpha=0.8,fill=NA)+
          annotate("rect", xmin=800, xmax=Inf, ymin=800, ymax=Inf,linetype="dashed",color="white",size=0.8,alpha=0.8,fill=NA)+
          annotate("pointrange", x =-300, y = -300, ymin = 0, ymax = 0,linetype="dashed",colour = "gray40", size =1,alpha=0.8)+
          annotate("pointrange", x =0, y = 0, ymin = 0, ymax = 0,linetype="dashed",colour = "gray40", size =1,alpha=0.8)+
          annotate("pointrange", x =800, y = 800, ymin = 0, ymax = 0,linetype="dashed",colour = "gray40", size =1,alpha=0.8)+
          #Aesthetics!-------------------------
          theme(panel.border=element_rect(size=1.5,fill=NA),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.margin=unit(rep(1,1,4),"line"),
              axis.line=element_line(size=1),
              legend.key = element_blank(),
              #legend.key.size=unit(1,"lines"),
              legend.position = "right",
              legend.title=element_blank(),
              legend.text =element_text(family="Helvetica",size=16),
              axis.ticks.length = unit(0.25, "lines"),
              axis.title  = element_text(family="Helvetica",vjust=1.8),
              axis.text   = element_text(size=18),
              axis.text.x = element_blank(),
              plot.title = element_text(hjust = 0, size = 20),
              strip.background=element_blank())
  )
  # Save figure
  ggsave(file = paste0(s,"_",modelV,"_",crossV,"_",sessN,".png"), width = 6.5, height = 5)
}
