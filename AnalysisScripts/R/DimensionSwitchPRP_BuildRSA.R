# DimensionSwitchPRP_BuildRSA
# load pred file (contains single-trial accuracy & confidence)
# fit RSA models to decoding resulf of the conjunction variable
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))

#Setting up directory
fsep<-.Platform$file.sep;
userName <- "Ashwin Srinivasan"; # Atsushi Kikumoto/Ashwin Srinivasan
if (userName == "Atsushi Kikumoto") {Dir_R<-path.expand("~/Dropbox/w_ONGOINGRFILES/w_OTHERS")
}else if(userName == "Ashwin Srinivasan"){Dir_R<-file.path("/users/asrini17/data/asrini17/R_Helper_Scripts")
}else{print("uh oh! check for typos")}
Dir_HERE<-dirname(rstudioapi::getSourceEditorContext()$path)
Dir_BDATA<-file.path(dirname(Dir_HERE),"Matlab","BEH","z_ALLGRAND(EEG,6sessions)")
Dir_EDATA<-file.path(dirname(Dir_HERE),"Matlab","EEG")
Dir_GRAND<-paste0(Dir_EDATA,"/w_ALLGRAND")
Dir_RSAMODEL<-paste0(Dir_EDATA,"/w_RSAMODELS")

# Load libraries
library(data.table)
library(binhf)
library(RSQLite)
library(lazyeval)
library(stringr)
library(psych)

# Source Files
setwd(Dir_R)
source('basic_lib.R')

# Load behavioral template
ds<-fread(file=file.path(Dir_BDATA,"DimensionSwitchPRP_BehPP.txt"))


 #===============================================================================
# RSA using S1ID/S2ID
#===============================================================================

modelV<-"TASKSET_CORSIDE_2"
grpV <- c("TASKSET_2", "CORSIDE_2") #c("F_NUMBER_1","FV_OE_1","FV_LH_1","F_COLOR_1","FV_RB_1","FV_BF_1")
pattern <-ds[,lapply(.SD,unique),.SDcols = grpV,by=modelV] %>% arrange(TASKSET_CORSIDE_2)

# Export
write.table(pattern,paste0(Dir_RSAMODEL,'/DimensionSwitchPRP_DesignSpace_TST_CS_12.txt'), sep="\t",row.names=FALSE)

#===============================================================================
# RSA using submatrix of the design 
# - full design space: 4 rules x 4 color feat x 4 number feat = 64
# - submatrix design space: 4 rule x 2 composite feature x 2 response = 16
#===============================================================================

# # Bivariant stimuli sets (4 four each feature dimensions)
# # Task sets 1, 2 (Numbers): High vs. Low / Odd vs. Even 
# # Task sets 3, 4 (Colors): Red vs. Blue / Faint vs. Bold 
# stim_f1 <- c(1:4) # Numbers: 3,4,6,7
# stim_f2 <- c(1:4) # Colors: B,Bf,Rf,R
# task_sets <- c(1:4) 
# stimspace <- expand.grid(stim_f1, stim_f2)
# sets <- expand.grid(stim_f1, stim_f2, task_sets)
# sets <- sets[, ncol(sets):1]
# names(sets) <- c("TASK", "STIM_F1", "STIM_F2")
# 
# # PRP task requires: two consecutive tasks (F1 - F2 or F2 - F1)
# sets <- data.table(sets)
# sets[,`:=` (RESP=NA,　EVENT = ifelse(TASK==1|TASK==2, "F1", "F2"))]
# sets[, RESP:=ifelse(STIM_F1==1 | STIM_F1==2, 1, 2)]
# sets[EVENT == "F1" & TASK==1, RESP:=ifelse(STIM_F1==1 | STIM_F1==2, 1, 2)]#Event_F1(Number) Task_1(HighLow)
# sets[EVENT == "F1" & TASK==2, RESP:=ifelse(STIM_F1==2 | STIM_F1==3, 1, 2)]#Event_F1(Number) Task_2(OddEven)
# sets[EVENT == "F2" & TASK==3, RESP:=ifelse(STIM_F2==1 | STIM_F2==2, 2, 1)]#Event_F2(Color) Task_3(RedBlue)
# sets[EVENT == "F2" & TASK==4, RESP:=ifelse(STIM_F2==2 | STIM_F2==3, 2, 1)]#Event_F2(Color) Task_2(FaintBold)
# 
# # Task x StimF1 x StimF2 x Resp
# rsapairs <- t(combn(4,2))
# pairID <- rsapairs[1,]
# sets_t_s1_s2_resp <- list()
# sets_t_s1_s2_resp[1] <- list(sets[STIM_F1 ==pairID[1] & STIM_F2 == pairID[1]])
# sets_t_s1_s2_resp[2] <- list(sets[STIM_F1 ==pairID[1] & STIM_F2 == pairID[2]])
# sets_t_s1_s2_resp[3] <- list(sets[STIM_F1 ==pairID[2] & STIM_F2 == pairID[1]])
# sets_t_s1_s2_resp[4] <- list(sets[STIM_F1 ==pairID[2] & STIM_F2 == pairID[2]])
# sets_i <- rbindlist(sets_t_s1_s2_resp) %>% print()


# # Generate models and save
# modelV <- "STIMSPACE(small)"
# pattern <- data.table(expand.grid(c(1,2),c(1,2)))
# nclass <- n_distinct(pattern)
# dir.create(paste0(Dir_RSAMODEL,"/",modelV))
# setwd(paste0(Dir_RSAMODEL,"/",modelV))
# 
# for (v in colnames(pattern)){
#   m <- matrix(nrow=nclass,ncol=nclass);diag(m)<-1;
#   for (p in unique(pattern[[v]])){i<-pattern[[v]]==p;m[i,i]<-1/length(which(i))} # m[i,i]<-1/length(which(i)
#   m[is.na(m)]<-0;assign(paste0(v,"_M"),as.table(m));
#   write.table(m, file=paste0(v,"_M.txt"),sep="\t",row.names=F,col.names=F)
# }
# 
# 
# # Save design?
# ds_p <- as.data.table(pattern)
# ds_p[,COND:=.I]
# setcolorder(ds_p, c("COND", names(ds_p)[1:(length(ds_p)-1)]))
# write.table(ds_p, file=paste0("Task_Design.txt"),sep="\t", row.names=F,col.names=T)

#===============================================================================
#Generate & Save models
#===============================================================================

# Generate models and save
dir.create(paste0(Dir_RSAMODEL,"/",modelV))
setwd(paste0(Dir_RSAMODEL,"/",modelV))
nclass<-n_distinct(ds[[modelV]])

for (v in colnames(pattern)){
  m <- matrix(nrow=nclass,ncol=nclass);diag(m)<-1;
  for (p in unique(pattern[[v]])){i<-pattern[[v]]==p;m[i,i]<-1/length(which(i))}
  m[is.na(m)]<-0;assign(paste0(v,"_M"),as.table(m));
  write.table(m, file=paste0(v,"_M.txt"),sep="\t",row.names=F,col.names=F)
}
