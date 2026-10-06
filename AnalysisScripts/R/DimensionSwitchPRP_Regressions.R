# DimensionSwitchPRP
#Ashwin's script for buidling/testing regressions for decoding accuracy on behavior

#========================================================
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
Dir_REPO <- path.expand("/oscar/home/asrini17/data/asrini17/DimensionSwitch_PRP/PRP_DimensionSwitch/")


# Load libraries
library(data.table)
#library(tidyverse)
library(tidyr)
library(dplyr)
library(lme4)
library(mgcv)
library(car)
library(RSQLite)
library(sjPlot)
library(lemon)
library(knitr)
library(readr)
#library(permuco)

# Source Files
setwd(Dir_R)
source('basic_theme.R')
source('basic_lib.R')

#Load behavioral data-----
setwd(Dir_BDATA)
ds_b<-fread("DimensionSwitchPRP_BehPP.txt") 
ds_b[is.nanM(ds_b)]<-NA;


########-Behavior models
# ds_m1<- ds_b %>%
#   filter(ds_b$RT_2<ds_b$cutRT2 & ds_b$RT_2>0 & NUMTASK == "double")
# 
# ds_m1$RESPCONG <- as.factor(ds_m1$RESPCONG)
# ds_m1$NUMTASK <- as.factor(ds_m1$NUMTASK)
# ds_m1$SOA <- as.factor(ds_m1$SOA)
# ds_m1$SESSION <- as.factor(ds_m1$SESSION)
# contrasts(ds_m1$SESSION, how.many = 2) <- contr.poly(5)
# contrasts(ds_m1$SOA, how.many = 2) <- contr.poly(5)
# 
# #response congruency
# #ds_m1[RT_2==0,RT_2:=1e-10];ds_m1[RT_2==0,RT_2:=1e-10]
# #ds_m1 <- ds_m1[-117382] #will have to fiiix
# m1 <- lmer(log(RT_2) ~SESSION*RESPCONG*SOA + (1|SUBID_S), data = ds_m1)
# summary(m1)
########-

#settings
sub_nums <-c(601:624)
sub_ids <-paste0('A',sub_nums)
modelV <- "TASKSET_CORSIDE_2"
sessN <- "nb_S2"
balanceV <- "balance"
evLock <- "RESP"
#modelIV <- "CORSIDE_1" #pick variable of choice
#modelDV <- "CORSIDE_2" #pick variable of choice

#load in decoding .RDS files
r_pred <- vector("list", length(sub_ids))

#r_predIV<- vector("list", length(sub_ids))
#r_predDV<- vector("list", length(sub_ids))

for (s in sub_ids){
  Dir_Data_i<-paste0(Dir_EDATA,fsep,s,fsep,"DATASETS") 
  setwd(Dir_Data_i);
  files <- list.files(getwd())
  
  if (exists("modelIV") & exists("modelDV")){
    iv_rds_object <- files[grepl(modelIV, files) & grepl("merged", files)]
    
    iv_pred_results <- data.frame(readRDS(iv_rds_object))
    
    item<- match(s,sub_ids)
    r_predIV[[item]] <- iv_pred_results
  
    dv_rds_object <- files[grepl(modelDV, files)]# & grepl("merged", files)]
    
    dv_pred_results <- data.frame(readRDS(dv_rds_object))
    r_predDV[[item]] <- dv_pred_results
    
    rm(iv_pred_results)
    rm(dv_pred_results)
  } else {
    #rds_object <- files[grepl(modelV, files)]# & grepl("merged", files)]
    # rds_object <- files[grepl(modelV, files) & grepl("merged", files) & grepl(sessN, files)]
    rds_object <- paste0(s, "_merged_",modelV,"_",sessN,"_",evLock, "_",balanceV,"_pred.rds")
    pred_results <- data.frame(readRDS(rds_object))
    
    item<- match(s,sub_ids)
    r_pred[[item]] <- pred_results
  }
  
}
#r_pred<- do.call(rbind, r_pred) #combine into one BIIIG df
r_pred<-rbindlist(r_pred)
#setDT(ds_b)

#filter behavioral trials
ds_b<-ds_b[RESPOMIT==0] # Exclude response omission trials?
ds_b<-ds_b[RESPORDER_K==1]# Exclude response flipping trials?
ds_b<-ds_b[!(c(ds_b$RT_1 > ds_b$cutRT1)),]; ds_b<-ds_b[!(c(ds_b$RT_1 < 0)),];
ds_b<-ds_b[!(c(ds_b$RT_2 > ds_b$cutRT2)),];ds_b<-ds_b[!(c(ds_b$RT_2 < 0)),];

for (v in modelV){ #converts model V levels into classifier labels (from preprop_decode() function)
  ALPHABETS<-LETTERS_ex(length(c(unique(ds_b[[v]]))));#comes from other function
  if (is.character(ds_b[[v]])){ds_b[[v]]<-as.numeric(as.factor(ds_b[[v]]))}# convert strings to number
  if (min(ds_b[[v]])==0){ds_b[[v]]<-ds_b[[v]]+1} # adjust numbers starting from 0
#  ds_b[[v]]<-ALPHABETS[as.matrix(ds_b[[v]])];
  ds_b[[v]]<-ds_b[[v]]%>%unlist%>%factor;
  ds_b[[v]]<-droplevels(ds_b[[v]])
  
}

#cols_to_add <- c("NUMTASK", "SESSION", "SOA","RT_2", modelV)

#r_pred <- ds_b[, .SD, .SDcols = c("SUBID_S", "BLOCK", "TRIAL", "NUMTASK", "SESSION", "SOA","RT_1" ,"RT_2", modelV)][r_pred, on = .(SUBID_S, BLOCK, TRIAL)]
r_pred <- r_pred %>% #merge behavioral variables and decoding accuracy
   left_join(ds_b %>% dplyr::select(SUBID_S, BLOCK, TRIAL,NUMTASK, SESSION, SOA, all_of(modelV), RT_2),# RT_2, # pick variable of choice
            by = c("SUBID_S","BLOCK", "TRIAL"))


##### this approach is not great but....
setDT(r_pred)
r_pred[, acc :=NULL]

#setDT(r_pred)
r_pred[, prob := as.matrix(.SD)[cbind(.I, TASKSET_CORSIDE_2)],  #grab probability of decoding relevant conjunction
       .SDcols = ALPHABETS]
r_pred[, (ALPHABETS):= NULL]

#####


#r_pred$prob <- ifelse(r_pred[[modelV]] ==1, r_pred$A, r_pred$B) #keep prob of correct side
#r_pred$prob <- r_pred[[modelV]]
#prob <- r_pred[cbind(seq_len(nrow(r_pred)), r_pred[[modelV]])]  

#r_predDV<- do.call(rbind, r_predDV) #combine into one BIIIG df
#r_predIV<- do.call(rbind, r_predIV) #combine into one BIIIG df

#r_predDV <- r_predDV %>% #merge behavioral variables and decoding accuracy
#  left_join(ds_b %>% dplyr::select( SUBID_S,BLOCK, TRIAL, NUMTASK, SESSION, ACC_2, RT_2, all_of(modelDV)),# pick variable of choice
#            by = c("SUBID_S","BLOCK", "TRIAL"))
#probIdx <- r_predDV[[modelDV]]+1
#prob <- mapply(function(i, j) r_predDV[i, j], 
#               i = seq_len(nrow(r_predDV)), 
#               j = probIdx)
#r_predDV$prob <- r_predGRAND[cbind(seq_len(nrow(r_predGRAND)), r_predGRAND$modelV + 1)]
#r_predGRAND$prob <- r_predGRAND[,(r_predGRAND[[modelV]])+1] # grab the column that corresponds to value in modelV column + 1, since acc is always first
#r_predDV$prob <- ifelse(r_predDV[[modelDV]] ==1, r_predDV$A, r_predDV$B) #keep prob of correct side

#r_predIV <- r_predIV %>% #merge behavioral variables and decoding accuracy
#  left_join(ds_b %>% dplyr::select( SUBID_S,BLOCK, TRIAL, NUMTASK, SESSION, ACC_1, RT_1, all_of(modelIV)),# pick variable of choice
#            by = c("SUBID_S","BLOCK", "TRIAL"))
#r_predIV$prob <- ifelse(r_predIV[[modelIV]] ==1, r_predIV$A, r_predIV$B) #keep prob of correct side


#cut data to time of interest and avg:
  # - mean RT_1 (sessions 2,3,6) is 777.6 --> time (index) range is 50-244
  # - mean RT_2 (sessions 2,3,6) is 830.7 --> time (index) range is 125-332 (but time index only goes to 301)
  # for response lock data (evLock == "RESP"), both time intervals are -700 to 200 relative to response
  # so for resp locked data cut time range to 175 (no sense in using decoding after response)
r_predTrial <- r_pred %>%
    filter(time <=175 & time >=125) %>%
    mutate(time_ms = (time*4)-700) %>% #convert back to MS 
    group_by(SUBID_S,BLOCK, TRIAL, SESSION, NUMTASK, SUBID, SOA) %>% #average across timepoints
    summarise(mean_probDV = mean(prob),
#             task_acc = mean(ACC_2),
              task_rt = mean(RT_2),
              log_rt = log(task_rt)) %>% setDT()
r_predTrial$SESSION <- as.factor(r_predTrial$SESSION)
r_predTrial$SUBID <- as.factor(r_predTrial$SUBID)
r_predTrial$SUBID_S <- as.factor(r_predTrial$SUBID_S)
r_predTrial$SOA <- as.factor(r_predTrial$SOA)
r_predTrial$NUMTASK <- factor(r_predTrial$NUMTASK, levels = c("single","double"))
contrasts(r_predTrial$SESSION, how.many = 1) <- contr.poly(3)
contrasts(r_predTrial$SOA, how.many = 1) <-  contr.poly(5)


# same as above but avg for each time point
r_predTime <- r_pred %>%
  filter(time <=175 & time >=125) %>%
  mutate(time_ms = (time*4)-700) %>% #convert back to MS 
  group_by(SUBID_S,time_ms, SESSION, NUMTASK, SUBID, SOA) %>% #average across timepoints
  summarise(mean_probDV = mean(prob),
            #             task_acc = mean(ACC_2),
            task_rt = mean(RT_2),
            log_rt = log(task_rt)) %>% setDT()
r_predTime$SESSION <- as.factor(r_predTime$SESSION)
r_predTime$SUBID <- as.factor(r_predTime$SUBID)
r_predTime$SUBID_S <- as.factor(r_predTime$SUBID_S)
r_predTime$SOA <- as.factor(r_predTime$SOA)
r_predTime$NUMTASK <- factor(r_predTime$NUMTASK, levels = c("single","double"))

  
#r_predDT <- r_predDV %>%
#  filter(time <=225 & time >125) %>%
#  mutate(time_ms = time*4) %>% #convert back to MS 
#  group_by(SUBID_S,BLOCK, TRIAL, SESSION, NUMTASK, SUBID) %>% #average across timepoints
#  summarise(mean_probDV = mean(prob),
#            task_acc = mean(ACC_2),
#            task_rt = mean(RT_2),
#            log_rt = log(task_rt)) 
#r_predIT <- r_predIV %>%
#  filter(time <=150 & time >50) %>% #for S1
#  mutate(time_ms = time*4) %>% #convert back to MS 
#  group_by(SUBID_S,BLOCK, TRIAL, SESSION, NUMTASK, SUBID) %>% #average across timepoints
#  summarise(mean_probIV = mean(prob),
#            task_acc = mean(ACC_1),
#            task_rt = mean(RT_1),
#            log_rt = log(task_rt)) 
#r_predGRANDT <- r_predDT %>% #merge behavioral variables and decoding accuracy
#  left_join(r_predIT %>% dplyr::select(SUBID_S,BLOCK, TRIAL,mean_probIV),# pick variable of choice
#           by = c("SUBID_S","BLOCK", "TRIAL", "SESSION","NUMTASK"))


#Vizualizing Relationships --------------------------------
theme_set(theme_bw(base_size = 20))#32/28
CSCALE_YlOrRd = rev(brewer.pal(9,"YlOrRd"));

#Plot 
r_predTime[SOA==500] %>%
  group_by(SESSION, NUMTASK, time_ms) %>% #average across timepoints
  summarise(prob = mean(mean_probDV),
            mean_rt = mean(task_rt),
            sem= sd(mean_probDV)/sqrt(n()),) %>% 
  ggplot(aes(x = time_ms, y = prob, group = NUMTASK)) +
#  geom_point(shape = 1, alpha = 0.2) +
  geom_ribbon(aes(ymin = prob-sem, ymax = prob+sem,fill = NUMTASK), alpha = 0.2) +
  geom_point(aes(color=NUMTASK), shape = 1, alpha = 0.4) +
  geom_line(aes(color = NUMTASK)) +
  facet_grid(~SESSION)

r_predTrial[SOA==250] %>%
  group_by(SESSION, NUMTASK) %>% #average across timepoints
  summarise(prob = mean(mean_probDV),
            mean_rt = mean(task_rt)) %>%
  ggplot(aes(x = SESSION, y = prob, group = NUMTASK)) +
  #  geom_point(shape = 1, alpha = 0.2) +
  geom_jitter(aes(color=NUMTASK), shape = 1, alpha = 0.4) + 
  geom_smooth(method = "lm", aes(color = NUMTASK))


# Models --------------------------------------------------
contrasts(r_predTime$SESSION, how.many = 2) <- contr.poly(3)
contrasts(r_predTime$SOA, how.many = 2) <- contr.poly(5)

# Model 0.5 - looking at decoding by numtask over time windows
#m5 <- lmer(formula = mean_probDV~SESSION*NUMTASK*SOA +(1|SUBID_S), data = r_predT, control=lmerControl(optimizer="bobyqa"))

# m5 <- aovperm( this will crash! takes a HUUUGe amnt of memory
#   formula = mean_probDV ~ NUMTASK*SESSION*time_ms + Error(SUBID/time_ms),
#   data = r_predTime,
#   np = 5000,
#   cluster = "time"
# )

# Model 1 - predicting decoding by session, numtask
setDT(r_predTrial)
r_predDouble <- r_predTrial[NUMTASK == 'double']
m1 <- lmer(formula = mean_probDV~SESSION*SOA +(1|SUBID_S), data = r_predDouble)
summary(m1)
plot(allEffects(m1), main = deparse(formula(m1)))
plot_model(m1, type = "int")
plot_model(m1, type="std", show.p=T, show.values=T, value.size=2,value.offset=.25) + ylim(c(-.05,.05))+scale_y_continuous(labels=dropLeadingZero) + ggtitle(paste(modelV, "decoding in", sessN))

#Model 1.5 - predicting decoding of CORSIDE_2 by decoding of CORSIDE_1
m15 <- lmer(formula = mean_probDV~SESSION*NUMTASK*mean_probIV +(1|SUBID_S), data = r_predGRANDT, control=lmerControl(optimizer="bobyqa"))
summary(m15)
plot_model(m15, type = "int")

# Model 2 - predicting RT by decoding 
#ds_m <- r_predT %>% filter(NUMTASK == "double")
#m <- lmer(formula = log_rt~SOA*RESPCONG*SESSION*mean_probDV + (1 |SUBID_S), data = ds_m, control=lmerControl(optimizer="bobyqa")) #just behavior
summary(m)

m_rt=lmer(log_rt~1+SESSION*SOA*NUMTASK*mean_probDV+(1|SUBID_S),data=r_predTrial);#std_beta_MLM(m_rt)@fixef
summary(m_rt)
plot_model(m_rt, type="std", show.p=T, show.values=T, value.size=4,value.offset=.25, title = "2nd Item RT") + ylim(c(-.05,.05))+scale_y_continuous(labels=dropLeadingZero)

#adopting the approach from behavior models -- 1 SOA at a time
setDT(r_predTrial)
ds_100 <- r_predTrial[SOA == 100]
ds_250 <- r_predTrial[SOA == 250]
ds_500 <- r_predTrial[SOA == 500]
ds_750 <- r_predTrial[SOA == 750]
ds_1000 <- r_predTrial[SOA == 1000]


#predicting change in decodability
m500<- lmer(formula = mean_probDV~SESSION*NUMTASK +(1|SUBID_S), data = ds_500)
plot(allEffects(m500), main = deparse(formula(m500)))
summary(m500)
setwd(paste0(Dir_REPO,"figures/TasksetXCorsideRSA/RESP-lock/decodingStats"))
#for saving summary output as .txt file
fileN <- paste0(modelV,"_bySOA_simple_summaries.txt")
sink(fileN)
print(summary(m1000))
sink()


#rt models
m100=lmer(log_rt~1+SESSION*NUMTASK*mean_probDV+(1|SUBID_S),data=ds_100);#std_beta_MLM(m_rt)@fixef
summary(m100)
plot(allEffects(m100), main = deparse(formula(m100)))
m250=lmer(log_rt~1+SESSION*NUMTASK*mean_probDV+(1|SUBID_S),data=ds_250);#std_beta_MLM(m_rt)@fixef
summary(m250)
plot(allEffects(m250), main = deparse(formula(m250)))
m500=lmer(log_rt~1+SESSION*NUMTASK*mean_probDV+(1|SUBID_S),data=ds_500);#std_beta_MLM(m_rt)@fixef
summary(m500)
plot(allEffects(m500), main = deparse(formula(m500)))
m750=lmer(log_rt~1+SESSION*NUMTASK*mean_probDV+(1|SUBID_S),data=ds_750);#std_beta_MLM(m_rt)@fixef
summary(m750)
plot(allEffects(m750), main = deparse(formula(m750)))
m1000=lmer(log_rt~1+SESSION*NUMTASK*mean_probDV+(1|SUBID_S),data=ds_1000);#std_beta_MLM(m_rt)@fixef
summary(m1000)
plot(allEffects(m1000), main = deparse(formula(m1000)))
ggplot(data = r_predTrial[SOA %in% c("100","500","1000")], aes(x=mean_probDV, y = log_rt))+
  geom_smooth(method = "glm",aes(color = NUMTASK, fill = NUMTASK)) +
  facet_grid(SOA~SESSION)+ ggtitle("Conjunction decoding and RT_2")+
  scale_color_manual(values = c("cornflowerblue", "orange")) + scale_fill_manual(values = c("cornflowerblue", "orange"))


setwd(paste0(Dir_REPO,"Figures&Stats/BEH/BEH(6sessions)/plots&models/"))
#for saving summary output as .txt file
fileN <- paste0(glue("RT_2_{modelV}_bySOA_simple_summaries.txt")) #"VI" for variable interaction (m1); "FI" for fixed interaction (m2); "simple" for no random slopes
sink(fileN)
print(summary(m1000))
sink()
saveRDS(m750, file = "RT_2_CONJ_SOA_750_simple_model.rds")

#Model 2.5 - predicting RT_2 by decoding of CORSIDE_2 and CORSIDE_1
m25 <- lmer(formula = log_rt~SESSION*NUMTASK*mean_probDV*mean_probIV + (1+SESSION|SUBID), data = r_predGRANDT, control=lmerControl(optimizer="bobyqa"))
summary(m25)
plot_model(m25, type = "pred", terms = c("SESSION", "NUMTASK", "mean_prob", "mean_probIV")) 
plot_model(m25, type = "int")


#=For RSA regressions===========================================================
setwd(Dir_GRAND)

fileN<- "DimensionSwitchPRP"
modelVs <- "TASKSET_CORSIDE_2"
sessN <- "2"
regs <- "BIS" #RTACC
value <- "betas"
test<- "byTrial" #timeSERIES
if(!is_empty(regs)){dName <- paste(fileN,modelVs,regs,value,test, sep="_")} else {dName <- paste(fileN,modelVs,value,test, sep="_")}

intV <- c("SUBID", "SESSION","NUMTASK","SOA")
dVs <- c("TASKSET_RSA", "CORSIDE_RSA", "CONJ_RSA") # list of values that will be DV

dsG <- readRDS(paste0(dName, ".rds")) #load in data
dsG[,SOA := as.factor(SOA)]
dsG[,NUMTASK := factor(NUMTASK, levels = c("single","double"))]
dsG[,SESSION := as.factor(SESSION)]
contrasts(dsG$SESSION, how.many = 2) <- contr.poly(3)
contrasts(dsG$SOA, how.many = 2) <- contr.poly(5)

if(any(grepl("early", names(dsG)))){dVs <- c(paste0("early",dVs), paste0("late",dVs))}
models <- vector("list", length(dVs)) #for early & late!
for (i in 1:length(dVs)){
  
  dv <- dVs[i]
  cols<- c(intV, dv)
  dsM <- dsG[, ..cols] #slice the big df 
  
  
  models[i] <- lmer(formula = get(dv) ~ SESSION*NUMTASK*SOA + (1 + SESSION|SUBID), data = dsM, control= lmerControl(optimizer= "bobyqa"))
  names(models)[i] <- paste0(dv, "_model")
  rm(dsM) #get rid of sliced df to reduce memory use
}

list2env(models, envir = .GlobalEnv);modelNs <- names(models); rm(models);

for (model in modelNs){
  if(!is_empty(regs)){fileN <- paste0(model,"_", regs,"_" ,value, "_output.txt")} else {fileN <- paste0(model,"_" ,value, "_output.txt")}
  #fileN <- paste0(model,"_", regs,"_" ,value, "_output.txt")
  sink(fileN)
  print(summary(get(model)))
  sink()
  
  #print(tab_model(get(model), file = paste0(fileN, ".doc")))
  #kable_out <- kable(sum, 'html')
  #write_file(sum, paste0(model,"_", regs,"_" ,value, "_output.html")[1])
  #webshot(paste0(file, ".html"), paste0(file, "png"))
}


# Regression for S1ID RSA predicting S2ID RSA ----------------------------------
dsS2ID <- fread("/oscar/data/dbadre/asrini17/DimensionSwitch_PRP/DataAnalysis/Matlab/EEG/w_ALLGRAND/DimensionSwitchPRP_S2ID_SOA_byTrial.txt") 
  #readRDS("/oscar/data/dbadre/asrini17/DimensionSwitch_PRP/DataAnalysis/Matlab/EEG/w_ALLGRAND/DimensionSwitchPRP_S2ID_test_timeSERIES.rds")
dsS1ID <- readRDS("/oscar/data/dbadre/asrini17/DimensionSwitch_PRP/DataAnalysis/Matlab/EEG/w_ALLGRAND/DimensionSwitchPRP_S1ID_test_timeSERIES.rds")

#code for relevant features
dsS1ID$relFeature <- case_when(
  dsS1ID$TASKSET_1 == "Low-High" ~ dsS1ID$LH_RSA,
  dsS1ID$TASKSET_1 == "Odd-Even" ~ dsS1ID$OE_RSA,
  dsS1ID$TASKSET_1 == "Red-Blue" ~ dsS1ID$RB_RSA,
  dsS1ID$TASKSET_1 == "Bold-Faint" ~ dsS1ID$BF_RSA)
#Irrelevant features
# Low-High
dsS1ID[TASKSET_1 == "Low-High", 
       irrelFeature := rowMeans(.SD, na.rm = TRUE), 
       .SDcols = c("OE_RSA", "RB_RSA", "BF_RSA")]
# Odd-Even
dsS1ID[TASKSET_1 == "Odd-Even", 
       irrelFeature := rowMeans(.SD, na.rm = TRUE), 
       .SDcols = c("LH_RSA", "RB_RSA", "BF_RSA")]
# Red-Blue
dsS1ID[TASKSET_1 == "Red-Blue", 
       irrelFeature := rowMeans(.SD, na.rm = TRUE), 
       .SDcols = c("OE_RSA", "LH_RSA", "BF_RSA")]
# Bold-Faint
dsS1ID[TASKSET_1 == "Bold-Faint", 
       irrelFeature := rowMeans(.SD, na.rm = TRUE), 
       .SDcols = c("OE_RSA", "RB_RSA", "LH_RSA")]

dsS2ID$relFeature <- case_when(
  dsS2ID$TASKSET_2 == "Low-High" ~ dsS2ID$LH_RSA,
  dsS2ID$TASKSET_2 == "Odd-Even" ~ dsS2ID$OE_RSA,
  dsS2ID$TASKSET_2 == "Red-Blue" ~ dsS2ID$RB_RSA,
  dsS2ID$TASKSET_2 == "Bold-Faint" ~ dsS2ID$BF_RSA)
#Irrelevant features
# Low-High
dsS2ID[TASKSET_2 == "Low-High", 
       irrelFeature := rowMeans(.SD, na.rm = TRUE), 
       .SDcols = c("OE_RSA", "RB_RSA", "BF_RSA")]
# Odd-Even
dsS2ID[TASKSET_2 == "Odd-Even", 
       irrelFeature := rowMeans(.SD, na.rm = TRUE), 
       .SDcols = c("LH_RSA", "RB_RSA", "BF_RSA")]
# Red-Blue
dsS2ID[TASKSET_2 == "Red-Blue", 
       irrelFeature := rowMeans(.SD, na.rm = TRUE), 
       .SDcols = c("OE_RSA", "LH_RSA", "BF_RSA")]
# Bold-Faint
dsS2ID[TASKSET_2 == "Bold-Faint", 
       irrelFeature := rowMeans(.SD, na.rm = TRUE), 
       .SDcols = c("OE_RSA", "RB_RSA", "LH_RSA")]

dsS1ID[, diff := (relFeature-irrelFeature)];dsS2ID[, diff := (relFeature-irrelFeature)]
dsS1ID[ds_b, SOA := i.SOA, on = .(SUBID, SESSION, NUMTASK, BLOCK, TRIAL)]; dsS2ID[ds_b, SOA := i.SOA, on = .(SUBID, SESSION, NUMTASK, BLOCK, TRIAL)]# add SOA info

r_predTrial <- dsS2ID %>%
  filter(time <=790 & time >=0) %>%
  group_by(SUBID, SESSION, BLOCK, TRIAL, NUMTASK, SOA) %>% #average across timepoints
  summarise(mean_probRel = mean(relFeature))
#             task_acc = mean(ACC_2),
#              task_rt = mean(RT_1),
#              log_rt = log(task_rt))
r_predTrial$SESSION <- as.factor(r_predTrial$SESSION)
r_predTrial$SUBID <- as.factor(r_predTrial$SUBID)
r_predTrial$SOA <- as.factor(r_predTrial$SOA)
r_predTrial$NUMTASK <- factor(r_predTrial$NUMTASK, levels = c("single","double"))
contrasts(r_predTrial$SESSION, how.many = 1) <- contr.poly(3)
contrasts(r_predTrial$SOA, how.many = 1) <-  contr.poly(5)

#quick plot if need be----
r_predTime <- dsS2ID %>%
  filter(time <=790 & time >=0) %>%
  group_by(SESSION, time, NUMTASK, SOA) %>% #average across timepoints
  summarise(mean_probRel = mean(relFeature))

ggplot(data = r_predTime, aes(x=time,y=mean_probRel))+ #,group=rel,color=rel)) +
  #scale_y_continuous(breaks = scales::pretty_breaks(n=2))+
  geom_line(aes(color=NUMTASK)) + 
  #geom_point(aes(group=rel))+
  #coord_cartesian(ylim=c(-0.05,0.1)) +#  c(-0.05,2.3) + #stat_summary(size = 0.5, shape = 1)+ 
#  geom_ribbon(alpha=0.2,linetype=0, aes(fill = rel, ymin = mean_RSA - sem_RSA, ymax = mean_RSA + sem_RSA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  facet_grid(SOA~SESSION)# +
#  scale_color_manual(values=c("cornflowerblue","orange", "darkolivegreen4")) + scale_fill_manual(values=c("cornflowerblue","orange", "darkolivegreen4"))

#----
m1 <- lmer(formula = mean_probRel~SESSION*NUMTASK*SOA +(1+ SESSION|SUBID), data = r_predTrial, control=lmerControl(optimizer="bobyqa"))
summary(m1)
plot_model(m1, type = "int")
plot_model(m1, type="std", show.p=T, show.values=T, value.size=2,value.offset=.25) +scale_y_continuous(labels=dropLeadingZero) + ggtitle(paste("Relevant Feature RSA in S2"))

