# DimensionSwitchPRP
#Aggregates single-session/participant data to plot group-level decoding

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
Dir_REPO <- path.expand("/oscar/home/asrini17/data/asrini17/DimensionSwitch_PRP/PRP_DimensionSwitch/figures")

# Load libraries
library(data.table)
#library(tidyverse)
library(dplyr)
library(tidyr)
library(broom)
library(RColorBrewer)
library(scales)
# library(forcats)
library(stringr)
library(plotrix)
library(lemon)
library(glue)

# Source Files
setwd(Dir_R)
source('basic_theme.R')
source('basic_lib.R')

#Load behavioral data
setwd(Dir_BDATA)
ds_b<-fread("DimensionSwitchPRP_BehPP.txt") 
#ds_b<-fread("DimensionSwitchPRP_BehPP_LONG.txt") #for cross-condition decoding!
ds_b[is.nanM(ds_b)]<-NA;

# Response congruency effects in behavior
# ds_b$RESPCONG <- as.factor(ds_b$RESPCONG)
# ds_b$SESSION <- as.factor(ds_b$SESSION)
# ds_b$SOA <- as.factor(ds_b$SOA)
# ds_b %>%
#   filter(ds_b$RT_2<ds_b$cutRT2 & ds_b$RT_2 >0 & NUMTASK == "double") %>%
#   group_by(SOA, SESSION,RESPCONG ) %>% 
#   summarise(rt = mean(log(RT_2)),
#             sem = sd(log(RT_2))/sqrt(sum(!is.na(log(RT_2)))))%>%
# ggplot( aes(x = SOA, y = rt)) + geom_point(size = 0.5,aes(color = RESPCONG)) + facet_wrap(~SESSION) + geom_linerange(aes(ymin = rt-sem, ymax = rt+sem, color = RESPCONG))
# 
#

#settings
#sub_nums <-c(601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624)
#sub_nums <- c(201:224,601:624)
#sub_nums <- c(601:618,620:624)
sub_nums <- c(601:624)
sub_ids <-paste0('A',sub_nums)
modelV <- "TASKSET_1" #pick variable of choice
sessN <- "nb_S1"#"nppr_S2", "nb_CROSS"
crossV <- "STIMNUM" #"STIMNUM"
balanceV <- "NUMTASK" #"_NUMTASK"
evLock <- "STIM" # "STIM"
decodingType <- "merged" #merged

#load in decoding .RDS files
r_predGRAND <- vector("list", length(sub_ids))

for (s in sub_ids){
  
  
  if(decodingType == 'merged') {
    Dir_Data_i<-paste0(Dir_EDATA,fsep,s,fsep,"DATASETS") 
    setwd(Dir_Data_i);
    files <- list.files(getwd())
  #if (){} fill this in later
  rds_object <- files[grepl(modelV, files) & grepl("merged", files) & grepl(sessN, files)][1] #grabs merged RDS specified by modelV and merged across sessions
  #rds_object <- paste0(s, "_merged_",modelV,"_",sessN,"_",evLock,"_",balanceV,"_pred.rds") #_TAKE2
  #rds_object <- files[grepl(modelV, files) & grepl(crossV, files) & grepl(sessN, files) & grepl("pred",files)]
  pred_results <- data.frame(readRDS(rds_object))
  }
  if(decodingType == 'unmerged') { #my plan here is to load classifiers data independently and then merge
    #load session 2 data
    ses2ID <- sub("[0-9]", "2", s)
    
    Dir_Data_s2<-paste0(Dir_EDATA,fsep,ses2ID,fsep,"DATASETS")
    Dir_Data_s6<-paste0(Dir_EDATA,fsep,s,fsep,"DATASETS") 
    setwd(Dir_Data_s2);
    ses2 <- data.frame(readRDS(paste0(ses2ID, "_unmerged_",modelV,"_",sessN,"_",evLock, "_",balanceV,"_pred.rds")))
    setwd(Dir_Data_s6);
    ses6 <- data.frame(readRDS(paste0(s, "_unmerged_",modelV,"_",sessN,"_",evLock, "_",balanceV,"_pred.rds")))
    pred_results <- rbind(ses2, ses6)
  }
  
  item<- match(s,sub_ids)
  r_predGRAND[[item]] <- pred_results
  
}

r_predGRAND<- rbindlist(r_predGRAND) #combine into one BIIIG df
if (!("SUBID_S" %in% names(r_predGRAND))) {
setnames(r_predGRAND, "SUBID", "SUBID_S")
}

if (sessN == "nb_CROSS"){
 r_predGRAND[, SUBID_S := as.integer(case_when(
   SESSION == 2 ~ sub("[0-9]", "2", SUBID_S),
   SESSION == 3 ~ sub("[0-9]", "3", SUBID_S),
   SESSION == 6 ~ sub("[0-9]", "6", SUBID_S),
 ))]
}
#faster?
#r_predGRAND[SESSION %in% c(2, 3, 6), 
# SUBID_S := sub("[0-9]", as.character(SESSION), SUBID_S)]

#filter behavioral trials
ds_b<-ds_b[RESPOMIT==0] # Exclude response omission trials?
ds_b<-ds_b[RESPORDER_K==1]# Exclude response flipping trials?
if(sessN == "nb_CROSS"){
  ds_b<-ds_b[!(c(ds_b$RT > ds_b$cutRT)),]
} else{  ds_b<-ds_b[!(c(ds_b$RT_1 > ds_b$cutRT1)),]; ds_b<-ds_b[!(c(ds_b$RT_1 < 0)),];
ds_b<-ds_b[!(c(ds_b$RT_2 > ds_b$cutRT2)),];ds_b<-ds_b[!(c(ds_b$RT_2 < 0)),];}

#relevel variables
ds_b[,NUMTASK := factor(NUMTASK, levels = c("single", "double"))]

# for looking at how decoding relates to within-subject variability in RT

subidRTs <- ds_b %>%
  filter(SESSION %in% c(2,3,6)) %>%
  group_by(SUBID, SESSION, NUMTASK) %>%
  summarise(median_RT_1 = median(RT_1),
            median_RT_2 = median(RT_2))
ds_b <- ds_b %>%
  left_join(subidRTs %>% dplyr::select(SUBID,SESSION, NUMTASK, median_RT_1, median_RT_2), #STIMNUM, pick variable of choice
            by = c("SUBID", "SESSION", "NUMTASK")) #,"STIMNUM"
ds_b$RT_1_s <- ifelse(ds_b$RT_1 < ds_b$median_RT_1, "fast", "slow")
ds_b$RT_2_s <- ifelse(ds_b$RT_2 < ds_b$median_RT_2, "fast", "slow")
####


if (sessN == "nb_CROSS"){
r_predGRAND <- r_predGRAND %>%
  left_join(ds_b %>% dplyr::select(SUBID_S, BLOCK,TRIAL, NUMTASK, SESSION, SOA, STIMNUM), #STIMNUM, pick variable of choice
            by = c("SUBID_S","SESSION","BLOCK", "TRIAL","STIMNUM", "NUMTASK", "SESSION")) #,"STIMNUM"
} else{
  r_predGRAND <- r_predGRAND %>%
    left_join(ds_b %>% dplyr::select(SUBID_S, BLOCK,TRIAL, SESSION, NUMTASK, SOA, RT_1_s, RT_2_s), #STIMNUM, pick variable of choice
              by = c("SUBID_S", "BLOCK", "TRIAL")) #,"STIMNUM"
}



r_predG <- r_predGRAND %>%
  group_by(NUMTASK,time) %>%
  summarise(mean_acc = mean(acc),
            sd_acc = sd(acc),
            sem = sd(acc)/sqrt(sum(!is.na(acc))))
# Plotting
theme_set(theme_classic(base_size = 20))#32/28
CSCALE_YlOrRd = rev(brewer.pal(9,"YlOrRd"));
chance<-1/length(unique(ds_b[[modelV]]))


p2<-ggplot(data=r_predG,aes(x=time,y=mean_acc, group = NUMTASK)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(126),linetype=1,linewidth=1)+annotate("text",x=150, y=chance-0.05,label="Stimulus + Cue")+
  geom_hline(yintercept=chance,linetype=1,linewidth=1)+
  geom_ribbon(alpha=0.2,linetype=0, aes(fill=NUMTASK, ymin = mean_acc - sem, ymax = mean_acc + sem))+
  ggtitle(modelV)+
  geom_line(linewidth=1.5,aes(color=NUMTASK, linetype = NUMTASK))+
  labs(x = "Time", y = "Decoding Accuracy")
p2

#broken down by session
r_predSESSION <- r_predGRAND %>%
  group_by(SESSION, NUMTASK, time) %>% #WITHIN
  summarise(mean_acc = mean(acc),
            sd_acc = sd(acc),
            sem = sd(acc)/sqrt(sum(!is.na(acc))))
r_predSESSION$time_ms <- r_predSESSION$time*4-200
p3<-ggplot(data=r_predSESSION,aes(x=time_ms,y=mean_acc, group = NUMTASK)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(0),linetype=1,size=1)+
  geom_hline(yintercept=chance,linetype=1,size=1)+
  geom_ribbon(alpha=0.2,linetype = 0, aes(fill=NUMTASK, ymin = mean_acc - sem, ymax = mean_acc + sem)) +
  ggtitle(paste(modelV,balanceV, sep = " bal "))+
  geom_line(size=1.5,aes(color=NUMTASK)) + 
  #ylim(c(0.247,0.3))+
  facet_wrap(~SESSION) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab(glue("Time relative to {evLock} {sessN}")) + ylab("Decoding Accuracy") +
  scale_color_manual(values = c("cornflowerblue", "orange")) + scale_fill_manual(values = c("cornflowerblue", "orange"))
p3


#broken down by session,SOA
r_predSOA <- r_predGRAND[SOA %in% c("100","500","1000")] %>% #for ease of visualizing, let's see what Atsushi thinks
  group_by(SESSION, SOA, NUMTASK, time) %>% #NUMTASK,
  summarise(mean_acc = mean(acc),
            sd_acc = sd(acc),
            sem = sd(acc)/sqrt(sum(!is.na(acc))))
r_predSOA$time_ms <- r_predSOA$time*4-700
r_predSOA$SOA<-factor(r_predSOA$SOA, levels = c("100","500","1000")) #"250","750"
p4<-ggplot(data=r_predSOA,aes(x=time_ms,y=mean_acc, group = NUMTASK)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(0),linetype=1,size=1)+#annotate("text",x=150, y=chance-0.05,label="Stimulus + Cue")+
  geom_hline(yintercept=chance,linetype=1,size=1)+
  geom_ribbon(alpha=0.2,linetype=0, aes(fill = NUMTASK, ymin = mean_acc - sem, ymax = mean_acc + sem)) +
  ggtitle(paste0(modelV, " bal ", balanceV, " dual-task"))+
  geom_line(size=1.5,aes(color=NUMTASK, group = NUMTASK)) + 
  #ylim(c(0.247,0.3))+
  facet_rep_grid(SOA~SESSION) +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
                                   legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  scale_color_manual(values = c("cornflowerblue", "orange")) + scale_fill_manual(values = c("cornflowerblue", "orange"))+
xlab("Time relative to Response 2") + ylab("Decoding Accuracy")
p4

p3 + p4

#looking at decoding and slow/fast RT
r_predSF <- r_predGRAND %>%
  group_by(SESSION, NUMTASK, time, RT_2_s) %>% #WITHIN
  summarise(mean_acc = mean(acc),
            sd_acc = sd(acc),
            sem = sd(acc)/sqrt(sum(!is.na(acc)))) %>%setDT
r_predSF$time_ms <- r_predSF$time*4-700
r_predSF[,taskXspeed:= paste(NUMTASK, RT_2_s, sep = "_")]
pp3<-ggplot(data=r_predSF,aes(x=time_ms,y=mean_acc, group = taskXspeed)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(0),linetype=1,size=1)+
  geom_hline(yintercept=chance,linetype=1,size=1)+
  geom_ribbon(alpha=0.2,linetype = 0, aes(fill=taskXspeed, ymin = mean_acc - sem, ymax = mean_acc + sem)) +
  ggtitle(paste(modelV,balanceV, sep = " bal "))+
  geom_line(size=1.5,aes(color=taskXspeed)) + 
  #ylim(c(0.247,0.3))+
  facet_grid(~SESSION) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab(glue("Time relative to {evLock} {sessN}")) + ylab("Decoding Accuracy") +
  scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired")
pp3


#save to Repo (but need to commit/push directly after!!)
outfile <- file.path(paste0(Dir_REPO), paste0(modelV,"_merged_bal_", balanceV,"_",evLock, "_combined_SOA_plot.png"))
print(outfile)
png(outfile, width = 2220, height = 1080, units = "px")
print(p3 + p4)
dev.off()


#RESPCONG
r_predGRAND[,RESPCONG := as.factor(RESPCONG)]
r_predGRAND[,NUMSW := as.factor(NUMSW)]
r_predCONG <- r_predGRAND %>%
  group_by( NUMTASK, RESPCONG, time, SOA) %>% #WITHIN, SESSION, 
  summarise(mean_acc = mean(acc),
            sd_acc = sd(acc),
            sem = sd(acc)/sqrt(sum(!is.na(acc)))) %>% setDT()
r_predCONG$time_ms <- r_predCONG$time*4-700
pCong<-ggplot(data=r_predCONG[NUMTASK == 'double'],aes(x=time_ms,y=mean_acc, group = RESPCONG)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(0),linetype=1,size=1)+
  geom_hline(yintercept=chance,linetype=1,size=1)+
  geom_ribbon(alpha=0.2,linetype = 0, aes(fill=RESPCONG, ymin = mean_acc - sem, ymax = mean_acc + sem)) +
  # ggtitle("single-task")+
  geom_line(size=1.5,aes(color=RESPCONG)) + 
  #ylim(c(0.247,0.3))+
  facet_grid(~SOA) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab("Time relative to Response 2 ") + ylab("Decoding Accuracy") 
  #scale_color_manual(values = c("cornflowerblue", "orange")) + scale_fill_manual(values = c("cornflowerblue", "orange"))

pCong + labs(caption = paste0(modelV, " bal ", balanceV, " plotted by response congruency"))

#feature switches...
r_predSW <- r_predGRAND %>%
  group_by(SESSION, NUMTASK, NUMSW, time) %>% #WITHIN, SESSION, 
  summarise(mean_acc = mean(acc),
            sd_acc = sd(acc),
            sem = sd(acc)/sqrt(sum(!is.na(acc)))) %>% setDT()
r_predSW$time_ms <- r_predSW$time*4-700
pSW<-ggplot(data=r_predSW,aes(x=time_ms,y=mean_acc, group = NUMSW)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(0),linetype=1,size=1)+
  geom_hline(yintercept=chance,linetype=1,size=1)+
  geom_ribbon(alpha=0.2,linetype = 0, aes(fill=NUMSW, ymin = mean_acc - sem, ymax = mean_acc + sem)) +
  # ggtitle("single-task")+
  geom_line(size=1.5,aes(color=NUMSW)) + 
  #ylim(c(0.247,0.3))+
  facet_grid(NUMTASK~SESSION) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab("Time relative to Response 2 ") + ylab("Decoding Accuracy") 
#scale_color_manual(values = c("cornflowerblue", "orange")) + scale_fill_manual(values = c("cornflowerblue", "orange"))

pSW + labs(caption = paste0(modelV, " bal ", balanceV, " plotted by no. feature switches"))



#for tasksets - look at decoding over the course of blocks
r_predTrial <- r_predGRAND[NUMTASK == 'double'] %>%
  group_by(SESSION, TRIAL) %>% #NUMTASK,
  summarise(mean_acc = mean(acc),
            sd_acc = sd(acc),
            sem = sd(acc)/sqrt(sum(!is.na(acc))))
#r_predTrial$SOA<-factor(r_predTrial$SOA, levels = c("100","250","500","750","1000"))
p5<-ggplot(data=r_predTrial,aes(x=TRIAL,y=mean_acc)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  #geom_vline(xintercept=c(126),linetype=1,size=1)+#annotate("text",x=150, y=chance-0.05,label="Stimulus + Cue")+
  geom_hline(yintercept=chance,linetype=1,size=1)+
  geom_ribbon(alpha=0.2,linetype=0, aes(fill = SESSION, ymin = mean_acc - sem, ymax = mean_acc + sem)) +
  ggtitle(paste(modelV, "Dual-task"))+
  geom_line(size=1.5,aes(color=SESSION, group = SESSION)) + 
  #   ylim(c(0.48,0.575))+
  facet_rep_grid(~SESSION) +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab("Time relative to Stimulus 2 onset") + ylab("Decoding Accuracy")
p5

#=For RSA ===========================================================
theme_set(theme_bw(base_size = 20))#32/28
CSCALE_YlOrRd = rev(brewer.pal(9,"YlOrRd"));

setwd(Dir_GRAND)

fileN<- "DimensionSwitchPRP"
modelVs <- "TASKSET_CORSIDE_2"
sessN <- "2"
value <- "betas"
test<- "timeSERIES" #timeSERIES or byTrial
regs <- "BIS" #RTACC
evLock <- "RESP"
decodingType <- "unmerged"
if(!is_empty(regs)){dName <- paste(fileN,modelVs,evLock,regs,value,test, sep="_")} else {dName <- paste(fileN,modelVs,value,test, sep="_")}
#dName<-"DimensionSwitchPRP_unmerged_TASKSET_CORSIDE_2_RESP__betas_timeSERIES"

intV <- c("SUBID", "SESSION","NUMTASK","SOA")
#dVs <- c("TASKSET_RSA", "CORSIDE_RSA", "CONJ_RSA") # list of values that will be DV
measureVs <- c("TASKSET", "CORSIDE", "TASKSET_CORSIDE")

dsG <- readRDS(paste0(dName, ".rds")) #load in data
dsG[,SOA := as.factor(SOA)]
dsG[,NUMTASK := factor(NUMTASK, levels = c("single","double"))]
dsG[,SESSION := as.factor(SESSION)]

##for test == "timeSERIES" plot by Session
dsM<- dsG %>%
  group_by(SESSION, NUMTASK, time) %>% #NUMTASK,SOA
  summarise(val_TASKSET = mean(TASKSET_RSA),
            val_CORSIDE = mean(CORSIDE_RSA),
            val_CONJ = mean(CONJ_RSA),
            sem_TASKSET = sd(TASKSET_RSA)/sqrt(sum(!is.na(TASKSET_RSA))),
            sem_CORSIDE = sd(CORSIDE_RSA)/sqrt(sum(!is.na(CORSIDE_RSA))),
            sem_CONJ = sd(CONJ_RSA)/sqrt(sum(!is.na(CONJ_RSA)))) %>% setDT() %>%
  melt(measure.vars = patterns("val_", "sem_"), variable.name = "feature", value.name = c("value", "sem"))
dsM[, feature := forcats::lvls_revalue(feature, c("TASKSET", "CORSIDE", "TASKSET_CORSIDE"))]


for (measure in measureVs){
print(dsM[feature == measure] %>% # & SOA %in% c("100", "500","1000")
  ggplot(aes(x=time,y=value, group = NUMTASK)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(0),linetype=1,size=1)+#annotate("text",x=150, y=chance-0.05,label="Stimulus + Cue")+
  geom_ribbon(alpha=0.2,linetype=0, aes(fill = NUMTASK, ymin = value - sem, ymax = value + sem)) +
  ggtitle(paste(measure, "RSA"))+
  geom_hline(yintercept = 0.00)+
  geom_line(size=1.5,aes(color=NUMTASK, group = NUMTASK)) + 
  #   ylim(c(0.48,0.575))+
  facet_grid(~SESSION, labeller = label_both) +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab(glue("Time relative to {evLock} {sessN}")) + ylab("beta") +
  scale_color_manual(values = c("cornflowerblue", "orange")) + scale_fill_manual(values = c("cornflowerblue", "orange"))
)
}

for(ttype in levels(dsM$NUMTASK)) {
print(dsM[NUMTASK == ttype & SOA %in% c("100", "500","1000")] %>%
  ggplot(aes(x=time,y=value, group = feature)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(0),linetype=1,size=1)+#annotate("text",x=150, y=chance-0.05,label="Stimulus + Cue")+
  geom_ribbon(alpha=0.2,linetype=0, aes(fill = feature, ymin = value - sem, ymax = value + sem)) +
  ggtitle(paste(ttype, "RSA"))+
  geom_line(size=1.5,aes(color=feature, group = feature)) + 
  #   ylim(c(0.48,0.575))+
  facet_grid(SOA~SESSION, labeller = label_both) +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab(glue("Time relative to Response {sessN}")) + ylab("beta") )
}


#plot by SOA
if(any(grepl("early", names(dsG)))){dVs <- c(paste0("early",dVs), paste0("late",dVs))}
dsG[NUMTASK == 'double'] %>%
  group_by(SESSION, SOA, time) %>% #NUMTASK,
  summarise(TASKSET = mean(TASKSET_RSA),
            CORSIDE = mean(CORSIDE_RSA),
            CONJ = mean(CONJ_RSA),
            sem_TST = sd(TASKSET_RSA)/sqrt(sum(!is.na(TASKSET_RSA))),
            sem_COR = sd(CORSIDE_RSA)/sqrt(sum(!is.na(CORSIDE_RSA))),
            sem_CONJ = sd(CONJ_RSA)/sqrt(sum(!is.na(CONJ_RSA)))) %>%
ggplot(aes(x=time,y=CONJ, group = SOA)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(0),linetype=1,size=1)+#annotate("text",x=150, y=chance-0.05,label="Stimulus + Cue")+
  geom_ribbon(alpha=0.2,linetype=0, aes(fill = SOA, ymin = CONJ - sem_CONJ, ymax = CONJ + sem_CONJ)) +
  ggtitle(paste(modelVs, regs, "dual-task"))+
  geom_line(size=1.5,aes(color=SOA, group = SOA)) + 
  #   ylim(c(0.48,0.575))+
  facet_rep_grid(~SESSION) +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab(glue("Time relative to Stimulus {sessN} onset")) + ylab("beta")
p1 <- dsG %>%
  pivot_longer(cols = 6:8,
               names_to = "value", values_to = "beta", names_transform = list(value = as.factor)) %>%
  group_by(time, SESSION, NUMTASK, value) %>% 
  summarise(beta = mean(beta))%>% 
  ggplot(aes(x = time, y = beta)) +
  geom_point(aes(color = value)) +
  geom_line(aes(color = value))+
  facet_grid(SESSION~NUMTASK)+ xlab(glue("Time relative to Response {sessN}")) +
  scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))  

dsG[NUMTASK == 'single'] %>%
  pivot_longer(cols = 6:8,
               names_to = "value", values_to = "beta") %>%
  group_by(time, SESSION, SOA, measure) %>%
  summarise(beta = mean(beta))%>%
  ggplot(aes(x = time, y = beta)) +
  geom_point(aes(color = value)) +
  geom_line(aes(color = value))+
  facet_grid(SESSION~SOA)


####
##for test == "byTrial"
dsL <- melt(dsG, 
            id.vars = intV,
            measure.vars = dVs,
            variable.name = "predictor", value.name = "beta")
dsL %>%
  group_by(SESSION, NUMTASK, SOA, predictor)%>%
  summarise(beta = mean(beta),
            sem = sd(beta)/sqrt(sum(!is.na(beta)))) %>%
  ggplot(aes(x=predictor,y=beta, group = NUMTASK)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_col(aes(fill = NUMTASK), position = "dodge")+
  ggtitle(paste(modelVs, regs, "dual-task"))+
  #   ylim(c(0.48,0.575))+
  facet_rep_grid(SOA~SESSION) +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab("predictor") + ylab("beta")
            



#plot stimulus RSA stuff if possible
#Load behavioral data
setwd(Dir_GRAND)
ds_rsa<-fread(paste("DimensionSwitchPRP_",modelV,"_SOA_byTrial.txt", sep = ""))
ds_rsa$TS1 <- case_when(
  ds_rsa$TASKSET_1 == "Low-High" ~ "LH",
  ds_rsa$TASKSET_1 == "Odd-Even" ~ "OE",
  ds_rsa$TASKSET_1 == "Red-Blue" ~ "RB",
  ds_rsa$TASKSET_1 == "Bold-Faint" ~ "BF")
ds_rsa$TS2 <- case_when(
  ds_rsa$TASKSET_2 == "Low-High" ~ "LH",
  ds_rsa$TASKSET_2 == "Odd-Even" ~ "OE",
  ds_rsa$TASKSET_2 == "Red-Blue" ~ "RB",
  ds_rsa$TASKSET_2 == "Bold-Faint" ~ "BF")

#dsAGGB<-ds_rsa[, lapply(.SD, mean, na.rm=TRUE), by=c("time", "SESSION", "NUMTASK", "TASKSET_2", "TS2"), .SDcols=c("LH_RSA","OE_RSA","BF_RSA","RB_RSA")]
#dsAGGB<-melt(dsAGGB,id.vars = c("time", "SESSION", "NUMTASK", "TASKSET_2", "TS2"))
#dsAGGB$rel <-  ifelse(mapply(function(pat, txt) grepl(pat, txt), dsAGGB$TS2, dsAGGB$variable), "Relevant Feature", "Irrelevant Feature")
#dsAGG<-dsAGGB[, lapply(.SD, mean, na.rm=TRUE), by=c(paste(c(grpV[-1], "rel"))), .SDcols="value"]


rsa_cols <- grep("_RSA$", names(ds_rsa), value = TRUE)

ds_long <- melt(
  ds_rsa,
  measure.vars = rsa_cols,
  variable.name = "RSA_type",
  value.name    = "RSA_value",
  variable.factor = FALSE
)
#ds_long$rel <- ifelse(mapply(function(pat, txt) grepl(pat, txt), ds_long$TS2, ds_long$RSA_type), "Relevant Feature", "Irrelevant Feature")

ds_long[, rel := "Irrelevant"]
ds_long[
  mapply(grepl, TS2, RSA_type),
  rel := "Relevant"
]

ds_long[
  mapply(grepl, TS1, RSA_type) & rel == "Irrelevant",
  rel := "PreviouslyRelevant"
]


ds_plot <- ds_long[,.(mean_RSA = mean(RSA_value, na.rm = TRUE),
                      sem_RSA  = sd(RSA_value, na.rm = TRUE) / sqrt(sum(!is.na(RSA_value)))),
  by = .(time, SESSION, NUMTASK,  rel) #SOA,
]

ggplot(data = ds_plot[NUMTASK == 'single'], aes(x=time,y=mean_RSA, color=rel))+#,group=rel,color=rel)) +
  #scale_y_continuous(breaks = scales::pretty_breaks(n=2))+
  geom_line(aes(group=rel),size=1) + 
  #geom_point(aes(group=rel))+
  #coord_cartesian(ylim=c(-0.05,0.1)) +#  c(-0.05,2.3) + #stat_summary(size = 0.5, shape = 1)+ 
  geom_ribbon(alpha=0.2,linetype=0, aes(fill = rel, ymin = mean_RSA - sem_RSA, ymax = mean_RSA + sem_RSA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  facet_grid(~SESSION) +
  scale_color_manual(values=c("cornflowerblue","orange", "darkolivegreen4")) + scale_fill_manual(values=c("cornflowerblue","orange", "darkolivegreen4")) +
  xlab("Time after Stimulus 2 onset (ms)") + ylab("t-value")




#=For Cross-condition============================================================
r_predP <- r_predGRAND %>%
  group_by(NUMTASK, time, WITHIN, SESSION) %>% #WITHIN
  summarise(mean_acc = mean(acc),
            sem = sd(acc)/sqrt(sum(!is.na(acc))))
p6<-ggplot(data=r_predP,aes(x=time,y=mean_acc, group = NUMTASK)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(0),linetype=1,size=1)+annotate("text",x=125, y=chance-0.01,label="Stimulus + Cue")+
  geom_hline(yintercept=chance,linetype=1,size=1)+
  geom_ribbon(alpha=0.2,linetype = 0, aes(fill=NUMTASK, ymin = mean_acc - sem, ymax = mean_acc + sem)) +
  ggtitle(modelV)+
  geom_line(size=1.5,aes(color=NUMTASK)) + 
  # ylim(c(0.05,0.08))+
  facet_grid(WITHIN~SESSION, labeller = label_both) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab("Time (ms)") + ylab("Decoding Accuracy") +
  scale_color_manual(values = c("cornflowerblue", "orange")) + scale_fill_manual(values = c("cornflowerblue", "orange"))

p6

#broken down by SOA
r_predSOA <- r_predGRAND[SOA %in% c("100","500","1000")] %>%
  group_by(SOA, time, WITHIN, NUMTASK) %>% #NUMTASK,
  summarise(mean_acc = mean(acc),
            sem = sd(acc)/sqrt(sum(!is.na(acc))))
r_predSOA$SOA<-factor(r_predSOA$SOA, levels = c("100","500","1000")) #"250","750",
p7<-ggplot(data=r_predSOA,aes(x=time,y=mean_acc, group = NUMTASK)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(0),linetype=1,size=1)+#annotate("text",x=150, y=chance-0.05,label="Stimulus + Cue")+
  geom_hline(yintercept=chance,linetype=1,size=1)+
  geom_ribbon(alpha=0.2,linetype=0, aes(fill = NUMTASK, ymin = mean_acc - sem, ymax = mean_acc + sem)) +
  ggtitle(paste(modelV))+
  geom_line(size=1.5,aes(color=NUMTASK, group = NUMTASK)) + 
  #   ylim(c(0.48,0.575))+
  facet_grid(SOA~WITHIN, labeller = label_both) +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab("Time (ms)") + ylab("Decoding Accuracy") +
  scale_color_manual(values = c("cornflowerblue", "orange")) + scale_fill_manual(values = c("cornflowerblue", "orange"))

p7


#Looking at feature dimension switch:

r_predGRAND <- r_predGRAND %>%
  left_join(ds_b %>% dplyr::select(SUBID_S, BLOCK,TRIAL, NUMTASK, SESSION, SOA, STIMNUM, FDIM), #STIMNUM, pick variable of choice
            by = c("SUBID_S","SESSION", "BLOCK", "TRIAL","STIMNUM", "NUMTASK", "SOA"))

r_predGRAND[, FDIM_SW := as.integer(uniqueN(FDIM) > 1), by = .(SUBID_S, SESSION, BLOCK, TRIAL)]

r_predF <- r_predGRAND %>%
  group_by(NUMTASK, time, WITHIN, FDIM_SW, STIMNUM, SESSION) %>% #WITHIN SESSION
  summarise(mean_acc = mean(acc),
            sem = sd(acc)/sqrt(sum(!is.na(acc)))) %>%setDT()
r_predF[, FDIM_SW := as.factor(FDIM_SW)]
p8<-ggplot(data=r_predF[NUMTASK == "double" & STIMNUM == 2],aes(x=time,y=mean_acc, group = FDIM_SW)) + #ymin=mean_acc-sd_acc,ymax=mean_acc+sd_acc, 
  geom_vline(xintercept=c(0),linetype=1,size=1)+annotate("text",x=125, y=chance-0.01,label="Stimulus + Cue")+
  geom_hline(yintercept=chance,linetype=1,size=1)+
  geom_ribbon(alpha=0.2,linetype = 0, aes(fill=FDIM_SW, ymin = mean_acc - sem, ymax = mean_acc + sem)) +
  ggtitle(paste0(modelV, " decoding by FDIM_SW"))+
  geom_line(size=1.5,aes(color=FDIM_SW)) + 
  # ylim(c(0.05,0.08))+
  facet_grid(WITHIN~SESSION, labeller = label_both) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  #scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))+
  xlab("Time (ms)") + ylab("Decoding Accuracy") 
p8

