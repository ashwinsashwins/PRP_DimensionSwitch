# DimensionSwitchPRP(6 session, Study2)
# To do 
# 1. think of more analyses for switching of features!
# 2. Re-code information in meaningful manner
# 3. Exclude flipped order responses...(how to detect this??)! <- cool thing to look into -- why might people respond to task 1 vs task 2 first (assumes a situation where participant waited for both to appear before responding instead of doing one after the other. Is this an individual differences thing or is there something else at play?)
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))

# Load libraries
library(data.table)
library(tidyverse)
library(broom)
library(effects)
library(afex)
library(psych)
library(stringr)
library(forcats)
library(lme4)
library(RColorBrewer)
library(ggplot2)
library(zoo)
library(emmeans)
library(lemon)
library(sjPlot)
library(optimx)
library(pagedown)

#Setting up directory
fsep<-.Platform$file.sep;
Dir_R<-path.expand("~/Brown Dropbox/Ashwin Srinivasan/")
Dir_BDATA<-path.expand("~/Brown Dropbox/Ashwin Srinivasan/DimensionSwitch_PRP/DataAnalysis/Matlab/BEH/z_ALLGRAND(EEG,6sessions)")
Dir_EDATA<-path.expand("~/Brown Dropbox/Atsushi Kikumoto/w_SCRIPTS_H/DimensionSwitch_PRP/DataAnalysis/Matlab/EEG")

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
ds <- fread('DimensionSwitchPRP_BehP.txt', header=TRUE, fill=TRUE) 
#fread('~/Brown Dropbox/Ashwin Srinivasan/DimensionSwitch_PRP/DataAnalysis/Matlab/BEH/z_ALLGRAND(EEG,6sessions)/DimensionSwitchPRP_BehP.txt', header=TRUE, fill=TRUE)

# Change NaNs to NAs
ds[is.nanM(ds)] <- NA;
ds <- ds[, !names(ds) %like%  "EV_", with=F]

# Recode variables
ds[,SUBID_S:=SUBID]
ds[,SUBID:=SUBID_S - 100*(SESSION-1)]
ds[,NUMTASK := factor(NUMTASK, labels = c("single","double"))]
ds[,RESPCONG :=factor(RESPCONG,labels = c("0","1"))]
ds[,EXACT_REP:= factor(EXACT_SW, labels = c("exact","diff"))] # flip of EXACT_SW
ds[,FDIM_1:= factor(FDIM_1, labels = c("Number","Color"))]
ds[,FDIM_2:= factor(FDIM_2, labels = c("Number","Color"))]
ds[,TASKSET_1:= factor(TASKSET_1, labels = c("Low-High","Odd-Even","Red-Blue","Bold-Faint"))]
ds[,TASKSET_2:= factor(TASKSET_2, labels = c("Low-High","Odd-Even","Red-Blue","Bold-Faint"))]
ds[, `:=` (SOAn = SOA, SOA = factor(SOA, labels = c("100","250","500","750","1000")))]
contrasts(ds$SOA, how.many=2) <- contr.poly(5)

# Recode stimulus variables
s_number <- c("3","4","6","7")
s_color <- c("red","red_faint",'blue_faint','blue')
s_set <- expand.grid(s_color, s_number)
ds[,F_COLOR_1:=factor(F_COLOR_1, label=s_color)]
ds[,F_COLOR_2:=factor(F_COLOR_2, label=s_color)]
ds[,F_NUMBER_1:=factor(F_NUMBER_1, label=s_number)]
ds[,F_NUMBER_2:=factor(F_NUMBER_2, label=s_number)]
ds[,S1_1:=factor(S1_1, label=s_number)] # first column of S1 = number dimension
ds[,S1_2:=factor(S1_2, label=s_color)] # second column of S1 = color dimension
ds[,S2_1:=factor(S2_1, label=s_number)] # first column of S2 = number dimension
ds[,S2_2:=factor(S2_2, label=s_color)] # second column of S2 = color dimension
setDT(ds)[, S1IDc := Reduce(function(...) paste(..., sep = "-"), .SD), .SDcols = c("S1_1","S1_2")]
setDT(ds)[, S2IDc := Reduce(function(...) paste(..., sep = "-"), .SD), .SDcols = c("S2_1","S2_2")]

# FDIM_MATCH (only switch blocks): 
# Match = switches between same features
# Mismatch = switches between different features
ds[,FDIM_MATCH:=NA]
ds[,FDIM_MATCH:=ifelse(NUMTASK=="double" & FDIM_1 == FDIM_2, "Match",FDIM_MATCH)]
ds[,FDIM_MATCH:=ifelse(NUMTASK=="double" & FDIM_1 != FDIM_2, "Mismatch",FDIM_MATCH)]

# RT cutoff (exclude outlier like RTs > 10 sec)
cutRT1i <- ds[ACC_1 == 1 & ds$RT_1 < 10000,list(cutRT1 = quantile(RT_1, 0.99)),by = c("SUBID","SESSION")]
cutRT2i <- ds[ACC_2 == 1 & ds$RT_2 < 10000,list(cutRT2 = quantile(RT_2, 0.99)),by = c("SUBID","SESSION")]
ds <- merge(ds, cutRT1i,by=c("SESSION","SUBID"))
ds <- merge(ds, cutRT2i,by=c("SESSION","SUBID"))
print(sprintf("Response cutoff %.3f for Response 1",mean(cutRT1i$cutRT)));
print(sprintf("Response cutoff %.3f for Response 2",mean(cutRT2i$cutRT)));

# Checking
rt<-ds[ACC_1==1 & ACC_2==1,]  %>% group_by(SUBID,SESSION) %>% summarise(RT_1=mean(RT_1,na.rm=TRUE),RT_2=mean(RT_2,na.rm=TRUE));as.data.frame(rt)
acc<-ds %>% group_by(SUBID,SESSION) %>% summarise(ACC_1=mean(ACC_1,na.rm=TRUE), ACC_2=mean(ACC_2,na.rm=TRUE));as.data.frame(acc)
numT<-ds %>% group_by(SUBID,SESSION) %>% summarise(T=max(TRIAL));as.data.frame(numT) #aha! because number of trials is variable -- good to keep in mind
numdata<-ds %>% group_by(SUBID,SESSION) %>% summarise(T=n());as.data.frame(numdata) 
participants<-ds %>% group_by(SESSION) %>% summarise(T=list(unique(SUBID)));as.data.frame(participants)
resporder<-ds %>% group_by(SUBID,SESSION, SOA) %>% summarise(RORDER=mean(RESPORDER_K,na.rm=T));as.data.frame(resporder)

# Export for updated data (not filtered at all!)
#write.table(ds,paste0(Dir_BDATA,'/DimensionSwitchPRP_BehPP.txt'), sep="\t",row.names=FALSE)

# ASHWIN EDIT: plotting # of response flipping trials by participants across sessions
# bumps seem to correspond to behavioral sessions -- dips correspond to EEG sessions

ds[, SESSION := factor(SESSION)]
ds <- ds %>%
  mutate(eeg = if_else(SESSION == "1" | SESSION == "4" | SESSION == "5", "beh", "eeg"))

ds_response_switched <- ds[RESPORDER_K == 0,]
ds_response_switched[, SUBID := factor(SUBID)]
ds_response_switched[, RESPORDER_K := factor(RESPORDER_K)]
ds_response_switched[, SOA := factor(SOA)] 
switched_by_session <- ds_response_switched %>% # will want to normalize this for # of trials per session
  group_by(SUBID, SESSION) %>%
  summarise(T = sum(RESPORDER_K == 0))
ggplot(dat = switched_by_session, aes(x = SESSION, y = T, color = SUBID)) + geom_line() + theme_bw()


switched_by_SOA <- ds_response_switched %>%
  group_by(SUBID, SOA) %>%
  summarise(T = sum(RESPORDER_K == 0))
ggplot(dat = switched_by_SOA, aes(x = SOA, y = T, group = SUBID)) + geom_line(aes(color = SUBID)) + theme_bw()

#how does switches by session change when removing 100ms SOA trials?
ds_response_switched %>%
  filter(SOA != 100) %>%
  group_by(SUBID, SESSION) %>%
  summarise(T = sum(RESPORDER_K == 0)) %>%
  ggplot(aes(x = SESSION, y = T, color = SUBID)) + geom_line() + theme_bw()
#seems to be driven by 100ms trials -- delve into this more
#seems to be not learning-dependent, rather strategy, self-pacing, tiring over time, not-rushing
# but maybe it is learning-dependent -- look at ordered vs. flipped trials @100ms over sessions -- does response quality (particularly accuracy) is tied to 

ds[, RESPORDER_K := factor(RESPORDER_K)]
fast_switches <- ds %>%
  filter(SOA == 100) %>%
  group_by(SESSION, RESPORDER_K) %>%
  summarise(ACC_1 = mean(ACC_1), ACC_2 = mean(ACC_2), T = length(RESPORDER_K)) # T corresponds to # of trials 
ggplot(data = fast_switches, aes(x = SESSION, y = ACC_1, color = RESPORDER_K)) + geom_line() + theme_bw()
ggplot(data = fast_switches, aes(x = SESSION, y = ACC_2, color = RESPORDER_K)) + geom_line() + theme_bw()
# performance is worse on flipped trials @ short SOAs -- why? can this be explained by features of the task? (i.e., response flip & poor performance are both outcomes of specific task setup)
# note that accuracy for flipped trials drops on day 2 because day 1 is single task only

fast_trials_switched <- ds_response_switched %>%
  filter(SOA == 100)

# flipping for fast trials over sessions
fast_trials_switched[, SUBID := factor(SUBID)]
fast_trials_switched[, SESSION := factor(SESSION)]
fast_trials_session <- fast_trials_switched %>%
  group_by(SUBID, SESSION, eeg) %>%
  summarise(switches = sum(RESPORDER_K == 0)) %>%
  merge(numdata,by=c("SESSION","SUBID"))
fast_trials_session$propSwitch <- fast_trials_session$switches/fast_trials_session$T

mean_fast_trials_session <- fast_trials_session %>%
  group_by(SESSION) %>%
  summarise(meanSwitches = mean(propSwitch))
fast_switch_plot <- ggplot(data = fast_trials_session, aes(x = SESSION, y = propSwitch)) +
  geom_point(aes(group = SUBID, color = SUBID), alpha = 0.6) +
  #geom_line(data = mean_fast_trials_session, aes(x = SESSION, y = meanSwitches), group = 1)+
  facet_wrap(~eeg) +
  theme_bw()
fast_switch_plot

#let's look at single vs. double trial trajectory over time
ds_response_switched %>%
  group_by(SESSION, NUMTASK) %>%
  summarise(switches = sum(RESPORDER_K == 0)) %>%
  ggplot(aes(x = SESSION, y = switches, color = NUMTASK)) +
  geom_line(aes(group = NUMTASK, color = NUMTASK)) +
  geom_point() +
  theme_bw()
#ok, so people switch more on double trials than single trials, but the SESSION-trajectory is the same across both
#but also switching this to proportion of trials switched changes things a bit: (proportion of trials of that type switched)
ds %>%
  group_by(SESSION, NUMTASK) %>%
  summarise(switches = (sum(RESPORDER_K == 0)/length(RESPORDER_K))) %>%
  ggplot(aes(x = SESSION, y = switches, color = NUMTASK)) +
  geom_line(aes(group = NUMTASK, color = NUMTASK)) +
  geom_point() +
  theme_bw()
#in that people switch double trials at a lower rate than single trials

fast_trials <- ds %>%
  filter(SOA == 100)

fast_trials %>%
  group_by(NUMTASK, RESPORDER_K) %>%
  summarise(ACC_1 = mean(ACC_1), ACC_2 = mean(ACC_2), T = length(RESPORDER_K)) %>%
ggplot(aes(x = NUMTASK, y = ACC_1, fill = RESPORDER_K)) + geom_col(position = "dodge") + theme_bw()  
#accuracy is very low on switched trials in the double task -- task-switching cost?
#flipped response is NOT associated with poor performance on the single task

#tracking performance on flipped double trials over time
fast_trials %>%
  group_by(SESSION, RESPORDER_K, NUMTASK) %>%
  summarise(ACC_1 = mean(ACC_1), ACC_2 = mean(ACC_2), T = length(RESPORDER_K)) %>%
ggplot(aes(x = SESSION, y = ACC_1, color = RESPORDER_K)) + geom_line() + geom_point() + facet_wrap(~NUMTASK) + theme_bw()

fast_trials %>%
  group_by(SESSION, RESPORDER_K, NUMTASK) %>%
  summarise(Acc_1 = mean(ACC_1), ACC_2 = mean(ACC_2), sd_1 = sd(ACC_1)) %>%
  ggplot(aes(x = SESSION, y = Acc_1, color = RESPORDER_K)) + 
  geom_line() + 
  geom_point() + 
  geom_linerange(aes(ymin = Acc_1-sd_1, ymax = Acc_1+sd_1), linewidth = 0.4) + 
  facet_wrap(~NUMTASK) + 
  theme_bw()
# ACC_1 of flipped, double trials gets better over time! ACC_2 does NOT. ACC_2 refers to accuracy of Item 2, which is the FIRST item responded to on flipped trials -- why is the accuracy so low (chance)?
# learning effect? -- flipped trials suggest a different strategy -- i.e., accumulate everything --> decide as opposed to 

#let's see if poor ACC_2 on flipped trials is because participants are flipping the item-hand mappings in their heads
ds[, SUBID := factor(SUBID)]

ds %>% # number of flipped trials by participant
  filter(RESPORDER_K == 0 & ACC_1 == 0 & ACC_2 == 0) %>% #flipped, incorrect trials
  group_by(SUBID) %>%
  summarise(T = length(RESPORDER_K)) %>%
  ggplot(aes(x = SUBID, y = T, color = SUBID, fill = SUBID)) + geom_col(position = "dodge") + theme_bw() #
 
swappedFlips <- ds %>%
  mutate(swap1 = ds$RESP_2, swap2 = ds$RESP_1, .after = RESP_2) %>% #this is incorrect. Keys are coded differently on opposite hands. 
  mutate(acc_1 = 0, acc_2 = 0, .after = ACC_2) %>%
  filter(RESPORDER_K == 0 & ACC_1 == 0 & ACC_2 == 0) %>%
  mutate(acc_1 = if_else(swap1 == CORRESP_1, 1, 0)) %>%
  mutate(acc_2 = if_else(swap2 == CORRESP_2, 1, 0))
ifelse(swappedFlips$swap1 == swappedFlips$CORRESP_1, acc_1 = 1, acc_1 = 0)

#similarly look at RT -- only for correct trials perhaps?
fast_trials %>%
  filter(ACC_1 == 1 & ACC_2 == 1) %>%
  group_by(SESSION, RESPORDER_K, NUMTASK) %>%
  summarise(RT_1 = mean(RT_1), RT_2 = mean(RT_2), T = length(RESPORDER_K)) %>%
  ggplot(aes(x = SESSION, y = RT_1, color = RESPORDER_K)) + geom_line() + geom_point() + facet_wrap(~NUMTASK) + theme_bw()
fast_trials %>%
  filter(ACC_1 == 1 & ACC_2 == 1) %>%
  group_by(SESSION, RESPORDER_K, NUMTASK) %>%
  summarise(RT_1 = mean(RT_1), RT_2 = mean(RT_2), T = length(RESPORDER_K)) %>%
  ggplot(aes(x = SESSION, y = RT_2, color = RESPORDER_K)) + geom_line() + geom_point() + facet_wrap(~NUMTASK) + theme_bw()


fast_trial_double <- fast_trials %>%
  filter(NUMTASK == 'double')


# ok. so participants flip responses when:
#    - SOA = 100ms
#    - double-task (also single-task, but accuracy is only affected in double-task condition)
# and this is reflected in poor accuracy, especially ACC_2, whereas ACC_1 shows some learning trajectory(stats needed). What next?
#    - are there more specific conditions that participants are flipping responses? (i.e., specific feature switches)
#    - what, if anything, can explain the drop in flipped responses between sessions 1 & 2? presumably combo of learning + EEG effect
#    - why (and under what conditions) do participants flip on the single-task? is there something in the data that explains this? especially since accuracy is not affected, and SOA is consistent (100ms)

fast_trial_double %>%
  group_by(SESSION, S2IDc) %>%
  summarise(switches = (sum(RESPORDER_K == 0)/length(RESPORDER_K))) %>%
  ggplot(aes(x = SESSION, y = switches,  color = S2IDc)) +
  geom_line(aes(group = S2IDc, color = S2IDc)) +
  geom_point()+
  theme_bw()
fast_trial_double %>%
  group_by(SESSION, S1IDc) %>%
  summarise(switches = (sum(RESPORDER_K == 0)/length(RESPORDER_K))) %>%
  ggplot(aes(x = SESSION, y = switches, color = S1IDc)) +
  geom_line(aes(group = S1IDc, color = S1IDc)) +
  geom_point() +
  theme_bw()
#all stim features show the same shape -- but certainly variability in number -- have to control for number of trials however

fast_trials[, TASKPAIR := factor(TASKPAIR)]
fast_trials %>%
  group_by(SESSION, TASKPAIR) %>%
  summarise(switches = (sum(RESPORDER_K == 0)/length(RESPORDER_K))) %>%
  ggplot(aes(x = SESSION, y = switches)) +
  geom_line(aes(group = TASKPAIR, color = TASKPAIR)) +
  theme_bw() 
  #geom_label(aes(label = TASKPAIR))
#interestingly, the taskpairs that have the most switches in terms of number are 1,11,6,16 -- all single-task
# also worth noting that taskpairs are independent of specific features of the stimulus on a given trial, so:

#number of features switching each trial
fast_trial_double[, NUMSW := factor(NUMSW)]
fast_trial_double %>%
  group_by(SESSION, NUMSW) %>%
  summarise(switches = (sum(RESPORDER_K == 0)/length(RESPORDER_K))) %>%
  ggplot(aes(x = SESSION, y = switches)) +
  geom_line(aes(group = NUMSW, color = NUMSW)) +
  theme_bw() +
  geom_label(aes(label = NUMSW))
#not a particularly telling relationship, if anything ppl switch more when NO features switch
fast_trial_double %>%
  group_by(SESSION, NUMSW) %>%
  summarise(ACC_1 = (mean(ACC_1))) %>%
  ggplot(aes(x = SESSION, y = ACC_1)) +
  geom_line(aes(group = NUMSW, color = NUMSW)) +
  theme_bw() +
  geom_label(aes(label = NUMSW))
#accuracy is lower when no features switch? idk abt significance 
#also why is accuracy so low on days 4 &5??
ggplot(data = ds, aes(x = SESSION, y = mean(ACC_1)))+ geom_point() + theme_bw() +
  facet_wrap(~SOA)
tab <- ds %>%
  group_by(SESSION, NUMTASK, SOA) %>%
  summarise(T = mean(ACC_1))
#what about relevant features
fast_trial_double[, REL_SW := factor(REL_SW)]
fast_trial_double %>%
  group_by(SESSION, REL_SW) %>%
  summarise(switches = (sum(RESPORDER_K == 0)/length(RESPORDER_K))) %>%
  ggplot(aes(x = SESSION, y = switches)) +
  geom_line(aes(group = REL_SW, color = REL_SW)) +
  theme_bw() +
  geom_label(aes(label = REL_SW))
#probs not a significant diff here

#what about response congruence (hand response)
fast_trial_double[, RESPCONG := factor(RESPCONG)]
fast_trial_double %>%
  group_by(SESSION, RESPCONG) %>%
  summarise(switches = (sum(RESPORDER_K == 0)/length(RESPORDER_K))) %>%
  ggplot(aes(x = SESSION, y = switches)) +
  geom_line(aes(group = RESPCONG, color = RESPCONG)) +
  theme_bw() +
  geom_label(aes(label = RESPCONG))
#most likely not an effect here, but interesting if it is

#what about just where things appear on the screen??? STIMPOS_2 = 1 means the second item appeared in the top position
fast_trial_double[, STIMPOS_1 := factor(STIMPOS_1)]
fast_trial_double %>%
  filter(ACC_1 == "1" & ACC_2 == "1") %>%
  group_by(SESSION, STIMPOS_1) %>%
  summarise(switches = (sum(RESPORDER_K == 0)/length(RESPORDER_K))) %>%
  ggplot(aes(x = SESSION, y = switches)) +
  geom_line(aes(group = STIMPOS_1, color = STIMPOS_1)) +
  theme_bw() +
  geom_label(aes(label = STIMPOS_1))

#ok this sort of stands to reason, especially since most switches happen on 100 ms SOA trials -- people wait for both and respond top -> bottom
#yet still some switched responses when the 2nd item presented at the bottom

ds[, STIMPOS_2 := factor(STIMPOS_2)]
ds %>%
  group_by(SESSION, NUMTASK, STIMPOS_2) %>%
  summarise(propSwitches = (sum(RESPORDER_K == 0)/length(RESPORDER_K)),
            error = sqrt((propSwitches * (1-propSwitches))/sum(RESPORDER_K == 0))) %>%
  ggplot(aes(x = SESSION, y = propSwitches)) +
  geom_line(aes(group = STIMPOS_2, color = STIMPOS_2)) + 
  geom_linerange(aes(ymin = propSwitches - error, ymax = propSwitches + error, color = STIMPOS_2), linewidth = 0.4)+
  facet_wrap(~NUMTASK) +
  theme_bw() 

#lets look at this and accuracy and RT
ds %>%
  filter(STIMPOS_2 == '1') %>%
  group_by(SESSION,NUMTASK, RESPORDER_K) %>%
  summarise(acc_2 = mean(ACC_2)) %>%
  ggplot(aes(x = SESSION, y = acc_2)) +
  geom_line(aes(group = RESPORDER_K, color = RESPORDER_K)) + 
#  geom_linerange(aes(ymin = propSwitches - error, ymax = propSwitches + error, color = STIMPOS_2), linewidth = 0.4)+
  facet_wrap(~NUMTASK) +
  theme_bw() 

#look at where flips happen within each session (could break down into block also)
test <- ds %>%
  group_by(SESSION,BLOCK) %>%
  summarise(flips = as.numeric(sum(RESPORDER_K == 0)))
test$BLOCK <- as.numeric(test$BLOCK)
ggplot(data = test, aes(x = BLOCK, y = flips,fill = SESSION)) + #plot of proportion of trials that are flipped in each block
  geom_col() +
  facet_wrap(~SESSION) +
  theme_bw()
# also makes sense (probably) to break up into days 1,4,5 (beh) and 2,3,6 (EEG)


# Filtering
# ds<-ds[BLOCK>0,]
ds<-ds[RESPOMIT==0] # Exclude response omission trials?
ds_roall <-ds # keep this for later analysis
ds<-ds[RESPORDER_K==1]# Exclude response flipping trials?
ds<-ds[!(c(ds$RT_1 > ds$cutRT1)),];
ds<-ds[!(c(ds$RT_2 > ds$cutRT2)),];
# ASHWIN EDIT: might be good to remove trials with Negative RT_2 (3 observations)

# Make few separate datasets
ds_f<-ds[ACC_1==1 & ACC_2==1,] 
ds_fdouble<-ds_f[NUMTASK=='double',c("SUBID","BLOCK","TRIAL","TASKSET_1","TASKSET_2",
                                     "S1IDc","S2IDc","FDIM_1","FDIM_2","F_NUMBER_1","F_NUMBER_2",
                                     "F_COLOR_1","F_COLOR_2","REL_SW"), with=F]

#======================================================================================================
#Process2-Data Aggregation
#======================================================================================================
# Number of observations / Grand average RT & ACC
# ds_i <- ds_f[,.N,by=c("SUBID")]
# ds_i$RT<-ds_f[,mean(RT),by=c("SUBID")]$V1
# ds_i$ACC<-ds[,mean(ACC),by=c("SUBID")]$V1

# calculation of within subject standard

# RESPORDER x SOA
RORDER_SOAi <- ds_roall %>% group_by(SUBID,SESSION,NUMTASK,SOA) %>% summarise(DV=mean(RESPORDER_K, na.rm=TRUE) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
RORDER_SOA <- RORDER_SOAi %>% group_by(SESSION,NUMTASK,SOA) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# RT_2:SESSION x SOA
RT2_SOAi <- ds_f %>% group_by(SUBID,SESSION,NUMTASK,SOA) %>% summarise(DV=mean(RT_2, na.rm=TRUE) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
RT2_SOA <- RT2_SOAi %>% group_by(SESSION,NUMTASK,SOA) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# ACC_2:SESSION x SOA
ERR2_SOAi <- ds %>% group_by(SUBID,SESSION,NUMTASK,SOA) %>% summarise(DV=100*(1-mean(ACC_2,na.rm=TRUE)) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
ERR2_SOA <- ERR2_SOAi %>% group_by(SESSION,NUMTASK,SOA) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# RT_2:SESSION x SOA x FDIMATCH (double task)
RT2_FDIMMi <- ds_f %>% group_by(SUBID,SESSION,NUMTASK,SOA,FDIM_MATCH) %>% summarise(DV=mean(RT_2, na.rm=TRUE) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
RT2_FDIMM <- RT2_FDIMMi %>% group_by(SESSION,NUMTASK,SOA,FDIM_MATCH) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# ACC_2:SESSION x SOA
ERR2_FDIMMi <- ds %>% group_by(SUBID,SESSION,NUMTASK,SOA,FDIM_MATCH) %>% summarise(DV=100*(1-mean(ACC_2,na.rm=TRUE)) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
ERR2_FDIMM <- ERR2_FDIMMi %>% group_by(SESSION,NUMTASK,SOA,FDIM_MATCH) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# RT_1:SESSION x SOA
RT1_SOAi <- ds_f %>% group_by(SUBID,SESSION,NUMTASK,SOA) %>% summarise(DV=mean(RT_1, na.rm=TRUE) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
RT1_SOA <- RT1_SOAi %>% group_by(SESSION,NUMTASK,SOA) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# ACC_1:SESSION x SOA
ERR1_SOAi <- ds %>% group_by(SUBID,SESSION,NUMTASK,SOA) %>% summarise(DV=100*(1-mean(ACC_1,na.rm=TRUE)) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
ERR1_SOA <- ERR1_SOAi %>% group_by(SESSION,NUMTASK,SOA) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# *************
# Output (response interference)
# *************
# RT_2:SESSION x SOA X RESPCONG
RT2_RCONGi <- ds_f %>% group_by(SUBID,SESSION,NUMTASK,SOA,RESPCONG) %>% summarise(DV=mean(RT_2, na.rm=TRUE) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
RT2_RCONG <- RT2_RCONGi %>% group_by(SESSION,NUMTASK,SOA,RESPCONG) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# ACC_2:SESSION x SOA X RESPCONG
ERR2_RCONGi <- ds %>% group_by(SUBID,SESSION,NUMTASK,SOA,RESPCONG) %>% summarise(DV=100*(1-mean(ACC_2,na.rm=TRUE)) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
ERR2_RCONG <- ERR2_RCONGi %>% group_by(SESSION,NUMTASK,SOA,RESPCONG) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()


# *************
# Task separation (on exact inputs)
# *************
RT2_EXACTi <- ds_f %>% group_by(SUBID,SESSION,NUMTASK,SOA) %>% summarise(DV=mean(RT_2, na.rm=TRUE) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
RT2_EXACT <- RT2_EXACTi %>% group_by(SESSION,NUMTASK,SOA) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

ERR2_EXACTi <- ds %>% group_by(SUBID,SESSION,NUMTASK,SOA) %>% summarise(DV=100*(1-mean(ACC_1,na.rm=TRUE)) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
ERR2_EXACT <- ERR2_EXACTi %>% group_by(SESSION,NUMTASK,SOA) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()


# *************
# TESTING STUFF
# *************

# RT_2:SESSION x N FACTOR
RT_TESTi <- ds_f %>% group_by(SUBID_S,NUMSW,NUMTASK,SOA) %>% summarise(DV=mean(RT_2, na.rm=TRUE) ,nn=n()) %>% normWS("SUBID_S","DV") %>% print()
RT_TEST <- RT_TESTi %>% group_by(NUMSW,NUMTASK,SOA) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# ACC_2:SESSION x N FACTOR
ERR_TESTi <- ds %>% group_by(SUBID_S,SESSION,NUMTASK,SOA,NUMSW) %>% summarise(DV=100*(1-mean(ACC_2,na.rm=TRUE)) ,nn=n()) %>% normWS("SUBID_S","DV") %>% print()
ERR_TEST <- ERR_TESTi %>% group_by(SESSION,NUMTASK,SOA,NUMSW) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()


#======================================================================================================
#Conjunction RTs for control analysis of RSA
#======================================================================================================
# # Minimum number of trials
# write.table(numdata,paste0(Dir_BDATA,'/HAL2017_PUSHMORE_NUMTRIAL.txt'), sep="\t",row.names=FALSE)
# 
# # PTSRCONJ (GroupSpecific)
# cj_DVi<-ds[!is.na(RT) & LACC==1 & ACC==1,.(RT=mean(RT),RTlog=log(mean(RT))),by=c("SUBID_S","PTSRGROUP","PTSRCONJ")]
# cj_DVi[,ACC:=ds[,.(ACC=mean(ACC)),by=c("SUBID_S","PTSRGROUP","PTSRCONJ")]$ACC] 
# cj_DVi<-cj_DVi[order(SUBID_S,PTSRGROUP,PTSRCONJ),]
# write.table(cj_DVi,paste0(Dir_BDATA,'/HAL2017_PUSHMORE_PTSRCONJ_CONTROL_GRP.txt'), sep="\t",row.names=FALSE)

#======================================================================================================
#Process2-Data Vincentising
#======================================================================================================
# # Vincentising?
# bin_RT<-ds_f %>% group_by(SUBID,RSI) %>% 
#   do(data.frame(t(quantile(.$RT,probs=seq(0,1,0.1)))))%>%
#   gather(prob, DV, X10.:X90.)%>% normWS("SUBID","DV")%>%
#   group_by(RSI,prob)%>%
#   summarize(DVm=mean(DV),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n());
# 
# bin_RT$prob<-as.numeric(str_extract(bin_RT$prob,"[[:digit:]]+"))

#======================================================================================================
#Data Aggregation for individual difference!
#======================================================================================================

# # Save results
# write.table(ds_i,paste0(Dir_BDATA,'/DimensionSwitchPRP_BEH_I.txt'), sep="\t",row.names=FALSE)
# #======================================================================================================
# #Process-Statistical Analysis 
# #======================================================================================================
# SESSION*SOA*NUMTASK effect
m_rt=lmer(RT_2~1+SESSION*SOA*NUMTASK+(1|SUBID),data=ds_f[SESSION!=1,]);std_beta_MLM(m_rt)@fixef
# m_acc=glmer(ACC_2~1+SOA+(1+SOA|SUBID),family=binomial,data=ds[SESSION!=1,]);std_beta_MLM(m_acc)@fixef

# SESSION*SOA*NUMTASK effect (group level)
aggDV<-ds_f[SESSION!=1,.(DV=mean(RT_2,na.rm=T)),by=c("SUBID","SESSION","SOA","NUMTASK")]
# aggDV<-ds[SESSION!=1,.(DV=mean(ACC_2)),by=c("SUBID","SESSION","SOA","NUMTASK)]
anova_r<-aov_ez("SUBID","DV",aggDV,within=c("SESSION","SOA","NUMTASK"),anova_table=list(es="pes"));print(anova_r)

# Testing stuff
dsACT <-ds_f[SESSION!=1,]
m_rt=lmer(RT_2~1+SESSION*SOA*NUMTASK*RESPCONG+(1|SUBID),data=dsACT);std_beta_MLM(m_rt)@fixef

# Plot the model
quartz(width=4,height=5) # RT
plot_model(m_rt, type="std", show.p=T, show.values=T, value.size=2,value.offset=.25) + ylim(c(-.05,.05))+scale_y_continuous(labels=dropLeadingZero)

# Make a table of the model
fname <- "DimensionSwitchPRP"
tab_model(m_rt, show.std=T,show.est=F, show.stat=T,
          string.std = "Beta", string.std_ci = "CI",
          file = paste0(fname, ".html"), use.viewer =F,
          col.order = c("std.est","std.ci","stat","p"))

# chrome_print(paste0(fname, ".html"), output = paste0(fname, ".pdf"))

#======================================================================================================
#Process3-Plotting-ggplot2
#======================================================================================================
theme_set(theme_bwL(base_size = 16))#32/28
CSCALE_PURD = rev(brewer.pal(9,"PuRd"));
CSCALE_BLUE = rev(brewer.pal(9,"Blues"));
CSCALE_PiYG = rev(brewer.pal(11,"PiYG"));
CSCALE_RdBu = rev(brewer.pal(11,"RdBu"));
CSCALE_PAIRED = rev(brewer.pal(12,"Paired"));
CSCALE_YlGnBu = rev(brewer.pal(9,"YlGnBu"));
CSCALE_BrBG = rev(brewer.pal(9,"BrBG"));
CSCALE_Set1 = (brewer.pal(9,"Set1"));
CSCALE_Pair = (brewer.pal(9,"Paired"));
CSCALE_Grey = (brewer.pal(9,"Greys"));
CSCALE_Accent = (brewer.pal(8,"Accent"));


# Individual's improvement
# Basic results
quartz(width=10,height=4.70) # RT
# quartz(width=3.5,height=2) # ERR
# <Multitasking costs>
# rw_DV<-RT_SOA;rw_DVi<-RT_SOAi;
# rw_DV<-ERR_SOA;rw_DVi<-ERR_SOAi;

# <Response congruency effect>
rw_DV<-RT2_RCONG;rw_DVi<-RT2_RCONGi;
# rw_DV<-ERR2_RCONG;rw_DVi<-ERR2_RCONGi;

# <Exact repeat of stimulus>
# rw_DV<-RT2_EXACT;rw_DVi<-RT2_EXACTi;
# rw_DV<-ERR2_EXACT;rw_DVi<-ERR2_EXACTi;

ggplot() + 
  #Panels ----------------------------
  facet_grid(~SESSION)+
  #geom_hline(yintercept = rw_DV[NUMTASK=="double" & SOA==100 & SESSION==2]$DVm, linetype=3, alpha = .8, color="red")+
  #geom_hline(yintercept = rw_DV[NUMTASK=="single" & SOA==100 & SESSION==1]$DVm, linetype=3, alpha = .8)+
  #Main layer-------------------------
  geom_errorbar(data=rw_DV,aes(x=SOA,ymin=DVm-wseb, ymax=DVm+wseb,group=NUMTASK,color=interaction(NUMTASK)),width=.0,size=.75)+#,color="black"
  #geom_point(data=rw_DVi,aes(x=SOA,y=DV,group=interaction(SUBID,NUMTASK),color=NUMTASK),size=1.5,alpha=0.10)+
  #geom_line(data=rw_DVi,aes(x=SOA,y=DV,group=interaction(SUBID,NUMTASK),color=NUMTASK),size=0.5,alpha=0.10) +
  geom_line(data=rw_DV,aes(x=SOA,y=DVm,group=interaction(NUMTASK,RESPCONG),color=interaction(NUMTASK),linetype=RESPCONG),size=1) +
  geom_point(data=rw_DV,aes(x=SOA,y=DVm,group=interaction(NUMTASK,RESPCONG),color=interaction(NUMTASK)),size=3.5) +
  #Axis & label-----------------------
  #coord_cartesian(ylim=c(375,2000))+
  #scale_y_continuous(breaks=seq(0,2000,250))+
  #coord_cartesian(ylim=c(0,10))+scale_y_continuous(breaks=c(0,10))+
  ylab("RT of 2nd response (ms)") + xlab("SOA (ms)")+
  #ylab("Error of 2nd response (ms)") + xlab("SOA (ms)")+
  #Colors ----------------------------
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  scale_color_manual(values=c("black","red","gray","red"))


