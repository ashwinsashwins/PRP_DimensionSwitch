# DimensionSwitchPRP
# goal of this script is to visualize and statistically test behavioral effects of interest
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
Dir_REPO <- path.expand("/oscar/home/asrini17/data/asrini17/DimensionSwitch_PRP/PRP_DimensionSwitch/")


#Load libraries
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
library(afex)
library(effects)

#Source Files
setwd(Dir_R)
source('basic_theme.R')
source('basic_lib.R')

#Load behavioral data
setwd(Dir_BDATA)
ds_b<-fread("DimensionSwitchPRP_BehPP.txt") 
ds_b[is.nanM(ds_b)]<-NA

#Filtering
ds_b <- ds_b[SUBID != 125] #we don't want this participant in here, righty?
ds_b<-ds_b[RESPOMIT==0] # Exclude response omission trials?
ds_b<-ds_b[RESPORDER_K==1]# Exclude response flipping trials?
ds_b<-ds_b[!(c(ds_b$RT_1 > ds_b$cutRT1)),]; ds_b<-ds_b[!(c(ds_b$RT_1 < 0)),]
ds_b<-ds_b[!(c(ds_b$RT_2 > ds_b$cutRT2)),];ds_b<-ds_b[!(c(ds_b$RT_2 < 0)),]

#recoding variables
ds_b[,SOA := as.factor(SOA)]
ds_b[,NUMTASK := factor(NUMTASK, levels = c("single","double"))]
ds_b[,SESSION := as.factor(SESSION)]

#for looking at rt!
ds_f<-ds_b[ACC_1==1 & ACC_2==1,] 

#for looking at eegtrials
ds_eeg <- ds_b %>%
  filter(SESSION %in% c(2, 3, 6))
ds_f_eeg <- ds_f %>%
  filter(SESSION %in% c(2, 3, 6))

#Setup for plotting
theme_set(theme_bw(base_size = 20))#32/28
CSCALE_YlOrRd = rev(brewer.pal(9,"YlOrRd"));
#chance<-1/length(unique(ds_b[[modelV]]))

#Settings
saveOn <- T # for whether plots/models, etc. should be saved to REPO/elsewhere


#=RT and Error Rate=============================================================

#Overall
paste("the overall RT for task one was", round(mean(ds_f$RT_1),2), "ms and", round(mean(ds_f$RT_2),2),"ms for task two")
paste("the overall ACC for task one was", round(mean(ds_b$ACC_1),2), "and", round(mean(ds_b$ACC_2),2),"for task two")

#EEG sessions
paste("for EEG sessions, the overall ACC for task one was", round(mean(ds_eeg$ACC_1),2), "and", round(mean(ds_eeg$ACC_2),2),"for task two")
paste("for EEG sessions, the overall RT for task one was", round(mean(ds_f_eeg$RT_1),2), "ms and", round(mean(ds_f_eeg$RT_2),2),"ms for task two")

#plot these!
ds_f %>%
  mutate(eeg_ses = as.factor(ifelse(SESSION %in% c(2, 3, 6), 1, 0))) %>%
  group_by(SUBID, eeg_ses) %>%
  summarise(RT_1 = mean(RT_1),
            RT_2 = mean(RT_2)) %>%
  pivot_longer(cols = c(3,4),
               names_to = "measure",
               values_to = "value") %>%
ggplot( aes(x= measure, y = value))+
  geom_violin(aes(color = eeg_ses))+
  geom_jitter(width = 0.2, aes(color = eeg_ses, group = SUBID)) +
  facet_grid(~eeg_ses) +
  ggtitle("RT_1 and RT_2")

ds_b %>%
  mutate(eeg_ses = as.factor(ifelse(SESSION %in% c(2, 3, 6), 1, 0))) %>%
  group_by(SUBID, eeg_ses) %>%
  summarise(ACC_1 = mean(ACC_1),
            ACC_2 = mean(ACC_2)) %>%
  pivot_longer(cols = c(3,4),
               names_to = "measure",
               values_to = "value") %>%
ggplot(aes(x= measure, y = value))+
  geom_violin(aes(color = eeg_ses))+
  geom_jitter(width = 0.2, aes(color = eeg_ses, group = SUBID)) +
  facet_grid(~eeg_ses) +
  ggtitle("ACC_1 and ACC_2")

#=Performance by task features=====================================
ds_b[,TASKSET_1 := factor(TASKSET_1, levels = c("Bold-Faint","Red-Blue", "Low-High", "Odd-Even"))]
ds_b[,TASKSET_2 := factor(TASKSET_2, levels = c("Bold-Faint","Red-Blue", "Low-High", "Odd-Even"))]
ds_b[,TASKPAIR := as.factor(TASKPAIR)]
ds_b[,TASKPAIRc := paste(TASKSET_1, TASKSET_2, sep = "_")]
ds_f<-ds_b[ACC_1==1 & ACC_2==1,] 

#Is any one task harder than another? RTs first, then ACC
ds_f %>%
  group_by(NUMTASK, TASKSET_1, SESSION) %>%
  summarise(rt_1 = mean(log(RT_1)),
            sem = sd(log(RT_1))/sqrt(sum(!is.na(RT_1)))) %>%
  ggplot(aes(x= TASKSET_1, y = rt_1))+
  geom_line(aes(color = NUMTASK, group = NUMTASK))+
  geom_point(aes(color = NUMTASK, group = NUMTASK), size = 0.5)+
  geom_linerange(aes(ymin = rt_1 - sem, ymax = rt_1 + sem, color= NUMTASK))+
  facet_grid(~SESSION)+
  ggtitle("RTs for TASKSET_1") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
ds_f %>%
    group_by(NUMTASK, TASKSET_2, SESSION) %>%
    summarise(rt_2 = mean(log(RT_2)),
              sem = sd(log(RT_2))/sqrt(sum(!is.na(RT_2)))) %>%
    ggplot(aes(x= TASKSET_2, y = rt_2))+
    geom_point(aes(color = NUMTASK, group = NUMTASK), size = 0.5)+
    geom_line(aes(color = NUMTASK, group = NUMTASK))+
    geom_linerange(aes(ymin = rt_2 - sem, ymax = rt_2 + sem, color= NUMTASK))+
    facet_grid(~SESSION)+
    ggtitle("RTs for TASKSET_2")+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ds_b %>%
  group_by(NUMTASK, TASKSET_1, SESSION) %>%
  summarise(acc_1 = mean(ACC_1),
            sem = sd(ACC_1)/sqrt(sum(!is.na(ACC_1)))) %>%
  ggplot(aes(x= TASKSET_1, y = acc_1))+
  geom_line(aes(color = NUMTASK, group = NUMTASK))+
  geom_point(aes(color = NUMTASK, group = NUMTASK), size = 0.5)+
  geom_linerange(aes(ymin = acc_1 - sem, ymax = acc_1 + sem, color= NUMTASK))+
  facet_grid(~SESSION)+
  ggtitle("ACCs for TASKSET_1") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
ds_b %>%
  group_by(NUMTASK, TASKSET_2, SESSION) %>%
  summarise(acc_2 = mean(ACC_2),
            sem = sd(ACC_2)/sqrt(sum(!is.na(ACC_2)))) %>%
  ggplot(aes(x= TASKSET_2, y = acc_2))+
  geom_point(aes(color = NUMTASK, group = NUMTASK), size = 0.5)+
  geom_line(aes(color = NUMTASK, group = NUMTASK))+
  geom_linerange(aes(ymin = acc_2 - sem, ymax = acc_2 + sem, color= NUMTASK))+
  facet_grid(~SESSION)+
  ggtitle("ACCs for TASKSET_2")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
# looks like in both S1 & S2, Odd-Even is 'hardest' task. Interestingly, different pattern of RTs/ACC depening on if trial is Single or Double task!!

# Now let's look at task combinations (TASKPAIR)
ds_f %>%
  group_by(NUMTASK, TASKPAIRc) %>%
  summarise(rt_2 = mean(log(RT_2)),
            rt_1 = mean(log(RT_1))) %>%
  pivot_longer(cols = c(3,4),
               names_to = "stimnum",
               values_to = "rt") %>%
  ggplot(aes(x= TASKPAIRc, y = rt))+
  geom_point(aes(color = NUMTASK, group = NUMTASK))+
#  geom_linerange(aes(ymin = rt - sem, ymax = rt + sem, color= NUMTASK))+
  facet_grid(~stimnum)+
  ggtitle("RTs for TASKPAIR")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = .5))



#look at number of observations across TASKSETxCORSIDE conj and split across SOAs (any people have very imbalanced?)

ds_f_eeg[NUMTASK == "single"] %>%
  group_by(SUBID,TASKSET_CORSIDE_2c, NUMTASK, SOA) %>% #, SUBID_S
  summarise(T = n(), .groups = "drop") %>%
  group_by(SUBID) %>%
  slice_min(order_by=T, with_ties = TRUE) %>% # returns column with minimum value of variable by specified grouping variable!
  arrange(T)%>%
print(n=32) #%>%
  ggplot(aes(x = TASKSET_CORSIDE_2c, y = SOA)) +
  geom_tile(aes(fill = T)) +
  facet_rep_wrap(~SUBID)
ds_f_eeg %>%
  group_by(TASKSET_CORSIDE_2c, SOA) %>% #, SUBID_S
  summarise(T = n()) %>%
  arrange(T)%>%
  print() %>%
  ggplot(aes(x = TASKSET_CORSIDE_2c, y = SOA)) +
  geom_tile(aes(fill = T))

#=Similar to previous section, but stimulus features============================
ds_b[,S1IDc := as.factor(S1IDc)]
ds_b[,S2IDc := as.factor(S2IDc)]
ds_f[,S1IDc := as.factor(S1IDc)]
ds_f[,S2IDc := as.factor(S2IDc)]


#Is any one stimulus harder than another? RTs first, then ACC 
ds_f %>%
  group_by(NUMTASK, S1IDc, SESSION) %>%
  summarise(rt_1 = mean(log(RT_1)),
            sem = sd(log(RT_1))/sqrt(sum(!is.na(RT_1)))) %>%
  ggplot(aes(x= S1IDc, y = rt_1))+
  geom_line(aes(color = NUMTASK, group = NUMTASK))+
  geom_point(aes(color = NUMTASK, group = NUMTASK), size = 0.5)+
  geom_linerange(aes(ymin = rt_1 - sem, ymax = rt_1 + sem, color= NUMTASK))+
  facet_grid(~SESSION)+
  ggtitle("RTs for S1IDc") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ds_f %>%
  group_by(NUMTASK, S2IDc, SESSION) %>%
  summarise(rt_2 = mean(log(RT_2)),
            sem = sd(log(RT_2))/sqrt(sum(!is.na(RT_2)))) %>%
  ggplot(aes(x= S2IDc, y = rt_2))+
  geom_point(aes(color = NUMTASK, group = NUMTASK), size = 0.5)+
  geom_line(aes(color = NUMTASK, group = NUMTASK))+
  geom_linerange(aes(ymin = rt_2 - sem, ymax = rt_2 + sem, color= NUMTASK))+
  facet_grid(~SESSION)+
  ggtitle("RTs for S2IDc")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ds_b %>%
  group_by(NUMTASK, S1IDc, SESSION) %>%
  summarise(acc_1 = mean(ACC_1),
            sem = sd(ACC_1)/sqrt(sum(!is.na(ACC_1)))) %>%
  ggplot(aes(x= S1IDc, y = acc_1))+
  #geom_line(aes(color = NUMTASK, group = NUMTASK))+
  geom_point(aes(color = NUMTASK, group = NUMTASK), size = 0.5)+
  geom_linerange(aes(ymin = acc_1 - sem, ymax = acc_1 + sem, color= NUMTASK))+
  facet_grid(~SESSION)+
  ggtitle("ACCs for S1IDc") +
  ylim(0.9,1)
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ds_b %>%
  group_by(NUMTASK, S2IDc, SESSION) %>%
  summarise(acc_2 = mean(ACC_2),
            sem = sd(ACC_2)/sqrt(sum(!is.na(ACC_2)))) %>%
  ggplot(aes(x= S2IDc, y = acc_2))+
  geom_point(aes(color = NUMTASK, group = NUMTASK), size = 0.5)+
  #geom_line(aes(color = NUMTASK, group = NUMTASK))+
  geom_linerange(aes(ymin = acc_2 - sem, ymax = acc_2 + sem, color= NUMTASK))+
  facet_grid(~SESSION)+
  ggtitle("ACCs for S2IDc")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

#These plots are super messy, but think that while there are 1 or 2 outliers, acc is not SUUUUper different across stim



#=Multitasking Costs===========================================================


#plotting multitasking costs (RT & errors) across SOAs
ds_f %>% #ds_m <- 
  group_by(NUMTASK, SOA)%>%
  summarise(rt_1 = mean(RT_1),
            #sem1 = sd(log(RT_1))/sqrt(sum(!is.na(RT_1))),
            rt_2 = mean(RT_2)
            #sem2= sd(log(RT_2))/sqrt(sum(!is.na(RT_2)))
            ) %>%
  pivot_longer(cols = c(3,4),
               names_to = "stimnum",
               values_to = "rt") %>%
  ggplot(aes(x = SOA, y=rt)) +
    geom_line(aes(color = NUMTASK)) +
    geom_point(aes(color = NUMTASK)) +
    facet_wrap(~stimnum)+
  ggtitle("RT x SOAs")

ds_eeg %>% #ds_m <- 
  group_by(NUMTASK, SOA)%>%
  summarise(acc_1 = mean(ACC_1),
            #sem1 = sd(log(RT_1))/sqrt(sum(!is.na(RT_1))),
            acc_2 = mean(ACC_2)
            #sem2= sd(log(RT_2))/sqrt(sum(!is.na(RT_2)))
  ) %>%
  pivot_longer(cols = c(3,4),
               names_to = "stimnum",
               values_to = "acc") %>%
  ggplot(aes(x = SOA, y=acc)) +
  geom_line(aes(color = NUMTASK)) +
  geom_point(aes(color = NUMTASK)) +
  facet_wrap(~stimnum) +
  ggtitle("ACC x SOAs")

# interesting pattern of results here:
# - SOA effects at shortest SOA (100ms)
# - participants are *faster* in responding the second task & seem to plateau RT_1
# - looking at ACC_1, they are also reaching ceiling at long SOAs, so maybe no need to 'try harder'
# - suggests strategic deployment of attention/resources? interesting to think how this may play out in S1 decoding
  
#looking at multitasking costs over sessions
ds_sub <- ds_f %>% #ds_m <- 
  group_by(NUMTASK, SESSION, SOA, SUBID)%>%
  summarise(sub_rt = mean(RT_2),.groups = "drop")

if(saveOn == T){outfile <- file.path(paste0(Dir_REPO,"Figures&Stats/BEH/BEH(6sessions)/", "RT_2_SOA_plot.tiff"));tiff(outfile, width = 12, height = 6, units = "in", res = 300)}
ds_sub %>% #ds_m <- 
  group_by(NUMTASK, SESSION, SOA)%>%
  summarise(rt = mean(sub_rt), #rt_1 = mean(RT_1), sem1 = sd(log(RT_1))/sqrt(sum(!is.na(RT_1))),
            sem= sd(sub_rt)/sqrt(n()), .groups = "drop"
            ) %>%
  ggplot(aes(x = SOA, y=rt, group = NUMTASK)) +
  geom_point(data = ds_sub, aes(x = SOA, y=sub_rt, group = SUBID, color = NUMTASK), alpha = 0.2, size = 1.5)+
  geom_line(data = ds_sub, aes(x = SOA, y=sub_rt, color = NUMTASK, group = interaction(SUBID, NUMTASK)), alpha = 0.2, linewidth = 0.5)+
  geom_line(aes(color = NUMTASK), linewidth = 1.5) +
  geom_point(aes(color = NUMTASK), size = 3) +
  facet_grid(~SESSION)+
  geom_linerange(aes(ymin = rt-sem, ymax = rt+sem, color= NUMTASK), linewidth = 1)+
  #ggtitle("Task 2 RT over Sessions") +
  xlab("SOA (ms)") + ylab("RT of 2nd response (ms)") +
  scale_color_manual(values=c("cornflowerblue", "orange")) + 
  theme_classic() +
  theme(legend.title = element_blank())
dev.off();rm(outfile)

#error rate!
ds_err <- ds_b %>% #ds_m <- 
  group_by(NUMTASK, SESSION, SOA, SUBID)%>%
  summarise(sub_err = mean(1- ACC_2))
if(saveOn == T){outfile <- file.path(paste0(Dir_REPO,"Figures&Stats/BEH/BEH(6sessions)/", "ERR_2_SOA_plot.tiff"));tiff(outfile, width = 12, height = 6, units = "in", res = 300)}
ds_err %>% #ds_m <- 
  group_by(NUMTASK, SESSION, SOA)%>%
  summarise(err = mean(sub_err),
            sem= sd(sub_err)/sqrt(n())
  ) %>%
  ggplot(aes(x = SOA, y=err, group = NUMTASK)) +
  geom_point(data = ds_err, aes(x = SOA, y=sub_err, group = SUBID, color = NUMTASK), alpha = 0.2, size = 1.5)+
  geom_line(data = ds_err, aes(x = SOA, y=sub_err, color = NUMTASK, group = interaction(SUBID, NUMTASK)), alpha = 0.2, linewidth = 0.5)+
  geom_line(aes(color = NUMTASK), linewidth = 1.5) +
  geom_point(aes(color = NUMTASK), size = 4) +
  facet_grid(~SESSION)+
  geom_linerange(aes(ymin = err-sem, ymax = err+sem, color= NUMTASK), linewidth = 1)+
  #ggtitle("Task 2 RT over Sessions") +
  xlab("SOA (ms)") + ylab("Error rate") +
  scale_color_manual(values=c("cornflowerblue", "orange")) + 
  theme_classic()+
  theme(legend.title = element_blank())
dev.off();rm(outfile)



#model RT relationships
ds_m <- ds_f %>% #create data for model
  filter(SESSION %in% c(2, 3, 6))%>% #only look at EEG days... is this right?
  #filter(SOA %in% c(250, 500, 750, 1000)) %>% #remove shortest SOA
setDT
ds_m[,SOA := as.factor(SOA)]
ds_m[,NUMTASK := factor(NUMTASK, levels = c("single","double"))]
ds_m[,SESSION := as.factor(SESSION)]
ds_m[,SESSION:= droplevels(ds_m$SESSION)]
contrasts(ds_m$SESSION, how.many = 2) <- contr.poly(3)
contrasts(ds_m$SOA, how.many = 2) <- contr.poly(5)

m1 <- lmer(rt_1 ~ SESSION*NUMTASK*SOA + (1|SUBID_S), data = ds_m)
summary(m1)
m2 <- lmer(RT_2 ~ SESSION*SOA + (1|SUBID_S), data = ds_m[NUMTASK=="double"])
summary(m2)

ds_100 <- ds_m[SOA == 100]
ds_250 <- ds_m[SOA == 250]
ds_500 <- ds_m[SOA == 500]
ds_750 <- ds_m[SOA == 750]
ds_1000 <- ds_m[SOA == 1000]



  #m1 <- lmer(RT_2 ~ SESSION*NUMTASK + (1+ SESSION*NUMTASK|SUBID_S), data = ds_1000) #maximally complex model
  #m2 <- lmer(RT_2 ~ SESSION*NUMTASK + (1+ SESSION+NUMTASK|SUBID_S), data = ds_1000) #less
  m3 <- lmer(RT_2 ~ SESSION*NUMTASK + (1|SUBID_S), data = ds_1000) #least
  
  setwd(paste0(Dir_REPO,"Figures&Stats/BEH/BEH(6sessions)/plots&models/"))
 #for saving summary output as .txt file
  fileN <- "RT_2_SOA_1000_simple_summary.txt" #"VI" for variable interaction (m1); "FI" for fixed interaction (m2); "simple" for no random slopes
  sink(fileN)
  print(summary(m3))
  sink()
  
#for saving rds object of model 
  saveRDS(m3, file = "RT_2_SOA_1000_simple_model.rds")
  
  
#for saving plots of effects
  
  outfile <- file.path(paste0(Dir_REPO,"Figures&Stats/BEH/BEH(6sessions)/plots&models/", "RT_2_SOA_1000_simple_plot.png"))
  png(outfile, width = 1920, height = 1080, units = "px")
  # plot(allEffects(m1), main = deparse(formula(m1)))
  # plot(allEffects(m2), main = deparse(formula(m2)))
   plot(allEffects(m3), main = deparse(formula(m3)))
  dev.off()
  

# m_500 <- lmer(RT_2 ~ SESSION*NUMTASK + (1+ SESSION+NUMTASK|SUBID_S), data = ds_500) 
# summary(m_SOA)
# plot(allEffects(m_500))

#plot_model(m1, type= "int")

#=Let's look at the behavior relating to geometry questions ====================
# Basically, partial overlap in representations (e.g.,response congruency, stim feature, task repeat (which is just single task trials)) should tell us
# something about whether people are reusing representations (expressive geometry) or rapidly switching geometries
# 
# Take, for example, a multitasking trial:
# - cued tasks are High/Low & Red/Blue
# - stimulus is a Red "7", so responses are congruent
# - as soon as 1st response is made, the 'High/Low' feature is no longer relevant
  
# Looking at RT_2 in multitask trials where responses are congruent
 
ds_f_eeg[,RESPCONG := factor(RESPCONG)]
  
#subject means, for SEM error bars
ds_sub_cong <- ds_f_eeg %>% #plot single-task trials, for reference. 
  group_by(RESPCONG, SESSION, SUBID, NUMTASK)%>%
  summarise(rt_2 = mean(RT_2)) %>% setDT()
    
ds_sub_cong[NUMTASK == 'single'] %>% #plot single-task trials, for reference. 
    group_by(RESPCONG, SESSION) %>%
    summarise(rt = mean(rt_2),
              sem = sd(rt_2)/sqrt(n())) %>%
  ggplot(aes(x= RESPCONG, y = rt))+
  geom_line(aes(color = SESSION, group = SESSION))+
   geom_point(aes(color = SESSION, group = SESSION), size = 0.5)+
   geom_linerange(aes(ymin = rt - sem, ymax = rt + sem, color= SESSION))+
  ggtitle("RT_2 by response congruency in single and dualtask trials")+ 
ds_sub_cong[NUMTASK == 'double'] %>% 
  group_by(RESPCONG, SESSION) %>%
  summarise(rt = mean(rt_2),
            sem = sd(rt_2)/sqrt(n())) %>%
  ggplot(aes(x= RESPCONG, y = rt))+
  geom_line(aes(color = SESSION, group = SESSION))+
  geom_point(aes(color = SESSION, group = SESSION), size = 0.5)+
  geom_linerange(aes(ymin = rt - sem, ymax = rt + sem, color= SESSION)) 

#[insert code here to save plots to repo]

#interesting pattern of results here... evidence of congruency-induced interference specific to multitasking
# next is to look into what task features are related to this. stimulus features? tasksets? then look @ decoding

#let's look at how stimulus feature switches affect response congruency effects on RT_2 in single task v. dual task trials
# start with NUMSW, which is # of stimulus feature switches from S1 to S2 (min 0 max 4)
ds_sub_cong_sw <- ds_f_eeg %>% #plot single-task trials, for reference. 
  group_by(RESPCONG, SESSION, SUBID, NUMTASK, NUMSW)%>%
  summarise(rt_2 = mean(RT_2)) %>% setDT()

ds_sub_cong_sw[NUMTASK == 'single'] %>% #plot single-task trials, for reference. 
  group_by(RESPCONG, SESSION, NUMSW) %>%
  summarise(rt = mean(rt_2),
            sem = sd(rt_2)/sqrt(n())) %>%
  ggplot(aes(x= RESPCONG, y = rt))+
  geom_line(aes(color = SESSION, group = SESSION))+
  geom_point(aes(color = SESSION, group = SESSION), size = 0.5)+
  geom_linerange(aes(ymin = rt - sem, ymax = rt + sem, color= SESSION))+
  facet_grid(~NUMSW) +
  ggtitle("RT_2 by response congruency*number of feature switches (single-task)")
ds_sub_cong_sw[NUMTASK == 'double'] %>% 
  group_by(RESPCONG, SESSION, NUMSW) %>%
  summarise(rt = mean(rt_2),
            sem = sd(rt_2)/sqrt(n())) %>%
  ggplot(aes(x= RESPCONG, y = rt))+
  geom_line(aes(color = SESSION, group = SESSION))+
  geom_point(aes(color = SESSION, group = SESSION), size = 0.5)+
  geom_linerange(aes(ymin = rt - sem, ymax = rt + sem, color= SESSION)) +
  facet_grid(~NUMSW) +
  ggtitle("RT_2 by response congruency*number of feature switches (dual-task)")
# now this second plot is super interesting!!! interaction of stimulus feature switches and response congruency on RT_2 in dual task trials!
# such that behavior is 'best' when: 0 feature switches (same simuli), response is repeated AND 4 feature switches (no feature overlap), response is different!!
# this is an example of partial overlap costs in RT_2? evidence for conjunctive representations?
# though ALL of these include a task switch. 


#Shared stimulus feature lane:
  # one hypothesis (M&C) is that early on, shared representations facilitate learning but also promote costs/interference when both representations must be leveraged at once
  # Presumably, stimulus features (number & color) are represented similarly between tasks. Let's look at behavior relating to this

  #RT costs when there is a feature dimension switch (e.g., number -> color or color -> number) vs. not (color -> color)

ds_f<- ds_f %>%
  mutate(FDIM_SW = as.factor(ifelse(FDIM_1 == FDIM_2, 0, 1))) # new column to code if there is a feature dimension switch between T1 and T2
ds_f[NUMTASK == "double"] %>%
  group_by(FDIM_SW, SESSION) %>%
  summarise(rt_2 = mean(RT_2))%>%
  ggplot(aes(x=SESSION, y = rt_2, group = FDIM_SW)) + geom_line(aes(color = FDIM_SW))

#cool let's model it
ds_m <- ds_f %>% #create data for model
  filter(SESSION %in% c(2, 3, 6) & NUMTASK == 'double')%>%
  setDT
ds_m[,SOA := as.factor(SOA)]
ds_m[,SESSION := as.factor(SESSION)]
ds_m[,SESSION:= droplevels(ds_m$SESSION)]
contrasts(ds_m$SESSION, how.many = 2) <- contr.poly(3)
contrasts(ds_m$SOA, how.many = 2) <- contr.poly(5)
m <- lmer(RT_2 ~ SOA*SESSION*FDIM_SW + (1|SUBID_S), data = ds_m)
summary(m)
ds_m %>%
  group_by(SUBID,SESSION, SOA, FDIM_SW)%>%
  summarize(rt = mean(RT_2),
            sd = sd(RT_2))%>% ungroup()%>%
  group_by(SESSION, SOA, FDIM_SW)%>%
  summarize(RT_2 = mean(rt),
            sem = sd(sd)/sqrt(sum(!is.na(rt))))%>% ungroup()%>%
ggplot(aes(x = SESSION, y = RT_2, group = FDIM_SW)) + 
  geom_line(aes(color = FDIM_SW)) + geom_point(aes(color = FDIM_SW))+
  geom_linerange(aes(ymin = RT_2 -sem, ymax = RT_2+sem, color =FDIM_SW ))+
  facet_grid(~SOA)

#checking that response-locked data is created correctly------------
  #should be (-700, 200) relative to response
  # EV RESP
setwd(Dir_BDATA)
ds_bb<-fread("DimensionSwitchPRP_BehP.txt") 




