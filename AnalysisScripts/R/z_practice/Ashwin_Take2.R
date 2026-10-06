#ASHWIN PRP Dimension Switch

# =============================================================

#clear environment
graphics.off()
rm(list = ls(all = TRUE))

#Setting up directories
fsep<-.Platform$file.sep;
Dir_Support <- path.expand("~/Brown Dropbox/Ashwin Srinivasan/R_support_scripts") 
Dir_Experiment <- path.expand("~/Brown Dropbox/Ashwin Srinivasan/DimensionSwitch_PRP")
Dir_BDATA<-path.expand("DataAnalysis/Matlab/BEH")


#load libraries
library(data.table)
library(tidyverse)
library(ggplot2)


#load data
setwd(Dir_Experiment)
setwd(Dir_BDATA)
df <- fread('z_ALLGRAND(EEG,6sessions)/DimensionSwitchPRP_BehP.txt', header=TRUE, fill=TRUE) 


#Clean and Recode Data
taskLabels <- c("Low-High", "Odd-Even", "Red-Blue", "Bold-Faint")
stimNumbers <- c("3", "4", "6", "7")
stimColor <- c("red", "red-faint", "blue-faint", "blue")

df[,SOA:= factor(SOA, labels = c("100ms","250ms", "500ms", "750ms", "1000ms"))]
df[,TASKSET_1:= factor(TASKSET_1, labels = taskLabels)]
df[,TASKSET_2:= factor(TASKSET_2, labels = taskLabels)]
df[,STIMPOS_1:= factor(STIMPOS_1, labels = c("Top", "Bottom"))]
df[,STIMPOS_2:= factor(STIMPOS_2, labels = c("Top", "Bottom"))]
df[,FDIM_1:= factor(FDIM_1, labels = c("Number", "Color"))]
df[,FDIM_2:= factor(FDIM_2, labels = c("Number", "Bottom"))]

df[,F_NUMBER_1:= factor(F_NUMBER_1, labels = stimNumbers)]
df[,F_NUMBER_2:= factor(F_NUMBER_2, labels = stimNumbers)]
df[,S1_1:= factor(S1_1, labels = stimNumbers)]
df[,S2_1:= factor(S2_1, labels = stimNumbers)]

df[,F_COLOR_1:= factor(F_COLOR_1, labels = stimColor)]
df[,F_COLOR_2:= factor(F_COLOR_2, labels = stimColor)]
df[,S1_2:= factor(S1_2, labels = stimColor)]
df[,S2_2:= factor(S2_2, labels = stimColor)]

df[,NUMTASK:= factor(NUMTASK, labels = c("single","double"))]
df[,RESPCONG:= factor(RESPCONG)]

df$SUBID_S <- df$SUBID
df <- transform(df, SUBID = paste0("1",substr(SUBID, 2, 3)))
df[,EXACT_REP:= factor(EXACT_SW, labels = c("exact","diff"))]

df$S1IDc <- paste(df$S1_1, df$S1_2, sep = "_" )
df$S2IDc <- paste(df$S2_1, df$S2_2, sep = "_" )


#remove bad trials (flips, long RT, etc.)
cutRT1i <- df[ACC_1 == 1 & df$RT_1 < 10000,list(cutRT1 = quantile(RT_1, 0.99)),by = c("SUBID","SESSION")]
cutRT2i <- df[ACC_2 == 1 & df$RT_2 < 10000,list(cutRT2 = quantile(RT_2, 0.99)),by = c("SUBID","SESSION")]
df <- merge(df, cutRT1i,by=c("SESSION","SUBID"))
df <- merge(df, cutRT2i,by=c("SESSION","SUBID"))
print(sprintf("Response cutoff %.3f for Response 1",mean(cutRT1i$cutRT)))
print(sprintf("Response cutoff %.3f for Response 2",mean(cutRT2i$cutRT)))

# Filtering
# ds<-ds[BLOCK>0,]
df<-df[RESPOMIT==0] # Exclude response omission trials?
df_roall <-df # keep this for later analysis
df<-df[RESPORDER_K==1]# Exclude response flipping trials?
df<-df[!(c(df$RT_1 > df$cutRT1)),];
df<-df[!(c(df$RT_2 > df$cutRT2)),];

#recreate Atsushi's lab meeting presentation

#1 - Practice induced changes in multitasking costs: comparing dual task and single task performance (Task 2 RT) across SOAs across sessions

RT_2s <- df%>%
  filter(ACC_1 ==1 & ACC_2==1) %>%
  group_by(NUMTASK, SOA, SESSION) %>%
  summarise(rt2 = mean(RT_2), sd = sd(RT_2)) %>%
  rename(task_condition = NUMTASK)

RT2_plot <- df %>%
  filter(ACC_1 ==1 & ACC_2==1) %>%
  group_by(NUMTASK, SOA, SESSION,SUBID) %>%
  summarise(rt2 = mean(RT_2), sd = sd(RT_2)) %>%
  ggplot(aes(x = SOA, y = rt2)) + geom_point(aes(color = NUMTASK), alpha = 0.1)+ facet_grid(~SESSION) + geom_line(alpha = 0.1, aes(color = NUMTASK, group=interaction(SUBID, NUMTASK)))
RT2_plot + geom_point(data = RT_2s, aes(x = SOA, y = rt2, color = task_condition), size = 2) + 
  geom_line(data = RT_2s, aes(x = SOA, y = rt2, group = task_condition, color = task_condition)) + 
  scale_color_manual(values = c("black", "red")) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1))   

#2 Practice induced changes in output interference (i.e., looking at response (side) repetition benfits)
#why? what does this tell us
# - i guess the idea is that task switching causes output interference, in that repetition benefits present in the single-task condition are absent in the dual-task condition
# - one explanation for this is that switch costs negate the repetition benefit

RespRep_RT2_plot <- df %>%
  filter(ACC_1 ==1 & ACC_2==1) %>%
  group_by(NUMTASK, SOA, SESSION,RESPCONG) %>%
  summarise(rt2 = mean(RT_2), sd = sd(RT_2)) %>%
  ggplot(aes(x = SOA, y = rt2)) + geom_point(aes(color = NUMTASK), alpha = 1)+ facet_grid(~SESSION) + geom_line(alpha = 1, aes(group = interaction(RESPCONG,NUMTASK), linetype = RESPCONG, color = NUMTASK))
RespRep_RT2_plot +
  scale_color_brewer(values = c("black", "red")) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1))   
