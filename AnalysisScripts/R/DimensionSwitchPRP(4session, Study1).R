# DimensionSwitchPRP(4 session, Study1)
# To do 
# 1. think of more analyses for switching of features!
# 2. Re-code information in meaningful manner
# 3. Exclude flipped order responses...(how to detect this??)!
# =============================================================
#Removing evrything from workspace
# graphics.off()
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
Dir_R<-path.expand("~/Dropbox/w_ONGOINGRFILES/w_OTHERS")
Dir_BDATA<-path.expand("~/Dropbox/w_SCRIPTS/P.HAL2017/DimensionSwitch_PRP/DataAnalysis/Matlab/BEH/z_ALLGRAND(4sessions)")
Dir_EDATA<-path.expand("~/Dropbox/w_SCRIPTS/P.HAL2017/DimensionSwitch_PRP/DataAnalysis/Matlab/EEG")

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

# Change NaNs to NAs
ds[is.nanM(ds)] <- NA;
ds <- ds[, !names(ds) %like%  "EV_", with=F]

# Recode variables
ds[,SUBID_S:=SUBID]
ds[,SUBID:=SUBID_S - 100*(SESSION-1)]
ds[,NUMTASK := factor(NUMTASK, labels = c("single","double"))]
ds[,RESPCONG :=factor(RESPCONG,labels = c("incong","cong"))]
ds[,EXACT_REP:= factor(EXACT_SW, labels = c("exact","diff"))] # flip of EXACT_SW
ds[,FDIM_1:= factor(FDIM_1, labels = c("Number","Color"))]
ds[,FDIM_2:= factor(FDIM_2, labels = c("Number","Color"))]
ds[,TASKSET_1:= factor(TASKSET_1, labels = c("Low-High","Odd-Even","Red-Blue","Bold-Faint"))]
ds[,TASKSET_2:= factor(TASKSET_2, labels = c("Low-High","Odd-Even","Red-Blue","Bold-Faint"))]
ds[, `:=` (SOAn = SOA, SOA = factor(SOA, labels = c("100","250","500","750","1000")))]
contrasts(ds$SOA, how.many=2) <- contr.poly(5)

# Add task variables
# FDIM_MATCH (only switch blocks): 
# Match = switches between same features
# Mismatch = switches between different features
ds[,FDIM_MATCH:=NA]
ds[,FDIM_MATCH:=ifelse(NUMTASK=="double" & FDIM_1 == FDIM_2, "Match",FDIM_MATCH)]
ds[,FDIM_MATCH:=ifelse(NUMTASK=="double" & FDIM_1 != FDIM_2, "Mismatch",FDIM_MATCH)]

# Add stimulus variables
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


# RT cutoff (exclude outlier like RTs > 10 sec)
cutRT1i <- ds[ACC_1 == 1 & ds$RT_1 < 10000,list(cutRT1 = quantile(RT_1, 0.99)),by = c("SUBID","SESSION")]
cutRT2i <- ds[ACC_2 == 1 & ds$RT_2 < 10000,list(cutRT2 = quantile(RT_2, 0.99)),by = c("SUBID","SESSION")]
ds <- merge(ds, cutRT1i,by=c("SESSION","SUBID"))
ds <- merge(ds, cutRT2i,by=c("SESSION","SUBID"))
print(sprintf("Response cutoff %.3f for Response 1",mean(cutRT1i$cutRT)));
print(sprintf("Response cutoff %.3f for Response 2",mean(cutRT2i$cutRT)));

# Checking
rt<-ds_f %>% group_by(SUBID,SESSION) %>% summarise(RT_1=mean(RT_1,na.rm=TRUE),RT_2=mean(RT_2,na.rm=TRUE));as.data.frame(rt)
acc<-ds %>% group_by(SUBID,SESSION) %>% summarise(ACC_1=mean(ACC_1,na.rm=TRUE), ACC_2=mean(ACC_2,na.rm=TRUE));as.data.frame(acc)
numT<-ds %>% group_by(SUBID,SESSION) %>% summarise(T=max(TRIAL));as.data.frame(numT)
numdata<-ds %>% group_by(SUBID,SESSION) %>% summarise(T=n());as.data.frame(numdata)
resporder<-ds %>% group_by(SUBID,SESSION, SOA) %>% summarise(RORDER=mean(RESPORDER_K,na.rm=T));as.data.frame(resporder)


# Filtering
ds<-ds[BLOCK>0  & TRIAL >1,]
ds<-ds[RESPOMIT==0] # Exclude response omission trials?
ds_roall <-ds # keep this for later analysis
ds<-ds[RESPORDER_K==1]# Exclude response flipping trials?
ds<-ds[!(c(ds$RT_1 > ds$cutRT1)),];
ds<-ds[!(c(ds$RT_2 > ds$cutRT2)),];
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

# RESPORDER x SOA
RORDER_SOAi <- ds_roall %>% group_by(SUBID,SESSION,NUMTASK,SOA) %>% summarise(DV=100*(1-mean(RESPORDER_K, na.rm=TRUE)) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
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
RT1_SOAi <- ds_f[RESPORDER_K==1] %>% group_by(SUBID,SESSION,NUMTASK,SOA) %>% summarise(DV=mean(RT_1, na.rm=TRUE) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
RT1_SOA <- RT1_SOAi %>% group_by(SESSION,NUMTASK,SOA) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# ACC_1:SESSION x SOA
ERR1_SOAi <- ds %>% group_by(SUBID,SESSION,NUMTASK,SOA) %>% summarise(DV=100*(1-mean(ACC_1,na.rm=TRUE)) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
ERR1_SOA <- ERR1_SOAi %>% group_by(SESSION,NUMTASK,SOA) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# *************
# Output (response interference)
# *************
# RT_2:SESSION x SOA X RESPCONG (mainly in dual task)
RT2_RCONGi <- ds_f %>% group_by(SUBID,SESSION,NUMTASK,SOA,RESPCONG) %>% summarise(DV=mean(RT_2, na.rm=TRUE) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
RT2_RCONG <- RT2_RCONGi %>% group_by(SESSION,NUMTASK,SOA,RESPCONG) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# ACC_2:SESSION x SOA X RESPCONG
ERR2_RCONGi <- ds %>% group_by(SUBID,SESSION,NUMTASK,SOA,RESPCONG) %>% summarise(DV=100*(1-mean(ACC_2,na.rm=TRUE)) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
ERR2_RCONG <- ERR2_RCONGi %>% group_by(SESSION,NUMTASK,SOA,RESPCONG) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()


# *************
# Feature overlap (color-color/num-num task pairs vs. others)
# *************
RT2_FDIMi <- ds_f[NUMTASK=="double"] %>% group_by(SUBID,SESSION,NUMTASK,SOA,FDIM_MATCH) %>% summarise(DV=mean(RT_2, na.rm=TRUE) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
RT2_FDIM <- RT2_FDIMi %>% group_by(SESSION,NUMTASK,SOA,FDIM_MATCH) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()

# ACC_2:SESSION x SOA X RESPCONG
ERR2_FDIMi <- ds[NUMTASK=="double"] %>% group_by(SUBID,SESSION,NUMTASK,SOA,FDIM_MATCH) %>% summarise(DV=100*(1-mean(ACC_2,na.rm=TRUE)) ,nn=n()) %>% normWS("SUBID","DV") %>% print()
ERR2_FDIM <- ERR2_FDIMi %>% group_by(SESSION,NUMTASK,SOA,FDIM_MATCH) %>% summarize(DVm=mean(DV),SD=sd(DV_n),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n()) %>% data.table() %>% print()


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
m_rt=lmer(RT_2~1+SESSION*SOA*FDIM_MATCH+(1|SUBID),data=ds_f[SESSION!=1 & NUMTASK=="double",]);std_beta_MLM(m_rt)@fixef

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
theme_set(theme_bwL(base_size = 14))#32/28
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


# <Multitasking costs>
# rw_DV<-RT2_SOA;rw_DVi<-RT2_SOAi;
# rw_DV<-ERR2_SOA;rw_DVi<-ERR2_SOAi;
# rw_DV<-RT1_SOA;rw_DVi<-RT1_SOAi;
# rw_DV<-ERR1_SOA;rw_DVi<-ERR1_SOAi;

# <Response congruency effect>
rw_DV<-RT2_RCONG;rw_DVi<-RT2_RCONGi;
# rw_DV<-ERR2_RCONG;rw_DVi<-ERR2_RCONGi;

# <Feature dimension overlap effect>
# rw_DV<-RT2_FDIM;rw_DVi<-RT2_FDIMi;
# rw_DV<-ERR2_FDIM;rw_DVi<-ERR2_FDIMi;

# <Exact repeat of stimulus>
# rw_DV<-RT2_EXACT;rw_DVi<-RT2_EXACTi;
# rw_DV<-ERR2_EXACT;rw_DVi<-ERR2_EXACTi;

# <Response order>
# rw_DV<-RORDER_SOA;rw_DVi<-RORDER_SOAi;

# Individual's improvement
# Basic results
quartz(width=8.5,height=4.70) # RT
# quartz(width=3.5,height=2) # ERR
ggplot() + 
  #Panels ----------------------------
  facet_grid(~SESSION)+
  #geom_hline(yintercept = rw_DV[NUMTASK=="double" & SOA==100 & SESSION==2]$DVm, linetype=3, alpha = .8, color="red")+
  #geom_hline(yintercept = rw_DV[NUMTASK=="single" & SOA==100 & SESSION==1]$DVm, linetype=3, alpha = .8)+
  #Main layer-------------------------
  geom_errorbar(data=rw_DV,aes(x=SOA,ymin=DVm-wseb, ymax=DVm+wseb,group=interaction(NUMTASK,RESPCONG),color=NUMTASK),width=.0,size=.75)+#,color="black"
  #geom_point(data=rw_DVi,aes(x=SOA,y=DV,group=interaction(SUBID,NUMTASK),color=NUMTASK),size=1.5,alpha=0.10)+
  #geom_line(data=rw_DVi,aes(x=SOA,y=DV,group=interaction(SUBID,NUMTASK),color=NUMTASK,linetype=RESPCONG),size=0.5,alpha=0.10) +
  geom_line(data=rw_DV,aes(x=SOA,y=DVm,group=interaction(NUMTASK,RESPCONG),color=NUMTASK,linetype=RESPCONG),size=1) +
  geom_point(data=rw_DV,aes(x=SOA,y=DVm,group=interaction(NUMTASK,RESPCONG),color=NUMTASK),size=3.5) +
  #Axis & label-----------------------
  coord_cartesian(ylim=c(375,2200))+scale_y_continuous(breaks=seq(0,3000,500))+
  #coord_cartesian(ylim=c(0,32))+scale_y_continuous(breaks=seq(0,30,10))+
  ylab("Task 2 RT (ms)") + xlab("SOA (ms)")+
  #ylab("Task 2 Error (p)") + xlab("SOA (ms)")+
  #ylab("Order Error (%)") + xlab("SOA (ms)")+
  #Colors ----------------------------
  theme(axis.text.x = element_text(angle = 45, vjust = 1.08, hjust = 1),
        legend.position = "none") + 
  scale_color_manual(values=c("black","red"))


# Make a graph of response order swapping!!


