# DimensionSwitchPRP_AggregateDB
# load pred file (contains single-trial accuracy & confidence)
# then produce aggregated data for specified conditions (with few attributes)
# SOA Realignment
# 1. add timeLK_SOA column, using correct timeL 
# 2. Use following settings
# timeL<-seq(-400,350,4);
# timeS<-seq(-400,350,20);
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))
gc()

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
Dir_RSAMODEL<-paste0(Dir_EDATA,"/w_RSAMODELS")
Dir_REPO <- path.expand("/oscar/home/asrini17/data/asrini17/DimensionSwitch_PRP/PRP_DimensionSwitch/figures")

# Load libraries
library(data.table)
#library(tidyverse)
library(tidyr)
library(dplyr)
library(broom)
library(rhdf5)
library(foreach)
library(zoo)
library(doMC)
library(binhf)
library(RColorBrewer)
library(RSQLite)
library(lemon)
library(DBI)
library(lme4)
library(sjPlot)
library(cluster)
#library(dbplyr)

# Source Files
setwd(Dir_R)
source('basic_lib.R')
source('basic_theme.R')

#Define basic info
saveOn<-T;# save results or not
fileN<-"test"
#timeL<-seq(-200,700,4); # labels of time for time indexes
#timeS<-seq(-200,700,12); # labels of time you want

#if sessN == "nb_S2":
#timeL<-seq(-500,700,4); # labels of time for time indexes
#timeS<-seq(-500,700,12);

#if evLock == "RESP":
timeL<-seq(-700,200,4); # labels of time for time indexes
timeS<-seq(-700,200,12);

refRLV <- NA;#NA, SOA, Realignment: re-epoch data referencing refRLV timing
srate<-4
if (!is.na(refRLV)){adj<-which(timeL==0)}else{adj<-0}
cutF <- function(t){as.numeric(cut(t,breaks=timeS,include.lowest=T))}
idxF<-function(i,adj) which.min(abs(timeL-i))-adj # normal

#load in behavioral data in order to include task relevance - test
setwd(Dir_BDATA)
ds_beh<-fread("DimensionSwitchPRP_BehPP.txt"); #loading in in order to get TASKSET_2/TASKSET_1
ds_beh[is.nanM(ds_beh)]<-NA


# Database= PP: REG / RSAFits: FIT_RSA_TSCONJ,FIT_RSA_L1TSCONJ,FIT_RSA_TSCONJ_RL,FIT_RSA_TCTSCONJ
setwd(Dir_GRAND)
expN<-"DimensionSwitchPRP";
modelV <- "TASKSET_CORSIDE_2"
sessN <- "nb_S2"
regs <- "" #include RT/ACC nuisance regressors? leave blank if no
tableN<-"FIT_RSA_TEST_CTRR";# 
rsaValueN <- "betas"
evLock <- "RESP"
decodingType <- "unmerged"

# Dependent variables
# depV<-c("ACCRESP_acc","RESP_ED_acc");# REG
#depV<-sprintf(c("LH_RSA%s","OE_RSA%s","BF_RSA%s","RB_RSA%s"),"");#
depV<-sprintf(c("TASKSET_RSA%s","CORSIDE_RSA%s", "CONJ_RSA%s"),"")

# Grouping variables (keep first two columns fixed!)
grpV<-c("SUBID","time","SESSION","NUMTASK", "SOA") #, "SOA"
taskV <- c("TASKSET_1","TASKSET_2") #"TASKSET_1", 

# Open connections
dbCon=dbConnect(dbDriver("SQLite"),paste(expN,decodingType,modelV,evLock,regs,tableN, sep = "_"))# Open database connection
dscheck=as.data.table(dbGetQuery(dbCon,paste0('SELECT * FROM ',tableN,'  LIMIT 10000')));
dbCL<-dbListFields(dbCon,dbListTables(dbCon))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Realignment: add timeLK_[event] variable, idx in same unit as time
# timeL will be arbitrary length, and 0 = onset of event
# it willadd timeLK_[ev] column in the dataset (when it does not exist)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# if (!is.na(refRLV)){
#   # Define RL eqn and update few variables for later process
#   timeRLV<-paste0("timeLK_",refRLV)
#   grpV<-gsub("time",timeRLV,grpV);
# 
#     # update original database (only when it does not exist before)
#   if (!any(dbCL %in% timeRLV)){
#     zeroIDX <- max(which(timeL<=0))
#     tleqn<-sprintf("-1 *(floor(%s/%d) + %d - time)",refRLV,srate,zeroIDX);
#     #eqn<-function(RT,time){-1 *(round(RT/srate) - time + adjIDX)}
#     print(paste0('Realignment is requested:',timeRLV," = ",tleqn))
#     dbExecute(dbCon, sprintf("ALTER TABLE %s ADD COLUMN %s NUMBER",tableN,timeRLV))
#     dbExecute(dbCon,sprintf("UPDATE %s SET %s = %s",tableN,timeRLV,tleqn))
#     }
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AGGREGATION: grouping SUBID x gvar in database with flexible query sentences:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Aggregation with SQ language
#tIDX<-sapply(range(timeS),idxF,adj);whereQ="";print(tIDX)
#if (sessN == "nb_S2"){tIDX<- c(1,301)} else if (sessN == "nb_S1"){tIDX <- c(1,226)}; print(tIDX) #cut to time windows of interest
tIDX <- c(1,226); print(tIDX) #response locked data
selectQ<-paste0('SELECT ',paste(unique(c(grpV)),collapse = ", "))
compQ<-paste(sprintf(", AVG(%s) AS %s",depV,depV),collapse = "")
if (!any(is.na(timeL))){whereQ <- sprintf(" WHERE %s >= %d AND %s <= %d",grpV[2],tIDX[1],grpV[2],tIDX[2])}
groupQ<-paste0(' GROUP BY ',paste(grpV,collapse = " , "))
querySTW<-paste0(selectQ,compQ,sprintf(" FROM %s",tableN),whereQ,groupQ) %>% print()
tic();dsA=as.data.table(dbGetQuery(dbCon,querySTW));toc()
#Aggregation with dbplyr package
#empty_name <- tbl(dbCon, tableN) #creates local reference for db
#querySTWb <- empty_name %>%
#  filter(time >=1& time<=250)%>%
#  group_by(SUBID, time, SESSION, NUMTASK) %>%
#  summarise(LH_RSA = mean(LH_RSA),
#            OE_RSA = mean(OE_RSA),
#            BF_RSA = mean(BF_RSA),
#            RB_RSA = mean(RB_RSA))
#querySTWb %>% show_query() #check query
#tic();dsB<- as.data.table(dbGetQuery(dbCon,querySTWb))();toc()

#Ashwin addition -- grouping RSA by relevant features
#tIDX<-sapply(range(timeS),idxF,adj);whereQ="";print(tIDX)
# if (modelV == "S2ID"){tIDX<- c(125,225)} else if (modelV == "S1ID"){tIDX <- c(50,150)}; print(tIDX)
# selectQT<-paste0('SELECT ',paste(unique(c(grpV, taskV)), collapse = ", "))
# compQT<-paste(sprintf(", AVG(%s) AS %s",depV,depV),collapse = "")
# if (!any(is.na(timeL))){whereQ <- sprintf(" WHERE %s >= %d AND %s <= %d",grpV[2],tIDX[1],grpV[2],tIDX[2])}
# groupQT<-paste0(' GROUP BY ',paste(c(grpV,taskV),collapse = ", "))
# querySTWb<-paste0(selectQT,compQT,sprintf(" FROM %s",tableN),whereQ,groupQT) %>% print()
# tic();dsB=as.data.table(dbGetQuery(dbCon,querySTWb));toc();

# Keep consistant time labels (this should be fixed in decoding step...)!
timeC<-colnames(dsA)[grepl("time",colnames(dsA))]
if (!any(is.na(timeL)) & diff(unique(dsA[,get(timeC)]))[1]!=srate){
  #tIDX<-sapply(range(timeS),idxF,0)# no adjustment needed here
  timeLL<- timeL[tIDX[1]:tIDX[2]] %>% print();
  if (is.na(refRLV)){dsA[,(timeC):=lapply(.SD,function(t)timeL[t]),.SDcols=timeC]}
  if (!is.na(refRLV)){dsA[,(timeC):=timeLL[as.numeric(as.factor(dsA[,get(timeC)]))]]}
}

# Now for dsB
# timeC<-colnames(dsB)[grepl("time",colnames(dsB))]
# if (!any(is.na(timeL)) & diff(unique(dsB[,get(timeC)]))[1]!=srate){
#   tIDX<-sapply(range(timeS),idxF,0)# no adjustment needed here
#   timeLL<- timeL[tIDX[1]:tIDX[2]] %>% print();
#   if (is.na(refRLV)){dsB[,(timeC):=lapply(.SD,function(t)timeL[t]),.SDcols=timeC]}
#   if (!is.na(refRLV)){dsB[,(timeC):=timeLL[as.numeric(as.factor(dsB[,get(timeC)]))]]}
# }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Segmentation: averaging data over time window (aggregation with time intervals)
# Need to be adjusted fro DYN, CROSS, JMB, etc...(those already labeled vs. not)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (diff(timeS[c(1,2)])>srate & !any(is.na(timeL))){
  # Define averaging window then aggregate
  dsA[,(timeC):=(lapply(.SD,cutF)),.SDcols=timeC]
  dsA<-dsA[,lapply(.SD, mean,na.rm=T),by=grpV]# time labels will be mid-point
  dsA[,(timeC):=lapply(.SD,function(t)rollmean(timeS,2)[t]),.SDcols=timeC]
}

#now for dsB
if (diff(timeS[c(1,2)])>srate & !any(is.na(timeL)) & !any(is.na(taskV))){
  # Define averaging window then aggregate
  dsB[,(timeC):=(lapply(.SD,cutF)),.SDcols=timeC]
  grpV2<- paste(c(grpV,taskV))
  dsB<-dsB[,lapply(.SD, mean,na.rm=T),by=grpV2]# time labels will be mid-point
  dsB[,(timeC):=lapply(.SD,function(t)rollmean(timeS,2)[t]),.SDcols=timeC]
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summary
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save session information as attributes
setwd(Dir_GRAND)
attributes(dsA)$timeL<-unique(dsA[,get(timeC)])
attributes(dsA)$timeS<-timeS
attributes(dsA)$grpV<-grpV
attributes(dsA)$depV<-depV
if(grepl("time",grpV[2])){code<-"timeSERIES"}else{code<-timeRLV}
if (saveOn){saveRDS(dsA,paste0(expN,"_",decodingType,"_",modelV,"_",evLock,"_",regs,"_",rsaValueN,"_",code,".rds"))}
#if (saveOn){write.table(dsA,paste0('DimensionSwitchPRP_S2ID_SOA.txt'), sep="\t",row.names=FALSE)} #be careful with naming 

# checking
print(dsA);

# Filtering --------------------------------------------------------------------

# Filtering (within subjects)
# dsA<-dsA[ACC_EDVALID==T,]


# Averaging
dsAGG<-dsA[, lapply(.SD, mean, na.rm=TRUE), by=c(grpV[2:5]), .SDcols=depV]
dsAGG<-melt(dsAGG,id.vars = grpV[2:5])
# dsAGG<- dsAGG[ACCRESP==0,]

#theme_set(theme_bwL(14));
theme_set(theme_bw(base_size = 20))#32/28
CSCALE_YlOrRd = rev(brewer.pal(9,"YlOrRd"));
#chance<-1/length(unique(dsA[[depV]]))


# quartz(width=10,height=4)
p<-ggplot(dsAGG[variable == 'CONJ_RSA'],aes(x=time,y=value,color=SOA)) +
  #annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1.75)+
  #annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1.75)+
  facet_rep_grid(SESSION ~ NUMTASK, repeat.tick.labels = F,scales = "fixed")+
  scale_y_continuous(breaks = scales::pretty_breaks(n=2))+
  ggtitle(paste(depV[3], sessN, value))+
  #geom_hline(yintercept = 0, linetype =2) + geom_vline(xintercept = c(0)) +
  geom_point(size=2) + geom_line(size=0.8) + #coord_cartesian(ylim=c(-0.05,1.85)) + # c(-0.05,2.3)
  theme(legend.position = "right", legend.text =element_text(size=10)) +
  theme(strip.text = element_text(angle = 0))
p

outfile <- file.path(paste0(Dir_REPO, "/TasksetXCorsideRSA/"), paste0(depV[3],"_", sessN, "_plot", ".png"))
png(outfile, width = 1920, height = 1080, units = "px")
print(p)
dev.off()

#=========================== for dsB
dsB$TS2 <- case_when(
  dsB$TASKSET_2 == "Low-High" ~ "LH",
  dsB$TASKSET_2 == "Odd-Even" ~ "OE",
  dsB$TASKSET_2 == "Red-Blue" ~ "RB",
  dsB$TASKSET_2 == "Bold-Faint" ~ "BF")
dsAGGB<-dsB[, lapply(.SD, mean, na.rm=TRUE), by=c(paste(c(grpV2[-1],"TS2"))), .SDcols=depV]
dsAGGB<-melt(dsAGGB,id.vars = paste(c(grpV2[-1],"TS2")))
dsAGGB$rel <-  ifelse(mapply(function(pat, txt) grepl(pat, txt), dsAGGB$TS2, dsAGGB$variable), "Relevant Feature", "Irrelevant Feature")
dsAGGB<-dsAGGB[, lapply(.SD, mean, na.rm=TRUE), by=c(paste(c(grpV[-1], "rel"))), .SDcols="value"]
# rel feature v irrel
ggplot(data = dsAGGB, aes(x=time,y=value, color=rel))+#,group=rel,color=rel)) +
  #scale_y_continuous(breaks = scales::pretty_breaks(n=2))+
  geom_line(aes(group=rel),size=1) + 
  #geom_point(aes(group=rel))+
  #coord_cartesian(ylim=c(-0.05,0.1)) +#  c(-0.05,2.3) + #stat_summary(size = 0.5, shape = 1)+ 
  #geom_ribbon(alpha=0.2,linetype=0, aes(fill = rel, ymin = m - sem, ymax = m + sem)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  facet_wrap(SESSION~NUMTASK)+
  scale_color_manual(values=c("cornflowerblue","orange")) + scale_fill_manual(values=c("cornflowerblue","orange")) 


# rel feature diff
dsAGGB<-dcast(dsAGGB, time + SESSION + NUMTASK  ~ rel, value.var = "value") #+ SOA
dsAGGB$diff <- dsAGGB$`Relevant Feature`-dsAGGB$`Irrelevant Feature`
theme_set(theme_bwL(14));


ggplot(dsAGGB,aes(x=time,y=diff, color=NUMTASK))+#,group=rel,color=rel)) +
  #annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1.75)+
  #annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1.75)+
  facet_rep_grid(~SESSION, repeat.tick.labels = F,scales = "free")+
  scale_y_continuous(breaks = scales::pretty_breaks(n=2))+
  
  #geom_hline(yintercept = 0, linetype =2) + geom_vline(xintercept = c(0)) +
  geom_point(shape = 1,size=2, alpha = 0.5) + geom_line(size=1)+# stat_summary( size = 0.5)+ # geom_line(aes(group=variable),size=0.8) + #coord_cartesian(ylim=c(-0.05,1.85)) + # c(-0.05,2.3)
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
        legend.position = c("right"),legend.key.size=unit(1,"lines"),legend.title=element_blank(),legend.text = element_text(size=12)) + 
  scale_color_manual(values=c("red","black")) + scale_fill_manual(values= c("red","black"))

#small test of regression
modelB<-melt(dsB,id.vars = paste(c(grpV2,"TS2")))
modelB$rel <-  ifelse(mapply(function(pat, txt) grepl(pat, txt), modelB$TS2, modelB$variable), "RelevantFeature", "IrrelevantFeature")
modelB <-modelB[, lapply(.SD, mean, na.rm=TRUE), by=c(paste(c(grpV[-2],"rel"))), .SDcols="value"]
modelB<-dcast(modelB, SUBID + SESSION + NUMTASK + SOA~ rel, value.var = "value")
modelB$SESSION <- as.factor(modelB$SESSION)
modelB$SUBID <- as.factor(modelB$SUBID)
modelB$SOA <- as.factor(modelB$SOA)
modelB$NUMTASK <- factor(modelB$NUMTASK, levels = c("single","double"))
contrasts(modelB$SESSION, how.many = 2) <- contr.poly(3)
contrasts(modelB$SOA, how.many = 2) <- contr.poly(5)

m <- lmer(formula = (RelevantFeature)~SESSION*NUMTASK*SOA + (1+ SESSION|SUBID), data = modelB, control=lmerControl(optimizer="bobyqa"))
summary(m)
plot_model(m,type="int")


# Predicting behavior performance from RSA value--------------------------------

#new query to access trial info
tIDX<-sapply(range(timeS),idxF,adj);whereQ="";print(tIDX)
selectQTC<-paste0('SELECT ',paste(unique(c(grpV, "BLOCK", "TRIAL", taskV)), collapse = ", "))
compQTC<-paste(sprintf(", AVG(%s) AS %s",depV,depV),collapse = "")
if (!any(is.na(timeL))){whereQC <- sprintf(" WHERE %s >= %d AND %s <= %d",grpV[2],tIDX[1],grpV[2],tIDX[2])} #change time indices
groupQTC<-paste0(' GROUP BY ',paste(c(grpV, "BLOCK", "TRIAL", taskV),collapse = ", "))
querySTWbC<-paste0(selectQTC,compQTC,sprintf(" FROM %s",tableN),whereQC,groupQTC) %>% print()
tic();dsC=as.data.table(dbGetQuery(dbCon,querySTWbC));toc();

#time stuff
timeC<-colnames(dsC)[grepl("time",colnames(dsC))]
if (!any(is.na(timeL)) & diff(unique(dsC[,get(timeC)]))[1]!=srate){
  tIDX<-sapply(range(timeS),idxF,0)# no adjustment needed here
  timeLL<- timeL[tIDX[1]:tIDX[2]] %>% print();
  if (is.na(refRLV)){dsC[,(timeC):=lapply(.SD,function(t)timeL[t]),.SDcols=timeC]}
  if (!is.na(refRLV)){dsC[,(timeC):=timeLL[as.numeric(as.factor(dsC[,get(timeC)]))]]}
}

if (diff(timeS[c(1,2)])>srate & !any(is.na(timeL)) & !any(is.na(taskV))){
  # Define averaging window then aggregate
  dsC[,(timeC):=(lapply(.SD,cutF)),.SDcols=timeC]
  grpV3<- paste(c(grpV, "BLOCK", "TRIAL", taskV))
  dsC<-dsC[,lapply(.SD, mean,na.rm=T),by=grpV3]# time labels will be mid-point
  dsC[,(timeC):=lapply(.SD,function(t)rollmean(timeS,2)[t]),.SDcols=timeC]
}

# Save session information as attributes
setwd(Dir_GRAND)
attributes(dsC)$timeL<-unique(dsC[,get(timeC)])
attributes(dsC)$timeS<-timeS
attributes(dsC)$grpV<-grpV
attributes(dsC)$depV<-depV
if(grepl("time",grpV[2])){code<-"timeSERIES"}else{code<-timeRLV}
if (saveOn){saveRDS(dsC,paste0(expN,"_",modelV,"_",fileN,"_",code,".rds"))}
if (saveOn){write.table(dsC,paste0('DimensionSwitchPRP_',modelV,'_SOA_byTrial.txt'), sep="\t",row.names=FALSE)}



dsC$relFeature <- case_when(
  dsC$TASKSET_2 == "Low-High" ~ dsC$LH_RSA,
  dsC$TASKSET_2 == "Odd-Even" ~ dsC$OE_RSA,
  dsC$TASKSET_2 == "Red-Blue" ~ dsC$RB_RSA,
  dsC$TASKSET_2 == "Bold-Faint" ~ dsC$BF_RSA)
#Irrelevant features
# Low-High
dsC[TASKSET_2 == "Low-High", 
    irrelFeature := rowMeans(.SD, na.rm = TRUE), 
    .SDcols = c("OE_RSA", "RB_RSA", "BF_RSA")]
# Odd-Even
dsC[TASKSET_2 == "Odd-Even", 
    irrelFeature := rowMeans(.SD, na.rm = TRUE), 
    .SDcols = c("LH_RSA", "RB_RSA", "BF_RSA")]
# Red-Blue
dsC[TASKSET_2 == "Red-Blue", 
    irrelFeature := rowMeans(.SD, na.rm = TRUE), 
    .SDcols = c("OE_RSA", "LH_RSA", "BF_RSA")]
# Bold-Faint
dsC[TASKSET_2 == "Bold-Faint", 
    irrelFeature := rowMeans(.SD, na.rm = TRUE), 
    .SDcols = c("OE_RSA", "RB_RSA", "LH_RSA")]
#merge behavioral data (RT_2)
dsC[ds_beh, on = .(SUBID, SESSION, BLOCK, TRIAL), RT_2 := i.RT_2] #add RT info

dsC$SESSION <- as.factor(dsC$SESSION)
dsC$SUBID <- as.factor(dsC$SUBID)
dsC$NUMTASK <- factor(dsC$NUMTASK, levels = c("single","double"))
contrasts(dsC$SESSION, how.many = 2) <- contr.poly(3)

#dsC <- dsC[time >= 300]


dsAGGC<-dsC[, lapply(.SD, mean, na.rm=TRUE), by=c(grpV3[c(1,3:7)]), .SDcols=c("relFeature", "irrelFeature", "RT_2")]
#dsAGGC<-dsC[, lapply(.SD, mean, na.rm=TRUE), by=c(grpV3[2:4]), .SDcols=c("relFeature", "irrelFeature")]
#dsAGGC<-melt(dsAGGC,id.vars = paste(c(grpV3[2:4])))

ggplot(dsAGGC,aes(x=relFeature,y=log(RT_2), color=NUMTASK))+#,group=rel,color=rel)) +
  #annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1.75)+
  #annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1.75)+
  facet_rep_grid(~SESSION, repeat.tick.labels = F,scales = "free")+
  scale_y_continuous(breaks = scales::pretty_breaks(n=2))+
#  
#  geom_hline(yintercept = 0, linetype =2) + geom_vline(xintercept = c(0)) +
  geom_point(shape = 1,size=2, alpha = 0.5) + geom_smooth(method = "lm")+ # geom_line(size=1)+# stat_summary( size = 0.5)+ # geom_line(aes(group=variable),size=0.8) + #coord_cartesian(ylim=c(-0.05,1.85)) + # c(-0.05,2.3)
  theme(legend.position = "right", legend.text =element_text(size=10)) +
  theme(strip.text = element_text(angle = 0))


m2<- lmer(relFeature~SESSION*NUMTASK*SOA + (1+SESSION|SUBID), data = dsAGGC, control=lmerControl(optimizer="bobyqa"))
summary(m2)
plot_model(m2, type = "pred", terms = c("SESSION", "NUMTASK", "SOA"))
