# Plot_RSAModel
# load pred file (contains single-trial accuracy & confidence)
# fit RSA models to decoding resulf of the conjunction variable
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))

#Setting up directory
fsep<-.Platform$file.sep;
Dir_R<-path.expand("~/Dropbox/w_ONGOINGRFILES/w_OTHERS")
Dir_EDATA<-path.expand("~/Dropbox/w_SCRIPTS/P.HAL2017/HAL2017_PUSHPRP/DataAnalysis/Matlab/EEG")
Dir_BDATA<-path.expand("~/Dropbox/w_SCRIPTS/P.HAL2017/HAL2017_PUSHPRP/DataAnalysis/Matlab/BEH/w_ALLGRAND")
Dir_RSAMODEL<-paste0(Dir_EDATA,"/w_RSAMODELS")


# Load libraries
library(data.table)
library(dplyr)
library(broom)
library(rhdf5)
library(caret)
library(foreach)
library(doMC)
library(tidyr)
library(binhf)
library(GGally)
library(RColorBrewer)
library(RSQLite)
library(lazyeval)
library(stringr)
library(psych)
library(lme4)
library(RcppEigen)
library(scales)


# Source Files
setwd(Dir_R)
source('basic_lib.R')


# Define all model (crossed and nested ones)
setwd(Dir_RSAMODEL)
modelV <- "STIMSPACE(small)"

# Load all models
setwd(paste(Dir_RSAMODEL,modelV,sep=fsep))
mlist <- list.files(pattern = ".txt")
models<-lapply(as.list(mlist),function(f){read.table(f,header = F)})
for (m in 1:length(mlist)){assign(gsub(".txt","",mlist)[m],models[[m]])}


# Transform model 
m <- Var1_M %>% data.table()
m <- m[nrow(m):1,]
# m <- matrix(0, nrow = 4, ncol = 4) %>% data.table()
m[,Prediction := paste0("V",nrow(m):1)]
m[,Prediction := factor(Prediction,levels=paste0("V",nrow(m):1))]
mcols <- names(m)[!names(m) %in% "Prediction"]
m = melt(m, is.vars = "Prediction", measure.vars = mcols, variable.name = "Reference")

# # White matrix
m[, value:=0]


# Plot
CSCALE_Grey = (brewer.pal(9,"Greys"));
CMAP_GRAY <- colorRampPalette(c(CSCALE_Grey[1],CSCALE_Grey[9]))

# Confusion matrix plot
theme_set(theme_bw(base_size = 20))#32/28
quartz(width=6,height=5.75);
# quartz(width=12.5,height=5.75);# diverse
ggplot(m, aes(Reference,Prediction, fill=value)) +
  geom_tile(color="black",size=.20) +
  scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+
  #scale_fill_gradientn(limits=c(chance-0.08,chance+0.08),colours=CMAP_GRAY(40),oob=squish)+
  scale_fill_gradientn(limits=c(0,1),colours=CMAP_GRAY(40),oob=squish)+# for Model
  #coord_cartesian(xlim = c(1.1, 11.9), ylim = c(1.1,11.9))+
  #Aesthetics!-------------------------
  theme(plot.margin=unit(rep(1,1,4),"line"),
      panel.border = element_rect(colour = "black", fill=NA, size=1.5),
      legend.key = element_blank(),
      legend.position = "off",
      legend.title=element_blank(),
      legend.text =element_blank(),
      axis.ticks.length = unit(0.5, "lines"),
      axis.title  = element_blank(),
      axis.text   = element_blank(),
      plot.title = element_text(hjust = 0, size = 20),
      strip.background=element_blank())


