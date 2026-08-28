###### R Script for publication by Croll et. al. 2026 about spawning timing of Atlantic herring ######
# This script contains the R-code to import and process the data form the model simulations and produce the figures from the manuscript

###### SETTINGS AND PACKAGES ######

# device specific settings
setwd("C:/Users/jacr0004/Documents/Local/EBTresults/evo_Rscript")

#  Load needed packages
library(ggplot2)
library(ggnewscale)
library(zoo)
library(patchwork)
library(dplyr)
library(viridis)
library(purrr)

# load function to import population structure data from model simulations in txt format converted from the model output csb using csb2txt in the EBT package
source("Rcsbtxtread.R")

#system settings
Sys.setlocale("LC_ALL", "English")


# create colour blind friendly pallet
evo_pallet <- c("#000000FF",viridis_pal(option="turbo", begin=0.3, end=0.9)(2))
eco_pallet <- viridis_pal(option="plasma", begin=0, end=0.8)(2)

# create semi square root scaling
abssqrt <- function(x){
  sign(x)*sqrt(abs(x))
}

iabssqrt <- function(x){
  sign(x)*x^2
}


###### GENERAL USER FUNCTIONS ######

### FUNCTION TO  CALCULATE SELECTION GRADIENT
calc_fitgrad <- function(mutfit){
  # loop dataset to prevent edge cases
  mutfit0 <- mutfit
  mutfit2 <- mutfit
  mutfit0$SSTART <- mutfit0$SSTART-365
  mutfit2$SSTART <- mutfit2$SSTART+365
  mutfit <- rbind(mutfit0, mutfit, mutfit2)
  
  # fill missing data to prevent edge cases
  mutfit$FitmutHappr <- na.approx(mutfit$FitmutH, na.rm = FALSE)
  mutfit$FitmutLappr <- na.approx(mutfit$FitmutL, na.rm = FALSE)
  mutfit$Fitresappr <- na.approx(mutfit$Fitres, na.rm = FALSE)
  
  # calculate running median to smooth data
  mutfit$filterH <- runmed(mutfit$FitmutHappr-mutfit$Fitresappr, k=48)
  mutfit$filterL <- runmed(mutfit$FitmutLappr-mutfit$Fitresappr, k=48)
  
  # calculate fitness graident
  mutfit$filterfitgrad <- mutfit$filterH - mutfit$filterL
  mutfit$secder <- mutfit$filterH + mutfit$filterL
  
  # restore missing data and looping
  mutfit$filterfitgrad[which(is.na(mutfit$Fitres))] <- NA
  mutfit$filterfitgrad[which(mutfit$SSTART <= 45)] <- NA
  mutfit$filterfitgrad[which(mutfit$SSTART >=410) ] <- NA
  
  return(mutfit[which(mutfit$SSTART > 45 & mutfit$SSTART < 410),])
}

### FUNCTION TO LOCATE ESS
detectESS <- function(mutfit){
  # detect ESS by change in sign of fitness gradient
  mutfit$detectESS <- NA
  mutfit$detectESS[1:(nrow(mutfit)-1)] <- mutfit$filterfitgrad[1:(nrow(mutfit)-1)] * mutfit$filterfitgrad[2:(nrow(mutfit))] 
  
  index <- which(mutfit$detectESS<=0) 
  
  # output results
  ESSdata <- data.frame(SSTART = (mutfit$SSTART[index]+mutfit$SSTART[index+1])/2,
                        MUF = mutfit$MUF[index],
                        ESSstab = - sign(mutfit$filterfitgrad[index]),
                        ESSsecder = sign((mutfit$secder[index]+mutfit$secder[index+1])/2))
  
}

### FUNCTION TO SUMMARIZE POPULATION STRUCTURE
summarisestruc  <- function(struc){
  
  envvar <- struc[[1]]
  popstruc <- as.data.frame(struc[[2]])
  
  popstruc$coho <- floor(popstruc$btime/365)
  popstruc$totage <- popstruc$age + popstruc$bage
  popstruc$matage <- popstruc$matage + popstruc$bage
  
  popstruc_summary <- summarise(.data = popstruc, 
                                totdens = sum(density),
                                maxlen=max(len), 
                                minlen=min(len),
                                meanlen = weighted.mean(len,density),
                                maxcond=max(cond), 
                                mincond=min(cond),
                                meancond=weighted.mean(cond,density),
                                maxage=max(totage), 
                                minage=min(totage),
                                meanage = weighted.mean(totage,density),
                                bdens = sum(bdens),
                                minmatage =min(matage),
                                maxmatage =max(matage),
                                meanmatage =weighted.mean(matage,density),
                                .by=c("coho"))
  
  
  popstruc_summary <- cbind( envvar, popstruc_summary )
  
  return(popstruc_summary)
  
} 

###### FIG 1: BIFURCATION OF ESS OVER FISHING INTENSITY ######

### EXTINCTION BOUNDARIES

# load invasion experiment in environment without herring
evo_MUFr0_list <- lapply(seq(100, 900, 25) ,  function(x){
  
  # read file
  name <- paste0("Model/evo_r0/herringEBT_MUF00",x,".out")
  
  r0 <- read.delim(name,
                  sep = "\t",
                  header=FALSE,
                  col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","LRO","R0","logR0","GT","SSTART"))
  
  # recalculate start of spawning date
  r0$SSTART <- (r0$SSTART-1)%%365+1
  allday <- data.frame(SSTART = 1:365) 
  r0 <- merge(r0, allday, all=TRUE)
  r0 <- r0[order(r0$SSTART),]
  
  # add fishing mortality
  r0$MUF <- x[1]/100000
  
  return(r0)
  
} )

# buid dataframe
evo_MUFr0 <- data.frame(do.call("rbind",evo_MUFr0_list))



#### DETECT EXTINCION BRANSHING POINTS

BP_MUF_list <- lapply(evo_MUFr0_list, function(r0){
 
  #loop dataset to prevent edge cages
  r0 <- rbind(r0, r0, r0)
  
  # find value at which R0 switches sign
  r0$detectBP <- NA
  r0$detectBP[1:(nrow(r0)-1)] <- r0$R0[1:(nrow(r0)-1)] * r0$R0[2:(nrow(r0))] 
  
  BPindex <- which(r0$detectBP<=0) 
  
  # calculate values of branshing points
  BPlist <- lapply(BPindex, function(index, r0){
    BP <- colMeans(r0[c(index,index+1),c("SSTART", "MUF")])
    return(BP)
  }, r0)
  
  # clean data from looped dataset and return
  BPs <- data.frame(do.call("rbind",BPlist))
  BPs <- BPs[!duplicated(BPs),]
  
  return(BPs)
} )

# build dataframe
BP_MUF <- data.frame(do.call("rbind",BP_MUF_list))
BP_MUF$SSTART[which(BP_MUF$SSTART<45)] <- BP_MUF$SSTART[which(BP_MUF$SSTART<45)]+365



### CALCULATE ESS

# load simulation of invation fitness over full range
evo_MUFmutfit_list <- lapply(seq(1, 8 ,1), function(x){
  
  # read file
  name_up <- paste0("Model/evo_MUFbif/herringEBT_MUF00",x,"up.out")
  name_down <- paste0("Model/evo_MUFbif/herringEBT_MUF00",x,"down.out")
  
  mutfit_up <- read.delim(name_up,
                          sep = "\t",
                          header=FALSE,
                          col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitmutL","FitmutH","LROres","LROmutL","LROmutH", "GTres", "GTL","GTH", "SSTART"))
  
  mutfit_down <- read.delim(name_down,
                            sep = "\t",
                            header=FALSE,
                            col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitmutL","FitmutH","LROres","LROmutL","LROmutH", "GTres", "GTL","GTH", "SSTART"))
  
  mutfit_up$dir ="up"
  mutfit_down$dir ="down"
  
  mutfit <- rbind(mutfit_up, mutfit_down)
  
  # recalculate start of spawning date
  mutfit$SSTART <- (mutfit$SSTART-1)%%365+1
  mutfit <- mutfit[!duplicated(mutfit$SSTART),]
  
  allday <- data.frame(SSTART = 1:365) 
  mutfit <- merge(mutfit, allday, all.y=TRUE, all.x = FALSE)
  mutfit <- mutfit[order(mutfit$SSTART),]
  
  # add fishing mortality
  mutfit$MUF <- x[1]/1000
  
  return(mutfit)
  
} )


# load additional data runs with invasion fitness over small rage
evo_MUFmutfit_add_list <- lapply(seq(100, 850 ,25), function(x){
  
  # read file
  name_atta <- paste0("Model/evo_MUFbif/herringEBT_MUF00",x,"att_a.out")
  name_attb <- paste0("Model/evo_MUFbif/herringEBT_MUF00",x,"att_b.out")
  name_rep <- paste0("Model/evo_MUFbif/herringEBT_MUF00",x,"rep.out")

  if(file.exists(name_atta)){
    mutfit_atta <- read.delim(name_atta,
                              sep = "\t",
                              header=FALSE,
                              col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitmutL","FitmutH","LROres","LROmutL","LROmutH", "GTres", "GTL","GTH", "SSTART"))
    mutfit_atta$dir ="up"
    mutfit <- mutfit_atta
  }else{
    mutfit <- data.frame()
  }
  
  if(file.exists(name_attb)){
    mutfit_attb <- read.delim(name_attb,
                              sep = "\t",
                              header=FALSE,
                              col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitmutL","FitmutH","LROres","LROmutL","LROmutH", "GTres", "GTL","GTH", "SSTART"))
    mutfit_attb$dir ="up"
    mutfit <- rbind(mutfit, mutfit_attb)
  }
  
  if(file.exists(name_rep)){
    mutfit_rep <- read.delim(name_rep,
                             sep = "\t",
                             header=FALSE,
                             col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitmutL","FitmutH","LROres","LROmutL","LROmutH", "GTres", "GTL","GTH", "SSTART"))
    mutfit_rep$dir ="up"
    mutfit <- rbind(mutfit, mutfit_rep)
  }
  
  if(nrow(mutfit) == 0){return(NA)}
  
  # recalculate start of spawning date
  mutfit$SSTART <- (mutfit$SSTART-1)%%365+1
  mutfit <- mutfit[!duplicated(mutfit$SSTART),]
  
  allday <- data.frame(SSTART = 1:365) 
  mutfit <- merge(mutfit, allday, all.y=TRUE, all.x = FALSE)
  mutfit <- mutfit[order(mutfit$SSTART),]
  
  # add fishing mortality
  mutfit$MUF <- x[1]/100000
  
  
  return(mutfit)
  
} )


evo_MUFmutfit_add_list <- evo_MUFmutfit_add_list[which(!is.na(evo_MUFmutfit_add_list))]

evo_MUFmutfit_all_list <- append(evo_MUFmutfit_list, evo_MUFmutfit_add_list)

# calculate fitness gradient
evo_MUFmutfit_all_list <- lapply(evo_MUFmutfit_all_list, calc_fitgrad)

# detect ESS
ESSdata_list <- lapply(evo_MUFmutfit_all_list, detectESS)
ESSdata <- data.frame(do.call("rbind",ESSdata_list))

# Extract and duplicate evolutionary limit point
ESSLP <- ESSdata[which(ESSdata$MUF==0.00525 & ESSdata$SSTART>265),]
ESSdata <- rbind(ESSdata, data.frame(SSTART=ESSLP$SSTART, MUF=ESSLP$MUF, ESSstab=-1, ESSsecder=0))
BP_MUF <- rbind(BP_MUF, ESSLP[,c("SSTART","MUF")])



### PLOT FIGURE 1

ESSdata <- ESSdata[order(1/ESSdata$MUF*ESSdata$ESSstab),]
BP_MUF <- BP_MUF[order(BP_MUF$MUF),]

ESS_bifplot <- ggplot(data=BP_MUF, aes(x=as.Date((SSTART + 365 - 45)%%365+45 + 38), y = 365*MUF))+
  geom_path(data=ESSdata[which(ESSdata$SSTART<145),], aes(x=as.Date((SSTART + 365 - 45 )%%365+45 + 38), colour=as.factor(-1), linetype = as.factor(-1)), linewidth=1)+
  geom_path(data=ESSdata[which(ESSdata$SSTART>200),], aes(x=as.Date((SSTART + 365 - 45)%%365+45 + 38), colour=as.factor(ESSstab), linetype = as.factor(ESSstab)), linewidth=1)+
  geom_line(aes(colour=as.factor(0), linetype=as.factor(0)), linewidth=1)+
  labs(x="Peak of spawning season", y="Fishing pressure (1/y)")+
  scale_x_date(date_labels = "%b", , limits=as.Date(c(45+ 38,410+ 38)))+
  scale_y_continuous( limits=c(365*0.001, 365*0.0085),breaks=seq(0,3,0.5) )+
  scale_linetype_manual(name="", limits=as.factor(c(0, -1, 1)), labels= c("Extinction boundary", "Evolutionary attractor", "Evolutionary repellor"), values=c("dotdash","solid","dashed"))+
  scale_colour_manual(name="",limits=as.factor(c(0, -1, 1)), labels= c("Extinction boundary", "Evolutionary attractor", "Evolutionary repellor"), values=evo_pallet)+
  scale_shape_manual(name="", breaks=as.factor(c(-1, 1, 0)), labels= c("Evolutionary attractor","Evolutionary repellor", "Extinction boundary"), values=c(19, 1, 0), guide="none")+
  theme_classic( )+
  theme(legend.position = c(0.75, 0.9))+
  new_scale_color()+
  geom_line( data=data.frame(SSTART=seq(240+ 38, 270+ 38, 0.1), MUF=0.365), aes(x=as.Date(SSTART), y =MUF, colour=SSTART-(240+38)), arrow=arrow(ends="last", length = unit(0.4, "cm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(360+ 38, 330+ 38, -0.1), MUF=0.365), aes(x=as.Date(SSTART), y =MUF, colour=(360+38)-SSTART), arrow=arrow(ends="first", length = unit(0.4, "cm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(70+ 38, 100+ 38, 0.1), MUF=0.365), aes(x=as.Date(SSTART), y =MUF, colour=SSTART-(70+38)), arrow=arrow(ends="last", length = unit(0.4, "cm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(175+ 38, 145+ 38, -0.1), MUF=0.365), aes(x=as.Date(SSTART), y =MUF, colour=(175+38)-SSTART), arrow=arrow(ends="first", length = unit(0.4, "cm")), linewidth=1 )+
  scale_colour_gradient(low=evo_pallet[3],high=evo_pallet[2], guide="none")
    
ggsave("Figures/Figure1R.pdf", ESS_bifplot , width=12.1, height=8, units="cm")



###### FIG 2: FITNESS GRADIENTS WITH SINGLE TEMPERATRE DEPENDENCE ######

# load simulation data into list and calculate invasionfitness
evo_single_list <- lapply(c("all","onlylarv", "onlyphys","onlyres"), function(x){
  
  # read file
  name_up <- paste0("Model/evo_single/herringEBT_",x[1],"_up",".out")
  name_down <- paste0("Model/evo_single/herringEBT_",x[1],"_down",".out")
  
  mutfit_up <- read.delim(name_up,
                            sep = "\t",
                            header=FALSE,
                            col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitmutL","FitmutH","LROres","LROmutL","LROmutH", "GTres", "GTL","GTH", "SSTART"))
  mutfit_up$dir <- "up"
  
  mutfit_down <- read.delim(name_down,
                          sep = "\t",
                          header=FALSE,
                          col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitmutL","FitmutH","LROres","LROmutL","LROmutH", "GTres", "GTL","GTH", "SSTART"))
  mutfit_down$dir <- "down"

  mutfit <- rbind(mutfit_up, mutfit_down)
  
  
  if(nrow(mutfit) == 0){return(NA)}
  
  # recalculate start of spawning date
  mutfit$SSTART <- (mutfit$SSTART-1)%%365+1
  mutfit <- mutfit[!duplicated(mutfit$SSTART),]
  
  allday <- data.frame(SSTART = 1:365) 
  mutfit <- merge(mutfit, allday, all.y=TRUE, all.x = FALSE)
  mutfit <- mutfit[order(mutfit$SSTART),]
  
  # add run info
  mutfit$MUF <- 0.001
  
  mutfit$run <- x
  
  return(mutfit)
  
} )

# calculate fitness gradient
evo_single_list <- lapply(evo_single_list, calc_fitgrad)
names(evo_single_list) <- c("all","onlylarv", "onlyphys","onlyres")

evo_single <- data.frame(do.call("rbind",evo_single_list))

# detect ESS
ESSdata_single_list <- lapply(evo_single_list, detectESS )

ESSdata_single <- data.frame(do.call("rbind",ESSdata_single_list))


### MAKE FIGURE 2

# temperature dependence in all
fitness_all <- ggplot(data=evo_single[which(evo_single$run=="all" & evo_single$SSTART>54.5 & evo_single$SSTART<401.5),], aes(x=as.Date((SSTART + 365 - 45)%%365+45+ 38)))+
  geom_hline(yintercept=0)+
  geom_path(aes(y=filterfitgrad), size=1)+
  geom_point(data=ESSdata_single_list$all, y=0, size=3, aes(colour=as.factor(ESSstab), shape = as.factor(ESSstab)))+
  theme_classic()+
  scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)) )+
  scale_y_continuous(trans=scales::trans_new("abssqrt",abssqrt, iabssqrt), limits=c(-0.05, 0.075))+
  labs(x="Peak of spawning season", y="Selection gradient")+
  scale_colour_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=evo_pallet[c(2,3)])+
  scale_shape_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=c(19, 15))+
  new_scale_color()+
  geom_line( data=data.frame(SSTART=seq(240+ 38, 270+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(240+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(360+ 38, 330+ 38, -0.1)), aes(x=as.Date(SSTART), colour=360+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(70+ 38, 100+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(70+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(175+ 38, 145+ 38, -0.1)), aes(x=as.Date(SSTART), colour=175+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
  scale_colour_gradient(low=evo_pallet[3],high=evo_pallet[2], guide="none")+
  theme(legend.position = c(0.8, 1))+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, , face="bold"), axis.text.x=element_blank() )

# temperature dependence in egg and larvea
fitness_onlylarv <- ggplot(, aes(x=as.Date((SSTART + 365 - 45)%%365+45+ 38)))+
  geom_hline(yintercept=0)+
  geom_path(data=evo_single[which( (evo_single$run=="all") & evo_single$SSTART>54.5 & evo_single$SSTART<401.5),] , aes(y=filterfitgrad), size=1, colour="grey80")+
  geom_path(data=evo_single[which( (evo_single$run=="onlylarv")),] , aes(y=filterfitgrad), size=1, colour="black")+
  geom_point(data=ESSdata_single_list$onlylarv, y=0, size=3, aes(colour=as.factor(ESSstab), shape = as.factor(ESSstab)))+
  theme_classic()+
  scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)) )+
  scale_y_continuous(trans=scales::trans_new("abssqrt",abssqrt, iabssqrt), limits=c(-0.05, 0.075))+
  labs(x="Peak of spawning season", y="Selection gradient")+
  scale_colour_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=evo_pallet[c(2,3)])+
  scale_shape_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=c(19, 15))+
  new_scale_color()+
  geom_line( data=data.frame(SSTART=seq(240+ 38, 270+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(240+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(360+ 38, 330+ 38, -0.1)), aes(x=as.Date(SSTART), colour=360+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(80+ 38, 110+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(80+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(195+ 38, 165+ 38, -0.1)), aes(x=as.Date(SSTART), colour=195+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
  scale_colour_gradient(low=evo_pallet[3],high=evo_pallet[2], guide="none")+
  theme(legend.position = "none")+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), , axis.text.x=element_blank())+
  annotate(geom="text",  x=as.Date(Inf), hjust=1, y=0.06, label="Only temperature dependent eggs and larvae")

# temperature dependence in physiology
fitness_onlyphys<- ggplot(, aes(x=as.Date((SSTART + 365 - 45)%%365+45+ 38)))+
  geom_hline(yintercept=0)+
  geom_path(data=evo_single[which( (evo_single$run=="all") & evo_single$SSTART>54.5 & evo_single$SSTART<401.5),] , aes(y=filterfitgrad), size=1, colour="grey80")+
  geom_path(data=evo_single[which( (evo_single$run=="onlyphys") ),] , aes(y=filterfitgrad), size=1, colour="black")+
  geom_point(data=ESSdata_single_list$onlyphys, y=0, size=3, aes(colour=as.factor(ESSstab),shape = as.factor(ESSstab)))+
  theme_classic()+
  scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)) )+
  scale_y_continuous(trans=scales::trans_new("abssqrt",abssqrt, iabssqrt), limits=c(-0.05, 0.075))+
  labs(x="Peak of spawning season", y="Selection gradient")+
  scale_colour_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=evo_pallet[c(2,3)])+
  scale_shape_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=c(19, 15))+
  new_scale_color()+
  geom_line( data=data.frame(SSTART=seq(310+ 38, 340+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(310+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(70+ 38, 100+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(70+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(200+ 38, 170+ 38, -0.1)), aes(x=as.Date(SSTART), colour=200+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
  scale_colour_gradient(low=evo_pallet[3],high=evo_pallet[2], guide="none")+
  theme(legend.position = "none")+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank()  )+
  annotate(geom="text",  x=as.Date(Inf), hjust=1, y=0.06, label="Only temperature dependent herring physiology")
    
# temperature dependence in prey
fitness_onlyres <- ggplot(, aes(x=as.Date((SSTART + 365 - 45)%%365+45+ 38)))+
  geom_hline(yintercept=0)+
  geom_path(data=evo_single[which( (evo_single$run=="all") & evo_single$SSTART>54.5 & evo_single$SSTART<401.5),] , aes(y=filterfitgrad), size=1, colour="grey80")+
  geom_path(data=evo_single[which( (evo_single$run=="onlyres") ),] , aes(y=filterfitgrad), size=1, colour="black")+
  geom_point(data=ESSdata_single_list$onlyres, y=0, size=3, aes(colour=as.factor(ESSstab), shape = as.factor(ESSstab)))+
  theme_classic()+
  scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)) )+
  scale_y_continuous(trans=scales::trans_new("abssqrt",abssqrt, iabssqrt), limits=c(-0.05, 0.075))+
  labs(x="Peak of spawning season", y="Selection gradient")+
  scale_colour_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=evo_pallet[c(2,3)])+
  scale_shape_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=c(19, 15))+
  new_scale_color()+
  geom_line( data=data.frame(SSTART=seq(230+ 38, 260+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(230+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(395+ 38, 365+ 38, -0.1)), aes(x=as.Date(SSTART), colour=395+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
  geom_line( data=data.frame(SSTART=seq(110+ 38, 80+ 38, -0.1)), aes(x=as.Date(SSTART), colour=110+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
  scale_colour_gradient(low=evo_pallet[3],high=evo_pallet[2], guide="none")+
  theme(legend.position = "none")+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold") )+
  annotate(geom="text", x=as.Date(Inf), hjust=1, y=0.06, label="Only temperature and light dependent prey")

fitness_figs<- wrap_plots(fitness_all, plot_spacer(), fitness_onlylarv, plot_spacer(), fitness_onlyphys, plot_spacer(), fitness_onlyres, ncol = 1 , axis_titles = "collect_x", tag_level = 'new', heights=c(1,-0.1, 1, -0.1, 1, -0.1, 1))+
 plot_annotation(tag_levels="A") 

 
ggsave("Figures/Figure2R.pdf", fitness_figs , width=12.1, height=22, units="cm")
  
    
###### FIG S6: FITNESS GRADIENT WITH DOUBLE TEMPERATRE DEPENDENCE ######
  
# load simulation data into list

# load simulation data into list and calculate invasionfitness
evo_excl_list <- lapply(c("all","exlarv", "exphys","exres"), function(x){
  
  # read file
  name_up <- paste0("Model/evo_double/herringEBT_",x[1],"_up",".out")
  name_down <- paste0("Model/evo_double/herringEBT_",x[1],"_down",".out")
  
  mutfit_up <- read.delim(name_up,
                          sep = "\t",
                          header=FALSE,
                          col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitmutL","FitmutH","LROres","LROmutL","LROmutH", "GTres", "GTL","GTH", "SSTART"))
  mutfit_up$dir <- "up"
  
  mutfit_down <- read.delim(name_down,
                            sep = "\t",
                            header=FALSE,
                            col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitmutL","FitmutH","LROres","LROmutL","LROmutH", "GTres", "GTL","GTH", "SSTART"))
  mutfit_down$dir <- "down"
  
  mutfit <- rbind(mutfit_up, mutfit_down)
  
  
  if(nrow(mutfit) == 0){return(NA)}
  
  # recalculate start of spawning date
  mutfit$SSTART <- (mutfit$SSTART-1)%%365+1
  mutfit <- mutfit[!duplicated(mutfit$SSTART),]
  
  allday <- data.frame(SSTART = 1:365) 
  mutfit <- merge(mutfit, allday, all.y=TRUE, all.x = FALSE)
  mutfit <- mutfit[order(mutfit$SSTART),]
  
 
  # add run info
  mutfit$MUF <- 0.001
  
  mutfit$run <- x
  
  return(mutfit)
  
} )

# calculate selection gradient
evo_excl_list <- lapply(evo_excl_list, calc_fitgrad)
names(evo_excl_list) <- c("all","exlarv", "exphys","exres")
evo_excl <- data.frame(do.call("rbind",evo_excl_list))

# detect ESS
ESSdata_excl_list <- lapply(evo_excl_list, detectESS)
ESSdata_excl <- data.frame(do.call("rbind",ESSdata_excl_list))


### MAKE FIGURE S6

fitness_all_excl <- ggplot(data=evo_excl[which(evo_excl$run=="all" & evo_excl$SSTART>54.5 & evo_excl$SSTART<401.5),], aes(x=as.Date((SSTART + 365 - 45)%%365+45+ 38)))+
      geom_hline(yintercept=0)+
      geom_path(aes(y=filterfitgrad), size=1)+
      geom_point(data=ESSdata[which(ESSdata$MUF==0.001),], y=0, size=3, aes(colour=as.factor(ESSstab), shape = as.factor(ESSstab)))+
      theme_classic()+
      scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)) )+
      scale_y_continuous(trans=scales::trans_new("abssqrt",abssqrt, iabssqrt), limits=c(-0.055, 0.1))+
      labs(x="Peak of spawning season", y="Selection gradient")+
      scale_colour_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=evo_pallet[c(2,3)])+
      scale_shape_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=c(19, 15))+
      new_scale_color()+
      geom_line( data=data.frame(SSTART=seq(240+ 38, 270+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(240+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
      geom_line( data=data.frame(SSTART=seq(360+ 38, 330+ 38, -0.1)), aes(x=as.Date(SSTART), colour=360+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
      geom_line( data=data.frame(SSTART=seq(70+ 38, 100+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(70+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
      geom_line( data=data.frame(SSTART=seq(175+ 38, 145+ 38, -0.1)), aes(x=as.Date(SSTART), colour=175+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
      scale_colour_gradient(low=evo_pallet[3],high=evo_pallet[2], guide="none")+
      theme(legend.position = c(0.8, 1))+
      theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, , face="bold"), axis.text.x=element_blank() )
    
fitness_exlarv <- ggplot(, aes(x=as.Date((SSTART + 365 - 45)%%365+45+ 38)))+
      geom_hline(yintercept=0)+
      geom_path(data=evo_excl[which( (evo_excl$run=="all") & evo_single$SSTART>54.5 & evo_single$SSTART<401.5),] , aes(y=filterfitgrad), size=1, colour="grey80")+
      geom_path(data=evo_excl[which( (evo_excl$run=="exlarv")),] , aes(y=filterfitgrad), size=1, colour="black")+
      geom_point(data=ESSdata_excl_list$exlarv, y=0, size=3, aes(colour=as.factor(ESSstab), shape = as.factor(ESSstab)))+
      theme_classic()+
      scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)) )+
      scale_y_continuous(trans=scales::trans_new("abssqrt",abssqrt, iabssqrt), limits=c(-0.05, 0.1))+
      labs(x="Peak of spawning season", y="Selection gradient")+
      scale_colour_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=evo_pallet[c(2,3)])+
      scale_shape_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=c(19, 15))+
      new_scale_color()+
      geom_line( data=data.frame(SSTART=seq(285+ 38, 315+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(285+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
      geom_line( data=data.frame(SSTART=seq(50+ 38, 80+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(50+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
      geom_line( data=data.frame(SSTART=seq(165+ 38, 135+ 38, -0.1)), aes(x=as.Date(SSTART), colour=165+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
      scale_colour_gradient(low=evo_pallet[3],high=evo_pallet[2], guide="none")+
      theme(legend.position = "none")+
      theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, , face="bold") )+
      annotate(geom="text",  x=as.Date(Inf), hjust=1, y=0.085, label="Temperature and light dependent herring physiology and prey")
    
    fitness_exphys <- ggplot(, aes(x=as.Date((SSTART + 365 - 45)%%365+45+ 38)))+
      geom_hline(yintercept=0)+
      geom_path(data=evo_excl[which( (evo_excl$run=="all") & evo_excl$SSTART>54.5 & evo_excl$SSTART<401.5),] , aes(y=filterfitgrad), size=1, colour="grey80")+
      geom_path(data=evo_excl[which( (evo_excl$run=="exphys") ),] , aes(y=filterfitgrad), size=1, colour="black")+
      geom_point(data=ESSdata_excl_list$exphys, y=0, size=3, aes(colour=as.factor(ESSstab),shape = as.factor(ESSstab)))+
      theme_classic()+
      scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)) )+
      scale_y_continuous(trans=scales::trans_new("abssqrt",abssqrt, iabssqrt), limits=c(-0.055, 0.1))+
      labs(x="Peak of spawning season", y="Selection gradient")+
      scale_colour_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=evo_pallet[c(2,3)])+
      scale_shape_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=c(19, 15))+
      new_scale_color()+
      geom_line( data=data.frame(SSTART=seq(210+ 38, 240+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(210+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
      geom_line( data=data.frame(SSTART=seq(360+ 38, 330+ 38, -0.1)), aes(x=as.Date(SSTART), colour=360+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
      geom_line( data=data.frame(SSTART=seq(70+ 38, 100+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(70+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
      geom_line( data=data.frame(SSTART=seq(160+ 38, 120+ 38, -0.1)), aes(x=as.Date(SSTART), colour=150+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
      scale_colour_gradient(low=evo_pallet[3],high=evo_pallet[2], guide="none")+
      theme(legend.position = "none")+
      theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, , face="bold"), axis.text.x=element_blank()  )+
      annotate(geom="text",  x=as.Date(Inf), hjust=1, y=0.085, label="Temperature and light dependent herring egg, larvea and prey")
    
    
    fitness_exres <- ggplot(, aes(x=as.Date((SSTART + 365 - 45)%%365+45+ 38)))+
      geom_hline(yintercept=0)+
      geom_path(data=evo_excl[which( (evo_excl$run=="all") & evo_excl$SSTART>54.5 & evo_excl$SSTART<401.5),] , aes(y=filterfitgrad), size=1, colour="grey80")+
      geom_path(data=evo_excl[which( (evo_excl$run=="exres") ),] , aes(y=filterfitgrad), size=1, colour="black")+
      geom_point(data=ESSdata_excl_list$exres, y=0, size=3, aes(colour=as.factor(ESSstab), shape = as.factor(ESSstab)))+
      theme_classic()+
      scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)) )+
      scale_y_continuous(trans=scales::trans_new("abssqrt",abssqrt, iabssqrt), limits=c(-0.055, 0.1))+
      labs(x="Peak of spawning season", y="Selection gradient")+
      scale_colour_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=evo_pallet[c(2,3)])+
      scale_shape_manual(name="", limits=as.factor(c( -1, 1)), labels= c("Evolutionary attractor","Evolutionary repellor"), values=c(19, 15))+
      new_scale_color()+
      geom_line( data=data.frame(SSTART=seq(280+ 38, 310+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(280+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
      geom_line( data=data.frame(SSTART=seq(380+ 38, 350+ 38, -0.1)), aes(x=as.Date(SSTART), colour=380+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
      geom_line( data=data.frame(SSTART=seq(70+ 38, 100+ 38, 0.1)), aes(x=as.Date(SSTART), colour=SSTART-(70+ 38)), y=0, arrow=arrow(ends="last", length = unit(4, "mm")), linewidth=1 )+
      geom_line( data=data.frame(SSTART=seq(215+ 38, 185+ 38, -0.1)), aes(x=as.Date(SSTART), colour=215+ 38-SSTART), y=0, arrow=arrow(ends="first", length = unit(4, "mm")), linewidth=1 )+
      scale_colour_gradient(low=evo_pallet[3],high=evo_pallet[2], guide="none")+
      theme(legend.position = "none")+
      theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank() )+
      annotate(geom="text", x=as.Date(Inf), hjust=1, y=0.085, label="Temperature dependent herring physiology, eggs and larvea")
    
    fitness_figs_excl <- wrap_plots(fitness_all_excl, plot_spacer(), fitness_exres, plot_spacer(), fitness_exphys, plot_spacer(), fitness_exlarv, ncol = 1 , axis_titles = "collect_x", tag_level = 'new', heights=c(1,-0.1, 1, -0.1, 1, -0.1, 1))+
      plot_annotation(tag_levels="A") 
    
    ggsave("Figures/FigureS6R.pdf", fitness_figs_excl , width=14, height=22, units="cm")

    
   
  
###### FIG S4: TIMESERIES WITH ONLY TEMPERATURE DEPENDENT RESOURCE #####
    
# read summary output data
onlyres_SSTART351_env <- read.delim(paste0("Model/evo_onepop/herringEBT_onlyres_SSTART351.out"),
                                   sep = "\t",
                                   header=FALSE,
                                   col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass", "adltnum", "adltmass","matage","maxlen","temp","light"))
    
onlyres_SSTART156_env <- read.delim(paste0("Model/evo_onepop/herringEBT_onlyres_SSTART156.out"),
                                        sep = "\t",
                                        header=FALSE,
                                        col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass", "adltnum", "adltmass","matage","maxlen","temp","light"))
    
    
# read full population structures
envnames <- c("time", "zoop","bent","eggbuff")
statenames <- c("xmass","ymass","age","mass","len","cond","btime","bage","matage","bdens")
    
onlyres_SSTART351_struclist <- csbtxtread("Model/evo_onepop/herringEBT_onlyres_SSTART351.csb.txt", envnames=envnames, statenames = statenames)
    
onlyres_SSTART156_struclist <- csbtxtread("Model/evo_onepop/herringEBT_onlyres_SSTART156.csb.txt", envnames=envnames, statenames = statenames)
    
# summarize population structures
onlyres_SSTART351_strucsum <- bind_rows(lapply(X=onlyres_SSTART351_struclist, FUN=summarisestruc))
onlyres_SSTART156_strucsum <- bind_rows(lapply(X=onlyres_SSTART156_struclist, FUN=summarisestruc))
    
# plot length
onlyres_dev_plot<- ggplot()+
      geom_ribbon(data=onlyres_SSTART156_strucsum, aes(x=(meanage+38.5)/365, ymin=minlen, ymax=maxlen, y=meanlen, group=coho, fill="Evolutionary repellor" ),  alpha=0.4)+
      geom_ribbon(data=onlyres_SSTART351_strucsum, aes(x=(meanage+38.5)/365, ymin=minlen, ymax=maxlen, y=meanlen, group=coho, fill="Evolutionary attractor"), alpha=0.4)+
      geom_hline(yintercept = 14, linetype="dashed")+
      theme_classic()+
      scale_y_continuous(limits=c(0, 17), name = "Herring length (cm)")+
      scale_x_continuous(limits=c(0, 5), name = element_blank())+
      scale_fill_manual(values=c(evo_pallet[2],evo_pallet[3]), name=element_blank())+
      theme( legend.position=c(0.75,0.4), plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank() )
      
# plot survival 
onlyres_surv_plot <- ggplot()+
      geom_line(data=onlyres_SSTART156_strucsum, aes(x=(meanage+38.5)/365, y=totdens/bdens, group=coho), colour=evo_pallet[3])+
      geom_line(data=onlyres_SSTART351_strucsum, aes(x=(meanage+38.5)/365,y=totdens/bdens, group=coho),  colour=evo_pallet[2])+
      theme_classic()+
      scale_y_continuous(trans="log10", limits=c(1E-7, 1E0), name = "Survival probability")+
      scale_x_continuous(limits=c(0, 5), name=element_blank())+
      theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank() )
    
# plot resource   
onlyres_env_plot <- ggplot()+
      theme_classic()+
      geom_line(data=onlyres_SSTART156_strucsum, aes(x=(meanage+38.5)/365, y=zoop, group=coho), colour=evo_pallet[3])+
      geom_line(data=onlyres_SSTART351_strucsum, aes(x=(meanage+38.5)/365,  y=zoop, group=coho), colour=evo_pallet[2])+
      scale_y_continuous(limits=c(0, 2.1), name="Zooplankton" )+
      scale_x_continuous(limits=c(0, 5), name=element_blank())+
      theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=5, face="bold"), axis.text.x=element_blank() )

# plot temperature   
onlyres_temp_plot <- ggplot()+
      theme_classic()+
      geom_line(data=onlyres_SSTART156_env, aes(x=(time)/365, y=temp), colour=evo_pallet[3])+
      geom_line(data=onlyres_SSTART351_env, aes(x=(time)/365, y=temp), colour=evo_pallet[2])+
      scale_y_continuous(limits=c(0, 15), name="Temperature" )+
      scale_x_continuous(limits=c(0, 5), name="Herring age (years)")+
      theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=5, face="bold") )
    
# combine plots 
onlyres_devenv_plot <- wrap_plots(onlyres_dev_plot, plot_spacer(), onlyres_surv_plot, plot_spacer(), onlyres_env_plot, plot_spacer(), onlyres_temp_plot, ncol = 1 , axis_titles = "collect_x", tag_level = 'new', heights=c(1,-0.1, 1, -0.1, 0.5, -0.1, 0.5))+
      plot_annotation(tag_levels="A")
    
ggsave("Figures/FigureS4.pdf", onlyres_devenv_plot , width=12.1, height=20, units="cm")
    


###### FIG S3: TIMESERIES WITH ONLY TEMPERATURE DEPENDENT PHYSIOLOGY #####
    
# read summary output
onlyphys_SSTART129_env <- read.delim(paste0("Model/evo_onepop/herringEBT_onlyphys_SSTART129.out"),
                                        sep = "\t",
                                        header=FALSE,
                                        col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass", "adltnum", "adltmass","matage","maxlen","temp","light"))
    
onlyphys_SSTART254_env <- read.delim(paste0("Model/evo_onepop/herringEBT_onlyphys_SSTART254.out"),
                                         sep = "\t",
                                         header=FALSE,
                                         col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass", "adltnum", "adltmass","matage","maxlen","temp","light"))
    
# read in full population structure
onlyphys_SSTART129_struclist <- csbtxtread("Model/evo_onepop/herringEBT_onlyphys_SSTART129.csb.txt", envnames=envnames, statenames = statenames)
    
onlyphys_SSTART254_struclist <- csbtxtread("Model/evo_onepop/herringEBT_onlyphys_SSTART254.csb.txt", envnames=envnames, statenames = statenames)
    
# summarize population structures
onlyphys_SSTART129_strucsum <- bind_rows(lapply(X=onlyphys_SSTART129_struclist, FUN=summarisestruc))
onlyphys_SSTART254_strucsum <- bind_rows(lapply(X=onlyphys_SSTART254_struclist, FUN=summarisestruc))


# plot length
onlyphys_dev_plot<- ggplot()+
      geom_ribbon(data=onlyphys_SSTART254_strucsum, aes(x=(meanage+38.5)/365, ymin=minlen, ymax=maxlen, y=meanlen, group=coho, fill="Evolutionary repellor"),  alpha=0.4)+
      geom_ribbon(data=onlyphys_SSTART129_strucsum, aes(x=(meanage+38.5)/365, ymin=minlen, ymax=maxlen, y=meanlen, group=coho, fill="Evolutionary attractor"), alpha=0.4)+
      geom_hline(yintercept = 14, linetype="dashed")+
      theme_classic()+
      scale_y_continuous(limits=c(0, 17), name = "Herring length (cm)")+
      scale_x_continuous(limits=c(0, 5), name = element_blank())+
      scale_fill_manual(values=c(evo_pallet[2],evo_pallet[3]), name=element_blank())+
      theme( legend.position=c(0.75,0.4), plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank() )

# plot survival
onlyphys_surv_plot <- ggplot()+
      geom_line(data=onlyphys_SSTART254_strucsum, aes(x=(meanage+38.5)/365, y=totdens/bdens, group=coho), colour=evo_pallet[3])+
      geom_line(data=onlyphys_SSTART129_strucsum, aes(x=(meanage+38.5)/365,y=totdens/bdens, group=coho),  colour=evo_pallet[2])+
      theme_classic()+
      scale_y_continuous(trans="log10", limits=c(1E-7, 1E0), name = "Survival probability")+
      scale_x_continuous(limits=c(0, 5), name=element_blank())+
      theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank() )

# plot prey
onlyphys_env_plot <- ggplot()+
      theme_classic()+
      geom_line(data=onlyphys_SSTART254_strucsum, aes(x=(meanage+38.5)/365, y=zoop, group=coho), colour=evo_pallet[3])+
      geom_line(data=onlyphys_SSTART129_strucsum, aes(x=(meanage+38.5)/365,  y=zoop, group=coho), colour=evo_pallet[2])+
      scale_y_continuous(limits=c(0, 0.2), name="Zooplankton" )+
      scale_x_continuous(limits=c(0, 5), name=element_blank())+
      theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=5, face="bold"), axis.text.x=element_blank() )

# plot temperature 
onlyphys_temp_plot <- ggplot()+
      theme_classic()+
      geom_line(data=onlyphys_SSTART254_env, aes(x=(time)/365, y=temp), colour=evo_pallet[3])+
      geom_line(data=onlyphys_SSTART129_env, aes(x=(time)/365, y=temp), colour=evo_pallet[2])+
      scale_y_continuous(limits=c(0, 15), name="Temperature" )+
      scale_x_continuous(limits=c(0, 5), name="Herring age (years)")+
      theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=5, face="bold") )
    

# combine plots
onlyphys_devenv_plot <- wrap_plots(onlyphys_dev_plot, plot_spacer(), onlyphys_surv_plot, plot_spacer(), onlyphys_env_plot, plot_spacer(), onlyphys_temp_plot, ncol = 1 , axis_titles = "collect_x", tag_level = 'new', heights=c(1,-0.1, 1, -0.1, 0.5, -0.1, 0.5))+
            plot_annotation(tag_levels="A")
    
ggsave("Figures/FigureS3.pdf", onlyphys_devenv_plot , width=12.1, height=20, units="cm")
    

###### FIG S2: TIMESERIES WITH ONLY TEMPERATURE DEPENDENT EGG SURVIVAL AND DURATION #####

# read summary data
onlylarv_SSTART151_env <- read.delim(paste0("Model/evo_onepop/herringEBT_onlylarv_SSTART151.out"),
                                     sep = "\t",
                                     header=FALSE,
                                     col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass", "adltnum", "adltmass","matage","maxlen","temp","light"))

onlylarv_SSTART206_env <- read.delim(paste0("Model/evo_onepop/herringEBT_onlylarv_SSTART206.out"),
                                     sep = "\t",
                                     header=FALSE,
                                     col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass", "adltnum", "adltmass","matage","maxlen","temp","light"))

onlylarv_SSTART294_env <- read.delim(paste0("Model/evo_onepop/herringEBT_onlylarv_SSTART294.out"),
                                     sep = "\t",
                                     header=FALSE,
                                     col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juvnum","juvmass", "adltnum", "adltmass","matage","maxlen","temp","light"))


# read population structures
onlylarv_SSTART151_struclist <- csbtxtread("Model/evo_onepop/herringEBT_onlylarv_SSTART151.csb.txt", envnames=envnames, statenames = statenames)

onlylarv_SSTART206_struclist <- csbtxtread("Model/evo_onepop/herringEBT_onlylarv_SSTART206.csb.txt", envnames=envnames, statenames = statenames)

onlylarv_SSTART294_struclist <- csbtxtread("Model/evo_onepop/herringEBT_onlylarv_SSTART294.csb.txt", envnames=envnames, statenames = statenames)

# summarize population structure
onlylarv_SSTART151_strucsum <- bind_rows(lapply(X=onlylarv_SSTART151_struclist, FUN=summarisestruc))
onlylarv_SSTART206_strucsum <- bind_rows(lapply(X=onlylarv_SSTART206_struclist, FUN=summarisestruc))
onlylarv_SSTART294_strucsum <- bind_rows(lapply(X=onlylarv_SSTART294_struclist, FUN=summarisestruc))


# plot size
onlylarv_dev_plot<- ggplot()+
  geom_ribbon(data=onlylarv_SSTART206_strucsum, aes(x=(meanage+38.5)/365, ymin=minlen, ymax=maxlen, y=meanlen, group=coho, fill="Evolutionary repellor"),  alpha=0.4)+
  geom_ribbon(data=onlylarv_SSTART294_strucsum, aes(x=(meanage+38.5)/365, ymin=minlen, ymax=maxlen, y=meanlen, group=coho, fill="Evolutionary attractor"), alpha=0.4)+
  geom_ribbon(data=onlylarv_SSTART151_strucsum, aes(x=(meanage+38.5)/365, ymin=minlen, ymax=maxlen, y=meanlen, group=coho, fill="Evolutionary attractor"), alpha=0.4)+
  geom_hline(yintercept = 14, linetype="dashed")+
  theme_classic()+
  scale_y_continuous(limits=c(0, 17), name = "Herring length (cm)")+
  scale_x_continuous(limits=c(0, 5), name = element_blank())+
  scale_fill_manual(values=c(evo_pallet[2],evo_pallet[3]), name=element_blank())+
  theme(legend.position=c(0.75,0.4),  plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank() )

onlylarv_surv_plot <- ggplot()+
  geom_line(data=onlylarv_SSTART206_strucsum, aes(x=(meanage+38.5)/365, y=totdens/bdens, group=coho), colour=evo_pallet[3], alpha=0.4)+
  geom_line(data=onlylarv_SSTART294_strucsum, aes(x=(meanage+38.5)/365, y=totdens/bdens, group=coho), colour=evo_pallet[2], alpha=0.4)+
  geom_line(data=onlylarv_SSTART151_strucsum, aes(x=(meanage+38.5)/365, y=totdens/bdens, group=coho),  colour=evo_pallet[2],  alpha=0.4)+
  theme_classic()+
  scale_y_continuous(trans="log10", limits=c(1E-7, 1E0), name = "Survival probability")+
  scale_x_continuous(limits=c(0, 5), name=element_blank())+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank() )

onlylarv_env_plot <- ggplot()+
  theme_classic()+
  geom_line(data=onlylarv_SSTART206_strucsum, aes(x=(meanage+38.5)/365, y=zoop, group=coho), colour=evo_pallet[3], alpha=0.4)+
  geom_line(data=onlylarv_SSTART294_strucsum, aes(x=(meanage+38.5)/365, y=zoop, group=coho), colour=evo_pallet[2], alpha=0.4)+
  geom_line(data=onlylarv_SSTART151_strucsum, aes(x=(meanage+38.5)/365, y=zoop, group=coho),  colour=evo_pallet[2],  alpha=0.4)+
  scale_y_continuous(limits=c(0, 0.3), name="Zooplankton" )+
  scale_x_continuous(limits=c(0, 5), name=element_blank())+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=5, face="bold"), axis.text.x=element_blank() )

onlylarv_temp_plot <- ggplot()+
  theme_classic()+
  geom_line(data=onlylarv_SSTART206_env, aes(x=(time)/365, y=temp), colour=evo_pallet[3])+
  geom_line(data=onlylarv_SSTART294_env, aes(x=(time)/365, y=temp), colour=evo_pallet[2])+
  geom_line(data=onlylarv_SSTART151_env, aes(x=(time)/365, y=temp), colour=evo_pallet[2])+
  scale_y_continuous(limits=c(0, 15), name="Temperature" )+
  scale_x_continuous(limits=c(0, 5), name="Herring age (years)")+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=5, face="bold") )


onlylarv_devenv_plot <- wrap_plots(onlylarv_dev_plot, plot_spacer(), onlylarv_surv_plot, plot_spacer(), onlylarv_env_plot, plot_spacer(), onlylarv_temp_plot, ncol = 1 , axis_titles = "collect_x", tag_level = 'new', heights=c(1,-0.1, 1, -0.1, 0.5, -0.1, 0.5))+
  plot_annotation(tag_levels="A")

ggsave("Figures/FigureS2.pdf", onlylarv_devenv_plot , width=12.1, height=20, units="cm")


###### FIG 3: BIFURCATION OF DENSITY OVER FISHING INTENSITY WITH TWO POPULATIONS #####

# load simulations
twopop_up <- read.delim(paste0("Model/evo_twopop/herringEBT_2pop_bbif_up.out"),
                         sep = "\t",
                         header=FALSE,
                         col.names=c("time","time2","zoop","bent","eggbuffSS","eggbuffAS",
                                     "egglarvSS", "egglarvAS",
                                     "juvnumSS","juvnumAS", 
                                     "juvmassSS","juvmassAS", 
                                     "adltnumSS", "adltnumAS",
                                     "adltmassSS", "adltmassAS",
                                     "matageSS","matageAS",
                                     "maxlenSS", "maxlenAS",
                                     "MUF"))

twopop_down <- read.delim(paste0("Model/evo_twopop/herringEBT_2pop_bbif_down.out"),
                         sep = "\t",
                         header=FALSE,
                         col.names=c("time","time2","zoop","bent","eggbuffSS","eggbuffAS",
                                     "egglarvSS", "egglarvAS",
                                     "juvnumSS","juvnumAS", 
                                     "juvmassSS","juvmassAS", 
                                     "adltnumSS", "adltnumAS",
                                     "adltmassSS", "adltmassAS",
                                     "matageSS","matageAS",
                                     "maxlenSS", "maxlenAS",
                                     "MUF"))


twopop_b <- read.delim(paste0("Model/evo_twopop/herringEBT_2pop_bbif_up_add.out"),
                             sep = "\t",
                             header=FALSE,
                             col.names=c("time","time2","zoop","bent","eggbuffSS","eggbuffAS",
                                         "egglarvSS", "egglarvAS",
                                         "juvnumSS","juvnumAS", 
                                         "juvmassSS","juvmassAS", 
                                         "adltnumSS", "adltnumAS",
                                         "adltmassSS", "adltmassAS",
                                         "matageSS","matageAS",
                                         "maxlenSS", "maxlenAS",
                                         "MUF"))

twopop <- rbind(twopop_down[which(twopop_down$MUF < 0.0033),],
                   twopop_b[which(twopop_b$MUF < 0.0037),],
                   twopop_up[which(twopop_up$MUF > 0.0037),])


# summarize simulations
twopop_summary <- summarise(.data = twopop, 
                            meanjuvmassSS = mean(juvmassSS),
                            meanjuvmassAS = mean(juvmassAS),
                            meanadltmassSS = mean(adltmassSS),
                            meanadltmassAS = mean(adltmassAS),
                            meaneggbuffSS = mean(eggbuffSS),
                            meaneggbuffAS = mean(eggbuffAS),
                            meantotmassSS = mean(juvmassSS + adltmassSS),
                            meantotmassAS = mean(juvmassAS + adltmassAS),
                            
                            minjuvmassSS = min(juvmassSS),
                            minjuvmassAS = min(juvmassAS),
                            minadltmassSS = min(adltmassSS),
                            minadltmassAS = min(adltmassAS),
                            mineggbuffSS = min(eggbuffSS),
                            mineggbuffAS = min(eggbuffAS),
                            mintotmassSS = min(juvmassSS + adltmassSS),
                            mintotmassAS = min(juvmassAS + adltmassAS),
                            
                            maxjuvmassSS = max(juvmassSS),
                            maxjuvmassAS = max(juvmassAS),
                            maxadltmassSS = max(adltmassSS),
                            maxadltmassAS = max(adltmassAS),
                            maxeggbuffSS = max(eggbuffSS),
                            maxeggbuffAS = max(eggbuffAS),
                            maxtotmassSS = max(juvmassSS + adltmassSS),
                            maxtotmassAS = max(juvmassAS + adltmassAS),
                            .by=c("MUF"))

# discard small values leading to extinction on long term
twopop_summary$mintotmassAS[which(twopop_summary$mintotmassAS<100)] <- 0
twopop_summary$meantotmassAS[which(twopop_summary$meantotmassAS<100)] <- 0
twopop_summary$maxtotmassAS[which(twopop_summary$maxtotmassAS<100)] <- 0
twopop_summary$mintotmassSS[which(twopop_summary$mintotmassSS<100)] <- 0
twopop_summary$meantotmassSS[which(twopop_summary$meantotmassSS<100)] <- 0
twopop_summary$maxtotmassSS[which(twopop_summary$maxtotmassSS<100)] <- 0

# make figure 3
popmass_MUFbif <- ggplot(data=twopop_summary, aes(x=365*MUF))+
  geom_ribbon(aes(ymin=mintotmassAS, ymax=maxtotmassAS, fill="Autumn spawners"), alpha=0.1)+
  geom_ribbon(aes(ymin=mintotmassSS, ymax=maxtotmassSS, fill="Spring spawners"), alpha=0.1)+
  geom_point(data=twopop, aes(y=adltmassAS+juvmassAS, colour="Autumn spawners"), alpha=0.05, size=0.5)+
  geom_point(data=twopop, aes(y=adltmassSS+juvmassSS, colour="Spring spawners"), alpha=0.05, size=0.5)+
  scale_y_continuous(trans="log10", limits=c(1E2,0.7E8))+
  scale_x_continuous(limits=c(0.365, 3.25))+
  scale_fill_manual(values = eco_pallet, guide="none")+
  scale_colour_manual(values = eco_pallet)+
  theme_classic()+
  theme(legend.position = c(0.8, 0.9), legend.title = element_blank())+
  labs(x="Fishing pressure (1/y)", y="Relative population biomass")

ggsave("Figures/Figure3R.pdf", popmass_MUFbif , width=12, height=8, units="cm")


###### FIG S8: MUTUAL INVASION FITNESS OF SPRING and AUTUMN SPAWNERS #####


# importat data
inv_AS_down <- read.delim("model/evo_mutinv/herringEBT_invAS_down.out",
                        sep = "\t",
                        header=FALSE,
                        col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitAS","LROres","LROAS", "GTewa", "GTAS","MUF"))

inv_AS_up <- read.delim("model/evo_mutinv/herringEBT_invAS_up.out",
                          sep = "\t",
                          header=FALSE,
                          col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitAS","LROres","LROAS", "GTewa", "GTAS","MUF"))

inv_AS <- rbind(inv_AS_down, inv_AS_up)

inv_SS_down <- read.delim("model/evo_mutinv/herringEBT_invSS_down.out",
                          sep = "\t",
                          header=FALSE,
                          col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitSS","LROres","LROSS", "GTewa", "GTSS","MUF"))

inv_SS_up <- read.delim("model/evo_mutinv/herringEBT_invSS_up.out",
                        sep = "\t",
                        header=FALSE,
                        col.names=c("time","time2","zoop","bent","eggbuff","egglarv","juv","adlt","matage","maxlen","temp","Fitres","FitSS","LROres","LROSS", "GTewa", "GTSS","MUF"))

inv_SS <- rbind(inv_SS_down, inv_SS_up)


invfit_plot <- ggplot()+
  geom_rect( aes(xmin=-Inf, xmax=365*0.00445, ymin=-Inf, ymax=Inf), fill="gray90")+
  geom_hline(yintercept=0)+
  geom_point(data=inv_SS, aes(x=365*MUF, y=FitSS, colour="Spring spawners"))+
  geom_point(data=inv_AS, aes(x=365*MUF, y=FitAS, colour="Autumn spawners"))+
  scale_y_continuous(limits=c(-1,1))+
  scale_x_continuous(limits=c(365*0.001, 365*0.009))+
  scale_fill_manual(values = eco_pallet, guide="none")+
  scale_colour_manual(values = eco_pallet)+
  theme_classic()+
  theme(legend.position = c(0.8, 0.9), legend.title = element_blank())+
  labs(x="Fishing pressure (1/y)", y="Invasion fitness")
 
ggsave("Figures/FigureS8.pdf", invfit_plot, width=12, height=8, units="cm")




###### FIG S5: EGG SURVIVAL AT ESS #####

# make seasonal temperature pattern
TC = 8.588

egglarv <- data.frame(doy = 1:730)

egglarv$TA <- 7.69 + 5.914 * sin( 2*pi/365 * (egglarv$doy - 273.1)  )
egglarv$TS <- 172.1 + 47.06 * sin( 2*pi/365 * (egglarv$doy - 4.6)  )
egglarv$temp <- TC+ egglarv$TA * sin( 2*pi/365 * (egglarv$doy - egglarv$TS) )

# calculate temperature dependent egg survival
egglarv$tempsurv <- exp(-(-0.5029)/(8.6173E-5)*(1/(egglarv$temp+273.15) - 1/(8.588+273.15)))/(1 + exp(-(-8.8534)/(8.6173E-5)*(1/(egglarv$temp+273.15) - 1/(4.5019+273.15))) )

# calculate duration dependent survival
egglarv$duration <- 25 + 4.2 + 9.08 * exp(-(-1.531)/(8.6173E-5)*(1/(egglarv$temp+273.15)- 1/(8.588+273.15)))

egglarv$dursurv <- exp(-0.03 * egglarv$duration )

# combine survival rates
egglarv$surv <- egglarv$tempsurv * egglarv$dursurv

# make figure
egglarv_surv_plot <- ggplot()+
  geom_rect(aes(xmin=as.Date(206), xmax=as.Date(206+76), ymin=-Inf, ymax=Inf),  fill=evo_pallet[3], alpha=0.2)+
  geom_rect(aes(xmin=as.Date(294), xmax=as.Date(294+76), ymin=-Inf, ymax=Inf), fill=evo_pallet[2], alpha=0.2)+
  geom_rect(aes(xmin=as.Date(151), xmax=as.Date(151+76), ymin=-Inf, ymax=Inf), fill=evo_pallet[2], alpha=0.2)+
  geom_line(data=egglarv, aes(x=as.Date(doy), y=surv), linewidth=1)+
  scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)), name=element_blank() )+
  scale_y_continuous(limits=c(0,0.35), name="Total egg and larval survival")+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=-0.4, vjust=2, face="bold") )

egglarv_dur_plot <-ggplot()+
  geom_rect(aes(xmin=as.Date(206), xmax=as.Date(206+76), ymin=-Inf, ymax=Inf),  fill=evo_pallet[3], alpha=0.2)+
  geom_rect(aes(xmin=as.Date(294), xmax=as.Date(294+76), ymin=-Inf, ymax=Inf), fill=evo_pallet[2], alpha=0.2)+
  geom_rect(aes(xmin=as.Date(151), xmax=as.Date(151+76), ymin=-Inf, ymax=Inf), fill=evo_pallet[2], alpha=0.2)+
  geom_line(data=egglarv, aes(x=as.Date(doy), y=duration), linewidth=1)+
  scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)), name=element_blank() )+
  scale_y_continuous(limits=c(30,70), name="Total egg and larval\nperiod duration (days)")+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=-0.4, vjust=2, face="bold") )

egglarv_temp_plot <-ggplot()+
  geom_rect(aes(xmin=as.Date(206), xmax=as.Date(206+76), ymin=-Inf, ymax=Inf),  fill=evo_pallet[3], alpha=0.2)+
  geom_rect(aes(xmin=as.Date(294), xmax=as.Date(294+76), ymin=-Inf, ymax=Inf), fill=evo_pallet[2], alpha=0.2)+
  geom_rect(aes(xmin=as.Date(151), xmax=as.Date(151+76), ymin=-Inf, ymax=Inf), fill=evo_pallet[2], alpha=0.2)+
  geom_line(data=egglarv, aes(x=as.Date(doy), y=temp), linewidth=1)+
  scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)), name="Day of birth" )+
  scale_y_continuous(limits=c(0,15), name="Temperature\n(Celcius)")+
  theme_classic()+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=-0.4, vjust=2, face="bold") )

egglarv_plot <- wrap_plots(egglarv_surv_plot, plot_spacer(), egglarv_dur_plot, plot_spacer(), egglarv_temp_plot,  ncol = 1 , axis_titles = "collect_x", tag_level = 'new', heights=c(1,-0.1, 1, -0.1, 0.5))+
  plot_annotation(tag_levels="A") 

ggsave("Figures/FigureS5R.pdf", egglarv_plot , width=12.1, height=15, units="cm")



###### FIG S9: TEMPERATURE DEPENDENT EGG SURVIVAL  #####

# make temperature pattern
egglarvtemp <- expand.grid(doy=1:730, TC=seq(5, 15, 1) )

egglarvtemp$TA <- 7.69 + 5.914 * sin( 2*pi/365 * (egglarvtemp$doy - 273.1)  )
egglarvtemp$TS <- 172.1 + 47.06 * sin( 2*pi/365 * (egglarvtemp$doy - 4.6)  )
egglarvtemp$temp <- egglarvtemp$TC+ egglarvtemp$TA * sin( 2*pi/365 * (egglarvtemp$doy - egglarvtemp$TS) )

# calculate temperature dependent survival
egglarvtemp$tempsurv <- exp(-(-0.5029)/(8.6173E-5)*(1/(egglarvtemp$temp+273.15) - 1/(8.588+273.15)))/(1 + exp(-(-8.8534)/(8.6173E-5)*(1/(egglarvtemp$temp+273.15) - 1/(4.5019+273.15))) )

# calculate duration dependent survival
egglarvtemp$duration <- 25 + 4.2 + 9.08 * exp(-(-1.531)/(8.6173E-5)*(1/(egglarvtemp$temp+273.15)- 1/(8.588+273.15)))

egglarvtemp$dursurv <- exp(-0.03 * egglarvtemp$duration )

# calculate total survival
egglarvtemp$surv <- egglarvtemp$tempsurv * egglarvtemp$dursurv

# make figure
egglarvtemp_surv_plot <- ggplot(data=egglarvtemp)+
  geom_line(aes(x=as.Date(doy), y=surv, colour=TC, group=TC), linewidth=1)+
  geom_line(data=egglarv, aes(x=as.Date(doy), y=surv), linewidth=1, linetype="dashed")+
  scale_x_date(date_labels = "%b", limits=as.Date(c(45+ 38,410+ 38)), name="Day of birth" )+
  scale_y_continuous(limits=c(0,0.35), name="Total egg and larval survival")+
  scale_colour_viridis(name="Central\ntemperature\n(Celcius)") +
  theme_classic()+
  theme()

ggsave("Figures/FigureS9.pdf", egglarvtemp_surv_plot , width=12.1, height=8, units="cm")




###### FIG S7: POPULATION DYNAMICS WITH TWO STRATEGIES #####

# load general output data
twopop_MUF001_env <- read.delim(paste0("Model/evo_twopop_struc/herringEBT_2pop_MUF001.out"),
                                sep = "\t",
                                header=FALSE,
                                col.names=c("time","time2","zoop","bent","eggbuff1","eggbuff2", "egglarv1","egglarv2", "juvnum1", "juvnum2","juvmass1","juvmass2", "adltnum1", "adltnum2", "adltmass1","adltmass2", "matage1","matage2","maxlen1","maxlen2"))


twopop_MUF004_env <- read.delim(paste0("Model/evo_twopop_struc/herringEBT_2pop_MUF004.out"),
                                sep = "\t",
                                header=FALSE,
                                col.names=c("time","time2","zoop","bent","eggbuff1","eggbuff2", "egglarv1","egglarv2", "juvnum1", "juvnum2","juvmass1","juvmass2", "adltnum1", "adltnum2", "adltmass1","adltmass2", "matage1","matage2","maxlen1","maxlen2"))


# load population structures
envnames <- c("time", "zoop","bent","eggbuff1","eggbuff2")
statenames <- c("xmass","ymass","age","mass","len","cond","btime","bage","matage","bdens")

twopop_MUF001_struclist <- csbtxtread("Model/evo_twopop_struc/herringEBT_2pop_MUF001.csb.txt", envnames=envnames, statenames = statenames)
twopop_MUF004_struclist <- csbtxtread("Model/evo_twopop_struc/herringEBT_2pop_MUF004.csb.txt", envnames=envnames, statenames = statenames)

# seperate out the spring and autumn spawning populations from structure list
twopop_MUF001_struclist_t <- transpose(twopop_MUF001_struclist)
twopop_MUF001_struclist_spring <- transpose(list(env=twopop_MUF001_struclist_t$env, pop1=twopop_MUF001_struclist_t$pop1))
twopop_MUF001_struclist_autumn <- transpose(list(env=twopop_MUF001_struclist_t$env, pop2=twopop_MUF001_struclist_t$pop2))


twopop_MUF004_struclist_t <- transpose(twopop_MUF004_struclist)
twopop_MUF004_struclist_spring <- transpose(list(env=twopop_MUF004_struclist_t$env, pop1=twopop_MUF004_struclist_t$pop1))
twopop_MUF004_struclist_autumn <- transpose(list(env=twopop_MUF004_struclist_t$env, pop2=twopop_MUF004_struclist_t$pop2))


# function to summazie spring structure
summarisestruc_spring  <- function(struc){
  
  envvar <- struc[[1]]
  popstruc <- as.data.frame(struc[[2]])
  
  popstruc$coho <- floor((popstruc$btime)/365)
  popstruc$totage <- popstruc$age + popstruc$bage
  popstruc$matage <- popstruc$matage + popstruc$bage
  
  popstruc$repinv <-  0.5*(popstruc$ymass - 0.7 * popstruc$xmass)/0.0007 
  popstruc$repinv <- ifelse(popstruc$len < 14, 0, popstruc$repinv)
  popstruc$relrepinv <- popstruc$repinv/popstruc$mass
  popstruc$adlt <- ifelse(popstruc$len < 14, 0, 1)
  
  popstruc_summary <- summarise(.data = popstruc, 
                                totdens = sum(density),
                                maxx=max(xmass), 
                                minx=min(xmass),
                                meanx = weighted.mean(xmass,density),
                                maxy=max(ymass), 
                                miny=min(ymass),
                                meany = weighted.mean(ymass,density),
                                maxlen=max(len), 
                                minlen=min(len),
                                meanlen = weighted.mean(len,density),
                                maxcond=max(cond), 
                                mincond=min(cond),
                                meancond=weighted.mean(cond,density),
                                maxage=max(totage), 
                                minage=min(totage),
                                meanage = weighted.mean(totage,density),
                                bdens = sum(bdens),
                                minmatage =min(matage),
                                maxmatage =max(matage),
                                meanmatage =weighted.mean(matage,density),
                                minrepinv =min(repinv ),
                                maxrepinv =max(repinv ),
                                meanrepinv =weighted.mean(repinv,density),
                                minrelrepinv =min(relrepinv),
                                maxrelrepinv =max(relrepinv),
                                meanrelrepinv =weighted.mean(relrepinv,density),
                                .by=c("coho"))
  
  
  popstruc_summary <- cbind( envvar, popstruc_summary )
  
  return(popstruc_summary)
  
} 

# function to summarize autumn struncture
summarisestruc_autumn  <- function(struc){
  
  envvar <- struc[[1]]
  popstruc <- as.data.frame(struc[[2]])
  
  popstruc$coho <- floor((popstruc$btime-301+118)/365)
  popstruc$totage <- popstruc$age + popstruc$bage
  popstruc$matage <- popstruc$matage + popstruc$bage
  
  popstruc$repinv <-  0.5*(popstruc$ymass - 0.7 * popstruc$xmass)/0.0007  
  popstruc$repinv <- ifelse(popstruc$len < 14, 0, popstruc$repinv)
  popstruc$relrepinv <- popstruc$repinv/popstruc$mass
  popstruc$adlt <- ifelse(popstruc$len < 14, 0, 1)
  
  popstruc_summary <- summarise(.data = popstruc, 
                                totdens = sum(density),
                                maxx=max(xmass), 
                                minx=min(xmass),
                                meanx = weighted.mean(xmass,density),
                                maxy=max(ymass), 
                                miny=min(ymass),
                                meany = weighted.mean(ymass,density),
                                maxlen=max(len), 
                                minlen=min(len),
                                meanlen = weighted.mean(len,density),
                                maxcond=max(cond), 
                                mincond=min(cond),
                                meancond=weighted.mean(cond,density),
                                maxage=max(totage), 
                                minage=min(totage),
                                meanage = weighted.mean(totage,density),
                                bdens = sum(bdens),
                                minmatage =min(matage),
                                maxmatage =max(matage),
                                meanmatage =weighted.mean(matage,density),
                                minrepinv =min(repinv),
                                maxrepinv =max(repinv ),
                                meanrepinv =weighted.mean(repinv,density),
                                minrelrepinv =min(relrepinv),
                                maxrelrepinv =max(relrepinv),
                                meanrelrepinv =weighted.mean(relrepinv,density),
                                .by=c("coho"))
  
  
  popstruc_summary <- cbind( envvar, popstruc_summary )
  
  return(popstruc_summary)
  
} 

# summarize structures
twopop_MUF001_strucsum_spring <- bind_rows(lapply(X=twopop_MUF001_struclist_spring, FUN=summarisestruc_spring))
twopop_MUF001_strucsum_autumn <- bind_rows(lapply(X=twopop_MUF001_struclist_autumn, FUN=summarisestruc_autumn))

twopop_MUF004_strucsum_spring <- bind_rows(lapply(X=twopop_MUF004_struclist_spring, FUN=summarisestruc_spring))
twopop_MUF004_strucsum_autumn <- bind_rows(lapply(X=twopop_MUF004_struclist_autumn, FUN=summarisestruc_autumn))


twopop_MUF001_strucsum_spring$timeage <- twopop_MUF001_strucsum_spring$time%%365 + 365*floor((twopop_MUF001_strucsum_spring$meanage+38.5)/365)
twopop_MUF001_strucsum_autumn$timeage <- twopop_MUF001_strucsum_autumn$time%%365 + 365*floor((twopop_MUF001_strucsum_autumn$meanage+38.5)/365)

twopop_MUF004_strucsum_spring$timeage <- twopop_MUF004_strucsum_spring$time%%365 + 365*floor((twopop_MUF004_strucsum_spring$meanage+38.5)/365)
twopop_MUF004_strucsum_autumn$timeage <- twopop_MUF004_strucsum_autumn$time%%365 + 365*floor((twopop_MUF004_strucsum_autumn$meanage+38.5)/365)

# size at age figures
twopop_MUF001_dev_plot <- ggplot()+
  geom_ribbon(data=twopop_MUF001_strucsum_spring, aes(x=(meanage+38.5)/365, ymin=minlen, ymax=maxlen, y=meanlen, group=coho), fill=eco_pallet[2],  alpha=0.3)+
  geom_ribbon(data=twopop_MUF001_strucsum_autumn, aes(x=(meanage+38.5)/365, ymin=minlen, ymax=maxlen, y=meanlen, group=coho), fill=eco_pallet[1],  alpha=0.3)+
  geom_hline(yintercept = 14, linetype="dashed")+
  theme_classic()+
  scale_y_continuous(limits=c(0, 19), name = "Herring length (cm)")+
  scale_x_continuous(limits=c(0, 5), name = element_blank())+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank() )+
  annotate(geom="text", x=2.5, hjust=0.5, y=19, label=bquote('Low fishing pressure ('~0.001~d^-1~')'))


twopop_MUF004_dev_plot <- ggplot()+
  geom_ribbon(data=twopop_MUF004_strucsum_spring, aes(x=(meanage+38.5)/365, ymin=minlen, ymax=maxlen, y=meanlen, group=coho, fill="Spring spawners" ),  alpha=0.3)+
  geom_ribbon(data=twopop_MUF004_strucsum_autumn, aes(x=(meanage+38.5)/365, ymin=minlen, ymax=maxlen, y=meanlen, group=coho, fill="Autumn spawners"),  alpha=0.3)+
  geom_hline(yintercept = 14, linetype="dashed")+
  theme_classic()+
  scale_y_continuous(limits=c(0, 19), name = element_blank())+
  scale_x_continuous(limits=c(0, 5), name = element_blank())+
  scale_fill_manual(values=eco_pallet, name=element_blank())+
  theme(legend.position=c(0.75, 0.4), plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank(), axis.text.y=element_blank() )+
  annotate(geom="text", x=2.5, hjust=0.5, y=19, label=bquote('Intermediate fishing pressure ('~0.004~d^-1~')'))

# reproductive investment figures
twopop_MUF001_rep_plot <- ggplot()+
  geom_point(data=twopop_MUF001_strucsum_autumn[which(twopop_MUF001_strucsum_autumn$time%%365==182 ),], aes(x=(meanage+38.5-10)/365, y=meanrepinv, group=coho), fill=eco_pallet[1], colour=eco_pallet[1])+
  geom_errorbar(data=twopop_MUF001_strucsum_autumn[which(twopop_MUF001_strucsum_autumn$time%%365==182 ),], aes(x=(meanage+38.5-10)/365, ymin=minrepinv, ymax=maxrepinv, y=meanrepinv, group=coho),  colour=eco_pallet[1])+
  geom_point(data=twopop_MUF001_strucsum_spring[which(twopop_MUF001_strucsum_spring$time%%365==364 ),], aes(x=(meanage+38.5+10)/365, y=meanrepinv, group=coho), fill=eco_pallet[2], colour=eco_pallet[2])+
  geom_errorbar(data=twopop_MUF001_strucsum_spring[which(twopop_MUF001_strucsum_spring$time%%365==364 ),], aes(x=(meanage+38.5+10)/365, ymin=minrepinv, ymax=maxrepinv, y=meanrepinv, group=coho),  colour=eco_pallet[2])+
  theme_classic()+
  scale_y_continuous(limits=c(0, 60000), name = "Number of eggs")+
  scale_x_continuous(limits=c(0-10/365, 5+10/365), name = element_blank())+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank(), legend.position="none" )

twopop_MUF004_rep_plot <- ggplot()+
  geom_point(data=twopop_MUF004_strucsum_autumn[which(twopop_MUF004_strucsum_autumn$time%%365==182 ),], aes(x=(meanage+38.5-10)/365, y=meanrepinv, group=coho), fill=eco_pallet[1], colour=eco_pallet[1])+
  geom_errorbar(data=twopop_MUF004_strucsum_autumn[which(twopop_MUF004_strucsum_autumn$time%%365==182 ),], aes(x=(meanage+38.5-10)/365, ymin=minrepinv, ymax=maxrepinv, y=meanrepinv, group=coho),  colour=eco_pallet[1])+
  geom_point(data=twopop_MUF004_strucsum_spring[which(twopop_MUF004_strucsum_spring$time%%365==364 ),], aes(x=(meanage+38.5+10)/365, y=meanrepinv, group=coho), fill=eco_pallet[2], colour=eco_pallet[2])+
  geom_errorbar(data=twopop_MUF004_strucsum_spring[which(twopop_MUF004_strucsum_spring$time%%365==364 ),], aes(x=(meanage+38.5+10)/365, ymin=minrepinv, ymax=maxrepinv, y=meanrepinv, group=coho),  colour=eco_pallet[2])+
  theme_classic()+
  scale_y_continuous(limits=c(0, 60000), name = element_blank())+
  scale_x_continuous(limits=c(0-10/365, 5+10/365), name = element_blank())+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank(),  axis.text.y=element_blank(), legend.position="none" )

# survival figures
twopop_MUF001_surv_plot <- ggplot()+
  geom_line(data=twopop_MUF001_strucsum_spring, aes(x=(meanage+38.5)/365, y=totdens/bdens, group=coho), colour=eco_pallet[2])+
  geom_line(data=twopop_MUF001_strucsum_autumn, aes(x=(meanage+38.5)/365,y=totdens/bdens, group=coho),  colour=eco_pallet[1])+
  theme_classic()+
  scale_y_continuous(trans="log10", limits=c(1E-10, 1E0), name = "Survival probability")+
  scale_x_continuous(limits=c(0, 5), name=element_blank())+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank() )

twopop_MUF004_surv_plot <- ggplot()+
  geom_line(data=twopop_MUF004_strucsum_spring, aes(x=(meanage+38.5)/365, y=totdens/bdens, group=coho), colour=eco_pallet[2])+
  geom_line(data=twopop_MUF004_strucsum_autumn, aes(x=(meanage+38.5)/365,y=totdens/bdens, group=coho),  colour=eco_pallet[1])+
  theme_classic()+
  scale_y_continuous(trans="log10", limits=c(1E-10, 1E0), name = element_blank())+
  scale_x_continuous(limits=c(0, 5), name=element_blank())+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=2, face="bold"), axis.text.x=element_blank(), axis.text.y=element_blank() )


# prey figures
twopop_MUF001_env_plot <- ggplot()+
  theme_classic()+
  geom_line(data=twopop_MUF001_strucsum_spring, aes(x=(meanage+38.5)/365, y=zoop, group=coho), colour=eco_pallet[2])+
  geom_line(data=twopop_MUF001_strucsum_autumn, aes(x=(meanage+38.5)/365,  y=zoop, group=coho), colour=eco_pallet[1])+
  scale_y_continuous(limits=c(0, 2), name="Zooplankton" )+
  scale_x_continuous(limits=c(0, 5), name="Herring age (years)")+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=5, face="bold") )

twopop_MUF004_env_plot <- ggplot()+
  theme_classic()+
  geom_line(data=twopop_MUF004_strucsum_spring, aes(x=(meanage+38.5)/365, y=zoop, group=coho), colour=eco_pallet[2])+
  geom_line(data=twopop_MUF004_strucsum_autumn, aes(x=(meanage+38.5)/365,  y=zoop, group=coho), colour=eco_pallet[1])+
  scale_y_continuous(limits=c(0, 2), name=element_blank() )+
  scale_x_continuous(limits=c(0, 5), name="Herring age (years)" )+
  theme( plot.tag.location="plot" ,  plot.tag = element_text(hjust=0.5, vjust=5, face="bold"), axis.text.y=element_blank() )


twopop_plot <- wrap_plots(twopop_MUF001_dev_plot, plot_spacer(), twopop_MUF004_dev_plot, 
                          plot_spacer(), plot_spacer(), plot_spacer(),
                          twopop_MUF001_rep_plot, plot_spacer(), twopop_MUF004_rep_plot, 
                          plot_spacer(), plot_spacer(), plot_spacer(),
                          twopop_MUF001_surv_plot, plot_spacer(), twopop_MUF004_surv_plot, 
                          plot_spacer(), plot_spacer(), plot_spacer(),
                          twopop_MUF001_env_plot,  plot_spacer(), twopop_MUF004_env_plot, 
                          ncol = 3,  tag_level = 'new', 
                          heights=c(1,-0.1, 1, -0.1, 1,  -0.1, 0.5),
                          widths =c(1, -0.1, 1) )+
  plot_annotation(tag_levels="A")



ggsave("Figures/FigureS7R.pdf", twopop_plot , width=18, height=20, units="cm")
