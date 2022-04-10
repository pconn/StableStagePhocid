### Calculate stable stage distributions for bearded, ribbon, spotted, and ringed seals in the 
### Bering and Chukchi Seas using Leslie matrix models and examine possible effects on aerial 
### survey correction factors

#  This code relies on TDR data and pre-calculated survival and reproductive schedules
#  The TDR data are from London et al. 2021 and codified in .RData objects
#  As with London et al. 2021, we rely on the "glmmLDTS" R package to fit alternative haul-out models,
#  including those that do and do not account for sex-age structure.  As of September 2021,
#  glmmLDTS can be downloaded from https://github.com/jayverhoef/glmmLDTS_package

#  We assume that the user is executing R from our repositories "home" directory (if, for example,
#  it is cloned from github).  If executing from another location, some of the path strings may need to
#  be changed.

library(glmmLDTS)

Maturity = read.csv("./StableStagePhocid/data/Maturity.csv")
Survival = read.csv("./StableStagePhocid/data/Survival_ests.csv")
Reprod = read.csv("./StableStagePhocid/data/Reproduction_table.csv")

remotes::install_github('jmlondon/berchukFits') #London et al. 2021 availability predictions
library(berchukFits)
data(ribbon_fit)
data(spotted_fit)

# set up leslie matrices - via an array (4 matrices, one for each species)
A = array(0,dim=c(4,40,40))  
for(isp in 1:4){
  for(iage in 1:39){
    A[isp,iage+1,iage]=Survival[iage,isp+1]  #assume post-breeding census
  }
}

#reproduction; nb: adults have to survive to next spring to reproduce
# nb: Leslie matrices are "female only" and assume a 50/50 sex ratio at birth
A[1,1,2:10]=0.5*Reprod$bearded*Survival[2:10,2]  
A[1,1,11:40]=A[1,1,10]
A[2,1,2:7]=0.5*Reprod$ribbon[1:6]*Survival[2:7,3]
A[2,1,8:40]=A[2,1,7]
A[3,1,2:10]=0.5*Reprod$ringed*Survival[2:10,4]
A[3,1,11:40]=A[3,1,10]
A[4,1,2:10]=0.5*Reprod$spotted*Survival[2:10,5]
A[4,1,11:40]=A[4,1,10]


#stable age distributions
Stable = matrix(0,4,40)
Stable[1,] = eigen(A[1,,])$vectors[,1]
Stable[1,] = Stable[1,]/sum(Stable[1,]) 
Stable[2,] = eigen(A[2,,])$vectors[,1]
Stable[2,] = Stable[2,]/sum(Stable[2,])
Stable[3,] = eigen(A[3,,])$vectors[,1]
Stable[3,] = Stable[3,]/sum(Stable[3,]) 
Stable[4,] = eigen(A[4,,])$vectors[,1]
Stable[4,] = Stable[4,]/sum(Stable[4,])
Stable= matrix(as.numeric(Stable),4,40)

### Maturity & stable stage distributions (converts stable age to expected proportions
#  of young-of-year, subadult, adult).  Note maturity vectors do not include YOY (first entry is for yearlings)
Stage = array(0,dim=c(4,2,3))  #species, sex (M/F), stage class

Stable_trunc = Stable[,1:10]
for(isp in 1:4)Stable_trunc[isp,10]=Stable_trunc[isp,10]+sum(Stable[isp,11:40])

Stage[1,,1] = Stable_trunc[1,1]
Stage[1,1,2]=sum((1-Maturity$Bearded.male)*Stable_trunc[1,2:10])
Stage[1,1,3]=sum(Maturity$Bearded.male*Stable_trunc[1,2:10])
Stage[1,2,2]=sum((1-Maturity$Bearded.fem)*Stable_trunc[1,2:10])
Stage[1,2,3]=sum(Maturity$Bearded.fem*Stable_trunc[1,2:10])

Stage[2,1,1]=Stable_trunc[2,1]
Stage[2,1,2]=sum((1-Maturity$Ribbon.male)*Stable_trunc[2,2:10])
Stage[2,1,3]=sum(Maturity$Ribbon.male*Stable_trunc[2,2:10])
Stage[2,2,1]=Stable_trunc[2,1]
Stage[2,2,2]=sum((1-Maturity$Ribbon.fem)*Stable_trunc[2,2:10])
Stage[2,2,3]=sum(Maturity$Ribbon.fem*Stable_trunc[2,2:10])

Stage[3,,1]=Stable_trunc[3,1]
Stage[3,1,2]=sum((1-Maturity$Ringed.male)*Stable_trunc[3,2:10])
Stage[3,1,3]=sum(Maturity$Ringed.male*Stable_trunc[3,2:10])
Stage[3,2,2]=sum((1-Maturity$Ringed.fem)*Stable_trunc[3,2:10])
Stage[3,2,3]=sum(Maturity$Ringed.fem*Stable_trunc[3,2:10])

Stage[4,,1]=Stable_trunc[4,1]
Stage[4,1,2]=sum((1-Maturity$Spotted.male)*Stable_trunc[4,2:10])
Stage[4,1,3]=sum(Maturity$Spotted.male*Stable_trunc[4,2:10])
Stage[4,2,2]=sum((1-Maturity$Spotted.fem)*Stable_trunc[4,2:10])
Stage[4,2,3]=sum(Maturity$Spotted.fem*Stable_trunc[4,2:10])

Stage = Stage*0.5 #so that age-sex proportions sum to 1.0

#plot stage classes using bar chart 
Plot_df = data.frame(Proportion=as.numeric(Stage))
Plot_df$Species = rep(c("Bearded","Ribbon","Ringed","Spotted"),6)
Plot_df$Sex = rep(c(rep("Male",4),rep("Female",4)),3)
Plot_df$Stage = c(rep("YOY",8),rep("Sub-adult",8),rep("Adult",8))
Plot_df$Stage = factor(Plot_df$Stage,levels=c("YOY","Sub-adult","Adult")) #change so order of x-axis isn't alphabetical
library(ggplot2)
SS_plot = ggplot(data=Plot_df,aes(x=Stage,y=Proportion,fill=Sex))+
  geom_bar(stat="identity",position = 'dodge') + 
  geom_text(aes(label=sprintf("%0.2f", round(Proportion, digits = 2))), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_wrap(~Species)+
  scale_fill_manual(values=c("black","gray"))+
  ylim(0,0.35)


####################
#  fit GLMPMs to spotted and ribbon seal data w/ (1) just DOY and TOD and (2) DOY, TOD, and age-sex class.
#  then produce predictions at solar noon and adjust with stable stage distributions to compare
####################
HO_ribbon = ribbon_fit$dataset
AS.ribbon <- glmmLDTS(fixed.formula = dry ~ age_sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                        sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                        sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2 +   
                        age_sex:day + age_sex:day2 + age_sex:day3,
                      random.formula = dry ~ speno,
                      data = HO_ribbon,
                      EstMeth="REML",
                      timecol = "time_vec", 
                      group.vec = "ar1_id")
Data_no_pup = HO_ribbon[-which(HO_ribbon$age=="YOUNG OF YEAR"),]
Data_no_pup$speno = factor(Data_no_pup$speno)  #reset factors 
Data_no_pup$ar1_id = factor(Data_no_pup$ar1_id) #reset factors
Const.ribbon <- glmmLDTS(fixed.formula = dry ~ sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                           sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                           sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2,
                         random.formula = dry ~ speno,
                         data = Data_no_pup,
                         EstMeth="ML",
                         timecol = "time_vec", 
                         group.vec = "ar1_id")

HO_spotted = spotted_fit$dataset

AS.spotted <- glmmLDTS(fixed.formula = dry ~ age_sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                         sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                         sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2 +   
                         age_sex:day + age_sex:day2 + age_sex:day3,
                       random.formula = dry ~ speno,
                       data = HO_spotted,
                       EstMeth="REML",
                       timecol = "time_vec", 
                       group.vec = "ar1_id")
Data_no_pup = HO_spotted[-which(HO_ribbon$age=="YOUNG OF YEAR"),]
Data_no_pup$speno = factor(Data_no_pup$speno)  #reset factors 
Data_no_pup$ar1_id = factor(Data_no_pup$ar1_id) #reset factors
Const.spotted <- glmmLDTS(fixed.formula = dry ~ sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                            sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                            sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2,
                          random.formula = dry ~ speno,
                          data = Data_no_pup,
                          EstMeth="ML",
                          timecol = "time_vec", 
                          group.vec = "ar1_id")
## ribbon seal predictions at solar noon
#  1) constant
data = Const.ribbon$dataset  
FE = Const.ribbon$fixed.effects
AS.vec = c("SUBADULT","ADULT.F","ADULT.M","YOUNG OF YEAR")
day.start = 91
day.end = 151
n.days=day.end-day.start+1
n.days.ribbon=n.days
Day = (c(day.start:day.end)-120)/10
Day2 = Day^2
Day3 = Day^3
hr = 12
Sin1 = sin(pi*hr/12)
Cos1 = cos(pi*hr/12)
Sin2 = sin(pi*hr/6)
Cos2 = cos(pi*hr/6)
Sin3 = sin(pi*hr/4)
Cos3 = cos(pi*hr/4)


const.eff = FE[FE[,"effect"]=="intercept","estimate"] +     
  FE[FE[,"effect"]=="sin1","estimate"]*Sin1 +
  FE[FE[,"effect"]=="cos1","estimate"]*Cos1 +
  FE[FE[,"effect"]=="sin2","estimate"]*Sin2 +
  FE[FE[,"effect"]=="cos2","estimate"]*Cos2 +
  FE[FE[,"effect"]=="sin3","estimate"]*Sin3 +
  FE[FE[,"effect"]=="cos3","estimate"]*Cos3
Pred = rep(const.eff,n.days)

for(iday in 1:n.days){
  Pred[iday]=Pred[iday]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1*Day[iday]+
    FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1*Day2[iday]+
    FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2*Day[iday]+
    FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2*Day2[iday]+
    FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3*Day[iday]+
    FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3*Day2[iday]+
    FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1*Day[iday]+
    FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1*Day2[iday]+
    FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2*Day[iday]+
    FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2*Day2[iday]+
    FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3*Day[iday]+
    FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3*Day2[iday]
}
Pred.ribbon.const = plogis(Pred)

# 1b) ribbon - age specific
FE = AS.ribbon$fixed.effects

const.eff = FE[FE[,"effect"]=="intercept","estimate"] +     
  FE[FE[,"effect"]=="sin1","estimate"]*Sin1 +
  FE[FE[,"effect"]=="cos1","estimate"]*Cos1 +
  FE[FE[,"effect"]=="sin2","estimate"]*Sin2 +
  FE[FE[,"effect"]=="cos2","estimate"]*Cos2 +
  FE[FE[,"effect"]=="sin3","estimate"]*Sin3 +
  FE[FE[,"effect"]=="cos3","estimate"]*Cos3
Pred = matrix(const.eff,4,n.days)

for(iday in 1:n.days){
  Pred[,iday]=Pred[,iday]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1*Day[iday]+
    FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1*Day2[iday]+
    FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2*Day[iday]+
    FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2*Day2[iday]+
    FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3*Day[iday]+
    FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3*Day2[iday]+
    FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1*Day[iday]+
    FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1*Day2[iday]+
    FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2*Day[iday]+
    FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2*Day2[iday]+
    FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3*Day[iday]+
    FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3*Day2[iday]
}

for(iage in 1:4){
  Pred[iage,]=Pred[iage,]+FE[FE[,"levels"]==AS.vec[iage],"estimate"]
}

for(iday in 1:n.days){
  for(iage in 1:4){
    Pred[iage,iday]=Pred[iage,iday]+
      FE[FE[,"effect"]=="age_sex:day" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day[iday]+
      FE[FE[,"effect"]=="age_sex:day2" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day2[iday]+
      FE[FE[,"effect"]=="age_sex:day3" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day3[iday]
  }
}

Pred = plogis(Pred)

Stage_no_pup = Stage[2,,2:3] #take out pups
Stage_no_pup = Stage_no_pup/sum(Stage_no_pup) #renormalize
Pred.ribbon.AS = Pred[1,]*sum(Stage_no_pup[,1]) + Pred[2,]*Stage_no_pup[2,2] + Pred[3,]*Stage_no_pup[1,2] 


#### spotted seals

## spotted seal predictions at solar noon
#  1) constant
data = Const.spotted$dataset  
FE = Const.spotted$fixed.effects
AS.vec = c("SUBADULT","ADULT.F","ADULT.M","YOUNG OF YEAR")
day.start = 91
day.end = 151
n.days=day.end-day.start+1
Day = (c(day.start:day.end)-120)/10
Day2 = Day^2
Day3 = Day^3
hr = 12
Sin1 = sin(pi*hr/12)
Cos1 = cos(pi*hr/12)
Sin2 = sin(pi*hr/6)
Cos2 = cos(pi*hr/6)
Sin3 = sin(pi*hr/4)
Cos3 = cos(pi*hr/4)


const.eff = FE[FE[,"effect"]=="intercept","estimate"] +     
  FE[FE[,"effect"]=="sin1","estimate"]*Sin1 +
  FE[FE[,"effect"]=="cos1","estimate"]*Cos1 +
  FE[FE[,"effect"]=="sin2","estimate"]*Sin2 +
  FE[FE[,"effect"]=="cos2","estimate"]*Cos2 +
  FE[FE[,"effect"]=="sin3","estimate"]*Sin3 +
  FE[FE[,"effect"]=="cos3","estimate"]*Cos3
Pred = rep(const.eff,n.days)

for(iday in 1:n.days){
  Pred[iday]=Pred[iday]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1*Day[iday]+
    FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1*Day2[iday]+
    FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2*Day[iday]+
    FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2*Day2[iday]+
    FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3*Day[iday]+
    FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3*Day2[iday]+
    FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1*Day[iday]+
    FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1*Day2[iday]+
    FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2*Day[iday]+
    FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2*Day2[iday]+
    FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3*Day[iday]+
    FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3*Day2[iday]
}
Pred.spotted.const = plogis(Pred)

# 1b) spotted - age specific
FE = AS.spotted$fixed.effects

const.eff = FE[FE[,"effect"]=="intercept","estimate"] +     
  FE[FE[,"effect"]=="sin1","estimate"]*Sin1 +
  FE[FE[,"effect"]=="cos1","estimate"]*Cos1 +
  FE[FE[,"effect"]=="sin2","estimate"]*Sin2 +
  FE[FE[,"effect"]=="cos2","estimate"]*Cos2 +
  FE[FE[,"effect"]=="sin3","estimate"]*Sin3 +
  FE[FE[,"effect"]=="cos3","estimate"]*Cos3
Pred = matrix(const.eff,4,n.days)

for(iday in 1:n.days){
  Pred[,iday]=Pred[,iday]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1*Day[iday]+
    FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1*Day2[iday]+
    FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2*Day[iday]+
    FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2*Day2[iday]+
    FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3*Day[iday]+
    FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3*Day2[iday]+
    FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1*Day[iday]+
    FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1*Day2[iday]+
    FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2*Day[iday]+
    FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2*Day2[iday]+
    FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3*Day[iday]+
    FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3*Day2[iday]
}

for(iage in 1:3){  #not including pups in weighted estimate
  Pred[iage,]=Pred[iage,]+FE[FE[,"levels"]==AS.vec[iage],"estimate"]
}

for(iday in 1:n.days){
  for(iage in 1:3){  #not including pups in estimate
    Pred[iage,iday]=Pred[iage,iday]+
      FE[FE[,"effect"]=="age_sex:day" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day[iday]+
      FE[FE[,"effect"]=="age_sex:day2" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day2[iday]+
      FE[FE[,"effect"]=="age_sex:day3" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day3[iday]
  }
}

Pred = plogis(Pred)
Stage_no_pup = Stage[4,,2:3] #take out pups
Stage_no_pup = Stage_no_pup/sum(Stage_no_pup) #renormalize
Pred.spotted.AS = Pred[1,]*sum(Stage_no_pup[,1]) + Pred[2,]*Stage_no_pup[2,2] + Pred[3,]*Stage_no_pup[1,2] 


#plot AS-specific vs. constant AS mean haulout behavior
library(ggplot2)
Plot.df = data.frame(Day=c(rep(c(day.start:day.end),2),rep(c(day.start:day.end),2)),
                     HO=c(Pred.ribbon.const,Pred.ribbon.AS,Pred.spotted.const,Pred.spotted.AS),
                     Species = c(rep("Ribbon",n.days.ribbon*2),rep("Spotted",n.days*2)),
                     Type = c(rep("No age-sex",n.days.ribbon),rep("Stable stage",n.days.ribbon),rep("No age-sex",n.days),rep("Stable stage",n.days)))
SS.plot = ggplot(Plot.df) + geom_line(size=1.2,aes(x=Day,y=HO,colour=Species,linetype=Type)) + 
  xlab("Julian day")+ylab("Haul-out probability")+ theme(text=element_text(size=14)) +
  scale_color_manual(values=c("darkgoldenrod3","darkorchid3"))+labs(title="A. Modeled availability")

Bias.df = data.frame(Day=c(c(day.start:day.end),c(day.start:day.end)),
                     Bias=c((Pred.ribbon.AS-Pred.ribbon.const)/Pred.ribbon.const,(Pred.spotted.AS-Pred.spotted.const)/Pred.spotted.const),  
                     Species = c(rep("Ribbon",n.days),rep("Spotted",n.days))
)
Bias.plot = ggplot(Bias.df) + geom_line(size=1.2,aes(x=Day,y=Bias,colour=Species)) + 
  xlab("Julian day")+ylab("Relative bias")+ theme(text=element_text(size=14)) +
  scale_color_manual(values=c("darkgoldenrod3","darkorchid3"))+labs(title="B. Expected bias in abundance")


jpeg("HO_stable_stage_bias.jpg")
gridExtra::grid.arrange(SS.plot,Bias.plot,nrow=2)
dev.off()

pdf("HO_stable_stage_bias.pdf")
gridExtra::grid.arrange(SS.plot,Bias.plot,nrow=2)
dev.off()


