# This is an R script for a linear mixed effects model for estimating differences between phase II and phase III clinical trial results

# It consists of three analyses:
# 1. An analysis of objective response rate (ORR)
# 2. An analysis of PFS
# 3. An analysis of OS

# 1. An analysis of ORR

# note that we have paired data: the ORR between a phase 2 and 3 trial is for the same drug for the same indication
# covariates that are in the database:
# 1.1
# 1.2
# 1.3
# 1.4
# 1.5
# 1.6
# 1.7
# 1.8
# 1.9

library(survival)
library(nlme)
install.packages("lme4")
library(lme4)
attach(namedata)
detach(namedata)
tapply()
interaction.plot(x.factor=Time,trace.factor=Cow,response=protein)
plot(protein~Time)
boxplot(split(protein,Time))
histogram(~protein | Time*Diet)

lme.fit<-lme(fixed=protein~ Diet*factor(Time), random=~1|Cow)
anova(lme.fit)
summary(lme.fit)
fit0<-lme(fixed=chol~factor(time),random=~1|id)
anova(fit0)
fit2ML<-lme(fixed=chol~factor(time)+gender,random=~1|id, method="ML")
anova(fit0ML, fit2ML, fit1ML)
cow.fit <-lme(protein~factor(Time)+factor(Time):Diet,random=~Time|Cow)
cow2.fit <- lme(protein~factor(Time)+factor(Time):Diet,random=~Time-1|Cow)
anova(cow.fit,cow2.fit)
cow4.fit <- lme(protein~factor(Time)+Diet,random=~Time|Cow,method="ML")

# To get the correlation between the observations on the 11 timepoints: first get the variance covariance matrix:
vc <- getVarCov(rat.fit,individual="1",type="marginal")
cov2cor(vc[[1]])

library(lme4)
fit.lgre <- glmer(oc~factor(ground)+I(height-mean(height))+(1|father),family=binomial)
summary(fit.lgre)


Twisk1<-read.table("Twisk1.dat",header=T)
epi <- read.table("episode.txt",header=TRUE)
dalmat <- read.table("dalmatian.csv",header=TRUE,sep=",")
d <- read.table(file="BeetleEggCrowding.txt", header = TRUE)

summary(Twisk1)
head(namedata)
Twisk1$id<-factor(Twisk1$id)

getwd()

# The command for the split plot ANOVA, with time as factor, is
summary(aov(chol~factor(time)+Error(id)))


## that was it for the base functions
# Now some data analysis


ORRfakev1 <- read.table("Fakedatav4.txt",header=TRUE)
head(ORRfakev1)

Id <- ORRfakev1$Id
ORR <- ORRfakev1$ORR
Phase <- ORRfakev1$Phase
Random <- ORRfakev1$Randomized

Random

fit0 <- lme(fixed=ORR~factor(Phase), random=~Phase|Id)


summary(fit0)
fit1 <- lme(fixed=ORR~factor(Phase)+Random, random=~Phase|Id)
fit2 <- lme(fixed=ORR~factor(Phase)+Random, random=~Phase+Random|Id)

anova(fit1, fit2)

