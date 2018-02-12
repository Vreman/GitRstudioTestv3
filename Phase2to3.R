# This is an R script for a linear mixed effects model for estimating differences between phase II and phase III clinical trial results

# It consists of four sections:
# 1. General R codes for reference
# 2. An analysis of objective response rate (ORR)
# 3. An analysis of PFS
# 4. An analysis of OS


# 1. General R codes for reference
library(survival)
library(nlme)
library(ggplot2)
library(reshape2)
install.packages("lme4")
library(lme4)
attach(namedata)
detach(namedata)
metafor #(package meta analyses and forest plots)

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


## End of base functions

# 2. An analysis of ORR
# note that we have paired data: the ORR between a phase 2 and 3 trial is for the same drug for the same indication

# Data analysis
getwd()
rm(list=ls())

ORRdata <- read.table("DataTryoutv4.txt",header=TRUE)
# ORRdata <- read.table("DataTryoutv3.csv",header=TRUE,sep=",")

head(ORRdata)
summary(ORRdata)

attach(ORRdata)

# plot(ORR~Phase)
boxplot(ORR~Phase)
boxplot(ORR~Dose)
boxplot(ORR~Combination)
boxplot(ORR~Randomization)

interaction.plot(x.factor = Phase,trace.factor = Pair,response = ORR, xlab="Phase", ylab="ORR", legend=F)


# Test the random effects
# Random slope + intercept
fit0 <- lme(fixed=ORR~Phase+factor(Dose)+Combination+Randomization, random=~Phase|Pair, method="REML")
summary(fit0)
# Only random intercept
fit1 <- lme(fixed=ORR~Phase+factor(Dose)+Combination+Randomization, random=~1|Pair, method="REML")
summary(fit1)
# Only random slope
fit2 <- lme(fixed=ORR~Phase+factor(Dose)+Combination+Randomization, random=~Phase-1|Pair, method="REML")
summary(fit2)

anova(fit0, fit1)
anova(fit0, fit2)

# With interaction:
fit3 <- lme(fixed=ORR~Phase+factor(Dose)+Combination+Randomization+Dose:Combination, random=~Phase|Pair, method="REML")
summary(fit3)
# Only random intercept
fit4 <- lme(fixed=ORR~Phase+factor(Dose)+Combination+Randomization+Dose:Combination, random=~1|Pair, method="REML")
summary(fit4)
# Only random slope
fit5 <- lme(fixed=ORR~Phase+factor(Dose)+Combination+Randomization+Dose:Combination, random=~Phase-1|Pair, method="REML")
summary(fit5)

anova(fit3, fit4)
anova(fit3, fit5)

# Conclusion based on this data (Datatryoutv4.txt): 
# fit 3 is way better than 5: we need the random intercept.
# fit 3 is not significantly better than fit4, so we use the simpler model without the random slope (fit4)

# Test the fixed effects
fit11 <- lme(fixed=ORR~Phase+factor(Dose)+Combination+Randomization+Dose:Combination, random=~1|Pair, method="ML")
summary(fit11)
anova(fit11)

# Anova: least significant = randomization: lose that one
fit12 <- lme(fixed=ORR~Phase+factor(Dose)+Combination+Dose:Combination, random=~1|Pair, method="ML")
summary(fit12)
anova(fit12)
anova(fit11, fit12)

# Anova: least significant = phase: lose that one
fit13 <- lme(fixed=ORR~factor(Dose)+Combination+Dose:Combination, random=~1|Pair, method="ML")
summary(fit13)
anova(fit13)
anova(fit12, fit13)

# Anova: least significant = the interaction: lose that one
fit14 <- lme(fixed=ORR~factor(Dose)+Combination, random=~1|Pair, method="ML")
summary(fit14)
anova(fit14)
anova(fit13, fit14)

# Anova: both remaining items are now significant = final model 

# Compare final model to model with Phase in it:
fit15 <- lme(fixed=ORR~Phase+factor(Dose)+Combination, random=~1|Pair, method="ML")
summary(fit15)
anova(fit15)
anova(fit14, fit15)

# Anova: phase is not significant & models do not differ significantly. 
# Thus, phase is not an explanatory factor for differences in ORR (but dose & combination are)

detach(ORRdata)
