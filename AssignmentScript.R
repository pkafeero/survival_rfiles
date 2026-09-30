#Loading packages
library(rio)
library(dplyr)
library(tableone)




library(janitor)
library(ggplot2)
library(survival)
library(survminer)
library(KMsurv)
library(mexhaz)
library(flexsurv)
library(aod)
library(cmprsk)
library(tableone)
library(pec)
library(adjustedCurves)
library(riskRegression)
library(survRM2)
library(riskRegression)
library(rms)


##Importing the data.
assignment <- rio::import("./assignment_data_2026.csv") #Ensure the dataset is in the same working directory
View(assignment)

##Labeling caterical values
assignment <- assignment %>% 
  mutate(
    status = factor(status, levels=c(0,1,2), labels = c("Censored","Disease relapse","Death without relapse")),
    trtgiven = factor(trtgiven, levels=c(0,1),labels=c("Radiotherapy only","Radiotherapy plus chemotherapy")),
    sex = factor(sex, levels = c(0,1), labels = c("Male","Female")),
    clinstg =factor(clinstg, levels = c(1,2), labels = c("Stage I","Stage II")),
    medwidsi = factor(medwidsi, levels = c(0,1,2), labels = c("None","Small","Large")),
    extranod = factor(extranod, levels = c(0,1), labels = c("No","Yes"))
  )
#Creating composite outcome variable
assignment$event <- ifelse(assignment$status==0,0,1)
assignment$event <- as.factor(assignment$event)

#-------------------------------------------------------------------------------PART I
#a
#Summarising characteristics
tab1<-CreateTableOne(
  vars=c("age","sex","trtgiven","medwidsi","extranod","clinstg", "time","status"),
  strata="trtgiven",
  data=assignment,
  factorVars=c("sex","trtgiven","status","medwidsi","extranod","clinstg"),
  includeNA = FALSE,
  test = F,
  addOverall = T
)
print(tab1,showAllLevels = T, format = "word")
###Age
assignment %>% group_by(trtgiven) %>% summarise(minage=min(age), medage=median(age),maxage=max(age))

###Distribution of time and age: Do we assume normality for them?
ggplot(assignment, aes(x=age)) + geom_density()
ggplot(assignment, aes(x=time)) + geom_density()

##Ttest for time
t.test(assignment$time~assignment$trtgiven)



#-----b
##Summarising time to composite outcome: generate a binary variable for  an event
summary(assignment$time)
ggplot(assignment, aes(x=time)) + geom_density()

##Creating a composite outcome
assignment$event <- ifelse(assignment$status==0,0,1)
#Plotting Kaplan-Meier estimates for survival
km_assign <- survfit(Surv(time, event)~trtgiven, data=assignment) #Fitting Kaplan Meier estimates for survival
summary(km_assign)

#Plot Kaplan Meier estimates for surival
plot(km_assign,lty = c(1,2), col = c("navyblue","maroon"),
     xlab="Time since diagnosis (months)",
     ylab="Survival probability"     )
title(main="A",font=2)
legend("topright",c("Radiotherapy alone","Radio + Chemo"),
       col=c("navyblue","maroon"), lwd = 2,lty = c(1,2))

#Log rank test to compare survival between the two treatment groups
survdiff(Surv(time, event)~trtgiven, data=assignment)

#-----c
#Summarise time to  disease relapse and the time to death without relapse.
summary(assignment$time[assignment$status==1]) ##Disease relapse
summary(assignment$time[assignment$status==2]) ##Death without relapse

#Association between treatment status and competing outcomes using non-parametric analysis
#Fitting the Aaelen-Johansen estimator for cumulative incidence
nonp_compet <- survfit(Surv(time, as.factor(status)) ~trtgiven, data=assignment) 
summary(nonp_compet)
#Plot the cumulative incidences from the Aalen-Johansen
plot(nonp_compet,lwd=2,lty=c(1,2,1,2),col=c("navyblue","navyblue","maroon","maroon"),
     xlab = "Time since diagnosis (months)",
     ylab = "Cumulative Incidence", main="B")
legend("bottomright",c("Radio, Disease relapse","Radio+Chemo, Disease relapse", "Radio, Death without relapse","Radio+Chemo, Death without relapse"), 
       lwd=2,lty=c(1,2,1,2),col=c("navyblue","navyblue","maroon","maroon"))
#Grays test for comparing differences in cumulative incidence
cuminc(ftime=assignment$time, fstatus=assignment$status, group=assignment$trtgiven)

#-------------------------------------------------------------------------------PART 2
#-----a
#Using the event variable from 1b
#Unadjusted model
cox_p2 <- coxph(Surv(time, event)~trtgiven, data = assignment)
summary(cox_p2)

#Adjusted
cox_p2_adj <- coxph(Surv(time, event)~factor(trtgiven)+age+sex+clinstg+factor(medwidsi)+factor(extranod), data = assignment, x=T)
summary(cox_p2_adj)

#Checking assumptions
#Proportional Hazards assumption: The hazard ratio is constant across the entire period
#Check using Schoenfield residuals
schoen_mod <- coxph(Surv(time, event)~trtgiven+age+sex+clinstg+medwidsi+extranod, data = assignment)
sch.resid_p2a<-cox.zph(schoen_mod, transform = 'identity')
sch.resid_p2a
par(mfrow=c(2,3))
plot(sch.resid_p2a,col="red",lwd=2, cex.lab=1.5, cex.axis=1.5)


#Martingale residuals model without age
p2.cox.noage<-coxph(Surv(time=time,event=event)~trtgiven+sex+clinstg+factor(medwidsi)+factor(extranod),data=assignment)
mgale.res.noage<-resid(p2.cox.noage,type="martingale")
##Residuals plot without age in the model
plot(assignment$age,mgale.res.noage,
     xlab="Age",
     ylab="Martingale residuals")
lines(lowess(assignment$age,mgale.res.noage),col="red",lwd=2)
abline(h=0,lwd=2,col="grey")

#Martingale residuals model including age in Coxmodel
p2.cox.age<-coxph(Surv(time=time,event=event)~trtgiven+age+sex+clinstg+factor(medwidsi)+factor(extranod), data=assignment)
mgale.res.age<-resid(p2.cox.age,type="martingale")
#Residuals plot with age in the model
plot(assignment$age, mgale.res.age, 
     xlab="Age", 
     ylab="Martingale residuals")
lines(lowess(assignment$age, mgale.res.age),
      col="red",lwd=2)
abline(h=0,lwd=2,col="grey")


#-----b
##Create a dataset with radio therapy only and then compare a new dataset
p2b_data_radio <- expand.grid(
  trtgiven = 0, ##Assigning Radiotherapy only to trt given
  clinstg=1,
  age=50,
  sex=1, ##1 for Female
  medwidsi=1,
  extranod=1
)
#Fitting survival model for ratiotherapy only
p2b_radio <- survfit(cox_p2_adj, newdata = p2b_data_radio)
summary(p2b_radio)

#Fitting survival model for radiotherapy plus chemotherapy
p2b_data_radiochem <- p2b_data_radio # Create dataset for radiotherapy plus chemotherapy
p2b_data_radiochem$trtgiven <- 1 ##Asssigning Radiotherapy + Chemotherapy to trt given
p2b_radiochemo <- survfit(cox_p2_adj, newdata = p2b_data_radiochem) ##model for radiotherapy plus chemotherapy

#Plot estimated survival curves for the treatment strategies
plot(p2b_radio, conf.int = F, col = "black",
     xlab = "Time since diagnosis (months)",
     ylab = "Survival probability", main="Estimated survival curves for the composite outcome")
lines(p2b_radiochemo, conf.int = F, col = "grey")
legend("bottomleft", c("Radio", "Radio+Chemo"), col=c("black","grey"), lwd=2, text.width = 5)

#Estimated risks of the outcome up to 12 months
p2b_riskradio <- 1-summary(p2b_radio, times=12)$surv
p2b_riskradiochemo <- 1-summary(p2b_radiochemo, times=12)$surv
cat("The estimated risk of the outcome up to 12 months for radiotherapy alone is (", round(p2b_riskradio*100, 2), ")% \n")
cat("The estimated risk of the outcome up to 12 months for radiotherapy plus chemotherapy is (", round(p2b_riskradiochemo*100, 2), ")% \n")


#-------------------------------------------------------------------------------PART 3
#-----a
#Without adjustment
cox_p3 <- coxph(Surv(time, as.factor(status))~trtgiven, data = assignment, id=id)
summary(cox_p3)

#With adjustment
cox_p3_adj <- coxph(Surv(time, as.factor(status))~trtgiven+age+sex+clinstg+factor(medwidsi)+extranod, data = assignment, id=id)
summary(cox_p3_adj)


#-----b
##Fitting the models
#Radiotherapy alone
p3b_radio <- survfit(cox_p3_adj, newdata = p2b_data_radio)
summary(p3b_radio)
#Radiotherapy plus chemotherapy
p3b_radiochemo <- survfit(cox_p3_adj, newdata = p2b_data_radiochem)
plot(p3b_radiochemo)

##Cumulative incidence plot
plot(p3b_radio, conf.int = F, col = c("black","grey"), lwd=2, lty = 1,
     xlab = "Time from diagnosis (months)",
     ylab = "Cumulative Incidence",
     main = "Estimated cumulative curves for competing risks ", ylim = c(0,0.6))
lines(p3b_radiochemo, conf.int = F, col = c("black","grey"), lwd=2, lty = 2)
legend("topleft",c("RO: DR","RC: DR", "RO: DWR","RC: DWR"), lwd=2,lty=c(1,2,1,2),col=c("black","black","grey","grey"), text.width = 6)

##Cumulative incidence at Month 12
summary(p3b_radio, times=c(1,12,24))
summary(p3b_radiochemo, times=c(1,12,24))


#-------------------------------------------------------------------------------PART 4
#-----a: Marginal Risks
#Based on the model in 2b: cox_p2_adj
assignment$trtgiven <- as.factor(assignment$trtgiven)
adjsurv<-adjustedsurv(data=assignment,
                      variable="trtgiven",
                      ev_time="time",
                      event="event",
                      method="direct",
                      outcome_model=cox_p2_adj,
                      conf_int=TRUE,
                      bootstrap=F,times=12)

#Risks at 12 months
adjsurv$ate_object$meanRisk

##risk diff at time 12
adjsurv$ate_object$diffRisk

#risk ratio at time 12
adjsurv$ate_object$ratioRisk

#-----b: Marginal treatement effects


#RMST Difference
#6 months
rmst2(assignment$time, assignment$event, assignment$trtgiven, tau=6)
#12 months
rmst2(assignment$time, assignment$event, assignment$trtgiven, tau=12)
#24 months
rmst2(assignment$time, assignment$event, assignment$trtgiven, tau=24)

#RISK RATIO
####6 MONTHS
risk_radio<-1-summary(km_p4,times=c(12))$surv[1]
risk_radiochem<-1-summary(km_p4,times=c(12))$surv[2]

varrisk_radio<-summary(km_p4,times=c(12))$std.err[1]^2

varrisk_radiochem<-summary(km_p4,times=c(12))$std.err[2]^2
var.riskdiff<-varrisk_radiochem+varrisk_radio #Variance for the two independent times

lower.ci.riskdiff<-(risk_radiochem/risk_radio)-1.96*sqrt(var.riskdiff) ##Risk Ratio Lower CI
upper.ci.riskdiff<-(risk_radiochem/risk_radio)+1.96*sqrt(var.riskdiff) ##Risk Ratio Upper CI

risk_radiochem/risk_radio ##Risk Ratio

lower.ci.riskdiff
upper.ci.riskdiff

2*pnorm((risk_radiochem-risk_radio)/sqrt(var.riskdiff))


####12 MONTHS
risk_radio<-1-summary(km_p4,times=c(12))$surv[1]
risk_radiochem<-1-summary(km_p4,times=c(12))$surv[2]

varrisk_radio<-summary(km_p4,times=c(12))$std.err[1]^2

varrisk_radiochem<-summary(km_p4,times=c(12))$std.err[2]^2
var.riskdiff<-varrisk_radiochem+varrisk_radio #Variance for the two independent times

lower.ci.riskdiff<-(risk_radiochem/risk_radio)-1.96*sqrt(var.riskdiff) ##Risk Ratio Lower CI
upper.ci.riskdiff<-(risk_radiochem/risk_radio)+1.96*sqrt(var.riskdiff) ##Risk Ratio Upper CI

risk_radiochem/risk_radio ##Risk Ratio

lower.ci.riskdiff
upper.ci.riskdiff

2*pnorm((risk_radiochem-risk_radio)/sqrt(var.riskdiff))

