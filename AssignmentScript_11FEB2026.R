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
#Re-import data
assignment <- rio::import("./assignment_data_2026.csv") #File should be in your working directory.

##Creating a composite outcome
assignment$event <- ifelse(assignment$status==0,0,1)


##Summarising time to composite outcome: generate a binary variable for  an event
km_overall <- survfit(Surv(time, event)~1, data=assignment) #Model for overall outcome for both treatments
summary(km_overall)
summary(km_overall)$table

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
compet_all <- survfit(Surv(time, as.factor(status)) ~1, data=assignment) ##Model for competing risks for both treatments
summary(assignment$time[assignment$status==1]) ##Disease relapse
summary(assignment$time[assignment$status==2]) ##Death without relapse

summary(compet_all)$table
summary(compet_all$time[assignment$status==1])
summary(compet_all$time[assignment$status==2])




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
cox_p2_adj <- coxph(Surv(time, event)~factor(trtgiven)+age+sex+factor(clinstg)+factor(medwidsi)+factor(extranod), data = assignment, x=T)
summary(cox_p2_adj)

#Checking assumptions
#Proportional Hazards assumption: The hazard ratio is constant across the entire period
#Check using Schoenfield residuals
schoen_mod <- coxph(Surv(time, event)~trtgiven+age+sex+clinstg+medwidsi+extranod, data = assignment)
sch.resid_p2a<-cox.zph(schoen_mod, transform = 'identity')
sch.resid_p2a
par(mfrow=c(2,3))
plot(sch.resid_p2a,col="red",lwd=2, cex.lab=1.5, cex.axis=1.5)
par(mfrow=c(1,1))

#####Martingale reiduals plots to check format for age
#Martingale residuals model without age
p2.cox.noage<-coxph(Surv(time=time,event=event)~trtgiven+sex+clinstg+factor(medwidsi)+factor(extranod),data=assignment)
mgale.res.noage<-resid(p2.cox.noage,type="martingale")
##Residuals plot without age in the model
plot(assignment$age,mgale.res.noage,
     xlab="Age",
     ylab="Martingale residuals")
lines(lowess(assignment$age,mgale.res.noage),col="red",lwd=2)
abline(h=0,lwd=2,col="grey")

#Martingale residuals model including age in Cox model
p2.cox.age<-coxph(Surv(time=time,event=event)~trtgiven+age+sex+clinstg+factor(medwidsi)+factor(extranod), data=assignment) ###CONSIDER USING age^2!!!!
mgale.res.age<-resid(p2.cox.age,type="martingale")
#Residuals plot with age in the model
plot(assignment$age, mgale.res.age, 
     xlab="Age", 
     ylab="Martingale residuals")
lines(lowess(assignment$age, mgale.res.age),
      col="red",lwd=2)
abline(h=0,lwd=2,col="grey")


##New model with stratified Age and Sex
##Create age strata based in the quantiles.
# assignment <- assignment %>% 
#   mutate(age_strata = case_when(age < 23 ~ 1, 
#                                 age >= 23 & age < 31 ~ 2,
#                                 age >= 31 & age < 44 ~ 3,
#                                 age >= 44 ~ 4))
# 
# cox_p2_adj_new <- coxph(Surv(time, event) ~ trtgiven + factor(clinstg) + factor(medwidsi) + 
#                           extranod + strata(sex) + strata(age_strata), data = assignment, x=T)
# 
# summary(cox_p2_adj_new)
# cox.zph(cox_p2_adj_new, transform = 'identity')
# par(mfrow=c(2,3))
# plot((cox.zph(cox_p2_adj_new, transform = 'identity')),col="red",lwd=2, cex.lab=1.5, cex.axis=1.5)
# par(mfrow=c(1,1))


#-----b
##Create a dataset with radio therapy only and then compare a new dataset
p2b_data_radio <- data.frame(
  trtgiven = 0, ##Assigning Radiotherapy only to trt given
  clinstg=1,
  age=50,
  sex=1, ##1 for Female
  medwidsi=2,
  extranod=1
)
#Fitting survival model for ratiotherapy only
p2b_radio <- survfit(cox_p2_adj, newdata = p2b_data_radio)


#Fitting survival model for radiotherapy plus chemotherapy
p2b_data_radiochem <- p2b_data_radio # Create dataset for radiotherapy plus chemotherapy
p2b_data_radiochem$trtgiven <- 1 ##Asssigning Radiotherapy + Chemotherapy to trt given
p2b_radiochemo <- survfit(cox_p2_adj, newdata = p2b_data_radiochem) ##model for radiotherapy plus chemotherapy

#Plot estimated survival curves for the treatment strategies
plot(p2b_radio, conf.int = F, col = "black",
     xlab = "Time since diagnosis (months)",
     ylab = "Survival probability", main="Estimated survival curves for the composite outcome")
lines(p2b_radiochemo, conf.int = F, col = "grey")
legend("bottomleft", c("Radio", "Radio+Chemo"), col=c("black","grey"), lwd=2, text.width = 9)

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
#coxph(Surv(time, status==3)~trtgiven+age+sex+clinstg+factor(medwidsi)+extranod, data = assignment) ##CHECK WITH ROBUST==t FOR THE SEPARATE MODELS.

#-----b
##Fitting the models
#Radiotherapy alone
p3b_radio <- survfit(cox_p3_adj, newdata = p2b_data_radio)

#Radiotherapy plus chemotherapy
p3b_radiochemo <- survfit(cox_p3_adj, newdata = p2b_data_radiochem)

##Cumulative incidence plot
plot(p3b_radio, conf.int = F, col = c("black","grey"), lwd=2, lty = 1,
     xlab = "Time from diagnosis (months)",
     ylab = "Cumulative Incidence",
     main = "Estimated cumulative curves for competing risks ", ylim = c(0,0.6))
lines(p3b_radiochemo, conf.int = F, col = c("black","grey"), lwd=2, lty = 2)
legend("bottomright",c("RO: DR","RC: DR", "RO: DWR","RC: DWR"), lwd=2,lty=c(1,2,1,2),col=c("black","black","grey","grey"), text.width = 5.9)



# plot(p3b_radio,lwd=2,lty=c(1,2,1,2),col=c("navyblue","maroon"),
#       xlab = "Time since diagnosis (months)",
#       ylab = "Cumulative Incidence", main="B")
# lines(p3b_radiochemo, conf.int = F, col = c("navyblue","maroon"), lwd=2, lty = 2)
# legend("bottomright",c("Radio, Disease relapse","Radio+Chemo, Disease relapse", "Radio, Death without relapse","Radio+Chemo, Death without relapse"),
#         lwd=2,lty=c(1,2,1,2),col=c("navyblue","navyblue","maroon","maroon"))



##Cumulative incidence at Month 12
summary(p3b_radio, times=c(1,12,24))
summary(p3b_radiochemo, times=c(1,12,24))


#-------------------------------------------------------------------------------PART 4
#-----a: Marginal Risks
adjsurv<-adjustedsurv(data=assignment,
                      variable="trtgiven",
                      ev_time="time",
                      event="event",
                      method="direct",
                      outcome_model=cox_p2_adj,
                      conf_int=TRUE,
                      bootstrap=F,times=12)
#risks at time 5
adjsurv$ate_object$meanRisk


#Based on the model in 2b: cox_p2_adj
#Using Standardisation
#--
#create dataset in which trtgiven is set to 1 for everyone
p3_data_radiochemo<-assignment
p3_data_radiochemo$trtgiven<-1
#create dataset in which trtgiven is set to 0 for everyone
p3_data_radio<-assignment
p3_data_radio$trtgiven<-0


t <- survfit(cox_p2_adj,newdata=p3_data_radiochemo)
mean(summary(t, times=12)$surv)

#Predicted risk up to time 12 for each individual using cox_p2_adj
risk.radiochemo<-1-predictSurvProb(cox_p2_adj,newdata=p3_data_radiochemo,times=12)
risk.radio<-1-predictSurvProb(cox_p2_adj,newdata=p3_data_radio,times=12)


#mean survival probability at each time,averaging overall individuals
riskmean.radiochemo<-mean(risk.radiochemo)
riskmean.radio<-mean(risk.radio)
cat("The marginal risk estimate up to 12 months for the radiotherapy alone group is", riskmean.radio, "\n")
cat("The marginal risk estimate up to 12 months for the radiotherapy plus chemotherapy group is", riskmean.radiochemo)


#-----b: Marginal treatment effects
#RMST Difference
#12 months
rmst2(assignment$time, assignment$event, assignment$trtgiven, tau=12)
#rmstd <- rmst2(assignment$time, assignment$event, assignment$trtgiven, tau=12, covariates = assignment[,c("age", "sex","clinstg","medwidsi","extranod")]) #-- based on the cox model
#OR rmst2(assignment$time, assignment$event, assignment$trtgiven, tau=12, covariates = c(assignment$age, "sex","clinstg","medwidsi","extranod")])
#rmstdiff <- rmstd$adjusted.result[1] #-- based on the cox model
cat("The RMST difference is", 1.211,"months.")
#cat("The RMST difference is", rmstdiff,"months.") -- based on the cox model: 1.34
assignment$trtgiven <- factor(assignment$trtgiven)

##This is RMST from practical 10 though. But the values are a bit different! 
#Designed to estimate covariate-adjusted survival curves (using g-computation, IPTW, AFT, Cox models, etc.) 
adjsurv_asign<-adjustedsurv(data=assignment,
                      variable="trtgiven",
                      ev_time="time",
                      event="event",
                      method="direct",
                      outcome_model=cox_p2_adj,
                      conf_int=TRUE,
                      bootstrap=T,n_boot=500)
#RMST under the two treatment strategies
adjrmst_assign<-adjusted_rmst(adjsurv_asign,
                       from=0,to=max(assignment$time[assignment$event==1]),
                       conf_int=TRUE,conf_level=0.95)
adjrmst_assign

#RMST difference
adjrmst.diff<-adjusted_rmst(adjsurv_asign,
                            from=0,to=12,
                            conf_int=TRUE,conf_level=0.95,
                            contrast= "diff")
adjrmst.diff














#RISK RATIO
#Kaplan-Meier by treatment
km_p4<-survfit(Surv(time,event)~trtgiven,data=assignment)
km_p4<-survfit(cox_p2_adj,data=assignment)
#km_p4<-survfit(Surv(time,event)~trtgiven,data=assignment) #--radio
#km_p4<-survfit(Surv(time,event)~trtgiven,data=assignment) #--radio+chemo
names(km_p4)
km_p4$surv


####12 MONTHS

##Are we using the risk: riskmean.radiochemo/riskmean.radio -- Joe suggests this!!


risk_radio<-1-summary(km_p4,times=c(12))$surv[1]


risk_radiochem<-1-summary(km_p4,times=c(12))$surv[2]




varrisk_radio<-summary(km_p4,times=c(12))$std.err[1]^2

varrisk_radiochem<-summary(km_p4,times=c(12))$std.err[2]^2
var.riskdiff<-varrisk_radiochem+varrisk_radio #Variance for the two independent times

lower.ci.riskdiff<-(risk_radiochem/risk_radio)-1.96*sqrt(var.riskdiff) ##Risk Ratio Lower CI
upper.ci.riskdiff<-(risk_radiochem/risk_radio)+1.96*sqrt(var.riskdiff) ##Risk Ratio Upper CI

risk_radiochem/risk_radio ##Risk Ratio
cat("The risk ratio for the composite outcome between the intervention and radiotherapy alone is",risk_radiochem/risk_radio, "\n")

lower.ci.riskdiff
upper.ci.riskdiff

2*pnorm((risk_radiochem-risk_radio)/sqrt(var.riskdiff))
cat("The p-value for the risk ratio is", 2*pnorm((risk_radiochem-risk_radio)/sqrt(var.riskdiff)))

