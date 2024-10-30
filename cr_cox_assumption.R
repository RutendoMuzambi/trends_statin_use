#--------------------------------------------------------------------------------------
# Project: Statins and CVD prevention
# Author: Rutendo Muzambi
# Description:  Cox proportional hazards assumption
#--------------------------------------------------------------------------------------
rm(list=ls())

sink(paste0(excel, "cox_assumptions.txt"))

####Discontinuation
##open parquet
prim_disc <- read_parquet(paste0(datafiles, "prim_discont_fact.parquet"))

prim_disc <- prim_disc[ethnicity!="missing" & deprivation!="missing"]

prim_disc[, c("age_cat", "ethnicity", "gender", "deprivation") := lapply(.SD, as.factor), .SDcols = c("age_cat", "ethnicity", "gender", "deprivation")]
prim_disc[ethnicity=="missing", ethnicity := NA]
prim_disc[deprivation=="missing", deprivation := NA]

prim_disc[, age_cat := factor(age_cat, levels = c("25-39", "40-49", "50-59", "60-69", "70+"))]
prim_disc[, gender := factor(gender, levels = c("Male", "Female"))]
prim_disc[, ethnicity := factor(ethnicity, levels = c("White", "Black", "South Asian", "Mixed", "Other"))]
prim_disc[, deprivation := factor(deprivation, levels = c("1 - least deprived", "2", "3", "4", "5 - most deprived"))] 

library("survival")
#Adjusted Cox model
prim_model <- coxph(Surv(time, discontinuation) ~ age_cat + gender + ethnicity + deprivation,
                    data=prim_disc)
prim_model

prim.ph <- cox.zph(prim_model)
prim.ph

# chisq df       p
# age_cat      80.112  4 < 2e-16
# gender        0.949  1    0.33
# ethnicity    69.247  4 3.3e-14
# deprivation   2.845  4    0.58
# GLOBAL      154.085 13 < 2e-16

#survival plots
disc.cox <- coxph(Surv(time, discontinuation)~as.factor(ethnicity),data=prim_disc)
summary(disc.cox)
disc.survfit=survfit(disc.cox,newdata=data.frame(ethnicity=c("Black","South Asian", "Mixed", "Other")))
summary(disc.survfit)
plot(disc.survfit,mark.time=F,col=c("black","grey"),xlab="Time",ylab="Estimated survivor function")

library(survival)

# Plot the log-log survival curve
# Log-Log Survival Curve with Kaplan-Meier

km=survfit(Surv(time,discontinuation)~ethnicity,data=prim_disc)

prim.km <- survfit(Surv(time,discontinuation)~as.factor(ethnicity),data=prim_disc)
ggsurvplot(prim.km, data = prim_disc,conf.int = T,censor=F,legend.title="",
           legend.labs = c("White","Black", "South Asian", "Mixed", "Other"))

###curve separates between 3months and 6 months
ggsurvplot(prim.km, data = prim_disc,conf.int = T,censor=F,legend.title="",
           legend.labs = c("White","Black", "South Asian", "Mixed", "Other") , xlim = c(0, 0.5))

prim_disc[, time2 := time*365.25]

ggsurvplot(prim.km, data = prim_disc,conf.int = T,censor=F,legend.title="",
           legend.labs = c("White","Black", "South Asian", "Mixed", "Other") , break.x.by=2, xlim = c(115, 180))
           
plot(prim.km,fun="cloglog",xlab="time (log scale)",ylab="log(-log S(t))",
     col=c("blue","red", "green", "orange","yellow"))
legend(0.02,0,c("White","Black", "South Asian", "Mixed", "Other"),col=c("blue","red", "green", "orange", "yellow"),lty=1,cex=0.5)

ggsurvplot(prim.km, data = prim_disc,conf.int = T,fun="cloglog",censor=F,legend.title="",
           legend.labs = c("White","Black", "South Asian", "Mixed", "Other"))

ggsurvplot(prim.km, data = prim_disc,conf.int = F,fun="cloglog",loglog =T, xlim = c(1, 12),
           legend.title="",
           legend.labs = c("White","Black", "South Asian", "Mixed", "Other"))


disc.cox <- coxph(Surv(time,discontinuation)~as.factor(ethnicity),data=prim_disc)

summary(disc.cox)

#schoenfelds level off after 6 months
sch.resid=cox.zph(disc.cox, transform = 'identity')
sch.resid
plot(sch.resid)

###split follow-up time for each individual into two time periods at 6 months

disc.split=survSplit(Surv(time,discontinuation)~., data=prim_disc, cut=0.5, end="time", event="discontinuation",
                     start="time0", episode="time_period")

#
prim.disc.new=coxph(Surv(time0,time,discontinuation)~as.factor(age_cat)*as.factor(time_period) + gender +
                     as.factor(ethnicity)*as.factor(time_period) + deprivation,data=disc.split)
summary(prim.disc.new)


prim.disc.ph <- cox.zph(prim.disc.new)
prim.disc.ph

prim.disc.df <- tidy(prim.disc.new, conf.int=T, exponentiate=T)

prim.disc.df$estimate <- round(prim.disc.df$estimate, 2)
prim.disc.df$conf.low <- round(prim.disc.df$conf.low, 2)
prim.disc.df$conf.high <- round(prim.disc.df$conf.high, 2)

# Merge into a single string
prim.disc.df$conf.interval <- paste(prim.disc.df$conf.low, prim.disc.df$conf.high, sep = "-")

prim.disc.df <- prim.disc.df[, c('term', 'estimate', 'conf.interval', 'p.value')]

#load into excel 
wb <- loadWorkbook(paste0(excel, "cox_assumptions.xlsx"))
names(wb)
writeData(wb, "disc_R", prim.disc.df, startCol=1, startRow=4, rowNames=T)
saveWorkbook(wb,paste0(excel, "cox_assumptions.xlsx"), overwrite=T)

#####################################################################################
##########################RISK SCORING###############################################
#####################################################################################
cox_risk_dt <- read_parquet(paste0(datafiles, "risk_scoring_fact.parquet"))
cox_risk_dt <- cox_risk_dt[ethnicity!="missing" & deprivation!="missing"]

risk_model <- coxph(Surv(time, risk_scoring) ~ age_cat + gender + ethnicity + deprivation,
                    data=cox_risk_dt)
summary(risk_model)

#schoenfelds residuals global
risk.ph <- cox.zph(risk_model)
risk.ph

#schoenfelds level off after 6 months
sch.resid=cox.zph(risk_model, transform = 'identity')
sch.resid
plot(sch.resid)

##km
risk.km <- survfit(Surv(time,risk_scoring)~as.factor(deprivation),data=cox_risk_dt)
ggsurvplot(risk.km, data = cox_risk_dt,conf.int = T,censor=F)

#loglog plot
ggsurvplot(risk.km, data = cox_risk_dt,conf.int = T,fun="cloglog",censor=F)

#Interaction term between ethnicity and time
risk.split=survSplit(Surv(time,risk_scoring)~., data=cox_risk_dt, cut=5, end="time", event="risk_scoring", start="time0", episode="time_period")

risk.cox.model=coxph(Surv(time0,time,risk_scoring)~as.factor(age_cat)*as.factor(time_period) +
                 as.factor(gender)*as.factor(time_period) + as.factor(ethnicity)*as.factor(time_period) +
                 as.factor(deprivation)*as.factor(time_period), data=risk.split)
summary(risk.cox.model)

cox_risk_ph <- cox.zph(risk.cox.model)
cox_risk_ph

risk.cox.df <- tidy(risk.cox.model, conf.int=T, exponentiate=T)

risk.cox.df$estimate <- round(risk.cox.df$estimate, 2)
risk.cox.df$conf.low <- round(risk.cox.df$conf.low, 2)
risk.cox.df$conf.high <- round(risk.cox.df$conf.high, 2)

# Merge into a single string
risk.cox.df$conf.interval <- paste(risk.cox.df$conf.low, risk.cox.df$conf.high, sep = "-")

risk.cox.df <- risk.cox.df[, c('term', 'estimate', 'conf.interval', 'p.value')]

#load into excel 
wb <- loadWorkbook(paste0(excel, "cox_assumptions.xlsx"))
names(wb)
writeData(wb, "risk_scoring_R", risk.cox.df, startCol=1, startRow=4, rowNames=T)
saveWorkbook(wb,paste0(excel, "cox_assumptions.xlsx"), overwrite=T)

######################RE-INITIATION###################################################
######################################################################################
##open parquet
prim_restart <- read_parquet(paste0(datafiles, "prim_restart_fact.parquet"))

prim_restart <- prim_restart[ethnicity!="missing" & deprivation!="missing"]


prim_restart[, ethnicity := factor(ethnicity, levels = c("White", "Black", "South Asian", "Mixed", "Other"))]

#Adjusted Cox model
restart_model <- coxph(Surv(time, restart) ~ age_cat + gender + ethnicity + deprivation,
                    data=prim_restart)
summary(restart_model)

#schoenfelds residuals global
restart.ph <- cox.zph(restart_model)
restart.ph

#km plot

#schoenfelds level off after 6 months
sch.resid=cox.zph(restart_model, transform = 'identity')
sch.resid
plot(sch.resid)

##km
restart.km <- survfit(Surv(time,restart)~as.factor(age_cat),data=prim_restart)
ggsurvplot(restart.km, data = prim_restart,conf.int = T,censor=F, break.x.by=1)

#loglog plot
ggsurvplot(restart.km, data = prim_restart,conf.int = T,fun="cloglog",censor=F)

restart.split=survSplit(Surv(time,restart)~., data=prim_restart, cut=8, end="time", event="restart", start="time0", episode="time_period")

restart.model.new=coxph(Surv(time0,time,restart)~as.factor(age_cat)*as.factor(time_period) +
                 gender + ethnicity + deprivation, data=restart.split)
summary(restart.model.new)

restart.ph <- cox.zph(restart.model.new)
restart.ph

restart.cox.df <- tidy(restart.model.new, conf.int=T, exponentiate=T)
restart.cox.df$estimate <- round(restart.cox.df$estimate, 2)
restart.cox.df$conf.low <- round(restart.cox.df$conf.low, 2)
restart.cox.df$conf.high <- round(restart.cox.df$conf.high, 2)

# Merge into a single string
restart.cox.df$conf.interval <- paste(restart.cox.df$conf.low, restart.cox.df$conf.high, sep = "-")

restart.cox.df <- restart.cox.df[, c('term', 'estimate', 'conf.interval', 'p.value')]

#load into excel 
wb <- loadWorkbook(paste0(excel, "cox_assumptions.xlsx"))
names(wb)
writeData(wb, "restart_R", restart.cox.df, startCol=1, startRow=4, rowNames=T)
saveWorkbook(wb,paste0(excel, "cox_assumptions.xlsx"), overwrite=T)

#########################SECONDARY PREVENTION###################################
################################################################################

sec_disc <- read_parquet(paste0(datafiles, "sec_discont_fact.parquet"))

sec_disc <- sec_disc[ethnicity!="missing" & deprivation!="missing"]

sec_disc[, c("age_cat", "ethnicity", "gender", "deprivation") := lapply(.SD, as.factor), .SDcols = c("age_cat", "ethnicity", "gender", "deprivation")]
sec_disc[ethnicity=="missing", ethnicity := NA]
sec_disc[deprivation=="missing", deprivation := NA]

sec_disc[, age_cat := factor(age_cat, levels = c("25-39", "40-49", "50-59", "60-69", "70+"))]
sec_disc[, gender := factor(gender, levels = c("Male", "Female"))]
sec_disc[, ethnicity := factor(ethnicity, levels = c("White", "Black", "South Asian", "Mixed", "Other"))]
sec_disc[, deprivation := factor(deprivation, levels = c("1 - least deprived", "2", "3", "4", "5 - most deprived"))] 

library("survival")
#Adjusted Cox model
sec_model <- coxph(Surv(time, discontinuation) ~ age_cat + gender + ethnicity + deprivation,
                    data=sec_disc)
sec_model

sec.ph <- cox.zph(sec_model)
sec.ph


#schoenfelds 
sch.resid=cox.zph(sec_model, transform = 'identity')
sch.resid
plot(sch.resid)

##km
sec.disc.km <- survfit(Surv(time,discontinuation)~as.factor(gender),data=sec_disc)
ggsurvplot(sec.disc.km, data = sec_disc, conf.int = T,censor=F)

#loglog plot
ggsurvplot(sec.disc.km, data = sec_disc,conf.int = T,fun="cloglog",censor=F)

sec.disc.split=survSplit(Surv(time,discontinuation)~., data=sec_disc, cut=6, end="time", event="discontinuation",
                     start="time0", episode="time_period")
#
sec.disc.new=coxph(Surv(time0,time,discontinuation)~as.factor(age_cat)*as.factor(time_period) + gender*as.factor(time_period) +
                     as.factor(ethnicity)*as.factor(time_period) + deprivation,data=sec.disc.split)
summary(sec.disc.new)

sec.disc.ph <- cox.zph(sec.disc.new)
sec.disc.ph

sec.disc.df <- tidy(sec.disc.new, conf.int=T, exponentiate=T)

sec.disc.df$estimate <- round(sec.disc.df$estimate, 2)
sec.disc.df$conf.low <- round(sec.disc.df$conf.low, 2)
sec.disc.df$conf.high <- round(sec.disc.df$conf.high, 2)

# Merge into a single string
sec.disc.df$conf.interval <- paste(sec.disc.df$conf.low, sec.disc.df$conf.high, sep = "-")

sec.disc.df <- sec.disc.df[, c('term', 'estimate', 'conf.interval', 'p.value')]

#load into excel 
wb <- loadWorkbook(paste0(excel, "cox_assumptions.xlsx"))
names(wb)
writeData(wb, "disc_R", sec.disc.df, startCol=7, startRow=4, rowNames=T)
saveWorkbook(wb,paste0(excel, "cox_assumptions.xlsx"), overwrite=T)


####Re-initiations secondary
##open parquet
sec_restart <- read_parquet(paste0(datafiles, "sec_restart_fact.parquet"))

sec_restart <- sec_restart[ethnicity!="missing" & deprivation!="missing"]

#Adjusted Cox model
sec_restart_model <- coxph(Surv(time, restart) ~ age_cat + gender + ethnicity + deprivation,
                       data=sec_restart)
summary(sec_restart_model)

#schoenfelds residuals global
sec.restart.ph <- cox.zph(sec_restart_model)
sec.restart.ph

####Not violated

sch.resid=cox.zph(sec_restart_model, transform = 'identity')
sch.resid
plot(sch.resid)

#log log plots
restart.km <- survfit(Surv(time,restart)~as.factor(age_cat),data=sec_restart)
ggsurvplot(restart.km, data = sec_restart,conf.int = T,fun="cloglog",censor=F)

sink()