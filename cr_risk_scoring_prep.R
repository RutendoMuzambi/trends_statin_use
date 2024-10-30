#######################################################################################-
#Project: Statins and CVD prevention
#Author: Rutendo Muzambi
#Description: study population for factors associated with risk scoring analysis
#######################################################################################
rm(list=ls())

##########################################################################################################
###################FACTORS ASSOCIATED WITH RISK SCORING##################################
##########################################################################################################
patients_dt <- read_parquet(paste0(datafiles, "clean_pt.prac.parquet"))
patients_dt <- patients_dt[,c("patid", "start_fup", "enddate", "dob", "age_startfup")]
uniqueN(patients_dt, by = "patid") #

patients_dt[, do40thbday := dob+40*365.25]
patients_dt[, do75thbday := dob+75*365.25]

#create new start date latest of 40th birthday or start of follow up
patients_dt[, newstart := pmax(start_fup,do40thbday, na.rm=T)]

#Remove those not in study by 40th birthday or starting study after 75th birthday
patients_elig <- patients_dt[enddate>newstart & newstart<do75thbday]
uniqueN(patients_elig, by = "patid") #

#merge with EXCLUSIONS for primary prevention
exclusions_all <- read_parquet(paste0(datafiles, "excl_init_all.parquet"))
uniqueN(exclusions_all, by = "patid") #

exclusions_all <- exclusions_all[,c("patid", "excl_flag", "excl_date1")]

patients_elig <- merge(patients_elig ,exclusions_all, by="patid", all.x=TRUE)
uniqueN(patients_elig , by = "patid") #

#######exclusion date  before new start of follow up
patients_elig [excl_date1 <= newstart & !is.na(excl_date1), exclude := 1]
patients_elig  <- patients_elig [is.na(exclude)]
uniqueN(patients_elig , by = "patid") #

#merge with outcome - risk scoring
cvdrisk_all <- read_parquet(paste0(datafiles, "cvdrisk_all.parquet"))
uniqueN(cvdrisk_all, by = "patid") #

cvdrisk_all <- cvdrisk_all[!is.na(obsdate) & !is.na(risk_score)]
uniqueN(cvdrisk_all, by = "patid") #1,812,916

#Merge in patients dataset to flag those with CVD risk score
cvdrisk_all <- merge(cvdrisk_all[, c("patid", "risk_score", "obsdate")], patients_dt[, c("patid", "dob", "age_startfup", "enddate")], by="patid", all.x=TRUE)

#risk record during study period
cvdrisk_all <- cvdrisk_all[obsdate>dob]
cvdrisk_all <- cvdrisk_all[obsdate<=enddate]
uniqueN(cvdrisk_all, by = "patid") #1,811,847

#Merge in patients dataset to flag those with CVD risk score
cvdrisk_dt <- merge(patients_elig,cvdrisk_all[, c("patid", "risk_score", "obsdate")], by="patid", all.x=TRUE)

cvdrisk_dt <- cvdrisk_dt[(obsdate>=newstart & !is.na(obsdate))|is.na(obsdate)]
cvdrisk_dt <- cvdrisk_dt[(obsdate<enddate & !is.na(obsdate))|is.na(obsdate)]
uniqueN(cvdrisk_dt, by = "patid") #2,355,037

cvdrisk_dt[is.na(obsdate), risk_scoring := 0]
cvdrisk_dt[!is.na(obsdate), risk_scoring := 1]

cvdrisk_dt <- cvdrisk_dt[,c("enddate", "dob") := NULL]

#Select first risk score
cvdrisk_dt[!is.na(obsdate), risk_date1 := lapply(.SD, min, na.rm=TRUE), by=patid, .SDcols="obsdate"]
cvdrisk_dt <- cvdrisk_dt[(obsdate==risk_date1 & !is.na(obsdate))|is.na(obsdate)]

#if multiple risk scores on same date keep highest score
cvdrisk_dt[!is.na(risk_date1), maxscore := lapply(.SD, max, na.rm=T), by = "patid", .SDcols = "risk_score"]

cvdrisk_dt <- cvdrisk_dt[maxscore==risk_score|is.na(risk_date1)]
cvdrisk_dt <- unique(cvdrisk_dt)
uniqueN(cvdrisk_dt, by = "patid") #2,355,037

studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))

analysis_risk <- merge(studypop_all[, c("patid", "dob", "enddate", "gender", "ethnicity", "deprivation")],cvdrisk_dt, by="patid")

analysis_risk[is.na(risk_date1), risk_scoring := 0]
uniqueN(analysis_risk, by = "patid") #

library(epiDisplay)
library(epitools)

#age from new start
analysis_risk[, age_newstart := round(as.numeric(difftime(newstart, dob, units = "days")) / 365.25)]


analysis_risk[, age_cat := cut(age_newstart, breaks = c(40, 50, 60,70, 75),
                               labels = c("40-49", "50-59", "60-69", "70-74"), include.lowest=T)]

analysis_risk[, .N, by=age_cat]

uniqueN(analysis_risk, by = "patid") #2,355,037

test <- analysis_risk[is.na(age_cat)]

###############################COX REGRESSION ANALYSIS#######################
#install.packages("survival")
library(survival)
library(survminer)

#install.packages("eha")
#library(eha)

#add 1 day to risk date so those with risk on same date of entry are included
analysis_risk[, risk_date1 := risk_date1 + days(1)]

#censor at enddate, cvd risk date, exclusion date and 75th birthday
analysis_risk[, timeout := pmin(enddate,risk_date1,excl_date1,do75thbday, na.rm=T)]
analysis_risk[, timein := newstart]

analysis_risk[, age := round(as.numeric(difftime(timeout, dob, units = "days") / 365.25))]

analysis_risk[, time := as.numeric(timeout-timein)/365.25]

#drop if exit before start of follow up
analysis_risk <- analysis_risk[timein<=timeout]

#
uniqueN(analysis_risk, by = "patid") #2,355,037

cox_risk_dt <- analysis_risk

#distribution of outcome and time

hist(cox_risk_dt$time)
hist(subset(cox_risk_dt, risk_scoring==1)$time)
hist(subset(cox_risk_dt, risk_scoring==0)$time)

cox_risk_dt[, .N, by="risk_scoring"]

##
km_fit <- survfit(Surv(time, risk_scoring) ~ 1,
                  data=cox_risk_dt,
                  type="kaplan-meier")
print(km_fit)

###format variables
cox_risk_dt[, .N, by=age_cat]
cox_risk_dt[, .N, by=deprivation]

cox_risk_dt[deprivation=="missing", deprivation := NA]
cox_risk_dt[ethnicity=="missing", ethnicity := NA]

cox_risk_dt[, age_cat := factor(age_cat, levels = c("40-49", "50-59", "60-69","70-74"))]
cox_risk_dt[, gender := factor(gender, levels = c("Male", "Female"))]
cox_risk_dt[, ethnicity := factor(ethnicity, levels = c("White", "Black", "South Asian", "Mixed", "Other"))]
cox_risk_dt[, deprivation := factor(deprivation, levels = c("1 - least deprived", "2", "3", "4", "5 - most deprived"))]

cox_risk_dt[, .N, by="deprivation"]


##save parquet
write_parquet(cox_risk_dt, paste0(datafiles, "risk_scoring_fact.parquet"))


# Fit the Cox proportional hazards model 
cox_model1 <- coxph(Surv(time, risk_scoring) ~ age_cat + gender + ethnicity + deprivation,
                    data=cox_risk_dt)
summary(cox_model1)
