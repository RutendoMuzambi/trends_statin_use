#######################################################################################-
#Project: Statins and CVD prevention
#Author: Rutendo Muzambi
#Description: factors associated with statin initiation secondary prevention
#######################################################################################

rm(list=ls())

###############################################################################################
statins_pts <- read_parquet(paste0(datafiles, "statins.parquet"))
cvd_incidence_dt <- read_parquet(paste0(datafiles, "cvd_incident_all.parquet"))

cvd_incidence_dt[, c("dob", "enddate") := NULL]


patients_dt <- read_parquet(paste0(datafiles, "clean_pt.prac.parquet"))
cvd_incidence_dt <- merge(patients_dt[,c("patid", "start_fup", "enddate")], cvd_incidence_dt, by = "patid")
uniqueN(cvd_incidence_dt, by = "patid") #433,152

#format CVD subtype
cvd_incidence_dt[, .N, by=cvd_event] #no missing cvd events
cvd_incidence_dt[mi==1, CVD := "MI"]
cvd_incidence_dt[stroke==1, CVD := "Stroke"]
cvd_incidence_dt[angina==1, CVD := "Angina"]
cvd_incidence_dt[tia==1, CVD := "TIA"]
cvd_incidence_dt[pad==1, CVD := "PAD"]
cvd_incidence_dt[acs==1 & is.na(mi) & is.na(angina), CVD := "Acute Coronary Syndrome (non-specific)"]
cvd_incidence_dt[chd==1 & is.na(mi) & is.na(angina) & is.na(acs), CVD := "Coronary Heart Disease (non-specific)"]
cvd_incidence_dt[cerebro_procs==1, CVD := "cerebrovascular procedures"]
cvd_incidence_dt[, .N, by=CVD]

#secondary prevention eligibiltiy
cvd_incidence_dt[, .N, by=cvd_event] #no missing cvd events
cvd_incidence_dt <- cvd_incidence_dt[, c("patid", "start_fup", "enddate", "cvd_event", "CVD", "incident_date1")]

#keep only CVD diagnosis on and after start of follow up
sec_elig_dt <- cvd_incidence_dt[incident_date1>=start_fup]
uniqueN(sec_elig_dt, by = "patid") #157,267
#keep cvd diagnosis before enddate
sec_elig_dt <- sec_elig_dt[incident_date1<enddate]
uniqueN(sec_elig_dt, by = "patid") #156,100

########If statins issued within 60 days of first cvd event then statins issued under secondary prevention
sec_elig_dt[, grace_period := incident_date1+60] #allowing for a grace period of 60 days

#Merge with statins
sec_elig_dt <- merge(sec_elig_dt,statins_pts[,c("patid", "issuedate")], by="patid", all.x=TRUE)
#keep statin issue date on the day or after cvd event date - this drops people on primary prevention before CVD
#sec_elig_dt <- sec_elig_dt[is.na(issuedate)|issuedate >=incident_date1]
#uniqueN(sec_elig_dt, by = "patid") #167,120

#label first earliest statin date
setorder(sec_elig_dt, "patid", "issuedate")
sec_elig_dt[!is.na(issuedate), first_statin_dt := lapply(.SD, min, na.rm=T), by = "patid", .SDcols = "issuedate"]

uniqueN(sec_elig_dt, by = "patid") #156,100

#keep those who remain alive during grace period
sec_elig_dt <- sec_elig_dt[enddate>grace_period]

uniqueN(sec_elig_dt, by = "patid") #146,811

sec_elig_dt[!is.na(first_statin_dt) & first_statin_dt<=grace_period, initiation := "yes"]
sec_elig_dt[is.na(first_statin_dt)|(first_statin_dt>grace_period & !is.na(first_statin_dt)), initiation := "no"]

#keep first statin row
sec_elig_dt <- sec_elig_dt[is.na(first_statin_dt)|first_statin_dt==issuedate]
sec_elig_dt <- sec_elig_dt[order(first_statin_dt), head(.SD, 1L), by="patid"]

uniqueN(sec_elig_dt, by = "patid") #146,811

sec_elig_dt[, .N, by=initiation]
# initiation      N
# 1:        yes 115968
# 2:         no  30843

sec_elig_dt[, year := year(incident_date1)]
sec_elig_dt[, c("start_fup", "enddate") := NULL]


#merge with covariates
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))
analysis_sec <- merge(sec_elig_dt,studypop_all, by="patid")

#drop missing levels
#analysis_sec <- analysis_sec[deprivation!="missing" & ethnicity!="missing" & smokingstatus!="missing" & bmi_cat!="missing"]
#uniqueN(analysis_sec, by = "patid") 

cols_to_delete <- grepl("date$", names(analysis_sec), ignore.case = TRUE)
analysis_sec[, !cols_to_delete, with = FALSE]

library(epiDisplay)
library(epitools)

# age at start of follow up
analysis_sec[, age_cat := cut(age_startfup, breaks = c(25, 40, 50, 60, 70, Inf),
                               labels = c("25-39", "40-49", "50-59", "60-69", "70+"), include.lowest=T)]

analysis_sec[, initiation := as.factor(initiation)]
analysis_sec[, age_cat := as.factor(age_cat)]
analysis_sec[, ethnicity := as.factor(ethnicity)]
analysis_sec[, deprivation := as.factor(deprivation)]
analysis_sec[, gender := as.factor(gender)]
analysis_sec[, bmi_cat := as.factor(bmi_cat)]
analysis_sec[, smokingstatus := as.factor(smokingstatus)]
analysis_sec[, af_cat := as.factor(af_cat)]
analysis_sec[, ckd_cat := as.factor(ckd_cat)]
analysis_sec[, hypertension_cat := as.factor(hypertension_cat)]
analysis_sec[, ra_cat := as.factor(ra_cat)]
analysis_sec[, t2dm_cat := as.factor(t2dm_cat)]

str(analysis_sec)

#save columns ending with _date as an object
#date_cols <- grep("_date$", names(analysis_sec), value = TRUE)

# Select all columns that do not end in _date
#analysis_sec <- analysis_sec[, !date_cols, with = FALSE]

#nrow(analysis_sec[is.na(analysis_sec$initiation)|is.na(analysis_sec$age_cat),
#                   is.na(analysis_sec$ethnicity)|is.na(analysis_sec$deprivation)])

#analysis_sec <- analysis_sec[ethnicity!="missing" & deprivation!="missing" & bmi_cat!="missing"]

uniqueN(analysis_sec, by = "patid") #146,811

xtabs(~ initiation + gender, data=analysis_sec)
xtabs(~ initiation + age_cat, data=analysis_sec)
xtabs(~ initiation + bmi_cat, data=analysis_sec)
xtabs(~ initiation + deprivation, data=analysis_sec)
xtabs(~ initiation + ethnicity, data=analysis_sec)

#logistic <- glm(initiation ~ gender1, data = analysis_sec, family = "binomial")
#summary(logistic)

analysis_sec[, bmi_cat := relevel(factor(bmi_cat), ref="Normal weight (18.5-24.9)")]
analysis_sec[, deprivation := relevel(factor(deprivation), ref="1 - least deprived")]
analysis_sec[, ethnicity := relevel(factor(ethnicity), ref="White")]
analysis_sec[, gender := relevel(factor(gender), ref="Male")]
analysis_sec[, smokingstatus := relevel(factor(smokingstatus), ref="non-smoker")]


write_parquet(analysis_sec, paste0(datafiles, "sec_initiators.parquet"))


pre_pandemic_cvd <- analysis_sec[incident_date1<"2020-03-01"]

#install.packages("broom")
library(broom)

pre_lg_cvd <- glm(initiation ~ age_cat + ethnicity + gender + deprivation + bmi_cat + smokingstatus +
                    hypertension_cat + af_cat + ckd_cat + ra_cat + t2dm_cat, data = pre_pandemic_cvd, family = binomial())
logistic.display(pre_lg_cvd)













