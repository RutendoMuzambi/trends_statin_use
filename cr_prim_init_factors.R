#######################################################################################-
#Project: Statins and CVD prevention
#Author: Rutendo Muzambi
#Description: Study population for factors associated with statin initiation
#######################################################################################

#####PRIMARY PREVENTION Initiation
#Those with a no previous statin prescription, no CVD event and a newly recorded cvd risk score eligible for 
#statin initiation
#NICE Guidelines 'Do not use a risk assessment tool for people who are at high risk of CVD,
#including people with type 1 diabetes, an estimated glomerular filtration rate less than 60 ml/min/1.73 m2 and/or albuminuria (CKD),
#familial hypercholesterolaemia

rm(list=ls())

############################Aim 3: Monthly statin initiators ###################################
####Primary prevention
#Include those without statins 

##Defining statin user
statins_dt <- read_parquet(paste0(datafiles, "statins.parquet"))
uniqueN(statins_dt, by = "patid") #

statins_dt <- statins_dt[, c("patid", "statin_type", "issuedate")]

###########MERGE with patient dataset
patients_dt <- read_parquet(paste0(datafiles, "clean_pt.prac.parquet"))
patients_dt <- patients_dt[,c("patid", "start_fup", "enddate", "dob")]
uniqueN(patients_dt, by = "patid") #

statins_pts <- merge(patients_dt,statins_dt, by="patid")
uniqueN(statins_pts, by = "patid") #

#keep statins after dob
setorder(statins_pts, issuedate)
statins_pts <- statins_pts[issuedate>dob]

#keep first statin Rx
setorder(statins_pts, patid, issuedate)
statins_pts <- statins_pts[order(issuedate), head(.SD, 1L), by="patid"]

everstatin_pts <- statins_pts

#keep statins after start of follow up and before enddate
statins_pts <- statins_pts[issuedate>=start_fup]
uniqueN(statins_pts, by = "patid") #
#keep first statin Rx before enddate
statins_pts<- statins_pts[issuedate<=enddate]
uniqueN(statins_pts, by = "patid") #

statins_pts <- statins_pts[, c("patid", "statin_type", "issuedate")]
setnames(statins_pts, "issuedate", "first_statin_dt")

write_parquet(statins_pts, paste0(datafiles, "incident_statins.parquet"))

####################################################################################################
##Initiation among those with a recorded risk score
cvdrisk_all <- read_parquet(paste0(datafiles, "cvdrisk_all.parquet"))
uniqueN(cvdrisk_all, by = "patid") #5,000,000

cvdrisk_all <- cvdrisk_all[!is.na(risk_score)]#681,241

#Match eligibility with updating guidelines
cvdrisk_all[, year := year(obsdate)]
cvdrisk_all[risk_score>=20 & year<2014, eligible := 1]
cvdrisk_all[risk_score<20 & year<2014, eligible := 0]

cvdrisk_all[risk_score>=10 & year>=2014, eligible := 1]
cvdrisk_all[risk_score<10 & year>=2014, eligible := 0]

##keep risk scores above threshold 
cvdrisk_all <- cvdrisk_all[eligible==1]

uniqueN(cvdrisk_all, by = "patid") #681,241

#keep only risk date on and after start of follow up
cvdrisk_all  <- cvdrisk_all[obsdate>=start_fup]
uniqueN(cvdrisk_all, by = "patid") #613,029

#first score date above threshold
setorder(cvdrisk_all, "patid", "obsdate")
cvdrisk_all[, risk_date1 := lapply(.SD, min, na.rm=TRUE), by=patid, .SDcols="obsdate"]

#keep first recorded risk score above threshold
cvdrisk_all <- cvdrisk_all[obsdate==risk_date1]
cvdrisk_all[, c("start_fup","enddate") := NULL]
uniqueN(cvdrisk_all, by = "patid") #613,029
###merge in statin information and patient information
prim_init_elig <- merge(patients_dt,statins_pts, by="patid", all.x=TRUE)

prim_init_elig <- merge(prim_init_elig, cvdrisk_all, by="patid")
uniqueN(prim_init_elig, by = "patid") ##613,029

##drop if age at risk date is not between NHS health check ages (40-74) 
prim_init_elig[, ageatrisk := (unclass(risk_date1) - unclass(dob)) / 365.25]
prim_init_elig[, ageatrisk := round(ageatrisk, digits=0)]
prim_init_elig <- prim_init_elig[ageatrisk>=40 & ageatrisk<75]

uniqueN(prim_init_elig, by = "patid") #537,219

####without exclusion codes for risk assessment before cvd risk score date
##Read and merge exclusions for secondary prevention population with data table with statins
#defining initiation with CVD risk score thus excluding those ineligible for risk scoring
exclusions_all <- read_parquet(paste0(datafiles, "excl_init_all.parquet"))
uniqueN(exclusions_all, by = "patid") #756,873

exclusions_all <- exclusions_all[,c("patid", "excl_flag", "excl_date1")]

prim_init_elig<- merge(prim_init_elig,exclusions_all, by="patid", all.x=TRUE)
uniqueN(prim_init_elig, by = "patid") #537,219

#######primary prevention eligibility before statins and start of follow 
#prim_init_dt[excl_date1 <= first_statin_dt & !is.na(first_statin_dt) & !is.na(excl_date1), exclude := 1]
#prim_init_elig[excl_date1 <= start_fup & !is.na(excl_date1), exclude := 2]
prim_init_elig[excl_date1 <= risk_date1 & !is.na(excl_date1), exclude := 1]
prim_init_elig <- prim_init_elig[is.na(exclude)]
uniqueN(prim_init_elig, by = "patid") #466,420

#28 day window following risk record
prim_init_elig[, grace_period := risk_date1+28]

#keep those who remain alive during grace period
prim_init_elig <- prim_init_elig[enddate>grace_period]

uniqueN(prim_init_elig, by = "patid") #462,711

#drop those on statins before risk date
prim_init_elig <- prim_init_elig[first_statin_dt>=risk_date1|is.na(first_statin_dt)]

uniqueN(prim_init_elig, by = "patid") #432,686

#remove unnecessary columns
prim_init_elig[, c("excl_flag","statin_type", "exclude", "eligible", "enddate", "dob", "start_fup", "obsdate") := NULL]

#eligible for initiation if statin Rx'ed within grace period
prim_init_elig[!is.na(first_statin_dt) & first_statin_dt<=grace_period, initiation := "yes"]
prim_init_elig[is.na(first_statin_dt)|(first_statin_dt>grace_period & !is.na(first_statin_dt)), initiation := "no"]
prim_init_elig[, .N, by=initiation]

######Merge with study population
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))
analysis_prim <- merge(prim_init_elig,studypop_all, by="patid")

#drop those with missing variables
#analysis_prim <- analysis_prim[deprivation!="missing" & ethnicity!="missing" & smokingstatus!="missing" & bmi_cat!="missing"]
#uniqueN(analysis_prim, by = "patid") #70,761


library(epiDisplay)
library(epitools)

# age at risk 
any(analysis_prim[, ageatrisk < 40]) #FALSE
analysis_prim[, age_cat := cut(ageatrisk, breaks = c(40, 50, 60, 70, Inf),
                              labels = c("40-49", "50-59", "60-69", "70+"), include.lowest=T)]


count(analysis_prim, 'age_cat') 

analysis_prim[, initiation := as.factor(initiation)]
analysis_prim[, age_cat := as.factor(age_cat)]
analysis_prim[, ethnicity := as.factor(ethnicity)]
analysis_prim[, deprivation := as.factor(deprivation)]
analysis_prim[, gender := as.factor(gender)]
analysis_prim[, bmi_cat := as.factor(bmi_cat)]
analysis_prim[, smokingstatus := as.factor(smokingstatus)]
analysis_prim[, t2dm_cat := as.factor(t2dm_cat)]
analysis_prim[, af_cat := as.factor(af_cat)]
analysis_prim[, hypertension_cat := as.factor(hypertension_cat)]
analysis_prim[, ra_cat := as.factor(ra_cat)]

str(analysis_prim)

#save columns ending with _date as an object
#date_cols <- grep("_date$", names(analysis_prim), value = TRUE)

# Select all columns that do not end in _date
#analysis_prim <- analysis_prim[, !date_cols, with = FALSE]

uniqueN(analysis_prim, by = "patid") #432,686

xtabs(~ initiation + gender, data=analysis_prim)
xtabs(~ initiation + age_cat, data=analysis_prim)
xtabs(~ initiation + bmi_cat, data=analysis_prim)
xtabs(~ initiation + deprivation, data=analysis_prim)
xtabs(~ initiation + ethnicity, data=analysis_prim)

analysis_prim[, age_cat := relevel(factor(age_cat), ref="40-49")]
analysis_prim[, bmi_cat := relevel(factor(bmi_cat), ref="Normal weight (18.5-24.9)")]
analysis_prim[, deprivation := relevel(factor(deprivation), ref="1 - least deprived")]
analysis_prim[, ethnicity := relevel(factor(ethnicity), ref="White")]
analysis_prim[, gender := relevel(factor(gender), ref="Male")]
analysis_prim[, smokingstatus := relevel(factor(smokingstatus), ref="non-smoker")]
analysis_prim[, hypertension_cat := relevel(factor(hypertension_cat), ref="0")]
analysis_prim[, ra_cat := relevel(factor(ra_cat), ref="0")]
analysis_prim[, t2dm_cat := relevel(factor(t2dm_cat), ref="0")]
analysis_prim[, treatedhyp := relevel(factor(treatedhyp), ref="0")]
analysis_prim[, af_cat := relevel(factor(af_cat), ref="0")]
analysis_prim[, ckd_cat := relevel(factor(ckd_cat), ref="0")]

analysis_prim[ethnicity=="missing", ethnicity := NA]
analysis_prim[smokingstatus=="missing", smokingstatus := NA]
analysis_prim[bmi_cat=="missing", bmi_cat := NA]
analysis_prim[deprivation=="missing", deprivation := NA]

test <- analysis_prim[!is.na(dm_date)]

test2 <- test[,c("patid", "start_fup", "risk_date1", "dm_date", "initiation", "first_statin_dt", "excl_date1")]
##save parquet
write_parquet(analysis_prim, paste0(datafiles, "prim_init_rec_risk.parquet"))
