#--------------------------------------------------------------------------------------
# Project: Statins and CVD prevention
# Author: Rutendo Muzambi
# Description: study populations  
#--------------------------------------------------------------------------------------

#####Study populations#######
rm(list=ls())  ##clear environment

######################################################################################################
#Merge all covariates into final study population
patients_dt <- read_parquet(paste0(datafiles, "clean_pt.prac.parquet"))

patients_dt <- patients_dt[,c("patid", "start_fup", "enddate", "dob", "gender", "age_startfup")]

#age groups
patients_dt[, age_cat := cut(age_startfup, breaks = c(25, 40, 50, 60, 70, 80, Inf),
                                  labels = c("25-39", "40-49", "50-59", "60-69", "70-79", "80+"), include.lowest=T)]

patients_dt[, .N, by=age_cat]

#read ethnicity
ethnicity_dt <- read_parquet(paste0(datafiles, "ethnicity_all.parquet"))
ethnicity_dt <-ethnicity_dt[,c("patid", "eth5")]
ethnicity_dt[eth5=="Not stated", eth5 := "missing"]
setnames(ethnicity_dt, "eth5", "ethnicity")
ethnicity_dt[, .N, by="ethnicity"]

#read townsend
townsend <- read_parquet(paste0(datafiles, "patient_townsend.parquet"))

townsend[, .N, by="e2011_townsend_20"] #16,607 missing town

townsend[e2011_townsend_20 <=4, townsend_cat := 1]
townsend[e2011_townsend_20>4 & e2011_townsend_20<=8, townsend_cat := 2]
townsend[e2011_townsend_20>8 & e2011_townsend_20<=12, townsend_cat := 3]
townsend[e2011_townsend_20>12 & e2011_townsend_20<=16, townsend_cat := 4]
townsend[e2011_townsend_20>16 & e2011_townsend_20<=20, townsend_cat := 5]

townsend[, .N, by="townsend_cat"]

townsend[townsend_cat==1, deprivation := "1 - least deprived"]
townsend[townsend_cat==2, deprivation := 2]
townsend[townsend_cat==3, deprivation := 3]
townsend[townsend_cat==4, deprivation := 4]
townsend[townsend_cat==5, deprivation := "5 - most deprived"]
townsend[is.na(townsend_cat), deprivation := "missing"]

townsend[, .N, by="deprivation"]

#read bmi, smoking
bmi <- read_parquet(paste0(datafiles, "bmi.parquet"))

smoke_dt <- read_parquet(paste0(datafiles, "smoking.parquet"))
smoke_dt[smokstatus==0, smokingstatus := "non-smoker"]
smoke_dt[smokstatus==1, smokingstatus := "current"]
smoke_dt[smokstatus==2, smokingstatus := "ex-smoking"]
smoke_dt[smokstatus==9, smokingstatus := "missing"]
smoke_dt[smokstatus==12, smokingstatus := "current"]

#Merge variables
studypop_all <- merge(patients_dt,ethnicity_dt, by="patid", all.x=T)
studypop_all <- merge(studypop_all,townsend[,c("patid", "townsend_cat", "deprivation")], by="patid", all.x=T)
studypop_all <- merge(studypop_all,bmi[,c("patid", "bmi", "bmi_cat")], by="patid", all.x=T)
studypop_all <- merge(studypop_all,smoke_dt[,c("patid", "smokstatus", "smokingstatus")], by="patid", all.x=T)

#Add covariates
af_dt <- read_parquet(paste0(datafiles, "af_all.parquet"))
bp_dt <- read_parquet(paste0(datafiles, "blood_pressure.parquet"))
cholesterol_dt <- read_parquet(paste0(datafiles, "cholesterol.parquet"))
hypertension_dt <- read_parquet(paste0(datafiles, "hypertension_all.parquet"))
ra_dt <- read_parquet(paste0(datafiles, "rheumatoid_arthritis_all.parquet"))
t2dm_dt <- read_parquet(paste0(datafiles, "t2dm_all.parquet"))
ckd_dt <- read_parquet(paste0(datafiles, "ckd_all.parquet"))
treatedhyp_dt <- read_parquet(paste0(datafiles, "treatedhyp.parquet"))
setnames(treatedhyp_dt, "htn_date", "treat_htn_date")
antihyperten_dt <- read_parquet(paste0(datafiles, "antihypertensives.parquet"))
cvdrisk_all <- read_parquet(paste0(datafiles, "cvdrisk_all.parquet"))
setnames(cvdrisk_all, "obsdate", "risk_date")
setorder(cvdrisk_all,"patid", "risk_date")
cvdrisk_all <- cvdrisk_all[, head(.SD, 1L), by="patid"]

#Merge remaining covariates 
studypop_all <- merge(studypop_all,af_dt[,c("patid", "af_cat", "af_date")], by="patid", all.x=T)
studypop_all <- merge(studypop_all,ckd_dt[,c("patid", "ckd_cat", "ckd_date")], by="patid", all.x=T)
studypop_all <- merge(studypop_all,bp_dt, by="patid", all.x=T)
studypop_all <- merge(studypop_all,cholesterol_dt, by="patid", all.x=T)
studypop_all <- merge(studypop_all,hypertension_dt[,c("patid", "hypertension_cat", "htn_date")], by="patid", all.x=T)
studypop_all <- merge(studypop_all,ra_dt[,c("patid", "ra_cat", "ra_date")], by="patid", all.x=T)
studypop_all <- merge(studypop_all,t2dm_dt[,c("patid", "t2dm_cat", "dm_date")], by="patid", all.x=T)
studypop_all <- merge(studypop_all,treatedhyp_dt[,c("patid", "treatedhyp", "treat_htn_date")], by="patid", all.x=T)
studypop_all <- merge(studypop_all,antihyperten_dt[,c("patid", "antihypertensives", "antihyp_date")], by="patid", all.x=T)
studypop_all <- merge(studypop_all,cvdrisk_all[,c("patid", "risk_cat", "risk_date")], by="patid", all.x=T)


#code missings categories
cat_cols <- c("af_cat", "ckd_cat", "hypertension_cat", "ra_cat", "t2dm_cat", "treatedhyp")
studypop_all[, (cat_cols) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = cat_cols]

studypop_all[is.na(deprivation), deprivation := "missing"]
studypop_all[is.na(bmi), bmi_cat := "missing"]
studypop_all[is.na(smokingstatus), smokingstatus := "missing"]
studypop_all[is.na(ethnicity), ethnicity := "missing"]

studypop_all <- studypop_all[, head(.SD, 1L), by="patid"]

###remove from environment
rm(list = c("af_dt", "bmi", "bp_dt", "cholesterol_dt", "ethnicity_dt", "hypertension_dt",
            "ckd_dt", "ra_dt", "smoke_dt", "patients_dt", "t2dm_dt", "townsend", "treatedhyp_dt", "antihyperten_dt"))

uniqueN(studypop_all, by = "patid") #5000000

###Save parquet
write_parquet(studypop_all, paste0(datafiles, "studypop_all.parquet"))


studypop_dt <- read_parquet(paste0(datafiles, "studypop_all.parquet"))

#generate categorical variables for blood pressure and cholesterol
cols_miss <- c("mean_sbp", "mean_dbp", "hdl", "ldl", "tchl", "thdl")
# 
 for (col in cols_miss) {
   new_col_name <- paste0(col, "_cat")
   studypop_dt[, (new_col_name) := ifelse(!is.na(get(col)), "yes", "no")]
 }
 
studypop_dt[, ethnicity_cat := ifelse(ethnicity!="missing", "yes", "no")]
studypop_dt[, depri_cat := ifelse(deprivation!="missing", "yes", "no")]
studypop_dt[, smoking_cat := ifelse(smokingstatus!="missing","yes", "no")]
studypop_dt[, bmi_group := ifelse(bmi_cat!="missing", "yes", "no")]


library(gtsummary)

# Overall
missing_table <- studypop_dt %>% 
  dplyr::select(ethnicity_cat, depri_cat, bmi_group, smoking_cat,
          mean_sbp_cat, mean_dbp_cat, hdl_cat, ldl_cat,
         tchl_cat, thdl_cat)  %>%
  tbl_summary(
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous()  ~ c(2, 2),
                  all_categorical() ~ c(0, 1)),
    type = list(ethnicity_cat   = "categorical",
                depri_cat    = "categorical",
                bmi_group   = "categorical",
                smoking_cat      = "categorical",
                   mean_sbp_cat   = "categorical",
                mean_dbp_cat    = "categorical",
                hdl_cat   = "categorical",
                ldl_cat  = "categorical",
                tchl_cat  = "categorical",
                thdl_cat = "categorical")
  ) %>%
  modify_header(label = "**Variable**") %>%
  modify_caption("Participant characteristics") %>%
  bold_labels()

studypop_dt[, .N, by="mean_sbp_cat"]

studypop_dt[, .N, by="depri_cat"]
str(studypop_dt)

test <- studypop_dt[, c("patid", "mean_sbp", "mean_sbp_cat")]

library(openxlsx)

wb <- loadWorkbook(paste0(excel, "baseline_descriptives.xlsx"))
names(wb)
writeData(wb, "missing_data_R", missing_table, startCol=2, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "baseline_descriptives.xlsx"), overwrite=T)


