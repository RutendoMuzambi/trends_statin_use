#--------------------------------------------------------------------------------------
# Project: Statins and CVD prevention
# Author: Rutendo Muzambi
# Description: Baseline baseline_descriptives
#--------------------------------------------------------------------------------------

########################
rm(list=ls())
library(finalfit)
########################################################################################################
################################Study population for descriptive analyses###############################
########################################################################################################
studypop_dt <- read_parquet(paste0(datafiles, "studypop_all.parquet"))
library(finalfit)
library(openxlsx)
###Primary prevention prevalence cohort
prev_pop_dt <- read_parquet(paste0(datafiles, "prim_prev_all.parquet"))
prev_pop_dt <- prev_pop_dt[,c("patid")]
#create ever statin user
# prev_pop_dt[, ever_statin := ifelse(any(statin_user == 1), 1, 0), by = patid]
# prev_pop_dt[ever_statin==0, prim_ever_statin := "Non-user"]
# prev_pop_dt[ever_statin==1, prim_ever_statin := "Ever statin user"]
prev_pop_dt[, prim_prev := 1]
# prev_pop_dt <- prev_pop_dt[,c("statin_user", "ever_statin") := NULL]
prev_pop_dt <- prev_pop_dt[, head(.SD, 1L), by="patid"]
uniqueN(prev_pop_dt, by = "patid") #922,813
prev_pop_dt <- merge(studypop_dt, prev_pop_dt, by="patid", all.x=T)

###initiation cohort
prim_init_dt <- read_parquet(paste0(datafiles, "prim_init_rec_risk.parquet"))
prim_init_dt <- prim_init_dt[,c("patid", "risk_date1", "initiation")]
prim_init_dt[, prim_initiation := 1]
prim_init_dt <- prim_init_dt[, head(.SD, 1L), by="patid"]
uniqueN(prim_init_dt, by = "patid") #83,599
prim_init_dt <- merge(studypop_dt, prim_init_dt, by="patid", all.x=T)

###discontinuation cohort
prim_disc_dt <- read_parquet(paste0(datafiles, "prim_discont_fact.parquet"))
prim_disc_dt <- prim_disc_dt[,c("patid", "first_statin_dt", "discontinuation")]
prim_disc_dt[, prim_discontinuation := 1]
prim_disc_dt <- prim_disc_dt[, head(.SD, 1L), by="patid"]
uniqueN(prim_disc_dt, by = "patid") #71,383
prim_disc_dt <- merge(studypop_dt, prim_disc_dt, by="patid", all.x=T)

###re-initiation cohort
prim_restart_dt <- read_parquet(paste0(datafiles, "prim_restart_fact.parquet"))
prim_restart_dt <- prim_restart_dt[,c("patid", "restart")]
prim_restart_dt[, .N, by="restart"]
# restart      N
# 1:       1 115043
# 2:       0  27541
prim_restart_dt[, prim_restart := 1]
prim_restart_dt <- prim_restart_dt[, head(.SD, 1L), by="patid"]
uniqueN(prim_restart_dt, by = "patid") #71,383
prim_restart_dt <- merge(studypop_dt, prim_restart_dt, by="patid", all.x=T)

##append all datasets together
prim_pop_dt <- rbind(prev_pop_dt, prim_init_dt, prim_disc_dt, prim_restart_dt, fill=TRUE)

prim_pop_dt[prim_prev==1, prim_subpop := "Prevalence"]
prim_pop_dt[prim_initiation==1, prim_subpop := "Initiation"]
prim_pop_dt[prim_discontinuation==1, prim_subpop := "Discontinuation"]
prim_pop_dt[prim_restart==1, prim_subpop := "Re-initiation"]
prim_pop_dt[, .N, by="prim_subpop"]

prim_total_dt <- prim_pop_dt[prim_prev==1|prim_initiation==1|prim_discontinuation==1|prim_restart==1]
prim_total_dt[, total_pop := 1]
prim_total_dt <- prim_total_dt[, head(.SD, 1L), by="patid"]
studypop_dt <- rbind(prim_pop_dt, prim_total_dt, fill=TRUE)
studypop_dt[total_pop==1, prim_subpop := "Total Cohort"]
studypop_dt[, .N, by="prim_subpop"]

###remove datasets from environment
rm(list = c("prim_pop_dt", "prim_init_dt", "prim_disc_dt", "prim_restart_dt"))

###Save parquet
##write_parquet(studypop_dt, paste0(datafiles, "studypop_descr.parquet"))
######table 1 baseline characteristics
explanatory =c("age_startfup", "age_cat", "gender", "ethnicity", "deprivation", "bmi_cat", "smokingstatus",
               "mean_sbp", "mean_dbp", "hdl", "ldl", "tchl", "thdl", "hypertension_cat", "treatedhyp","t2dm_cat",  "af_cat", "ckd_cat",
               "ra_cat", "initiation", "discontinuation", "restart")
dependent= "prim_subpop"
studypop_dt %>%
  summary_factorlist(dependent, explanatory, p= F, cont="mean",
                     add_col_totals = TRUE,
                     include_col_totals_percent = FALSE, 
                     column=T) -> prim_descript_dt
prim_descript_dt <- setDT(prim_descript_dt)
setcolorder(prim_descript_dt, c("label", "levels", "Total Cohort", "Prevalence", "Initiation", "Discontinuation", "Re-initiation"))

#load into excel 
wb <- loadWorkbook(paste0(excel, "baseline_descriptives.xlsx"))
names(wb)
writeData(wb, "table1_R", prim_descript_dt, startCol=1, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "baseline_descriptives.xlsx"), overwrite=T)

##median age
explanatory =c("age_startfup")
dependent= "prim_subpop"
studypop_dt %>%
  summary_factorlist(dependent, explanatory, p= F, cont="median",
                     add_col_totals = TRUE,
                     include_col_totals_percent = FALSE, 
                     column=T) -> prim_descript_dt
prim_descript_dt <- setDT(prim_descript_dt)
setcolorder(prim_descript_dt, c("label", "levels", "Total Cohort", "Prevalence", "Initiation", "Discontinuation", "Re-initiation"))

#load into excel 
wb <- loadWorkbook(paste0(excel, "baseline_descriptives.xlsx"))
names(wb)
writeData(wb, "table1_R", prim_descript_dt, startCol=1, startRow=2, rowNames=T)
saveWorkbook(wb,paste0(excel, "baseline_descriptives.xlsx"), overwrite=T)

###remove datasets from environment
rm(list = c("prim_descript_dt", "prim_restart_dt", "prim_disc_dt"))

studypop_dt <- read_parquet(paste0(datafiles, "studypop_all.parquet"))

########Secondary prevention
###Secondary prevalence
sec_pop_dt <- read_parquet(paste0(datafiles, "sec_prev_pop.parquet"))
sec_pop_dt <- sec_pop_dt[,c("patid", "CVD")]
sec_pop_dt[, .N, by="CVD"]
sec_pop_dt[, sec_prev := 1]
sec_pop_dt <- sec_pop_dt[, head(.SD, 1L), by="patid"]
uniqueN(sec_pop_dt, by = "patid") #99,082
sec_pop_dt <- merge(studypop_dt, sec_pop_dt, by="patid", all.x=T)

###Secondary initiation
sec_init_dt <- read_parquet(paste0(datafiles, "sec_initiators.parquet"))
sec_init_dt <- sec_init_dt[,c("patid", "CVD", "initiation")]
sec_init_dt[, sec_initiation := 1]
sec_init_dt <- sec_init_dt[, head(.SD, 1L), by="patid"]
uniqueN(sec_init_dt, by = "patid") #37,291
sec_init_dt <- merge(studypop_dt, sec_init_dt, by="patid", all.x=T)

####Secondary discontinuation
sec_disc_dt <- read_parquet(paste0(datafiles, "sec_discont_fact.parquet"))
sec_disc_dt <- sec_disc_dt[,c("patid", "CVD", "discontinuation")]
sec_disc_dt[, sec_discontinuation := 1]
sec_disc_dt <- sec_disc_dt[, head(.SD, 1L), by="patid"]
uniqueN(sec_disc_dt, by = "patid") #30,339
sec_disc_dt<- merge(studypop_dt, sec_disc_dt, by="patid", all.x=T)

####Secondary re-initiation
sec_restart_dt <- read_parquet(paste0(datafiles, "sec_restart_fact.parquet"))
sec_restart_dt <- sec_restart_dt[,c("patid", "CVD", "restart")]
sec_restart_dt[, .N, by="restart"]
# restart     N
# 1:       1 26578
# 2:       0  2576
sec_restart_dt[, sec_restart := 1]
sec_restart_dt <- sec_restart_dt[, head(.SD, 1L), by="patid"]
uniqueN(sec_restart_dt, by = "patid") #30,339
sec_restart_dt<- merge(studypop_dt, sec_restart_dt, by="patid", all.x=T)

##append all datasets together
studypop_dt <- rbind(sec_pop_dt, sec_init_dt, sec_disc_dt, sec_restart_dt, fill=TRUE)

###Formatting
studypop_dt[sec_prev==1, sec_subpop := "Prevalence"]
studypop_dt[sec_initiation==1, sec_subpop := "Initiation"]
studypop_dt[sec_discontinuation==1, sec_subpop := "Discontinuation"]
studypop_dt[sec_restart==1, sec_subpop := "Re-initiation"]

sec_total_dt <- studypop_dt[sec_prev==1|sec_initiation==1|sec_discontinuation==1|sec_restart==1]
sec_total_dt[, total_pop := 1]
sec_total_dt <- sec_total_dt[, head(.SD, 1L), by="patid"]
studypop_dt <- rbind(studypop_dt, sec_total_dt, fill=TRUE)
studypop_dt[total_pop==1, sec_subpop := "Total Cohort"]

studypop_dt[, .N, by="sec_subpop"]
studypop_dt[, .N, by="CVD"]

studypop_dt[CVD=="Myocardial Infarction", CVD := "MI"]
studypop_dt[CVD=="Peripheral Arterial Disease", CVD := "PAD"]
studypop_dt[CVD=="Transient Ischeamic Attack", CVD := "TIA"]
studypop_dt[CVD=="Cerebrovascular Procedures", CVD := "cerebrovascular procedures"]

studypop_dt[, .N, by="CVD"]

###remove datasets from environment
#rm(list=ls())
explanatory =c("age_startfup", "age_cat", "gender", "ethnicity", "deprivation", "bmi_cat", "smokingstatus",
               "mean_sbp", "mean_dbp", "hdl", "ldl", "tchl", "thdl", "hypertension_cat", "treatedhyp","t2dm_cat",  "af_cat", "ckd_cat",
               "ra_cat", "CVD","initiation", "discontinuation", "restart")
dependent= "sec_subpop"
studypop_dt %>%
  summary_factorlist(dependent, explanatory, p= F, cont="mean",
                     add_col_totals = TRUE,
                     include_col_totals_percent = FALSE, 
                     column=T) -> sec_descript_dt

sec_descript_dt <- setDT(sec_descript_dt)
setcolorder(sec_descript_dt, c("label", "levels", "Total Cohort", "Prevalence", "Initiation", "Discontinuation", "Re-initiation"))

wb <- loadWorkbook(paste0(excel, "baseline_descriptives.xlsx"))
names(wb)
writeData(wb, "table1_R", sec_descript_dt, startCol=10, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "baseline_descriptives.xlsx"), overwrite=T)

##median age
explanatory =c("age_startfup")
dependent= "sec_subpop"
studypop_dt %>%
  summary_factorlist(dependent, explanatory, p= F, cont="median",
                     add_col_totals = TRUE,
                     include_col_totals_percent = FALSE, 
                     column=T) -> sec_prev_age_dt
sec_prev_age_dt <- setDT(sec_prev_age_dt)

setcolorder(sec_prev_age_dt, c("label", "levels", "Total Cohort", "Prevalence", "Initiation", "Discontinuation", "Re-initiation"))

#load into excel 
wb <- loadWorkbook(paste0(excel, "baseline_descriptives.xlsx"))
names(wb)
writeData(wb, "table1_R", sec_prev_age_dt, startCol=10, startRow=2, rowNames=T)
saveWorkbook(wb,paste0(excel, "baseline_descriptives.xlsx"), overwrite=T)
