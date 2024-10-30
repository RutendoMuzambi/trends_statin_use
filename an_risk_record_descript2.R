rm(list=ls())

######################################################################################

rec_risk_dt <- read_parquet(paste0(datafiles, "risk_scoring_fact.parquet"))

studypop_dt <- read_parquet(paste0(datafiles, "studypop_all.parquet"))

studypop_dt <- studypop_dt[, c("age_startfup", "age_cat", "gender",
                               "ethnicity", "deprivation") := NULL]

rec_risk_dt <- merge(rec_risk_dt, studypop_dt, by="patid", all.x=T)

rec_risk_dt[risk_scoring==1, risk_record := "Risk Record"]
rec_risk_dt[risk_scoring==0, risk_record := "No Risk Record"]

str(rec_risk_dt)
###BASELINE CHARACTERISTICS TABLE of those with and without recorded risk record

library(finalfit)
explanatory =c("age_startfup", "age_cat", "gender", "ethnicity", "deprivation", "bmi_cat", "smokingstatus",
               "mean_sbp", "mean_dbp", "hdl", "ldl", "tchl", "thdl", "hypertension_cat", "treatedhyp","t2dm_cat",  "af_cat", "ckd_cat",
               "ra_cat")
dependent= "risk_record"
rec_risk_dt %>%
  summary_factorlist(dependent, explanatory, p= F, cont="mean",
                     add_col_totals = TRUE,
                     include_col_totals_percent = FALSE, 
                     column=T, total_col = TRUE) -> cvdrisk_rec_dt
cvdrisk_rec_dt <- setDT(cvdrisk_rec_dt)
setcolorder(cvdrisk_rec_dt, c("label", "levels", "Total", "Risk Record", "No Risk Record"))

#load into excel 
library(openxlsx)
wb <- loadWorkbook(paste0(excel, "baseline_descriptives.xlsx"))
names(wb)
writeData(wb, "rec_risk2_R", cvdrisk_rec_dt, startCol=1, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "baseline_descriptives.xlsx"), overwrite=T)

##median age
explanatory =c("age_startfup")
dependent= "risk_record"
rec_risk_dt %>%
  summary_factorlist(dependent, explanatory, p= F, cont="median",
                     add_col_totals = TRUE,
                     include_col_totals_percent = FALSE, 
                     column=T, total_col = TRUE) -> cvdrisk_rec_dt
cvdrisk_rec_dt <- setDT(cvdrisk_rec_dt)
setcolorder(cvdrisk_rec_dt, c("label", "levels", "Total", "Risk Record", "No Risk Record"))

#load into excel 
wb <- loadWorkbook(paste0(excel, "baseline_descriptives.xlsx"))
names(wb)
writeData(wb, "rec_risk2_R", cvdrisk_rec_dt, startCol=1, startRow=2, rowNames=T)
saveWorkbook(wb,paste0(excel, "baseline_descriptives.xlsx"), overwrite=T)


#############TABLE 1

library(gtsummary)

explanatory =c("age_startfup", "age_cat", "gender", "ethnicity", "deprivation", "bmi_cat", "smokingstatus",
               "mean_sbp", "mean_dbp", "hdl", "ldl", "tchl", "thdl", "hypertension_cat", "treatedhyp","t2dm_cat",  "af_cat", "ckd_cat",
               "ra_cat")

# Overall
missing_table <- studypop_dt %>% 
  select(ethnicity, deprivation, bmi_cat, smokingstatus,
         mean_sbp_cat, mean_dbp_cat, hdl_cat, ldl_cat,
         tchl_cat, thdl_cat)  %>%
  tbl_summary(
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous()  ~ c(2, 2),
                  all_categorical() ~ c(0, 1)),
    type = list(ethnicity   ~ "categorical",
                deprivation    ~ "categorical",
                bmi_cat   ~ "categorical",
                smokingstatus      ~ "categorical",
                mean_sbp_cat   ~ "categorical",
                mean_dbp_cat    ~ "categorical",
                hdl_cat   ~ "categorical",
                ldl_cat  ~ "categorical",
                tchl_cat  ~ "categorical",
                thdl_cat ~ "categorical")
  ) %>%
  modify_header(label = "**Variable**") %>%
  modify_caption("Participant characteristics") %>%
  bold_labels()

wb <- loadWorkbook(paste0(excel, "baseline_descriptives.xlsx"))
names(wb)
writeData(wb, "missing_data_R", missing_table, startCol=2, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "baseline_descriptives.xlsx"), overwrite=T)


