#######################################################################################-
#Project: Statins and CVD prevention
#Author: Rutendo Muzambi
#Description: Number of CVD events prevented 
#######################################################################################

#####Study populations#######
##remove environment
rm(list=ls())

############################################################################################
#####################################CENSUS DATA#######################################
#############################################################################################
#Import england census data to calculate total of England population aged 25 years and older
library(readxl)
eng_census_dt <- read_excel(paste0(datafiles, "health_files\\eng_census_2021.xlsx"))
sum(eng_census_dt$Observation) #56,490,045
eng_census_dt <- setDT(eng_census_dt)
setnames(eng_census_dt, "Age (101 categories) Code", "age")

eng_pop_25 <- eng_census_dt[age>=25, sum(Observation)]
print(eng_pop_25) #40,017,720

#Include only those aged 40 years and older
eng_census_dt <- eng_census_dt[age>=25]

############################Statin prevalence###################################
####Primary prevention
#Include those without statins 

##Defining statin user
statins_dt <- read_parquet(paste0(datafiles, "statins.parquet"))
uniqueN(statins_dt, by = "patid") #197,996

statins_dt <- statins_dt[, c("patid", "statin_type", "issuedate", "duration")]

###########MERGE with patient dataset
patients_dt <- read_parquet(paste0(datafiles, "clean_pt.prac.parquet"))
patients_dt <- patients_dt[,c("patid", "start_fup", "enddate", "dob")]
uniqueN(patients_dt, by = "patid") #1,000,000

statins_pts <- merge(patients_dt,statins_dt, by="patid")
uniqueN(statins_pts, by = "patid") #197,996

#keep statins after dob, start_fup and before enddate
setorder(statins_pts, issuedate)
statins_pts <- statins_pts[issuedate>dob & (issuedate>=start_fup & issuedate<=enddate)]

#run out date and grace period
statins_pts[,runoutdate := issuedate+duration]
statins_pts[, grace_period := runoutdate+28]

#label first statin Rx
# setorder(statins_pts, patid, issuedate)
# statins_pts[, first_statin := lapply(.SD, min, na.rm=T), by = "patid", .SDcols = "issuedate"]

#keep first statin Rx
setorder(statins_pts, patid, -issuedate)
statins_pts <- statins_pts[order(-issuedate), head(.SD, 1L), by="patid"]


####Primary prevention
#Include those without statins 

#######population with eligible RECORDED CVD RISK SCORES 
cvdrisk_all <- read_parquet(paste0(datafiles, "cvdrisk_all.parquet"))
uniqueN(cvdrisk_all, by = "patid") #5,000,000

cvdrisk_all <- cvdrisk_all[!is.na(risk_score)]
uniqueN(cvdrisk_all, by = "patid") #1,812,926

#Match eligibility with updating guidelines
cvdrisk_all[, year := year(obsdate)]
cvdrisk_all[risk_score>=20 & year<2014, eligible := 1]
cvdrisk_all[risk_score<20 & year<2014, eligible := 0]

cvdrisk_all[risk_score>=10 & year>=2014, eligible := 1]
cvdrisk_all[risk_score<10 & year>=2014, eligible := 0]

##keep risk scores above 10%
recordrisk_dt <- cvdrisk_all[risk_score>=10]
uniqueN(recordrisk_dt, by = "patid") #823,059

#keep only risk date on and after start of follow up and before enddate
recordrisk_dt  <- recordrisk_dt[obsdate>=start_fup]
recordrisk_dt  <- recordrisk_dt[obsdate<=enddate]
uniqueN(recordrisk_dt, by = "patid") #703,785

#first score date above threshold
recordrisk_dt[, risk_date1 := lapply(.SD, min, na.rm=TRUE), by=patid, .SDcols="obsdate"]

#keep first recorded risk score above threshold
setorder(recordrisk_dt, "risk_date1")
recordrisk_dt <- recordrisk_dt[, head(.SD, 1L), by="patid"]
recordrisk_dt[, c("start_fup","enddate") := NULL]
uniqueN(recordrisk_dt, by = "patid") #703,785

#Merge those with risk score below 10
cvdrisk_ineglib <- cvdrisk_all[risk_score<10]
setnames(cvdrisk_ineglib, "obsdate", "obsdate10")
setnames(cvdrisk_ineglib, "risk_cat", "risk_cat10")
setnames(cvdrisk_ineglib, "risk_score", "risk_score10")
recordrisk_dt <- merge(recordrisk_dt,cvdrisk_ineglib[, c("patid", "risk_cat10", "obsdate10", "risk_score10")], by="patid", all=TRUE)
recordrisk_dt <- recordrisk_dt[risk_score>=10|is.na(risk_cat), head(.SD, 1L), by="patid"]
uniqueN(recordrisk_dt, by = "patid") #1,725,933
recordrisk_dt[is.na(risk_date1), risk_date1 := obsdate10]
recordrisk_dt[is.na(risk_score), risk_score := risk_score10]
recordrisk_dt[is.na(risk_cat), risk_cat := "<10%"]

#merge with statins 
recordrisk_dt <- merge(recordrisk_dt,statins_pts[, c("patid", "issuedate")], by="patid", all.x=T)
recordrisk_dt[, .N, by="risk_cat"]

#######EXCLUSIIONS#######UPDATE TO EXCLUSIONS OF CVD, statin contraindications & risk score
exclusions_all <- read_parquet(paste0(datafiles, "exclusions_all.parquet"))
uniqueN(exclusions_all, by = "patid") #692,515
exclusions_all <- exclusions_all[,c("patid", "excl_date1")]

#merge CVD events and ineligibility for QRISK scoring
recordrisk_dt <- merge(recordrisk_dt,exclusions_all, by="patid", all.x=TRUE)

#######primary prevention eligibility before and on risk date
recordrisk_dt <- recordrisk_dt[(excl_date1 > risk_date1 & !is.na(excl_date1)) | is.na(excl_date1)]
uniqueN(recordrisk_dt, by = "patid") #1,598,352  # inc with prior statins

###calculate number of CVD events prevented
#probabilities so divide by 100
setnames(recordrisk_dt, "risk_score", "score")

recordrisk_dt[, risk_score := round(score/100, digits=2)]

recordrisk_dt[score<10, risk_cat := "<10%"]
recordrisk_dt[score>=10 & risk_score<20, risk_cat := "10-19%"]
recordrisk_dt[score>=20 & risk_score<30, risk_cat := "20-29%"]
recordrisk_dt[score>=30, risk_cat := "30%+"]

recordrisk_dt[, .N, by=risk_cat]

#total risk score
recordrisk_dt[, total_pop := .N]
recordrisk_dt[, total_risk := sum(risk_score, na.rm=T), by="risk_cat"]
recordrisk_dt[, risk_reduction := 0.25]

recordrisk_dt[, total_events_avoided := total_risk*risk_reduction]

##total population in each risk subgroup
recordrisk_dt[risk_cat=="<10%", subpop := .N]
recordrisk_dt[risk_cat=="10-19%", subpop := .N]
recordrisk_dt[risk_cat=="20-29%", subpop := .N]
recordrisk_dt[risk_cat=="30%+", subpop := .N]

#total statins
recordrisk_dt[!is.na(issuedate), statin_use := 1]
recordrisk_dt[, total_statin := sum(statin_use, na.rm=T), by="risk_cat"]
#number not on statins
recordrisk_dt[, no_statin_use := subpop-total_statin]
#statin prevalence
recordrisk_dt[, curr_prev := total_statin/subpop]

#total number of expected events based on risk score
recordrisk_dt[, sumrisk := sum(risk_score, na.rm=T)]
recordrisk_dt[, sumrisk0 := sum(risk_score[risk_cat=="<10%"], na.rm=T)]
recordrisk_dt[, sumrisk10 := sum(risk_score[risk_cat=="10-19%"], na.rm=T)]
recordrisk_dt[, sumrisk20 := sum(risk_score[risk_cat=="20-29%"], na.rm=T)]
recordrisk_dt[, sumrisk30 := sum(risk_score[risk_cat=="30%+"], na.rm=T)]

#number of events avoided if all eligible patients were treated (100% prevalence)
#recordrisk_dt[, events_avoided_all := sumrisk*risk_reduction]
recordrisk_dt[risk_cat=="<10%", events_avoided100 := sumrisk0*risk_reduction]
recordrisk_dt[risk_cat=="10-19%", events_avoided100 := sumrisk10*risk_reduction]
recordrisk_dt[risk_cat=="20-29%", events_avoided100 := sumrisk20*risk_reduction]
recordrisk_dt[risk_cat=="30%+", events_avoided100 := sumrisk30*risk_reduction]

#number of events avoided if patients were treated based on current prevalence
#recordrisk_dt[, events_avoided_all := sumrisk*risk_reduction]
recordrisk_dt[risk_cat=="<10%", events_avoided_curr := sumrisk0*risk_reduction*curr_prev]
recordrisk_dt[risk_cat=="10-19%", events_avoided_curr := sumrisk10*risk_reduction*curr_prev]
recordrisk_dt[risk_cat=="20-29%", events_avoided_curr := sumrisk20*risk_reduction*curr_prev]
recordrisk_dt[risk_cat=="30%+", events_avoided_curr := sumrisk30*risk_reduction*curr_prev]

# recordrisk_dt[, statin_prev_all:= round(statin_rx_all/total_pop*100, digits=1)]
# recordrisk_dt[, statin_prev:= round(statin_rx/subpop*100, digits=1)]
# 
# #number of events avoided based on current prevalence
# recordrisk_dt[, events_avoided_curr_all := events_avoided_all*(statin_prev_all/100)]
# recordrisk_dt[, events_avoided_curr := events_avoided*(statin_prev/100)]
# 
# #difference between events calculated from 100% prevalence to current prevalence
#recordrisk_dt[, diff_gain_all := events_avoided-events_avoided_curr]
recordrisk_dt[, diff_gain := events_avoided100-events_avoided_curr]

#scale to 2021 England population 25 and older
#recordrisk_dt[, scaled_diff_all := (eng_pop_25*diff_gain_all)/5000000]
recordrisk_dt[, scaled_diff := (eng_pop_25*diff_gain)/5000000]

#prepare for export to excel
rec_risk_table <- recordrisk_dt
rec_risk_table <- rec_risk_table[, c("risk_cat", "subpop", "total_risk", "no_statin_use", "curr_prev", "events_avoided100", "events_avoided_curr",
                   "diff_gain", "scaled_diff")]
rec_risk_table<- unique(rec_risk_table)
setorder(rec_risk_table, -"subpop")

library(openxlsx)

wb <- loadWorkbook(paste0(excel, "cvd_events_avoided.xlsx"))
names(wb)
writeData(wb, "cvd_events_R", rec_risk_table, startCol=2, startRow=2, rowNames=T)
saveWorkbook(wb,paste0(excel, "cvd_events_avoided.xlsx"), overwrite=T)

########################################################################################
###############################SECONDARY Prevention#####################################
########################################################################################
cvd_incidence_dt <- read_parquet(paste0(datafiles, "cvd_incident_all.parquet"))
uniqueN(cvd_incidence_dt, by = "patid") #433,152
#format CVD subtype
cvd_incidence_dt[, .N, by=cvd_event] #no missing cvd events
cvd_incidence_dt[mi==1, CVD := "MI"]
cvd_incidence_dt[stroke==1, CVD := "Stroke"]
cvd_incidence_dt[angina==1, CVD := "Angina"]
cvd_incidence_dt[tia==1, CVD := "TIA"]
#cvd_incidence_dt[pad==1, CVD := "PAD"]

cvd_incidence_dt <- cvd_incidence_dt[!is.na(CVD)]
uniqueN(cvd_incidence_dt, by = "patid") #316,439

#merge with statins to keep those without ever statins before first risk date
cvd_incidence_dt <- merge(cvd_incidence_dt[, c("patid", "CVD", "incident_date1")],
                       statins_pts[, c("patid", "issuedate")], by="patid", all.x=T)
uniqueN(cvd_incidence_dt, by = "patid") #316,439

#exclude statin contraindications
statin_contraindic_dt <- read_parquet(paste0(datafiles,"statin_contraindic.parquet"))
statin_contraindic_dt[, statin_excl := 1]
setnames(statin_contraindic_dt, "eventdate_gp", "excl_date1")
cvd_incidence_dt <- merge(cvd_incidence_dt, statin_contraindic_dt[,c("patid", "excl_date1", "statin_excl")], by = "patid", all.x=TRUE)
cvd_incidence_dt <- cvd_incidence_dt[is.na(excl_date1) | incident_date1<excl_date1]
uniqueN(cvd_incidence_dt, by = "patid") #310,460

#total number on statins
cvd_incidence_dt[!is.na(issuedate), statin_use := 1]
cvd_incidence_dt[, total_statin := sum(statin_use, na.rm=T), by="CVD"]

#subpopulation
cvd_incidence_dt[, subpop := .N, by = CVD]
cvd_incidence_dt[, curr_prev := total_statin/subpop]

####number of events avoided if those not on statins are treated
#number not on statins
cvd_incidence_dt[, no_statin_use := subpop-total_statin]
#absolute benefit (NNT)
cvd_incidence_dt[, NNT := 10]
cvd_incidence_dt[, events_avoided := round(no_statin_use/NNT)]

#prepare for export to excel
cvd_table <- cvd_incidence_dt
cvd_table <- cvd_table[, c("CVD", "subpop", "total_statin", "curr_prev", "no_statin_use", "NNT",
                                     "events_avoided")]
cvd_table<- unique(cvd_table)
setorder(cvd_table, -"curr_prev")


wb <- loadWorkbook(paste0(excel, "cvd_events_avoided.xlsx"))
names(wb)
writeData(wb, "cvd_events_R", cvd_table, startCol=2, startRow=9, rowNames=T)
saveWorkbook(wb,paste0(excel, "cvd_events_avoided.xlsx"), overwrite=T)


############################################################################################
#####################################Drug costs#############################################
############################################################################################
#import csv
library(readxl)
pca_costs_dt <- read_excel("J:\\EHR-Working\\Rutendo\\Postdoc_Statins_CVD\\Data\\datafiles\\health_files\\pca_summary_tables_2021_22_v001.xlsx",
                           sheet="Presentations", col_names=F)
pca_costs_dt <- setDT(pca_costs_dt)

######edit column names
column_names <- pca_costs_dt[4, ]
pca_costs_dt <- pca_costs_dt[-4, ]
colnames(pca_costs_dt) <- as.character(column_names)
colnames(pca_costs_dt) <- tolower(colnames(pca_costs_dt))
#add _ to spaces in column name
column_names <- colnames(pca_costs_dt)
updated_column_names <- gsub("\\s", "_", column_names)
colnames(pca_costs_dt) <- updated_column_names

#PCA shows Rx's dispensed and submitted to NHS Prescription Services
#20mg = BNF dose for primary prevention, 80mg=bnf dose for secondary prevention
atorv_costs <-pca_costs_dt[bnf_presentation_name=="Atorvastatin 20mg tablets"|bnf_presentation_name=="Atorvastatin 80mg tablets"]

atorv_costs <- atorv_costs[preparation_class=="01"]

atorv_costs <- atorv_costs[, c( "bnf_presentation_name",
                              "total_items", "total_quantity", "total_cost_(gbp)", "cost_per_item_(gbp)",
                               "cost_per_quantity_(gbp)", "quantity_per_item")]

wb <- loadWorkbook(paste0(excel, "cvd_events_avoided.xlsx"))
names(wb)
writeData(wb, "cvd_events_R", atorv_costs, startCol=2, startRow=17, rowNames=T)
saveWorkbook(wb,paste0(excel, "cvd_events_avoided.xlsx"), overwrite=T)


############################################################################################
#####################################CVD EVENTS COSTS#######################################
############################################################################################
cvd_costs_dt <- read_excel("J:\\EHR-Working\\Rutendo\\Postdoc_Statins_CVD\\Data\\datafiles\\health_files\\2_national_schedule_of_NHS_costs_FY21-22_v3.xlsx",
                           sheet="Total HRGs", col_names=F)
cvd_costs_dt <- setDT(cvd_costs_dt)

cvd_costs_dt <- cvd_costs_dt[-c(1:5)]

setnames(cvd_costs_dt, as.character(cvd_costs_dt[1,]))

# Remove the first row
cvd_costs_dt <- cvd_costs_dt[-1,]

# Convert column names to lowercase
setnames(cvd_costs_dt, tolower(names(cvd_costs_dt)))

cvd_costs_dt <- cvd_costs_dt[, c("currency", "currency description", "activity", "unit cost", "total cost")]

setnames(cvd_costs_dt, "currency description", "description")
setnames(cvd_costs_dt, "unit cost", "unit_cost")
setnames(cvd_costs_dt, "total cost", "total_cost")

#keep costs of CVD events
cvd_costs_dt[, description := tolower(description)]
cvd_costs <- cvd_costs_dt[str_detect(description, "myocardial infarction")|
                          str_detect(description, "stroke")| 
                            str_detect(description, "angina") |
                          str_detect(description, "transient ischaemic attack")|
                          str_detect(description, "peripheral arterial disease")]

cvd_costs[, mi := as.integer(str_detect(description, "myocardial infarction"))]
cvd_costs[, stroke := as.integer(str_detect(description, "stroke"))]
cvd_costs[, angina := as.integer(str_detect(description, "angina"))]
cvd_costs[, tia := as.integer(str_detect(description, "transient ischaemic attack"))]

#sum total costs of CVD event
cvd_costs[, total_cost := as.numeric(total_cost)]
cvd_costs[mi==1, mi_cost := sum(total_cost)]
cvd_costs[stroke == 1, stroke_cost := sum(total_cost)]
cvd_costs[angina == 1, angina_cost := sum(total_cost)]
cvd_costs[tia == 1, tia_cost := sum(total_cost)]

wb <- loadWorkbook(paste0(excel, "cvd_events_avoided.xlsx"))
names(wb)
writeData(wb, "cvd_costs_R", cvd_costs, startCol=2, startRow=2, rowNames=T)
saveWorkbook(wb,paste0(excel, "cvd_events_avoided.xlsx"), overwrite=T)
