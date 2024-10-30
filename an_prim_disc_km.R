#######################################################################################-
#Project: Statins and CVD prevention
#Author: Rutendo Muzambi
#Description: Statin discontinuation - time from initiation to cessation for primary prevention
#######################################################################################

###Definition of cessation
rm(list=ls())

#explore statin restart
statins_dt <- read_parquet(paste0(datafiles, "statins.parquet"))
uniqueN(statins_dt, by = "patid") #992,248

statins_dt[, .N, by="statin_type"]

##import dosages
dosage_dt <- haven::read_dta(file = stata_dosages_22)
dosage_dt <- as.data.table(dosage_dt)

statins_dt <- merge(statins_dt, dosage_dt, by ="dosageid", all.x = T)

statins_dt[, duration := as.numeric(duration)]

#Define duration of statins using quantity and daily dose
statins_dt[duration==0 & !is.na(quantity) & !is.na(daily_dose) & quantity!=0 & daily_dose!=0, duration := quantity/daily_dose]
statins_dt[duration==0 & (daily_dose==0 | is.na(daily_dose)) & quantity!=0 & !is.na(quantity), duration := quantity]
statins_dt[duration<=1, duration := 28]

#only keep the prescription with the longest duration for a given date if multiple Rx on same date for the same patient or runoutdate will not be created
setorder(statins_dt, patid, issuedate, -duration)
statins_dt[, dupl := rowid(patid, issuedate)]
statins_dt[, .N, by=dupl] 
statins_dt <- statins_dt[dupl==1]

#run out date based on issuedate and duration
statins_dt[,runoutdate := issuedate+duration]

#first statin date - statin initiation
#delete statins before dob
patients_dt <- read_parquet(paste0(datafiles, "clean_pt.prac.parquet"))
patients_dt <- patients_dt[,c("patid", "dob", "start_fup", "enddate")]

statins_dt <- merge(statins_dt,patients_dt, by="patid", all.x=T)
statins_dt <- statins_dt[issuedate>dob]

statins_dt <- statins_dt[,c("patid", "statin_type", "issuedate", "runoutdate", "duration", "start_fup", "enddate")]
uniqueN(statins_dt, by = "patid") #
write_parquet(statins_dt, paste0(datafiles, "statins_sec.parquet"))

setorder(statins_dt, "patid", "issuedate")
statins_dt[, first_statin_dt := lapply(.SD, min, na.rm=T), by = "patid", .SDcols = "issuedate"]

#get rid of initiators before study period or Rx after enddate
uniqueN(statins_dt, by = "patid") #992,248
statins_dt <- statins_dt[first_statin_dt >= start_fup]
statins_dt <- statins_dt[first_statin_dt < enddate]
uniqueN(statins_dt, by = "patid") #481,790 with first statin Rx during study period

#if duration more than 180 days change to the median (28 days)
#statins_dt[duration>180, duration := 28]

#statin restarting - time from grace period to first subsequent statin Rx
# if issuedate in subeseqent Rx is greater than runoutdate

statins_dt[, cessation_days := shift(issuedate, type = "lead") - runoutdate, by = list(patid)]

#keep the rows with negative values i.e those that are past runout date
cessation_dt <- statins_dt[cessation_days<0]
cessation_dt <- cessation_dt[, cessation_days := as.numeric(abs(cessation_days))]

cessation_dt[cessation_days>=200, cessation_days := 200]
ggplot(cessation_dt, aes(x=cessation_days)) + 
  geom_histogram(binwidth=10) +
  scale_x_continuous(breaks=seq(0,200,28)) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE, digits = 2))
ggsave("time_to_restart.png", plot=last_plot(), path = graphs)

summary(cessation_dt$cessation)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.50    3.00    7.00   24.55   21.00 4808.00 

#####Definition of discontinuation
#90 days without statin Rx - can do sensitivity analysis for 90 days which is cut-off used in previous studies
#create date of next statin prescription 
setorder(statins_dt, "patid", "issuedate")
statins_dt[, next_statin_dt := shift(issuedate, type = "lead"), by = patid]
statins_dt[(next_statin_dt-runoutdate)>90, discontinuation := 1] 

statins_dt[, count := .N, by=patid]
statins_dt[count==1 & (enddate-runoutdate)>90, discontinuation := 1]
statins_dt[discontinuation==1, flag_date := lapply(.SD, min, na.rm=T), by = "patid", .SDcols = "runoutdate"]
statins_dt[flag_date!=runoutdate, discontinuation := NA]
statins_dt[, disc_flag := discontinuation]
statins_dt[, disc_flag := disc_flag[!is.na(disc_flag)][1L], by = patid]

statins_dt[is.na(discontinuation), last_statin_date := lapply(.SD, max, na.rm=T), by = "patid", .SDcols = "issuedate"]
#statins_dt[last_statin_date==issuedate, discontinuation := 0]
statins_dt[is.na(discontinuation) & last_statin_date==issuedate, discontinuation := 0]
uniqueN(statins_dt, by = "patid") #373,266
statins_dt <- statins_dt[(discontinuation==1 & disc_flag==1)|(discontinuation==0 & is.na(disc_flag))]
statins_dt[discontinuation==1, discont_date := runoutdate+90]
statins_dt[, c("count", "flag_date", "duration", "disc_flag", "start_fup", "enddate", "cessation_days") := NULL]
uniqueN(statins_dt, by = "patid") #373,266

#exclude those on secondary prevention
cvd_dt <- read_parquet(paste0(datafiles, "cvd_incident_all.parquet"))
uniqueN(cvd_dt, by = "patid") #433,152

cvd_dt <- cvd_dt[,c("patid", "cvd", "incident_date1")]

primarypop_dt<- merge(statins_dt,cvd_dt, by="patid", all.x=TRUE)
uniqueN(primarypop_dt, by = "patid") #373,266

#######primary prevention eligibility
primarypop_dt[incident_date1 <= first_statin_dt & !is.na(incident_date1), exclude := 1]
primarypop_dt <- primarypop_dt[is.na(exclude)]
uniqueN(primarypop_dt, by = "patid") #303,292 without CVD event before statin initiation

#######2. Those with CVD risk scores
cvdrisk_all <- read_parquet(paste0(datafiles, "cvdrisk_all.parquet"))
uniqueN(cvdrisk_all, by = "patid") #5,000,000

cvdrisk_all <- cvdrisk_all[!is.na(risk_score)]#
cvdrisk_all[, recorded_risk := 1]

setnames(cvdrisk_all, "obsdate", "risk_date")
cvdrisk_all <- cvdrisk_all[, c("patid", "recorded_risk", "risk_date", "risk_cat")]

#######3. Those with high scores
cvdrisk_all <- read_parquet(paste0(datafiles, "cvdrisk_all.parquet"))
uniqueN(cvdrisk_all, by = "patid") #5,000,000
cvdrisk_all <- merge(cvdrisk_all,patients_dt[, c("patid", "dob")], by="patid", all.x=T)
cvdrisk_all <- cvdrisk_all[obsdate>dob]

cvdrisk_all <- cvdrisk_all[!is.na(risk_score)]#

#Match eligibility with updating guidelines. Some people with risk dates coded incorrectly - 1860s
cvdrisk_all[, year := year(obsdate)]
cvdrisk_all[risk_score>=20 & year<2014, eligible := 1]
cvdrisk_all[risk_score<20 & year<2014, eligible := 0]

cvdrisk_all[risk_score>=10 & year>=2014, eligible := 1]
cvdrisk_all[risk_score<10 & year>=2014, eligible := 0]

cvdrisk_all[,c("dob", "enddate") := NULL]

#merge with eligible for statins
primarypop_dt <- merge(primarypop_dt,patients_dt[, c("patid", "enddate")], by="patid", all.x=T)
primarypop_dt <- merge(primarypop_dt, cvdrisk_all, by="patid", all.x=T)

###initiators being those with most recent risk score and non being those without risk score or outdated
setnames(primarypop_dt, "obsdate", "risk_date")
primarypop_dt[(first_statin_dt-risk_date<=1826.25) & !is.na(risk_date), initiation_5y := 1]
primarypop_dt[, initiation_5y := initiation_5y[!is.na(initiation_5y)][1L], by = patid]

primarypop_dt[initiation_5y==1, risk_date1 := risk_date[which.min(abs(first_statin_dt-risk_date))], by = "patid"]

primarypop_dt[, count := .N, by=patid]

primarypop_dt <- primarypop_dt[(count>1 & initiation_5y==1 & risk_date1==risk_date)|(count==1)|is.na(initiation_5y)]
uniqueN(primarypop_dt, by="patid") #380,733

primarypop_dt[, dupl := rowid(patid)]
primarypop_dt[, .N, by=dupl]
primarypop_dt <- primarypop_dt[dupl==1]
uniqueN(primarypop_dt, by="patid") #380,733
primarypop_dt[, c("dupl", "count") := NULL]

primarypop_dt[(first_statin_dt-risk_date>1826.25|risk_date<start_fup), initiation_5y := NA]

######Kaplan Meier - time from initiation  to discontinuation
###
#primarypop_dt[discontinuation==1, discont_date := runoutdate]
primarypop_dt[, exit := pmin(enddate,discont_date, incident_date1, na.rm=T)]

primarypop_dt[, time := as.numeric(exit-first_statin_dt)/365.25]

str(primarypop_dt)
uniqueN(primarypop_dt, by="patid") #380,733

###Merge characteristics
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))
primarypop_all <- merge(primarypop_dt, studypop_all[,c("patid", "gender", "dob", "ethnicity", "deprivation", "bmi", "bmi_cat")], by="patid")

primarypop_all[, age_initiation := round(as.numeric(difftime(first_statin_dt, dob, units = "days")) / 365.25)]

primarypop_all[, age_cat := cut(age_initiation, breaks = c(24, 40, 50, 60, 70, Inf),
                                labels = c("25-39", "40-49", "50-59", "60-69", "70+"), include.lowest=T)]

primarypop_all[, .N, by=age_cat]

#####save parquet
write_parquet(primarypop_all, paste0(datafiles, "prim_discont_fact.parquet"))




#stratify by subpopulation - add low risk - who is in this group - one subgroup containing everyone
primarypop_all[!is.na(risk_date1) & eligible==1, subpop := "Above threshold"]
primarypop_all[!is.na(risk_date1) & eligible==0, subpop := "Below threshold"]
primarypop_all[is.na(initiation_5y), subpop := "No risk"]
primarypop_all[, .N, by=subpop]

uniqueN(primarypop_all, by = "patid") #63,764

subpop_plot1 <- survfit2(Surv(time, discontinuation) ~ subpop, data = primarypop_all) %>% 
  ggsurvfit() +
  labs(
    #x = "Time to statin discontinuation (years)",
    title="CVD risk assessment",
  ) + 
  theme(axis.title.x = element_blank()) +
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5) 
#ggsave("discont_subpop_km.png", plot=last_plot(), path = graphs)

##Age at initiation
primarypop_all[, age_initiation := round(as.numeric(difftime(first_statin_dt, dob, units = "days")) / 365.25)]

primarypop_all[, age_cat := cut(age_initiation, breaks = c(24, 40, 50, 60, 70, Inf),
                                  labels = c("25-39", "40-49", "50-59", "60-69", "70+"), include.lowest=T)]

primarypop_all[, .N, by=age_cat]

age_plot1 <- survfit2(Surv(time, discontinuation) ~ age_cat, data = primarypop_all) %>% 
  ggsurvfit() +
  labs(
    #x = "Time to statin discontinuation (years)",
    title="Age group",
  ) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)
#ggsave("discont_age_km.png", plot=last_plot(), path = graphs)

sex_plot1 <- survfit2(Surv(time, discontinuation) ~ gender, data = primarypop_all) %>% 
  ggsurvfit() +
  labs(
    #x = "Time to statin discontinuation (years)",
    title="Gender",
  ) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)
#ggsave("discont_sex_km.png", plot=last_plot(), path = graphs)

primarypop_min <- primarypop_all[ethnicity!="missing"]
uniqueN(primarypop_min, by="patid") #

primarypop_min[, .N, by="ethnicity"]

eth_plot1 <- survfit2(Surv(time, discontinuation) ~ ethnicity, data = primarypop_min) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin discontinuation (years)",
    title="Ethnicity",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  theme(axis.title.y = element_blank()) + 
  xlim(0, 5)
#ggsave("discont_ethnic_km.png", plot=last_plot(), path = graphs)

primarypop_min <- primarypop_all[deprivation!="missing"]
uniqueN(primarypop_min, by="patid") #

primarypop_min[, .N, by="deprivation"]

depr_plot1 <- survfit2(Surv(time, discontinuation) ~ deprivation, data = primarypop_min) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin discontinuation (years)",
    title="Deprivation",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)
#ggsave("discont_depriv_km.png", plot=last_plot(), path = graphs)

############COMBINE ALL PLOTS############
prim_disc_comb <- cowplot::plot_grid(subpop_plot1,
                                 age_plot1, 
                                 sex_plot1,
                                 eth_plot1,
                                 #NULL,
                                 depr_plot1,
                                 align="cc",
                                 axis="1",
                                 ncol = 2
)
#Statin initiation in those with recorded risk score
print(prim_disc_comb)
ggsave("disc_prim_km.pdf", plot=last_plot(),  width=10, height=10, path = final_graphs)

