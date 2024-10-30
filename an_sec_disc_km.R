#######################################################################################-
#Project: Statins and CVD prevention
#Author: Rutendo Muzambi
#Description: Statin discontinuation - time from initiation to cessation
#######################################################################################
rm(list=ls())

################################################################################################################
#######################################SECONDARY PREVENTION#####################################################
################################################################################################################
#load in first CVD event - secondary prevention
cvd_incidence_dt <- read_parquet(paste0(datafiles, "cvd_incident_all.parquet"))
uniqueN(cvd_incidence_dt , by = "patid") #433,152

cvd_incidence_dt <- cvd_incidence_dt[, c("patid", "incident_date1", "mi", "angina",
                                         "chd", "stroke", "tia", "pad", "acs", "cerebro_procs", "cvd_event")]
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

#merge with patient dataset
patients_dt <- read_parquet(paste0(datafiles, "clean_pt.prac.parquet"))
cvd_incidence_dt <- merge(cvd_incidence_dt,patients_dt[, c("patid", "start_fup", "enddate")], by="patid", all.x=T)
uniqueN(cvd_incidence_dt, by = "patid") #433,152

#keep those with incidence of CVD during study period
cvd_incidence_dt <- cvd_incidence_dt[incident_date1 >= start_fup]
cvd_incidence_dt <- cvd_incidence_dt[incident_date1 < enddate]
uniqueN(cvd_incidence_dt, by = "patid") #156,100

cvd_incidence_dt[, c("start_fup", "enddate") := NULL]

#merge statins and cvd codes ##keeping everyone with CVD event
statins_sec <- read_parquet(paste0(datafiles, "statins_sec.parquet"))
secondary_dt <- merge(cvd_incidence_dt,statins_sec, by="patid")
uniqueN(secondary_dt, by = "patid") #131,526

#keep those with statins on or after CVD event
secondary_dt <- secondary_dt[issuedate>=incident_date1]
uniqueN(secondary_dt, by = "patid") #122,769

#generate first statin
secondary_dt[, first_statin_dt := lapply(.SD, min, na.rm=T), by = "patid", .SDcols = "issuedate"]

#keep only those initated on or after start of follow up
secondary_dt <- secondary_dt[first_statin_dt >= start_fup]
secondary_dt <- secondary_dt[first_statin_dt < enddate]
uniqueN(secondary_dt, by = "patid") #122,743

#secondary_dt[, c("start_fup", "enddate") := NULL]

#generate initiator column
secondary_dt[first_statin_dt-incident_date1<=60, initiation_60 := 1]
secondary_dt[first_statin_dt-incident_date1>60, initiation_60 := 0]

#####Definition of discontinuation
#90 days without statin Rx - can do sensitivity analysis for 90 days which is cut-off used in previous studies
setorder(secondary_dt, "patid", "issuedate")
secondary_dt[, next_statin_dt := shift(issuedate, type = "lead"), by = patid]
secondary_dt[(next_statin_dt-runoutdate)>90, discontinuation := 1] 

secondary_dt[, count := .N, by=patid]
secondary_dt[, .N, by=count]
secondary_dt[count==1 & (enddate-runoutdate)>90, discontinuation := 1]

secondary_dt[discontinuation==1, flag_date := lapply(.SD, min, na.rm=T), by = "patid", .SDcols = "runoutdate"]

secondary_dt[, .N, by=discontinuation]
secondary_dt[flag_date!=runoutdate, discontinuation := NA]

secondary_dt[, disc_flag := discontinuation]
secondary_dt[, disc_flag := disc_flag[!is.na(disc_flag)][1L], by = patid]

#for those not without discontinuation
secondary_dt[, last_statin_date := lapply(.SD, max, na.rm=T), by = "patid", .SDcols = "issuedate"]
secondary_dt[is.na(discontinuation) & last_statin_date==issuedate, discontinuation := 0]
uniqueN(secondary_dt, by = "patid") #122,743

secondary_dt <- secondary_dt[(discontinuation==1 & disc_flag==1)|(discontinuation==0 & is.na(disc_flag))]
secondary_dt[discontinuation==1, discont_date := runoutdate+90]
secondary_dt[, c("count", "flag_date", "duration", "disc_flag") := NULL]

uniqueN(secondary_dt, by = "patid") #122,743

##Merge with study population dataset
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))
secondpop_all <- merge(secondary_dt, studypop_all[,c("patid", "gender", "dob", "age_startfup", "ethnicity", "deprivation", "bmi", "bmi_cat")], by="patid")
uniqueN(secondpop_all, by = "patid") #122,743

###format variables
secondpop_all[, c("ethnicity", "gender", "deprivation") := lapply(.SD, as.factor), .SDcols = c("ethnicity", "gender", "deprivation")]
secondpop_all[ethnicity=="missing", ethnicity := NA]
secondpop_all[deprivation=="missing", deprivation := NA]

secondpop_all[, gender := factor(gender, levels = c("Female", "Male"))]
secondpop_all[, ethnicity := factor(ethnicity, levels = c("White", "Black", "South Asian", "Mixed", "Other"))]
secondpop_all[, deprivation := factor(deprivation, levels = c("1 - least deprived", "2", "3", "4", "5 - most deprived"))]

###conditions for kaplan meier
secondpop_all[, exit := pmin(enddate,discont_date, na.rm=T)]
secondpop_all[, time := as.numeric(exit-first_statin_dt)/365.25]
str(secondpop_all)

#keep only those alive throughout run out period otherwise they won't be able discontinue
#immortal time bias?
#secondpop_all <- secondpop_all[enddate>runoutdate]
#uniqueN(secondpop_all, by = "patid") #18,465

secondpop_all[, age_initiation := round(as.numeric(difftime(first_statin_dt, dob, units = "days")) / 365.25)]

secondpop_all[, age_cat := cut(age_initiation, breaks = c(24, 40, 50, 60, 70, Inf),
                               labels = c("25-39", "40-49", "50-59", "60-69", "70+"), include.lowest=T)]

secondpop_all[, .N, by=age_cat]

write_parquet(secondpop_all, paste0(datafiles, "sec_discont_fact.parquet"))

secondpop_cvd <- secondpop_all

#not including cerebrovascular procedures, ACS incl unspecific codes, CHD including non specific codes

secondpop_cvd[mi==1, CVD_type := "MI"]
secondpop_cvd[stroke==1, CVD_type := "Stroke"]
secondpop_cvd[angina==1, CVD_type := "Angina"]
secondpop_cvd[tia==1, CVD_type := "TIA"]
secondpop_cvd[pad==1, CVD_type := "PAD"]

cvd_plot <- survfit2(Surv(time, discontinuation) ~ CVD_type, data = secondpop_cvd) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin discontinuation (years)",
    title="CVD subtype",
  ) + 
  add_confidence_interval() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  add_risktable() +
  xlim(0, 5)
ggsave("discont_cvdtype_km.png", plot=last_plot(), path = graphs)


age_plot <- survfit2(Surv(time, discontinuation) ~ age_cat, data = secondpop_all) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin discontinuation (years)",
    title="Age Group",
  ) + 
  add_confidence_interval() +
  theme(axis.title.x = element_blank()) +
  add_risktable() +
  xlim(0, 5)
ggsave("discont_cvd_age_km.png", plot=last_plot(), path = graphs)

secondpop_all[gender=="Male", gender := "Men"]
secondpop_all[gender=="Female", gender := "Women"]

sex_plot <- survfit2(Surv(time, discontinuation) ~ gender, data = secondpop_all) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin discontinuation (years)",
    title="Gender",
  ) + 
  add_confidence_interval() +
  theme(axis.title.x = element_blank()) +
  add_risktable() +
  xlim(0, 5)
ggsave("discont_cvd_sex_km.png", plot=last_plot(), path = graphs)

secondpop_eth <- secondpop_all[!is.na(ethnicity)]

eth_plot <- survfit2(Surv(time, discontinuation) ~ ethnicity, data = secondpop_eth) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin discontinuation (years)",
    title="Ethnicity",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)
ggsave("discont_cvd_ethnic_km.png", plot=last_plot(), path = graphs)

secondpop_depr <- secondpop_all[!is.na(deprivation)]


depr_plot <- survfit2(Surv(time, discontinuation) ~ deprivation, data = secondpop_depr) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin discontinuation (years)",
    title="Deprivation",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)
ggsave("discont_cvd_depriv_km.png", plot=last_plot(), path = graphs)

survfit2(Surv(time, discontinuation) ~ bmi_cat, data = secondpop_all) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin discontinuation (years)",
    title="Statin discontinuation following statin initiation, stratified by bmi category",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)
ggsave("discont_cvd_bmi_km.png", plot=last_plot(), path = graphs)


############COMBINE ALL PLOTS############
sec_disc_comb <- cowplot::plot_grid(cvd_plot,
                                 age_plot, 
                                 sex_plot,
                                 eth_plot,
                                 #NULL,
                                 depr_plot,
                                 align="cc",
                                 axis="1",
                                 ncol = 2
)
#Statin initiation in those with recorded risk score
print(sec_disc_comb)
ggsave("all_sec_discont_km.pdf", plot=last_plot(), width=10, height=10, path = final_graphs)


secondpop_all[first_statin_dt<"2020-03-01", pandemic := "Pre-pandemic"]
secondpop_all[first_statin_dt>="2020-03-01", pandemic := "Pandemic"]

survfit2(Surv(time, discontinuation) ~ pandemic, data = secondpop_all) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin discontinuation (years)",
    title="Statin discontinuation following statin initiation, stratified by pandemic period",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 2)
ggsave("discont_cvd_pand_km.png", plot=last_plot(), path = graphs)


survfit2(Surv(time, discontinuation) ~ statin_type, data = secondpop_all) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin discontinuation (years)",
    title="Statin discontinuation following statin initiation, stratified by statin type",
  ) + 
  add_confidence_interval() +
  add_risktable(times= NULL) +
  xlim(0, 5) 
ggsave("discont_cvd_statintype_km.png", plot=last_plot(), path = graphs)

