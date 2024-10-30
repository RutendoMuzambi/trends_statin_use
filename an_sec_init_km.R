#######################################################################################-
#Project: Statins and CVD prevention
#Author: Rutendo Muzambi
#Description: Statin initiation - time from CVD event to initiation
#######################################################################################
rm(list=ls())

##################CVD##################################
statins_pts <- read_parquet(paste0(datafiles, "statins.parquet"))
cvd_incidence_dt <- read_parquet(paste0(datafiles, "cvd_incident_all.parquet"))
cvd_incidence_dt <- cvd_incidence_dt[, c("enddate", "dob") := NULL]
uniqueN(cvd_incidence_dt, by = "patid") #433,152

##load statin contraindications
statin_contraindic_dt <- read_parquet(paste0(datafiles,"statin_contraindic.parquet"))
statin_contraindic_dt[, statin_excl := 1]
setnames(statin_contraindic_dt, "eventdate_gp", "excl_date1")
cvd_incidence_dt <- merge(cvd_incidence_dt, statin_contraindic_dt[,c("patid", "excl_date1", "statin_excl")], by = "patid", all.x=TRUE)
cvd_incidence_dt <- cvd_incidence_dt[is.na(excl_date1) | incident_date1<excl_date1]
uniqueN(cvd_incidence_dt, by = "patid") #425,112

patients_dt <- read_parquet(paste0(datafiles, "clean_pt.prac.parquet"))
cvd_incidence_dt <- merge(patients_dt[,c("patid", "dob", "start_fup", "enddate", "gender")], cvd_incidence_dt, by = "patid")
uniqueN(cvd_incidence_dt, by = "patid") #425,112

#secondary prevention eligibilty
uniqueN(cvd_incidence_dt, by = "patid") 
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

cvd_incidence_dt <- cvd_incidence_dt[, c("patid", "start_fup", "enddate", "dob", "gender", "CVD", "incident_date1")]

#keep only CVD diagnosis on and after start of follow up
sec_elig_dt <- cvd_incidence_dt[incident_date1>=start_fup]
uniqueN(sec_elig_dt, by = "patid") #49,796
#keep cvd diagnosis before enddate
sec_elig_dt <- sec_elig_dt[incident_date1<=enddate]
uniqueN(sec_elig_dt, by = "patid") #42,982

#merge statins and cvd codes ##keeping everyone with CVD event
sec_elig_dt <- merge(sec_elig_dt,statins_pts[,c("patid", "statin_type", "issuedate")], by="patid", all.x=TRUE)

#create column for first statin following CVD
#tag statin Rx before CVD diagnosis
sec_elig_dt[!is.na(issuedate) & issuedate<incident_date1, prev_statin := 1]
#make constant
sec_elig_dt[, prev_statin := prev_statin[!is.na(prev_statin)][1L], by = patid]

#if 

#keep only those with first statin following CVD or no statin following CVD
sec_elig_dt <- sec_elig_dt[issuedate>=incident_date1|is.na(issuedate)]
uniqueN(sec_elig_dt, by = "patid") #39,991

#keep first statin Rx following CVD or missing statin
setorder(sec_elig_dt, patid, issuedate)
sec_elig_dt[!is.na(issuedate), dupl := rowid(patid)]

cvd_km <- sec_elig_dt[dupl==1|is.na(dupl)]
uniqueN(cvd_km, by = "patid") #39,991

#generate statin flag column 
cvd_km[!is.na(issuedate), statins := 1]
cvd_km[is.na(issuedate), statins := 0]

setnames(cvd_km, "issuedate", "statin_cvd1")
cvd_km[is.na(prev_statin)& !is.na(statin_cvd1), first_ever_statin :=1]

test <- cvd_km[!is.na(statin_cvd1) & is.na(prev_statin)] ##13,857 ever statin after CVD
#what about cvd event before start of follow up?
cvd_km[, exit := pmin(enddate,statin_cvd1, na.rm=T)]

cvd_km[, time := as.numeric(exit-incident_date1)/365.25]
cvd_km[, time_days := as.numeric(exit-incident_date1)]

###merge deprivation and ethnicity
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))
cvd_km <- merge(cvd_km, studypop_all[,c("patid", "ethnicity", "deprivation")], by="patid", all.x=TRUE)


# Get rid of non specific codes
cvd_km_red <- cvd_km[CVD!="Acute Coronary Syndrome (non-specific)" & CVD!="Coronary Heart Disease (non-specific)" & CVD!="cerebrovascular procedures"]
uniqueN(cvd_km_red, by = "patid") #29,114

#among those with and without previous statin Rx
cvd_type_plot <- survfit2(Surv(time_days, statins) ~ CVD, data = cvd_km_red) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin initiation following CVD diagnosis (days)",
    title="CVD subtype",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  theme(axis.text.x = element_blank(), axis.title.x= element_blank()) +
  xlim(0, 365)
#ggsave("initiat_cvd_km.png", plot=last_plot(), path = graphs)


#among those without previous statin Rx
cvd_km1 <- cvd_km_red[first_ever_statin==1]
uniqueN(cvd_km1, by = "patid") #11,438

survfit2(Surv(time_days, statins) ~ CVD, data = cvd_km1) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin initiation following CVD diagnosis (days)",
    title="Statin initiation among those with a CVD diagnosis",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 365)

#ggsave("cvd_firstever_km.png", plot=last_plot(), path = graphs)


###################AGE GROUP########################
# Calculate the age at stroke
library(lubridate)
cvd_km[, age_at_cvd := as.numeric(difftime(incident_date1, dob, units = "days")) / 365.25]

cvd_km[, age_group := cut(age_at_cvd, breaks = c(25, 40, 50, 60, 70, 80, Inf),
                                  labels = c("25-39", "40-49", "50-59", "60-69", "70-79", "80+"), include.lowest=T)]

cvd_km[, .N, by=age_group]

age_cvd_plot <- survfit2(Surv(time_days, statins) ~ age_group, data = cvd_km) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin initiation following CVD diagnosis (days)",
    title="Age Group",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  theme(axis.text.x = element_blank(), axis.title.x= element_blank()) +
  xlim(0, 365)
#ggsave("cvd_age_km.png", plot=last_plot(), path = graphs)

cvd_km[gender=="Male", gender := "Men"]
cvd_km[gender=="Female", gender := "Women"]

sex_cvd_plot <- survfit2(Surv(time_days, statins) ~ gender, data = cvd_km) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin initiation following CVD diagnosis (days)",
    title="Gender",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 365)
#ggsave("cvd_sex_km.png", plot=last_plot(), path = graphs)

cvd_km_eth <- cvd_km[ethnicity!="missing"]

eth_cvd_plot <- survfit2(Surv(time_days, statins) ~ ethnicity, data =cvd_km_eth) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin initiation following CVD diagnosis (days)",
    title="Ethnicity",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 365)
ggsave("cvd_ethnic_km.png", plot=last_plot(), path = graphs)

cvd_km_depr <- cvd_km[deprivation!="missing"]

depr_cvd_plot <- survfit2(Surv(time_days, statins) ~ deprivation, data = cvd_km_depr) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin initiation following CVD diagnosis (days)",
    title="Deprivation",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 365)
#ggsave("cvd_depriv_km.png", plot=last_plot(), path = graphs)


############COMBINE ALL PLOTS############
km_comb_cvd <- cowplot::plot_grid(cvd_type_plot,
                                    age_cvd_plot, 
                                    sex_cvd_plot,
                                    eth_cvd_plot,
                                    depr_cvd_plot,
                                    ncol = 2)
#Statin initiation following calculated QRISK2 eligibility
print(km_comb_cvd)
ggsave("all_cvd_km.pdf", plot=last_plot(), width=15, height=15, units = "in", dpi = 400, path = final_graphs)






