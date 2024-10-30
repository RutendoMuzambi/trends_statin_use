#######################################################################################-
#Project: Statins and CVD prevention
#Author: Rutendo Muzambi
#Description: kaplan meier plots for time to statin initiation
#######################################################################################

######################################################################################################
rm(list=ls())

######################################################################################################

#INITIATION #Defining statin user
statins_dt <- read_parquet(paste0(datafiles, "statins.parquet"))
uniqueN(statins_dt, by = "patid") #992,248
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

#keep statins after start of follow up and before enddate
statins_pts <- statins_pts[issuedate>=start_fup]
uniqueN(statins_pts, by = "patid") #96,343
#keep first statin Rx before enddate
statins_pts<- statins_pts[issuedate<=enddate]
uniqueN(statins_pts, by = "patid") #96,324

statins_pts <- statins_pts[, c("patid", "statin_type", "issuedate")]
setnames(statins_pts, "issuedate", "first_statin_dt")

write_parquet(statins_pts, paste0(datafiles, "incident_statins.parquet"))


#########################################################################################
##################Time to statin prescription among those with recorded CVD risk
#########################################################################################
cvdrisk_all <- read_parquet(paste0(datafiles, "cvdrisk_all.parquet"))
uniqueN(cvdrisk_all, by = "patid") #1,000,000

#Match eligibility with updating guidelines
cvdrisk_all[, year := year(obsdate)]
cvdrisk_all[risk_score>=20 & year<2014, eligible := 1]
cvdrisk_all[risk_score<20 & year<2014, eligible := 0]

cvdrisk_all[risk_score>=10 & year>=2014, eligible := 1]
cvdrisk_all[risk_score<10 & year>=2014, eligible := 0]

setnames(cvdrisk_all, "obsdate", "risk_date")

cvdrisk_all <- cvdrisk_all[!is.na(risk_score) & !is.na(risk_date)]#362,240

#risk record during study period
cvdrisk_all <- cvdrisk_all[risk_date>=start_fup]
cvdrisk_all <- cvdrisk_all[risk_date<=enddate]

#keep first recorded risk
setorder(cvdrisk_all, patid, risk_date)
cvdrisk_all <- cvdrisk_all[order(risk_date), head(.SD, 1L), by="patid"]

uniqueN(cvdrisk_all, by = "patid") #1,577,405

#read and merge exclusions
exclusions_all <- read_parquet(paste0(datafiles, "exclusions_all.parquet"))
uniqueN(exclusions_all, by = "patid") #73,755
exclusions_all <- exclusions_all[,c("patid", "excl_date1", "excl_flag")]

cvdrisk_all <- merge(cvdrisk_all,exclusions_all, by="patid", all.x=TRUE)
uniqueN(cvdrisk_all, by = "patid") #1,000,000
cvdrisk_all <- cvdrisk_all[,c("patid", "risk_tool", "risk_cat", "risk_score", "risk_date", "excl_date1", "excl_flag", "start_fup", "enddate", "eligible")]

#######risk scoring eligibility
#exclude those with cvd event before start of follow up
cvdrisk_all[excl_date1 <= risk_date & !is.na(excl_date1), exclude := 1]
cvdrisk_all[is.na(excl_date1), exclude := 0]
cvdrisk_all[exclude==0|is.na(exclude), primary_prev := 1]
recordrisk_pop <- cvdrisk_all[exclude==0|is.na(exclude)]

uniqueN(recordrisk_pop, by = "patid") #290,010

#merge with incident cvd to censor when cvd occurs
#excl_dt <- read_parquet(paste0(datafiles, "excl_init_all.parquet"))
#uniqueN(excl_dt, by = "patid") #104,225
#excl_dt <- excl_dt[,c("patid", "excl_date1")]

#recordrisk_pop <- merge(recordrisk_pop,excl_dt, by="patid", all.x=TRUE)

##merge those with first ever statin prescription
incident_statins <- read_parquet(paste0(datafiles, "incident_statins.parquet"))
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))
recordrisk_pop <- merge(recordrisk_pop,incident_statins[, c("patid", "first_statin_dt")], by="patid", all.x=TRUE)

recordrisk_pop <- recordrisk_pop[first_statin_dt>=risk_date|is.na(first_statin_dt)]
recordrisk_pop[, exit := pmin(enddate,first_statin_dt,excl_date1, na.rm=T)]
recordrisk_pop[, time := as.numeric(exit-risk_date)/365.25]
recordrisk_pop[!is.na(first_statin_dt), statins := 1]
recordrisk_pop[is.na(first_statin_dt), statins := 0]

uniqueN(recordrisk_pop, by = "patid") #280,521

recordrisk_km <- merge(recordrisk_pop, studypop_all[,c("patid", "gender","age_startfup","dob", "ethnicity", "deprivation", "bmi", "bmi_cat")], by="patid")

recordrisk_km[, ageatrisk := (unclass(risk_date) - unclass(dob)) / 365.25]
recordrisk_km[, ageatrisk := round(ageatrisk, digits=0)]

#recordrisk_km[, age_group := cut(ageatrisk, breaks = c(24, 40, 75, 85, Inf),
#                                 labels = c("25-39", "40-74", "75-84", "85+"), include.lowest=T)]

recordrisk_km[, age_group := cut(ageatrisk, breaks = c(24, 40, 50, 60, 70, Inf),
                               labels = c("25-39", "40-49", "50-59", "60-69", "70+"), include.lowest=T)]


str(recordrisk_km)

risk_plot2 <- survfit2(Surv(time, statins) ~ risk_cat, data = recordrisk_km) %>% 
  ggsurvfit() +
     labs(
         #x = "Time to statin initiation (years)",
           title="A. Recorded CVD risk category",
      ) + 
     theme(axis.title.x = element_blank()) +
  add_confidence_interval() +
  add_risktable() +
  theme(legend.text = element_text(size = 10)) +
  xlim(0, 5)
#ggsave("initiation_recordrisk_km.png", plot=last_plot(), path = graphs)

recordrisk_km[, .N, by="age_group"]

recordrisk_age <- recordrisk_km[age_group!="85+"]

age_plot2 <- survfit2(Surv(time, statins) ~ age_group, data = recordrisk_age) %>% 
  ggsurvfit() +
  labs(
    #x = "Time to statin initiation (years)",
    title="B. Age group",
  ) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  add_confidence_interval() +
  add_risktable() +
  theme(legend.text = element_text(size = 10)) +
  xlim(0, 5)
#ggsave("recordrisk_age_km.png", plot=last_plot(), path = graphs)

recordrisk_km[, .N, by=age_group]

recordrisk_km[, gender := ifelse(gender == "Male", "Men", ifelse(gender == "Female", "Women", gender))]

sex_plot2 <- survfit2(Surv(time, statins) ~ gender, data = recordrisk_km) %>% 
  ggsurvfit() +
    labs(
        x = "Time to statin initiation (years)",
        title="C. Gender",
       ) + 
  theme(axis.title.x = element_blank()) +
  add_confidence_interval() +
  add_risktable() +
  theme(legend.text = element_text(size = 10)) +
  xlim(0, 5)
#ggsave("recordrisk_sex_km.png", plot=last_plot(), path = graphs)

recordrisk_eth <- recordrisk_km[ethnicity!="missing"]

eth_plot2 <- survfit2(Surv(time, statins) ~ ethnicity, data = recordrisk_eth) %>% 
  ggsurvfit() +
     labs(
         x = "Time to statin initiation (years)",
         title="D.Ethnicity",
       ) +  
  theme(axis.title.y = element_blank()) +
  add_confidence_interval() +
  add_risktable() +
  theme(legend.text = element_text(size = 10)) +
  xlim(0, 5)
#ggsave("recordrisk_ethnic_km.png", plot=last_plot(), path = graphs)

recordrisk_depriv <- recordrisk_km[deprivation!="missing"]

depr_plot2 <- survfit2(Surv(time, statins) ~ deprivation, data = recordrisk_depriv) %>% 
  ggsurvfit() +
    labs(
         x = "Time to statin initiation (years)",
        title="E. Deprivation",
      ) + 
  add_confidence_interval() +
  add_risktable() +
  #theme(legend.position = "right") +
  theme(legend.text = element_text(size = 10)) +
  xlim(0, 5)


############COMBINE ALL PLOTS############
km_comb_recrisk <- cowplot::plot_grid(risk_plot2,
                                    age_plot2, 
                                    sex_plot2,
                                    eth_plot2,
                                    #NULL,
                                    depr_plot2,
                                    align="cc",
                                    axis="1",
                                    ncol = 2
                                    )
#Statin initiation in those with recorded risk score
print(km_comb_recrisk)
ggsave("all_recrisk_km.pdf", plot=last_plot(), width=20, height=15, units = "in", dpi = 400, path = final_graphs)



