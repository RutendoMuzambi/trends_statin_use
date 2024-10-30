#######################################################################################-
#Project: Statins and CVD prevention
#Author: Rutendo Muzambi
#Description: Statin restarting - time from cessation to subsequent statin Rx for secondary
#######################################################################################

rm(list=ls())

############################################################
sec_disc_pop <- read_parquet(paste0(datafiles, "sec_discont_fact.parquet"))
uniqueN(sec_disc_pop, by = "patid") #122,743
sec_disc_pop <- sec_disc_pop[, c("patid", "discontinuation", "discont_date", "last_statin_date",
                         "incident_date1", "start_fup", "enddate", "CVD")]

#keep only those who discontinued
sec_disc_pop <- sec_disc_pop[discontinuation==1]
uniqueN(sec_disc_pop, by = "patid") #27,257

#merge with all statin Rx's
statins_dt <- read_parquet(paste0(datafiles, "statins.parquet"))
statins_dt <- statins_dt[, c("patid", "issuedate")]
sec_restart_pop <- merge(sec_disc_pop,statins_dt, by="patid", all.x=T)
uniqueN(sec_restart_pop, by = "patid") #27,257
rm(sec_disc_pop)
#If statin Rx after discontinuation date code restart
setorder(sec_restart_pop, "patid", "issuedate")
sec_restart_pop[issuedate>=discont_date, statin_date1 := lapply(.SD, min, na.rm=T), by = "patid", .SDcols = "issuedate"]
#apply statin_date1 to rows for same patid
sec_restart_pop[, statin_date1 := statin_date1[!is.na(statin_date1)][1L], by = patid]

sec_restart_pop <- sec_restart_pop[issuedate==statin_date1 | is.na(statin_date1)]
uniqueN(sec_restart_pop, by = "patid") #27,257
setorder(sec_restart_pop, "patid", -"issuedate")
sec_restart_pop <- unique(sec_restart_pop, by = c("patid", "statin_date1"))

sec_restart_pop[statin_date1>=discont_date, restart := 1]
sec_restart_pop[is.na(statin_date1), restart := 0]
sec_restart_pop[, .N, by="restart"]
# restart     N
# 1:       1 24876
# 2:       0  2381

sec_restart_pop[, restart_date := statin_date1]

######Kaplan Meier - time from discontinuation to re-initiation
sec_restart_pop[, exit := pmin(enddate,restart_date, na.rm=T)]
sec_restart_pop[, time := as.numeric(exit-discont_date)/365.25]

uniqueN(sec_restart_pop, by="patid") #29,154


###Merge characteristics
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))
sec_restart_pop <- merge(sec_restart_pop, studypop_all[,c("patid", "gender", "dob", "ethnicity", "deprivation", "bmi", "bmi_cat")], by="patid")

##age at discontinuation
sec_restart_pop[, age_disc := round(as.numeric(difftime(discont_date, dob, units = "days")) / 365.25)]
sec_restart_pop[, age_cat := cut(age_disc, breaks = c(24, 40, 50, 60, 70, Inf),
                                labels = c("25-39", "40-49", "50-59", "60-69", "70+"), include.lowest=T)]
sec_restart_pop[, .N, by=age_cat]


#####save parquet
write_parquet(sec_restart_pop, paste0(datafiles, "sec_restart_fact.parquet"))

# 
# ###CVD SUBTYPE
# cvd_incidence_dt <- read_parquet(paste0(datafiles, "cvd_incident_all.parquet"))
# uniqueN(cvd_incidence_dt , by = "patid") #
# 
# cvd_incidence_dt <- cvd_incidence_dt[, c("patid",  "mi", "angina",
#                                           "stroke", "tia", "pad")]
# 
# cvd_incidence_dt[mi==1, CVD := "MI"]
# cvd_incidence_dt[stroke==1, CVD := "Stroke"]
# cvd_incidence_dt[angina==1, CVD := "Angina"]
# cvd_incidence_dt[tia==1, CVD := "TIA"]
# cvd_incidence_dt[pad==1, CVD := "PAD"]

#MERGE WITH RESTART POPULATION
restartpop_cvd <- sec_restart_pop[CVD=="MI"|CVD=="Stroke"|CVD=="Angina"|CVD=="MI"|CVD=="PAD"|CVD=="TIA"]

restartpop_cvd [, .N, by="CVD"]

cvd_plot <- survfit2(Surv(time, restart) ~ CVD, data = restartpop_cvd) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin re-initiation (years)",
    title="CVD subtype",
  ) + 
  add_confidence_interval() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  add_risktable() +
  xlim(0, 5)
#ggsave("restart_cvdtype_km.png", plot=last_plot(), path = graphs)

##Age at which eligible for statins
sec_restart_pop[, age_disc := round(as.numeric(difftime(discont_date, dob, units = "days")) / 365.25)]

sec_restart_pop[, age_group := cut(age_disc, breaks = c(24, 40, 50, 60, 70, Inf),
                                  labels = c("25-39", "40-49", "50-59", "60-69", "70+"), include.lowest=T)]

sec_restart_pop[, .N, by=age_group]

age_plot <- survfit2(Surv(time, restart) ~ age_group, data = sec_restart_pop) %>% 
  ggsurvfit() +
  labs(
    #x = "Time to statin re-initiation (years)",
    title="Age group",
  ) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)
#ggsave("restart_sec_age_km.png", plot=last_plot(), path = graphs)

gender_plot <- survfit2(Surv(time, restart) ~ gender, data = sec_restart_pop) %>% 
  ggsurvfit() +
  labs(
    #x = "Time to statin re-initiation (years)",
    title="Gender",
  ) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)
#ggsave("restart_sec_sex_km.png", plot=last_plot(), path = graphs)

sec_restartpop_min <- sec_restart_pop[ethnicity!="missing"]
uniqueN(sec_restartpop_min, by="patid") #

sec_restartpop_min[, .N, by="ethnicity"]

eth_plot <-  survfit2(Surv(time, restart) ~ ethnicity, data = sec_restartpop_min) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin re-initiation (years)",
    title="Ethnicity",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  theme(axis.title.y = element_blank()) + 
  xlim(0, 5)
#ggsave("restart_sec_ethnic_km.png", plot=last_plot(), path = graphs)

sec_restartpop_min <- sec_restart_pop[deprivation!="missing"]
uniqueN(sec_restartpop_min, by="patid") #

sec_restartpop_min[, .N, by="deprivation"]
sec_restartpop_min[deprivation=="1 - least deprived", deprivation := "1"]

depr_plot <- survfit2(Surv(time, restart) ~ deprivation, data = sec_restartpop_min) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin re-initiation (years)",
    title="Deprivation",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)
#ggsave("restart_sec_depriv_km.png", plot=last_plot(), path = graphs)


############COMBINE ALL PLOTS############
km_restart <- cowplot::plot_grid( cvd_plot,
                                  age_plot, 
                                 gender_plot,
                                 eth_plot,
                                 depr_plot,
                                 align="cc",
                                 axis="1",
                                 ncol = 2
)
#Statin initiation in those with recorded risk score
print(km_restart)
ggsave("all_restart_sec_km.pdf", plot=last_plot(), width=10, height=10, path = final_graphs)




