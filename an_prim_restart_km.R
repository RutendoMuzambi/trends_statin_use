#######################################################################################-
#Project: Statins and CVD prevention
#Author: Rutendo Muzambi
#Description: Statin restarting - time from cessation to subsequent statin Rx for primary prevention
#######################################################################################
rm(list=ls())

############################################################
disc_pop <- read_parquet(paste0(datafiles, "prim_discont_fact.parquet"))
uniqueN(disc_pop, by = "patid") #303,292
disc_pop <- disc_pop[, c("patid", "discontinuation", "discont_date", "last_statin_date",
                         "incident_date1", "start_fup", "enddate", "first_statin_dt")]

#keep only those who discontinued
disc_pop <- disc_pop[discontinuation==1]
uniqueN(disc_pop, by = "patid") #119,853

#exclude those with CVD event prior to discontinuation as they'd be in secondary prevention
disc_pop <- disc_pop[is.na(incident_date1) | incident_date1>discont_date]
uniqueN(disc_pop, by = "patid") #116,468

#merge with all statin Rx's
statins_dt <- read_parquet(paste0(datafiles, "statins.parquet"))
statins_dt <- statins_dt[, c("patid", "issuedate")]

restart_pop <- merge(disc_pop,statins_dt, by="patid", all.x=T)
uniqueN(restart_pop, by = "patid") #116,468
rm(disc_pop)
#If statin Rx after discontinuation date code restart
setorder(restart_pop, "patid", "issuedate")
restart_pop[issuedate>=discont_date, statin_date1 := lapply(.SD, min, na.rm=T), by = "patid", .SDcols = "issuedate"]
restart_pop[issuedate>=discont_date, days_since := issuedate-discont_date]
#apply statin_date1 to rows for same patid
restart_pop[, statin_date1 := statin_date1[!is.na(statin_date1)][1L], by = patid]

restart_pop <- restart_pop[issuedate==statin_date1 | is.na(statin_date1)]
uniqueN(restart_pop, by = "patid") #116,468
restart_pop <- unique(restart_pop)
#restart_pop <- unique(sec_restart_pop, by = c("patid", "statin_date1"))


restart_pop[statin_date1>=discont_date, restart := 1]
restart_pop[is.na(restart), restart := 0]
restart_pop[, .N, by="restart"]
#update
# restart     N
# 1:       1 92583
# 2:       0 23885

restart_pop[, restart_date := statin_date1]

uniqueN(restart_pop, by="patid") #116,468

restart_pop[incident_date1>=discont_date, cvd_event := 1]
restart_pop[incident_date1<discont_date | is.na(incident_date1), cvd_event := 0]

###Merge characteristics
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))
restartpop_all <- merge(restart_pop, studypop_all[,c("patid", "gender", "dob", "start_fup", "ethnicity", "deprivation", "bmi", "bmi_cat")], by="patid")


restartpop_all[, age_disc := round(as.numeric(difftime(discont_date, dob, units = "days")) / 365.25)]

restartpop_all[, age_cat := cut(age_disc, breaks = c(24, 40, 50, 60, 70, Inf),
                                labels = c("25-39", "40-49", "50-59", "60-69", "70+"), include.lowest=T)]
restartpop_all[, .N, by=age_cat]

######Kaplan Meier - time from discontinuation to re-initiation
restartpop_all[, exit := pmin(enddate,restart_date,incident_date1, na.rm=T)]
restartpop_all[, time := as.numeric(exit-discont_date)/365.25]
restartpop_all[, days_since_exit := as.numeric(exit-discont_date)]

uniqueN(restartpop_all, by="patid") #142,584

restartpop_all[, .N, by="restart"]
# restart      N
# 1:       1 115043
# 2:       0  27541

#####save parquet
write_parquet(restartpop_all, paste0(datafiles, "prim_restart_fact.parquet"))



test <- restartpop_all[, c("patid", "discont_date", "incident_date1", "enddate",
                        "first_statin_dt", "statin_date1", "restart_date", "exit", 
                        "time", "days_since", "days_since_exit", "restart")]


summary(restartpop_all$days_since_exit)

str(restartpop_all)

#####################KAPLAN MEIER###############################################

survfit2(Surv(time, restart) ~ cvd_event, data = restartpop_all) %>% 
  ggsurvfit() +
  labs(
    #x = "Time to statin re-initiation (years)",
    title="CVD_event",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)


survfit2(Surv(time, restart) ~ 1, data = restartpop_all) %>% 
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


age_plot <- survfit2(Surv(time, restart) ~ age_cat, data = restartpop_all) %>% 
  ggsurvfit() +
  labs(
    #x = "Time to statin re-initiation (years)",
    title="Age group",
  ) + 
  #theme(axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank()) +
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)
#ggsave("restart_age_km.png", plot=last_plot(), path = graphs)

gender_plot <- survfit2(Surv(time, restart) ~ gender, data = restartpop_all) %>% 
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
#ggsave("restart_sex_km.png", plot=last_plot(), path = graphs)

restartpop_min <- restartpop_all[ethnicity!="missing"]
uniqueN(restartpop_min, by="patid") #

restartpop_min[, .N, by="ethnicity"]

eth_plot <-  survfit2(Surv(time, restart) ~ ethnicity, data = restartpop_min) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin re-initiation (years)",
    title="Ethnicity",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  theme(axis.title.y = element_blank()) + 
  xlim(0, 5)
#ggsave("restart_ethnic_km.png", plot=last_plot(), path = graphs)

restartpop_min <- restartpop_all[deprivation!="missing"]
uniqueN(restartpop_min, by="patid") #

restartpop_min[, .N, by="deprivation"]

depr_plot <- survfit2(Surv(time, restart) ~ deprivation, data = restartpop_min) %>% 
  ggsurvfit() +
  labs(
    x = "Time to statin re-initiation (years)",
    title="Deprivation",
  ) + 
  add_confidence_interval() +
  add_risktable() +
  xlim(0, 5)
#ggsave("restart_depriv_km.png", plot=last_plot(), path = graphs)


############COMBINE ALL PLOTS############
km_restart <- cowplot::plot_grid(age_plot, 
                                 gender_plot,
                                 eth_plot,
                                 depr_plot,
                                 align="cc",
                                 axis="1",
                                 ncol = 2
)
#Statin initiation in those with recorded risk score
print(km_restart)
ggsave("all_restart_prim_km.pdf", plot=last_plot(), width=15, height=15, path = final_graphs)






