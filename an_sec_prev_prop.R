#--------------------------------------------------------------------------------------
# Project: Statins and CVD prevention
# Author: Rutendo Muzambi
# Description: Secondary prevention prevalence monthly proportions of statin user
#--------------------------------------------------------------------------------------

######################################################################################################
#rm(list=ls())

######################################################################################
###############Defining statin user
statins_dt <- read_parquet(paste0(datafiles, "statins.parquet"))
uniqueN(statins_dt, by = "patid") #197,996

##drug type
statins_dt[, .N, by=productname]
statins_list <- c('atorvastatin', 'fluvastatin', 'pravastatin', 'rosuvastatin', 'simvastatin')

#dose
dose_list <- c('100MICROGRAM', '10MG', '200MICROGRAM',  '20MG', "20MG + 10MG", '20MG/5ML',
               '300MICROGRAM', '400MICROGRAM', '40MG', "40MG + 10MG", '5MG', '80MG', "80MG + 10MG")

##import dosages
dosage_dt <- haven::read_dta(file = stata_dosages_22)

dosage_dt <- as.data.table(dosage_dt)

statins_dt <- merge(statins_dt, dosage_dt, by ="dosageid", all.x = T)

summary(statins_dt$duration)

statins_dt[, duration := as.numeric(duration)]

#Define duration of statins using quantity and daily dose
statins_dt[duration==0 & !is.na(quantity) & !is.na(daily_dose) & quantity!=0 & daily_dose!=0, duration := quantity/daily_dose]
statins_dt[duration==0 & (daily_dose==0 | is.na(daily_dose)) & quantity!=0 & !is.na(quantity), duration := quantity]
statins_dt[duration<=1, duration := 28]

statins_dt[, .N, by=duration]

#only keep the prescription with the longest duration for a given date if multiple Rx on same date for the same patient or runoutdate will not be created
setorder(statins_dt, patid, issuedate, -duration)
statins_dt[, dupl := rowid(patid, issuedate)]
statins_dt[, .N, by=dupl]
statins_dt <- statins_dt[dupl==1]

#remove prescriptions in period not interested in
statins_dt[,runoutdate := issuedate+duration]
statins_dt[, year := year(issuedate)]
statins_dt <- statins_dt[year>="2008"] ##study starts in 2009

#if duration more than 180 days change to the median (28 days)
statins_dt[duration>180, duration := 28]

statins_dt <- statins_dt[,c("patid", "issuedate", "duration", "runoutdate", "drug_type", "statin_type")]

rm(dosage_dt)
# ##################################################Secondary prevention###############################################
# #load in first CVD event - secondary prevention
cvd_incidence_dt <- read_parquet(paste0(datafiles, "cvd_incident_all.parquet"))
uniqueN(cvd_incidence_dt, by = "patid") #433,152
##load statin contraindications
statin_contraindic_dt <- read_parquet(paste0(datafiles,"statin_contraindic.parquet"))
statin_contraindic_dt[, statin_excl := 1]
setnames(statin_contraindic_dt, "eventdate_gp", "excl_date1")
cvd_incidence_dt <- merge(cvd_incidence_dt, statin_contraindic_dt[,c("patid", "excl_date1", "statin_excl")], by = "patid", all.x=TRUE)
cvd_incidence_dt <- cvd_incidence_dt[is.na(excl_date1) | incident_date1<excl_date1]
uniqueN(cvd_incidence_dt, by = "patid") #425,112

#merge statins and cvd codes ##keeping everyone with CVD event
secondary_dt <- merge(cvd_incidence_dt,statins_dt, by="patid", all.x=TRUE)

#secondary prevention eligibility
uniqueN(secondary_dt, by = "patid") #425,112
secondary_dt[, .N, by=cvd_event] #no missing cvd events
secondary_dt <- secondary_dt[, c("patid", "cvd_event", "issuedate","runoutdate", "incident_date1", "excl_date1", "drug_type", "statin_type")]
#generate statin flag column
secondary_dt[!is.na(issuedate), statin := 1]
#if statins issued on the date of or after first cvd event then statins issued under secondary prevention
secondary_dt[issuedate>=incident_date1 & !is.na(issuedate) & !is.na(incident_date1),  second_prev := 1]
#if statins issued before first cvd event then primary prevention
secondary_dt[issuedate<incident_date1 & !is.na(issuedate) & !is.na(incident_date1),  second_prev := 0]

secondary_dt <- secondary_dt[second_prev==1|is.na(second_prev)]
uniqueN(secondary_dt, by = "patid") #416,849


#create matching prop dates column in secondary dataset
# secondary_dt[, prop_dates := issuedate]
secondary_dt[, prop_dates := ceiling_date(issuedate, "month") -1]

setorder(secondary_dt, patid, prop_dates)

sec_sample_dt <- secondary_dt[!is.na(prop_dates)]

#create data table with month and years of study period for those eligible for secondary preventon
month_year_dt <- secondary_dt[,c("patid", "incident_date1", "cvd_event")]
uniqueN(month_year_dt, by = "patid") #416,849
month_year_dt <- unique(month_year_dt) #delete repeating patids

month_year_dt <- month_year_dt[rep(seq_len(nrow(month_year_dt)), each = 153), ]
prop_dates <- seq(as.Date("2009/05/01"), as.Date("2022/01/01"), "month", format="%d/%m/%Y") -1
month_year_dt <- cbind(month_year_dt, prop_dates)
setorder(month_year_dt, patid, prop_dates)

#append
second_prev <- rbind(month_year_dt,sec_sample_dt, fill=T)
uniqueN(second_prev, by = "patid") #416,849
#setorder(second_prev, patid, incident_date1, prop_dates, na.last=T)
setorder(second_prev, patid, prop_dates)

#fill in constant variables
constant_vars <- c("issuedate", "runoutdate", "second_prev", "statin",  "drug_type", "statin_type", "excl_date1")

setorder(second_prev, patid, prop_dates)
id_change = second_prev[, c(TRUE, patid[-1] != patid[-.N])]
second_prev[, (constant_vars) := lapply(.SD, function(x) x[cummax(((!is.na(x)) | id_change) * .I)]), .SDcols=constant_vars]

#if more than one row in per person in same month keep latest date
# second_prev[, month_year := format(as.Date(prop_dates), "%Y-%m")]
setorder(second_prev, patid, -prop_dates, -issuedate)
second_prev[, dupl := seq_len(.N), by = list(patid, prop_dates)]
second_prev <- second_prev[dupl==1]
uniqueN(second_prev, by = "patid") #416,849


rm(id_change, sec_sample_dt)

#A patient will be defined as a current statin user if they have an ongoing statin prescription based on the
#number of days prescribed, allowing for a grace period of 30 days.

secondary_plot <- second_prev
secondary_plot[, grace_period := runoutdate+30] #allowing for a grace period of 30 days
secondary_plot[prop_dates>grace_period, statin_user := 0]
secondary_plot[prop_dates<=grace_period, statin_user := 1]
secondary_plot[is.na(issuedate), statin_user := 0]
uniqueN(secondary_plot, by = "patid") #416,849
#
#exclude those with first cvd event after proportion date as they would not have been eligible for secondary
#prevention at that exact point in time.
secondary_plot[(incident_date1>prop_dates), exclude := 1]
secondary_plot <- secondary_plot[is.na(exclude)]
uniqueN(secondary_plot, by = "patid") #416,849


#merge with study population dataset for grouped variables
incident_cvd_all <- read_parquet(paste0(datafiles, "cvd_incident_all.parquet"))
secondary_plot <- merge(secondary_plot, incident_cvd_all[,c("patid", "acs", "mi", "stroke", "angina", "tia", "pad", "chd", "cerebro_procs")], by="patid", all.x=TRUE)
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))
secondary_plot <- merge(secondary_plot, studypop_all[,c("patid", "gender", "dob", "ethnicity", "deprivation", "start_fup", "enddate")], by="patid", all.x=TRUE)

#cvd event on and after start of follow up - historical recording of CVD
# secondary_plot <- secondary_plot[incident_date1>=start_fup]
#
secondary_plot[mi==1, CVD := "Myocardial Infarction"]
secondary_plot[stroke==1, CVD := "Stroke"]
secondary_plot[angina==1, CVD := "Angina"]
secondary_plot[tia==1, CVD := "Transient Ischeamic Attack"]
secondary_plot[pad==1, CVD := "Peripheral Arterial Disease"]
secondary_plot[acs==1 & is.na(mi) & is.na(angina), CVD := "Acute Coronary Syndrome (non-specific)"]
secondary_plot[chd==1 & is.na(mi) & is.na(angina) & is.na(acs), CVD := "Coronary Heart Disease (non-specific)"]
secondary_plot[cerebro_procs==1, CVD := "Cerebrovascular Procedures"]
secondary_plot[, .N, by=CVD]
#
#
# #####
secondary_plot[, start := pmax(start_fup,as.Date("2009-05-01"), na.rm=T)]
secondary_plot[, exit := pmin(enddate,excl_date1, na.rm=T)]

#remove rows where start_fup occurs after prop_date
secondary_plot <- secondary_plot[start<prop_dates]
secondary_plot[, .N,]
secondary_plot[, .N, by=statin_user]
#remove rows in which enddate is before prop_date
#test3 <- secondary_plot[,c("patid", "start_fup", "enddate", "prop_dates", "issuedate", "incident_date1", "statin", "statin_user")]
secondary_plot <- secondary_plot[exit>prop_dates]
secondary_plot[, .N,]
secondary_plot[, .N, by=statin_user]

 #only keep May 2009 dates as per start of study
 secondary_plot <- secondary_plot[prop_dates>= as.Date("2009-05-01")]
 secondary_plot[, .N,]
 secondary_plot[, .N, by=statin_user]
 uniqueN(secondary_plot, by = "patid") #406,746

###Save parquet
write_parquet(secondary_plot, paste0(datafiles, "sec_prev_pop.parquet"))

########################################################################################

#open parquet
secondary_plot <- read_parquet(paste0(datafiles, "sec_prev_pop.parquet"))

#generate age group
secondary_plot[, age := as.numeric(prop_dates - dob) / 365.25]

secondary_plot[, age_group := cut(age, breaks = c(25, 40, 50, 60, 70, 130),
                                  labels = c("25-39", "40-49", "50-59", "60-69", "70+"))]
secondary_plot[is.na(age_group), age_group := "70+"]

#generate denominator
secondary_plot[, total_pop := 1]
uniqueN(secondary_plot, by = "patid") #406,746

summary(secondary_plot)

ethnicity_non_miss2 <- secondary_plot[ethnicity!="missing", uniqueN(patid)] 
print(ethnicity_non_miss2)
deprivation_non_miss2 <- secondary_plot[deprivation!="missing", uniqueN(patid)]
print(deprivation_non_miss2)


uniqueN(secondary_plot, by = "patid") ##406,746

secondary_plot[, .N, by="CVD"]


##################################SECONDARY main plot#####################################
main_sec_plot <- secondary_plot[, .(total_users = sum(statin_user), 
                                    total_population = sum(total_pop)), 
                                by = .(prop_dates)]

main_sec_plot[, proportion := (total_users / total_population)*100]

setorder(main_sec_plot, "prop_dates")

overall_sec_plot <- ggplot(main_sec_plot, aes(x = prop_dates, y = proportion, group=1)) + 
  geom_line() + 
  ylab("Monthly proportion of statin users %") + 
  xlab("Month and year") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("Overall") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  ylim(20, 85) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12)) +
  theme(axis.title.y= element_blank()) +
theme(axis.title.x= element_blank()) 

#("sec_prevalence_overall.png", plot=last_plot(), path = graphs)


##################################SECONDARY CVD plot#####################################
CVD_plot <- secondary_plot[, .(total_users = sum(statin_user), 
                               total_population = sum(total_pop)), 
                           by = .(prop_dates, CVD)]

CVD_plot[, proportion := (total_users / total_population)*100]

CVD_plot_red <- CVD_plot[CVD!="Acute Coronary Syndrome (non-specific)" & CVD!="Coronary Heart Disease (non-specific)" & CVD!="Cerebrovascular Procedures"]

cvd_plot1 <- ggplot(CVD_plot_red, aes(x = prop_dates, y = proportion, col=CVD)) + 
  geom_line() +
  ylab("Monthly proportion of statin users %") + 
  xlab("Calendar Month") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12)) +
  ggtitle("CVD subtype") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  ylim(20, 85) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
#ggsave("sec_prevalence_CVD.png", plot=last_plot(), path = graphs)


##################################secondary Age-group plot#####################################
age_cvd_plot <- secondary_plot[, .(total_users = sum(statin_user), 
                             total_population = sum(total_pop)), 
                         by = .(prop_dates, age_group)]

age_cvd_plot[, proportion := (total_users / total_population)*100]
setorder(age_cvd_plot, "prop_dates")

age_cvd_plot1 <- ggplot(age_cvd_plot, aes(x = prop_dates, y = proportion, col=age_group)) + 
  geom_line() +
  ylab("Monthly proportion of statin users %") + 
  xlab("Calendar Month") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  theme(axis.text.x = element_text(size=12, angle = 45)) +
  labs(color = "Age group") +
  ggtitle("Age group") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  ylim(20, 85) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12)) +
  theme(axis.title.x= element_blank())  +
  theme(axis.title.y= element_blank())  +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
  
#ggsave("sec_prevalence_age.png", plot=last_plot(), path = graphs)

##################################secondary gender plot#####################################

gender_cvd_plot <- secondary_plot[, .(total_users = sum(statin_user), 
                                total_population = sum(total_pop)), 
                            by = .(prop_dates, gender)]

gender_cvd_plot[, proportion := (total_users / total_population)*100]
gender_cvd_plot <- gender_cvd_plot[proportion<100]

gender_cvd_plot[, gender := ifelse(gender == "Male", "Men", ifelse(gender == "Female", "Women", gender))]

gender_cvd_plot1 <- ggplot(gender_cvd_plot, aes(x = prop_dates, y = proportion, col=gender)) + 
  geom_line() + 
  ylab("Monthly proportion of statin users %") + 
  xlab("Calendar Month") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  labs(color = "Gender") +
  ggtitle("Gender") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  ylim(20, 85) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12)) +
  theme(axis.title.y= element_blank()) +
  theme(axis.title.x= element_blank()) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )

#ggsave("sec_prevalence_gender.png", plot=last_plot(), path = graphs)

##################################secondary ethnicity plot#####################################

ethnicity_cvd_plot <- secondary_plot[, .(total_users = sum(statin_user), 
                                   total_population = sum(total_pop)), 
                               by = .(prop_dates, ethnicity)]

ethnicity_cvd_plot[, proportion := (total_users / total_population)*100]
ethnicity_cvd_plot <- ethnicity_cvd_plot[proportion<100]

setnames(ethnicity_cvd_plot, "ethnicity", "Ethnicity")

ethnicity_cvd_plot <- ethnicity_cvd_plot[Ethnicity!="missing"]

eth_cvd_plot1 <- ggplot(ethnicity_cvd_plot, aes(x = prop_dates, y = proportion, col=Ethnicity)) + 
  geom_line() + 
  ylab("Monthly proportion of statin users %") + 
  xlab("Calendar Month") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  labs(color = "Ethnicity") +
  ggtitle("Ethnicity") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  ylim(20, 85) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12)) +
  theme(axis.title.y= element_blank()) +
  theme(axis.title.x= element_blank()) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
#ggsave("sec_prevalence_ethnicity.png", plot=last_plot(), path = graphs)

##################################secondary deprivation plot#####################################
deprivation_cvd_plot <- secondary_plot[, .(total_users = sum(statin_user), 
                                     total_population = sum(total_pop)), 
                                 by = .(prop_dates, deprivation)]

deprivation_cvd_plot[, proportion := (total_users / total_population)*100]
deprivation_cvd_plot <- deprivation_cvd_plot[proportion<100]

setnames(deprivation_cvd_plot, "deprivation", "Deprivation")

deprivation_cvd_plot <- deprivation_cvd_plot[Deprivation!="missing"]

depr_cvd_plot1 <- ggplot(deprivation_cvd_plot, aes(x = prop_dates, y = proportion, col=Deprivation)) + 
  geom_line() + 
  ylab("Monthly proportion of statin users %") + 
  xlab("Calendar Month") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  labs(color = "Deprivation") +
  ggtitle("Deprivation") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  ylim(20, 85) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12)) +
  theme(axis.title.x= element_blank()) +
  theme(axis.title.y= element_blank()) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )

#############################################################################################

############COMBINE secondary PREVENTION PLOTS

library(cowplot)

# Combine the plots into one graph with a shared x and y axis

sec_comb <- cowplot::plot_grid( overall_sec_plot,
                                          age_cvd_plot1, 
                                       gender_cvd_plot1,
                                       eth_cvd_plot1,
                                       depr_cvd_plot1,
                                cvd_plot1,
                                       ncol = 1)
sec_comb <- plot_grid(sec_comb, labels=c('B. Secondary Prevention'), label_size = 14)
print(sec_comb)
ggsave("all_sec_prevalence_plots.pdf", plot=last_plot(), width=20, height=15, units = "in", dpi = 400, path = final_graphs)

 total_comb <- plot_grid(prim_comb, sec_comb, ncol=2)
 print(total_comb)
 ggsave("prim_sec_prev.pdf", plot=last_plot(), width=20, height=25, units = "in", dpi = 400, path = final_graphs)

