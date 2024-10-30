# #--------------------------------------------------------------------------------------
# # Project: Statins and CVD prevention
# # Program Name: cr_prim_prev_subpops
# # Author: Rutendo Muzambi
# # Date version created: 05/2023
# # Description: primary prevalence proportion graphs
# #--------------------------------------------------------------------------------------
# 
# ######################################################################################################
rm(list=ls())

# rm(extract_medcode, extract_prodcode)
# # # ############################Aim 1: Prevalence of statins over time ###################################
# # 
# # ##Defining statin user
# statins_dt <- read_parquet(paste0(datafiles, "statins.parquet"))
# uniqueN(statins_dt, by = "patid") #992,248
# 
# statins_dt[, c("statin_type", "drug_type", "estnhscost", "drugsubstancename",
#                "productname", "enterdate", "pracid", "prodcodeid") := NULL]
# 
# ##import dosages
# dosage_dt <- haven::read_dta(file = stata_dosages_22)
# 
# dosage_dt <- as.data.table(dosage_dt)
# 
# statins_dt <- merge(statins_dt, dosage_dt, by ="dosageid", all.x = T)
# 
# rm(dosage_dt)
# 
# summary(statins_dt$duration)
# 
# statins_dt[, duration := as.numeric(duration)]
# 
# #Define duration of statins using quantity and daily dose
# statins_dt[duration==0 & !is.na(quantity) & !is.na(daily_dose) & quantity!=0 & daily_dose!=0, duration := quantity/daily_dose]
# statins_dt[duration==0 & (daily_dose==0 | is.na(daily_dose)) & quantity!=0 & !is.na(quantity), duration := quantity]
# statins_dt[duration<=1, duration := 28]
# 
# statins_dt[, .N, by=duration]
# 
# #only keep the prescription with the longest duration for a given date if multiple Rx on same date for the same patient or runoutdate will not be created
# setorder(statins_dt, patid, issuedate, -duration)
# statins_dt[, dupl := rowid(patid, issuedate)]
# statins_dt[, .N, by=dupl]
# statins_dt <- statins_dt[dupl==1]
# 
# #remove prescriptions in period not interested in
# statins_dt[,runoutdate := issuedate+duration]
# statins_dt[, year := year(issuedate)]
# statins_dt <- statins_dt[year>="2008"] ##study starts in 2009
# 
# #if duration more than 180 days change to the median (28 days)
# statins_dt[duration>180, duration := 28]
# 
# statins_dt <- statins_dt[,c("patid", "issuedate", "duration", "runoutdate")]
# 
# ###########MERGE with patient dataset
# patients_dt <- read_parquet(paste0(datafiles, "clean_pt.prac.parquet"))
# patients_dt <- patients_dt[,c("patid", "start_fup", "enddate")]
# 
# uniqueN(patients_dt, by = "patid") #5,000,000
# 
# ##################################################################################
# preval_dt <- merge(patients_dt,statins_dt, by="patid", all.x=TRUE)
# 
# rm(patients_dt)
# 
# ##Read and merge statin exclusions - cvd event or statin contraindications
# excl_prev_dt <- read_parquet(paste0(datafiles, "excl_prev_all.parquet"))
# uniqueN(excl_prev_dt, by = "patid") #563,550
# excl_prev_dt <- excl_prev_dt[,c("patid", "excl_date1")]
# 
# primary_dt <- merge(preval_dt,excl_prev_dt, by="patid", all.x=TRUE)
# uniqueN(primary_dt, by = "patid") #5,000,000
# primary_dt <- primary_dt[,c("patid", "runoutdate", "duration", "issuedate", "excl_date1", "start_fup", "enddate")]
# 
# rm(excl_prev_dt)
# #######primary prevention eligibility
# #primary_dt[excl_date1 <= issuedate & !is.na(issuedate) & !is.na(excl_date1), exclude := 1]
# primary_dt[excl_date1 <= start_fup & !is.na(excl_date1), exclude := 1]
# primary_dt <- primary_dt[is.na(exclude)]
# uniqueN(primary_dt, by = "patid") #4,658,400
# primary_dt[, primary_prev := 1]
# 
# #statins
# primary_dt[!is.na(issuedate), statins := 1]
# primary_dt[is.na(issuedate), statins := 0]
# 
# #create matching prop dates column
# # primary_dt[, prop_dates := issuedate]
# primary_dt[, prop_dates := ceiling_date(issuedate, "month") -1]
# 
# primary_dt <- primary_dt[,c("patid", "statins", "issuedate", "duration",
#                               "runoutdate", "excl_date1", "primary_prev",
#                               "prop_dates", "start_fup", "enddate")]
# setorder(primary_dt, patid, prop_dates)
# 
# #create data table with month and years of study period using patients eligible for primary prevention
# month_year_dt <- primary_dt[,c("patid", "start_fup", "enddate")]
# month_year_dt <- unique(month_year_dt) #delete repeating patids
# uniqueN(month_year_dt, by = "patid") #4,658,400
# 
# month_year_dt <- month_year_dt[rep(seq_len(nrow(month_year_dt)), each = 153), ]
# prop_dates <- seq(as.Date("2009/05/01"), as.Date("2022/01/01"), "month", format="%d/%m/%Y") -1
# month_year_dt <- cbind(month_year_dt, prop_dates)
# setorder(month_year_dt, patid, prop_dates)
# month_year_dt[, .N] #136287810
# #remove irrelevant rows
# #remove where start_fup is after prop_date
# month_year_dt <- month_year_dt[start_fup<prop_dates]
# #remove rows in which enddate is before prop_date
# month_year_dt <- month_year_dt[enddate>prop_dates]
# month_year_dt[, .N] #367334405
# 
# #append proportion dates and primary prevention dataset
# prim_prev <- rbind(primary_dt,month_year_dt, fill=T)
# uniqueN(prim_prev, by = "patid") #4,658,400
# 
# rm(month_year_dt)
# 
# constant_vars <- c("issuedate", "duration",
#                    "runoutdate", "excl_date1", "primary_prev")
# 
# #fill in constant variables
# setorder(prim_prev, patid, prop_dates)
# id_change = prim_prev[, c(TRUE, patid[-1] != patid[-.N])]
# prim_prev[, (constant_vars) := lapply(.SD, function(x) x[cummax(((!is.na(x)) | id_change) * .I)]), .SDcols=constant_vars]
# 
# rm(id_change)
# 
# #if more than one row in per person in same month keep latest date
# setorder(prim_prev, patid, -prop_dates, -issuedate)
# prim_prev[, dupl := seq_len(.N), by = list(patid, prop_dates)]
# prim_prev <- prim_prev[dupl==1]
# 
# #A patient will be defined as a current statin user if they have an ongoing statin prescription based on the number of days
# #prescribed, allowing for a grace period of 30 days.
# prim_prev <- prim_prev[,c("patid", "prop_dates", "issuedate", "duration", "runoutdate", "excl_date1", "start_fup", "enddate")]
# prim_prev[, grace_period := runoutdate+30] #allowing for a grace period of 30 days
# prim_prev[prop_dates>grace_period, statin_user := 0]
# prim_prev[prop_dates<=grace_period, statin_user := 1]
# prim_prev[is.na(issuedate), statin_user := 0]
# 
# uniqueN(prim_prev, by = "patid") #4,658,400
# rm(preval_dt, statins_dt, primary_dt)
# #merge with study population dataset for grouped variables
# studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))
# 
# prim_prev <- merge(prim_prev, studypop_all[,c("patid", "gender", "dob", "ethnicity", "deprivation")], by="patid")
# uniqueN(prim_prev, by = "patid") #4,658,400
# 
# rm(studypop_all)
# 
# prim_prev[, prop_dates := floor_date(prop_dates, "day")]
# #####
# prim_prev[, start := pmax(start_fup,as.Date("2009-05-01"), na.rm=T)]
# prim_prev[, exit := pmin(enddate,excl_date1, na.rm=T)]
# 
# #remove where start_fup is after prop_date
# prim_prev[, .N]
# prim_prev <- prim_prev[start<prop_dates]
# uniqueN(prim_prev, by = "patid") #4,600,119
# 
# #remove rows in which enddate is before prop_date
# prim_prev <- prim_prev[exit>prop_dates]
# uniqueN(prim_prev, by = "patid") #4,597,500
# 
# prim_prev[, total_pop := 1]
# 
# #generate age group
# prim_prev[, age := round(as.numeric(difftime(prop_dates, dob, units = "days")) / 365.25)]
# 
# prim_prev[, age_group := cut(age, breaks = c(25, 40, 50, 60, 70, 80, Inf),
#                               labels = c("25-39", "40-49", "50-59", "60-69", "70-79", "80+"))]
# prim_prev[is.na(age_group), age_group := "80+"]
# prim_prev[, .N, age_group]
# 
# prim_prev <- prim_prev[, c("patid", "statin_user", "total_pop", "prop_dates", "ethnicity", "deprivation", "gender", "age_group")]
# 
# uniqueN(prim_prev, by = "patid")  #4,597,500
# ###Save parquet
# write_parquet(prim_prev, paste0(datafiles, "prim_prev_all.parquet"))

library(openxlsx)

###Open parquet
prim_prev <- read_parquet(paste0(datafiles, "prim_prev_all.parquet"))

#remove unnecessary datasets
#rm(primary_dt, primary_prev, studypop_all)
##################################PRIMARY main plot#####################################
uniqueN(prim_prev, by = "patid") #4,597,500

main_plot <- prim_prev[, .(total_users = sum(statin_user), 
                              total_population = sum(total_pop)), 
                          by = .(prop_dates)]

main_plot[, proportion := (total_users / total_population)*100]
main_plot <- main_plot[proportion<100]

setorder(main_plot, "prop_dates")

overall_plot <- ggplot(main_plot, aes(x = prop_dates, y = proportion, group=1)) + 
  geom_line() + 
  ggtitle("Overall") +
  ylab("Monthly proportion of statin users %") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlab("Calendar Month") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  ylim(0, 40) +
  theme(axis.title.x= element_blank())  +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
#ggsave("prim_prevalence_overall.png", plot=last_plot(), path = graphs)

#load into excel 
wb <- loadWorkbook(paste0(excel, "proportions.xlsx"))
names(wb)
writeData(wb, "primary_prevention", main_plot, startCol=1, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

rm(main_plot)
##################################PRIMARY Age-group plot#####################################
prim_prev[, age_cat := age_group]
prim_prev[age_group=="70-79"|age_group=="80+", age_cat := "70+"]
prim_prev[, .N, by=c("age_group", "age_cat")]
prim_prev[, .N, by=c("age_cat")]

age_plot <- prim_prev[, .(total_users = sum(statin_user), 
                             total_population = sum(total_pop)), 
                         by = .(prop_dates, age_cat)]

age_plot[, proportion := (total_users / total_population)*100]

age_plot1 <- ggplot(age_plot, aes(x = prop_dates, y = proportion, col=age_cat)) + 
  geom_line() +
  ylab("Monthly proportion of statin users %") + 
  xlab("Calendar Month") + 
  labs(color = "Age group") +
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("Age group") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  theme(axis.title.x= element_blank()) +
  ylim(0, 40) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
#ggsave("prim_prevalence_age.png", plot=last_plot(), path = graphs)

writeData(wb, "primary_prevention", age_plot, startCol=7, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

rm(age_plot)
##################################PRIMARY gender plot#####################################

gender_plot <- prim_prev[, .(total_users = sum(statin_user), 
                                total_population = sum(total_pop)), 
                            by = .(prop_dates, gender)]

gender_plot[, proportion := (total_users / total_population)*100]
gender_plot <- gender_plot[proportion<100]

gender_plot[, gender := ifelse(gender == "Male", "Men", ifelse(gender == "Female", "Women", gender))]

gender_plot1 <- ggplot(gender_plot, aes(x = prop_dates, y = proportion, col=gender)) + 
  geom_line() + 
  ylab("Monthly proportion of statin users %") + 
  xlab("Calendar Month") + 
  labs(color = "Gender") +
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("Gender") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  theme(axis.title.x= element_blank()) +
  ylim(0, 40) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )

#("prim_prevalence_gender.png", plot=last_plot(), path = graphs)

writeData(wb, "primary_prevention", gender_plot, startCol=14, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

rm(gender_plot)
##################################PRIMARY ethnicity plot#####################################

ethnicity_plot <- prim_prev[, .(total_users = sum(statin_user), 
                                   total_population = sum(total_pop)), 
                               by = .(prop_dates, ethnicity)]

ethnicity_plot[, proportion := (total_users / total_population)*100]
ethnicity_plot <- ethnicity_plot[proportion<100]

setnames(ethnicity_plot, "ethnicity", "Ethnicity")

ethnicity_plot <- ethnicity_plot[Ethnicity!="missing"]

eth_plot1 <- ggplot(ethnicity_plot, aes(x = prop_dates, y = proportion, col=Ethnicity)) + 
  geom_line() + 
  ylab("Monthly proportion of statin users %") + 
  xlab("Calendar Month") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("Ethnicity") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(axis.text.y = element_text(size = 9)) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  theme(axis.title.x= element_blank()) +
  ylim(0, 40) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )

#ggsave("prim_prevalence_ethnicity.png", plot=last_plot(), path = graphs)

writeData(wb, "primary_prevention", ethnicity_plot, startCol=21, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

rm(ethnicity_plot)
##################################PRIMARY deprivation plot#####################################


prim_prev <- prim_prev[, c("patid") := NULL]
deprivation_plot <- prim_prev[, .(total_users = sum(statin_user), 
                                     total_population = sum(total_pop)), 
                                 by = .(prop_dates, deprivation)]

deprivation_plot[, proportion := (total_users / total_population)*100]
deprivation_plot <- deprivation_plot[proportion<100]

setnames(deprivation_plot, "deprivation", "Deprivation")

class(deprivation_plot$Deprivation)

deprivation_plot <- deprivation_plot[Deprivation!="missing"]


depr_plot1 <- ggplot(deprivation_plot, aes(x = prop_dates, y = proportion, col=Deprivation)) + 
  geom_line() + 
  ylab("Monthly proportion of statin users %") + 
  xlab("Calendar Month") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("Deprivation") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12)) +
  ylim(0, 40) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )

#ggsave("prim_prevalence_depriv.png", plot=last_plot(), path = graphs)

writeData(wb, "primary_prevention", deprivation_plot, startCol=28, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

#############################################################################################

############COMBINE PRIMARY PREVENTION PLOTS

library(cowplot)

# Combine the plots into one graph with a shared x and y axis

prim_comb <- cowplot::plot_grid(overall_plot,
                                 age_plot1, 
                                       gender_plot1,
                                       eth_plot1,
                                       depr_plot1,
                                NULL,
                                       ncol = 1)
print(prim_comb)

prim_comb <- plot_grid(prim_comb, labels=c('A. Primary Prevention'), label_size = 14)

ggsave("all_prim_prev_plots.pdf", plot=last_plot(), width=10, height=30, units = "in", dpi = 400, path = final_graphs)

rm(age_plot1, deprivation_plot, depr_plot1, eth_plot1, gender_plot1, overall_plot, prim_prev)



# test <- age_plot[age_cat=="70+"]
# 
# setorder(test, "prop_dates")



