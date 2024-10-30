#---------------------------------------------------------------------------------------
# Project: Statins and CVD prevention
# Program Name: an_risk_prev_plots
# Author: Rutendo Muzambi
# Date version created: 20/02/2023
# Description: recorded risk score prevalence plots
#--------------------------------------------------------------------------------------
rm(list=ls())

####Append risk  plots
cvdrisk_pop_1 <- read_parquet(paste0(datafiles, "cvdrisk_pop_1.parquet"))
uniqueN(cvdrisk_pop_1, by = "patid") #1,460,551

#pop1 <- cvdrisk_pop_1[, c("patid")]

cvdrisk_pop_2 <- read_parquet(paste0(datafiles, "cvdrisk_pop_2.parquet"))
uniqueN(cvdrisk_pop_2, by = "patid") #1,458,261

#pop2 <- cvdrisk_pop_2[, c("patid")]

cvdrisk_pop_3 <- read_parquet(paste0(datafiles, "cvdrisk_pop_3.parquet"))
uniqueN(cvdrisk_pop_3, by = "patid") #1,460,970

#pop3 <- cvdrisk_pop_3[, c("patid")]

cvdrisk_pop_4 <- read_parquet(paste0(datafiles, "cvdrisk_pop_4.parquet"))
uniqueN(cvdrisk_pop_4, by = "patid") #1,459,167

#pop4 <- cvdrisk_pop_4[, c("patid")]

#pop_all <- rbind(pop1, pop2, pop3, pop4)
#uniqueN(pop_all, by = "patid") #4,470,870

#write_parquet(pop_all, paste0(datafiles, "cvdrisk_prevpop_all.parquet"))
#rm(pop1,pop2,pop3,pop4)
#rm(pop1,pop2,pop3,pop4, pop_all)
#Overall population

# Create a list of data tables
dt_list <- list(cvdrisk_pop_1, cvdrisk_pop_2, cvdrisk_pop_3, cvdrisk_pop_4)

# Create a list to store the results
mainrisk_list <- list()

# Loop through the data tables
for (i in 1:length(dt_list)) {
  # Calculate 'users' and 'population' for the current data table
  result <- dt_list[[i]][, .(users = sum(risk_record), population = sum(total_pop)), by = .(prop_dates)]
  
  # Store the result in the result list
  mainrisk_list[[i]] <- result
}
list2env(setNames(mainrisk_list, paste0("mainrisk_", 1:4)), envir = .GlobalEnv)

rm(result, dt_list, mainrisk_list)

mainrisk_final <- merge(mainrisk_1, mainrisk_2, by="prop_dates")

mainrisk_final[, users := users.x + users.y]
mainrisk_final[, population := population.x + population.y]
mainrisk_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

mainrisk_final <- merge(mainrisk_final, mainrisk_3, by="prop_dates")

mainrisk_final[, users := users.x + users.y]
mainrisk_final[, population := population.x + population.y]
mainrisk_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

mainrisk_final <- merge(mainrisk_final, mainrisk_4, by="prop_dates")
mainrisk_final[, users := users.x + users.y]
mainrisk_final[, population := population.x + population.y]
mainrisk_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]
mainrisk_final[, proportion := (users / population)*100]

rm(mainrisk_1, mainrisk_2, mainrisk_3, mainrisk_4)

#uniqueN(mainrisk_final, by = "patid") #1,811,847



###########################################################################################
mainrisk_plot <- ggplot(mainrisk_final, aes(x = prop_dates, y = proportion, group=1)) + 
  geom_line() + 
  ylab("Proportion with CVD risk score %") + 
  xlab("Calendar Month") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-11-30"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("Overall") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank()) +
  ylim(0, 100) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
#ggsave("main_risk_scoring.png", plot=last_plot(), path = graphs)

library(openxlsx)

#load into excel 
wb <- loadWorkbook(paste0(excel, "proportions.xlsx"))
names(wb)
writeData(wb, "cvd_risk", mainrisk_final, startCol=1, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

rm(mainrisk_final)
##################################PRIMARY risk type plot#####################################
##subpop 1
risktype_1 <- cvdrisk_pop_1[risk_record==1 & (risk_tool=="qrisk"|risk_tool=="framingham"|risk_tool=="jbs")]

risktool_1 <- risktype_1[, .(users = sum(risk_record)), 
                            by = .(prop_dates, risk_tool)]

risktool_2 <- risktype_1[, .(population = sum(total_pop)), 
                            by = .(prop_dates)]

setorder(risktool_1, prop_dates)
setorder(risktool_2, prop_dates)

risktype_1 <- merge(risktool_1, risktool_2, by="prop_dates")

##subpop 2
risktype_2 <- cvdrisk_pop_2[risk_record==1 & (risk_tool=="qrisk"|risk_tool=="framingham"|risk_tool=="jbs")]

risktool_1 <- risktype_2[, .(users = sum(risk_record)), 
                            by = .(prop_dates, risk_tool)]

risktool_2 <- risktype_2[, .(population = sum(total_pop)), 
                            by = .(prop_dates)]

setorder(risktool_1, prop_dates)
setorder(risktool_2, prop_dates)

risktype_2 <- merge(risktool_1, risktool_2, by="prop_dates")

##subpop 3
risktype_3 <- cvdrisk_pop_3[risk_record==1 & (risk_tool=="qrisk"|risk_tool=="framingham"|risk_tool=="jbs")]

risktool_1 <- risktype_3[, .(users = sum(risk_record)), 
                         by = .(prop_dates, risk_tool)]

risktool_2 <- risktype_3[, .(population = sum(total_pop)), 
                         by = .(prop_dates)]

setorder(risktool_1, prop_dates)
setorder(risktool_2, prop_dates)

risktype_3 <- merge(risktool_1, risktool_2, by="prop_dates")

##subpop 4
risktype_4 <- cvdrisk_pop_4[risk_record==1 & (risk_tool=="qrisk"|risk_tool=="framingham"|risk_tool=="jbs")]

risktool_1 <- risktype_4[, .(users = sum(risk_record)), 
                         by = .(prop_dates, risk_tool)]

risktool_2 <- risktype_4[, .(population = sum(total_pop)), 
                         by = .(prop_dates)]

setorder(risktool_1, prop_dates)
setorder(risktool_2, prop_dates)

risktype_4 <- merge(risktool_1, risktool_2, by="prop_dates")

rm(risktool_1, risktool_2)
################subpop 4

risktool_final <- merge(risktype_1, risktype_2, by=c("prop_dates", "risk_tool"), all = T)

risktool_final[, users := fcoalesce(users.x, 0) + fcoalesce(users.y, 0)]
risktool_final[, population := fcoalesce(population.x, 0) + fcoalesce(population.y, 0)]
risktool_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

risktool_final <- merge(risktool_final, risktype_3, by=c("prop_dates", "risk_tool"), all = T)

risktool_final[, users := fcoalesce(users.x, 0) + fcoalesce(users.y, 0)]
risktool_final[, population := fcoalesce(population.x, 0) + fcoalesce(population.y, 0)]
risktool_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

risktool_final <- merge(risktool_final, risktype_4, by=c("prop_dates", "risk_tool"), all = T)

risktool_final[, users := fcoalesce(users.x, 0) + fcoalesce(users.y, 0)]
risktool_final[, population := fcoalesce(population.x, 0) + fcoalesce(population.y, 0)]
risktool_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

risktool_final[, proportion := (users / population)*100]

risktool_final[risk_tool=="jbs", risk_tool := "JBS"]
risktool_final[risk_tool=="framingham", risk_tool := "Framingham"]
risktool_final[risk_tool=="qrisk", risk_tool := "QRISK"]

risktool_plot <- ggplot(risktool_final, aes(x = prop_dates, y = proportion, col=risk_tool)) + 
  geom_line() + 
  ylab("Proportion with CVD risk score %") +  
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("CVD risk tools") +
  labs(color = "Risk Tool") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y= element_blank()) +
  ylim(0, 100) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
#ggsave("risk_scoring_tool.png", plot=last_plot(), path = graphs)

writeData(wb, "cvd_risk", risktool_final, startCol=7, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

rm(risktool_final, risktype_1, risktype_2, risktype_3, risktype3, risktype_4)
##################################risk scoring category plot#####################################

##subpop 1
riskcat_pop1 <- cvdrisk_pop_1[risk_record==1 & !is.na(risk_cat)]

riskcat_1 <-riskcat_pop1[, .(users = sum(risk_record)), 
                          by = .(prop_dates, risk_cat)]

riskcat_2 <- riskcat_pop1[, .(population = sum(total_pop)), 
                           by = .(prop_dates)]

setorder(riskcat_1, prop_dates)
setorder(riskcat_2, prop_dates)

riskcat_pop1 <- merge(riskcat_1, riskcat_2, by="prop_dates")

##subpop 2
riskcat_pop2 <- cvdrisk_pop_2[risk_record==1 & !is.na(risk_cat)]

riskcat_1 <-riskcat_pop2[, .(users = sum(risk_record)), 
                          by = .(prop_dates, risk_cat)]

riskcat_2 <- riskcat_pop2[, .(population = sum(total_pop)), 
                          by = .(prop_dates)]

setorder(riskcat_1, prop_dates)
setorder(riskcat_2, prop_dates)

riskcat_pop2 <- merge(riskcat_1, riskcat_2, by="prop_dates")

##subpop 3
riskcat_pop3 <- cvdrisk_pop_3[risk_record==1 & !is.na(risk_cat)]

riskcat_1 <- riskcat_pop3[, .(users = sum(risk_record)), 
                          by = .(prop_dates, risk_cat)]

riskcat_2 <- riskcat_pop3[, .(population = sum(total_pop)), 
                           by = .(prop_dates)]

setorder(riskcat_1, prop_dates)
setorder(riskcat_2, prop_dates)

riskcat_pop3 <- merge(riskcat_1, riskcat_2, by="prop_dates")

##subpop 4
riskcat_pop4 <- cvdrisk_pop_4[risk_record==1 & !is.na(risk_cat)]

riskcat_1 <-riskcat_pop4[, .(users = sum(risk_record)), 
                          by = .(prop_dates, risk_cat)]

riskcat_2 <- riskcat_pop4[, .(population = sum(total_pop)), 
                           by = .(prop_dates)]

setorder(riskcat_1, prop_dates)
setorder(riskcat_2, prop_dates)

riskcat_pop4 <- merge(riskcat_1, riskcat_2, by="prop_dates")
rm(riskcat_1, riskcat_2)
################

riskcat_final <- merge(riskcat_pop1, riskcat_pop2, by=c("prop_dates", "risk_cat"), all = T)

riskcat_final[, users := fcoalesce(users.x, 0) + fcoalesce(users.y, 0)]
riskcat_final[, population := fcoalesce(population.x, 0) + fcoalesce(population.y, 0)]
riskcat_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

riskcat_final <- merge(riskcat_final, riskcat_pop3, by=c("prop_dates", "risk_cat"), all = T)

riskcat_final[, users := fcoalesce(users.x, 0) + fcoalesce(users.y, 0)]
riskcat_final[, population := fcoalesce(population.x, 0) + fcoalesce(population.y, 0)]
riskcat_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

riskcat_final <- merge(riskcat_final, riskcat_pop4, by=c("prop_dates", "risk_cat"), all = T)

riskcat_final[, total_users := fcoalesce(users.x, 0) + fcoalesce(users.y, 0)]
riskcat_final[, total_population := fcoalesce(population.x, 0) + fcoalesce(population.y, 0)]
riskcat_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

riskcat_final[, proportion := (total_users / total_population)*100]

riskcat_plot <- ggplot(riskcat_final, aes(x = prop_dates, y = proportion, col=risk_cat)) + 
  geom_line() + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("CVD Risk Category") +
  labs(color = "Risk Category") +
  ylab("Proportion with CVD risk score %") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank()) +
  ylim(0, 100) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
#ggsave("risk_scoring_cate.png", plot=last_plot(), path = graphs)

writeData(wb, "cvd_risk", riskcat_final, startCol=14, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

rm(riskcat_final, riskcat_pop1, riskcat_pop2, riskcat_pop3, riskcat_pop4)

##################################PRIMARY Age-group plot#####################################

cvdrisk_pop_1 <- cvdrisk_pop_1[is.na(age_group), age_group := "25-39"]
cvdrisk_pop_2 <- cvdrisk_pop_2[is.na(age_group), age_group := "25-39"]
cvdrisk_pop_3 <- cvdrisk_pop_3[is.na(age_group), age_group := "25-39"]
cvdrisk_pop_4 <- cvdrisk_pop_4[is.na(age_group), age_group := "25-39"]

age_1 <- cvdrisk_pop_1[, .(users = sum(risk_record), 
                           population = sum(total_pop)), 
                       by = .(prop_dates, age_group)]

age_2 <- cvdrisk_pop_2[, .(users = sum(risk_record), 
                           population = sum(total_pop)), 
                       by = .(prop_dates, age_group)]

age_3 <- cvdrisk_pop_3[, .(users = sum(risk_record), 
                           population = sum(total_pop)), 
                       by = .(prop_dates, age_group)]

age_4 <- cvdrisk_pop_4[, .(users = sum(risk_record), 
                           population = sum(total_pop)), 
                       by = .(prop_dates, age_group)]

age_final <- merge(age_1, age_2, by=c("prop_dates", "age_group"))

age_final[, users := users.x + users.y]
age_final[, population := population.x + population.y]
age_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

age_final <- merge(age_final, age_3, by=c("prop_dates", "age_group"))

age_final[, users := users.x + users.y]
age_final[, population := population.x + population.y]
age_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

age_final <- merge(age_final, age_4, by=c("prop_dates", "age_group"))
age_final[, total_users := users.x + users.y]
age_final[, total_population := population.x + population.y]
age_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

age_final[, proportion := (total_users / total_population)*100]

rm(age_1, age_2, age_3, age_4)

age_final[, .N, by="age_group"]

age_plot <- ggplot(age_final, aes(x = prop_dates, y = proportion, col=age_group)) + 
  geom_line() + 
  ylab("Proportion with CVD risk score %") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("Age group") +
  labs(color = "Age Group") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y= element_blank()) +
  ylim(0, 100) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
#ggsave("riskscoring_age.png", plot=last_plot(), path = graphs)

writeData(wb, "cvd_risk", age_final, startCol=21, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

rm(age_final)
##################################risk scoring gender plot#####################################

gender_1 <- cvdrisk_pop_1[, .(users = sum(risk_record), 
                                population = sum(total_pop)), 
                            by = .(prop_dates, gender)]

gender_2 <- cvdrisk_pop_2[, .(users = sum(risk_record), 
                                population = sum(total_pop)), 
                            by = .(prop_dates, gender)]

gender_3 <- cvdrisk_pop_3[, .(users = sum(risk_record), 
                                population = sum(total_pop)), 
                            by = .(prop_dates, gender)]

gender_4 <- cvdrisk_pop_4[, .(users = sum(risk_record), 
                                population = sum(total_pop)), 
                            by = .(prop_dates, gender)]

gender_final <- merge(gender_1, gender_2, by=c("prop_dates", "gender"))

gender_final[, users := users.x + users.y]
gender_final[, population := population.x + population.y]
gender_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

gender_final <- merge(gender_final, gender_3, by=c("prop_dates", "gender"))

gender_final[, users := users.x + users.y]
gender_final[, population := population.x + population.y]
gender_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

gender_final <- merge(gender_final, gender_4, by=c("prop_dates", "gender"))
gender_final[, total_users := users.x + users.y]
gender_final[, total_population := population.x + population.y]
gender_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

gender_final[, proportion := (total_users / total_population)*100]

setnames(gender_final, "gender", "Gender")

gender_final[Gender=="Male", Gender := "Men"]
gender_final[Gender=="Female", Gender := "Women"]

rm(gender_1, gender_2, gender_3, gender_4)

gender_plot <- ggplot(gender_final, aes(x = prop_dates, y = proportion, col=Gender)) + 
  geom_line() + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("Gender") +
  ylab("Proportion with CVD risk score %") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank()) +
  ylim(0, 100) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
#ggsave("risk_scoring_sex.png", plot=last_plot(), path = graphs)

writeData(wb, "cvd_risk", gender_final, startCol=28, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

rm(gender_final)

##################################PRIMARY ethnicity plot#####################################
eth_1 <- cvdrisk_pop_1[!ethnicity!="missing" & !is.na(ethnicity)]
eth_1 <- cvdrisk_pop_1[, .(users = sum(risk_record), 
                              population = sum(total_pop)), 
                          by = .(prop_dates, ethnicity)]

eth_2 <- cvdrisk_pop_2[!ethnicity!="missing" & !is.na(ethnicity)]
eth_2 <- cvdrisk_pop_2[, .(users = sum(risk_record), 
                              population = sum(total_pop)), 
                          by = .(prop_dates, ethnicity)]

eth_3 <- cvdrisk_pop_3[!ethnicity!="missing" & !is.na(ethnicity)]
eth_3 <- cvdrisk_pop_3[, .(users = sum(risk_record), 
                              population = sum(total_pop)), 
                          by = .(prop_dates, ethnicity)]

eth_4 <- cvdrisk_pop_4[!ethnicity!="missing" & !is.na(ethnicity)]
eth_4 <- cvdrisk_pop_4[, .(users = sum(risk_record), 
                              population = sum(total_pop)), 
                          by = .(prop_dates, ethnicity)]

ethnicity_final <- merge(eth_1, eth_2, by=c("prop_dates", "ethnicity"))

ethnicity_final[, users := users.x + users.y]
ethnicity_final[, population := population.x + population.y]
ethnicity_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

ethnicity_final <- merge(ethnicity_final, eth_3, by=c("prop_dates", "ethnicity"))

ethnicity_final[, users := users.x + users.y]
ethnicity_final[, population := population.x + population.y]
ethnicity_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

ethnicity_final <- merge(ethnicity_final, eth_4, by=c("prop_dates", "ethnicity"))
ethnicity_final[, total_users := users.x + users.y]
ethnicity_final[, total_population := population.x + population.y]
ethnicity_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

ethnicity_final[, proportion := (total_users / total_population)*100]

setnames(ethnicity_final, "ethnicity", "Ethnicity")

rm(eth_1, eth_2, eth_3, eth_4)

ethnicity_final <- ethnicity_final[Ethnicity!="missing" & !is.na(Ethnicity)]

eth_plot <- ggplot(ethnicity_final, aes(x = prop_dates, y = proportion, col=Ethnicity)) + 
  geom_line() + 
  ylab("Proportion with CVD risk score %") +  
  xlab("Calendar Month") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("Ethnicity") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y= element_blank()) +
  ylim(0, 100) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
#ggsave("riskscoring_ethnicity.png", plot=last_plot(), path = graphs)

writeData(wb, "cvd_risk", ethnicity_final, startCol=35, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

##################################PRIMARY deprivation plot#####################################
depr_1 <- cvdrisk_pop_1[deprivation!="missing" & !is.na(deprivation)]
depr_1 <- cvdrisk_pop_1[, .(users = sum(risk_record), 
                            population = sum(total_pop)), 
                        by = .(prop_dates, deprivation)]

depr_2 <- cvdrisk_pop_2[deprivation!="missing" & !is.na(deprivation)]
depr_2 <- cvdrisk_pop_2[, .(users = sum(risk_record), 
                            population = sum(total_pop)), 
                        by = .(prop_dates, deprivation)]

depr_3 <- cvdrisk_pop_3[deprivation!="missing" & !is.na(deprivation)]
depr_3 <- cvdrisk_pop_3[, .(users = sum(risk_record), 
                            population = sum(total_pop)), 
                        by = .(prop_dates, deprivation)]

depr_4 <- cvdrisk_pop_4[deprivation!="missing" & !is.na(deprivation)]
depr_4 <- cvdrisk_pop_4[, .(users = sum(risk_record), 
                            population = sum(total_pop)), 
                        by = .(prop_dates, deprivation)]

depr_final <- merge(depr_1, depr_2, by=c("prop_dates", "deprivation"))

depr_final[, users := users.x + users.y]
depr_final[, population := population.x + population.y]
depr_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

depr_final <- merge(depr_final, depr_3, by=c("prop_dates", "deprivation"))

depr_final[, users := users.x + users.y]
depr_final[, population := population.x + population.y]
depr_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

depr_final <- merge(depr_final, depr_4, by=c("prop_dates", "deprivation"))
depr_final[, total_users := users.x + users.y]
depr_final[, total_population := population.x + population.y]
depr_final[, c("users.x", "users.y", "population.x", "population.y") := NULL]

depr_final[, proportion := (total_users / total_population)*100]

setnames(depr_final, "deprivation", "Deprivation")

depr_final <- depr_final[Deprivation!="missing" & !is.na(Deprivation)]

depr_plot <- ggplot(depr_final, aes(x = prop_dates, y = proportion, col=Deprivation)) + 
  geom_line() + 
  xlab("Calendar Month") + 
  ylab("Proportion with CVD risk score %") + 
  scale_x_date(limits = c(as.Date("2009-05-31"), c(as.Date("2021-12-31"))), 
               breaks = date_breaks("1 year"), labels = date_format("%b %Y")) +
  ggtitle("Deprivation") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)) +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), 
             linetype="longdash", lwd=0.5) +
  theme(axis.text.x = element_text(size=12, angle = 45),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12)) +
  ylim(0, 100) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal", 
    legend.text = element_text(size = 12), legend.title = element_text(size=12)
  )
#ggsave("riskscoring_depriv.png", plot=last_plot(), path = graphs)

writeData(wb, "cvd_risk", depr_final, startCol=42, startRow=5, rowNames=T)
saveWorkbook(wb,paste0(excel, "proportions.xlsx"), overwrite=T)

rm(depr_final, depr_1, depr_2, depr_3, depr_4)

#############################################################################################

############COMBINE secondary PREVENTION PLOTS

library(cowplot)

# Combine the plots into one graph with a shared x and y axis

riskscore_plots <- cowplot::plot_grid( mainrisk_plot,
                                          risktool_plot,
                                          riskcat_plot,
                                          age_plot, 
                                          gender_plot,
                                          eth_plot,
                                          depr_plot,
                                          ncol = 2)
print(riskscore_plots)
ggsave("all_risk_score_plots.pdf", plot=last_plot(), width=15, height=15, units = "in", dpi = 400, path = final_graphs)


