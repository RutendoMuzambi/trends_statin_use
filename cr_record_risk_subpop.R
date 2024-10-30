#---------------------------------------------------------------------------------------
# Project: Statins and CVD prevention
# Author: Rutendo Muzambi
# Description: create sub-populations from recorded risk score
#--------------------------------------------------------------------------------------

rm(list=ls())

patients_dt <- read_parquet(paste0(datafiles, "clean_pt.prac.parquet"))

##Keep only patid and start of follow up
patients_dt <- patients_dt[,c("patid", "start_fup", "yob", "month")]

########RUN FROM HERE ####Monthly proportion of individuals with a CVD risk score
cvdrisk_all <- read_parquet(paste0(datafiles, "cvdrisk_all.parquet"))

uniqueN(cvdrisk_all, by = "patid") #5,000,000

setnames(cvdrisk_all, "obsdate", "risk_date")
##exclude duplicate same month and keep earliest
cvdrisk_all[, year := lubridate::year(risk_date)]
cvdrisk_all[, month := lubridate::month(risk_date)]

setorder(cvdrisk_all, patid, month, year, -risk_date)
cvdrisk_all[, dupl := seq_len(.N), by = list(patid, month, year)]
cvdrisk_all <- cvdrisk_all[dupl==1]
uniqueN(cvdrisk_all, by = "patid") #4,583,824

cvdrisk_all[, c("dupl", "exclude", "month", "year") := NULL]

# Set a random seed for reproducibility
set.seed(123)  # You can use any seed value

split_parts <- split(cvdrisk_all, sample(1:4, nrow(cvdrisk_all), replace = TRUE))

# Access the individual parts
cvdrisk_1 <- split_parts$`1`
cvdrisk_2 <- split_parts$`2`
cvdrisk_3 <- split_parts$`3`
cvdrisk_4 <- split_parts$`4`

rm(split_parts)
exclusions_all <- read_parquet(paste0(datafiles, "excl_init_all.parquet"))
uniqueN(exclusions_all, by = "patid") #73,755
exclusions_all <- exclusions_all[,c("patid", "excl_date1", "excl_flag")]

cvdrisk_1 <- merge(cvdrisk_1,exclusions_all, by="patid", all.x=TRUE)
uniqueN(cvdrisk_1, by = "patid") #5,000,000
cvdrisk_1 <- cvdrisk_1[,c("patid", "risk_tool", "risk_cat", "risk_score", "risk_date", "excl_date1", "excl_flag", "start_fup", "enddate")]

#######risk scoring eligibility
#exclude those with cvd event before start of follow up
cvdrisk_1[excl_date1 <= start_fup & excl_flag==1, exclude := 1]
cvdrisk_1[is.na(excl_date1), exclude := 0]
cvdrisk_1[exclude==0|is.na(exclude), primary_prev := 1]
riskscoring_pop <- cvdrisk_1[exclude==0|is.na(exclude)]
rm(cvdrisk_1, cvdrisk_all)

uniqueN(riskscoring_pop, by = "patid") #4,583,824

##merge those with an ever statin prescription
incident_statins <- read_parquet(paste0(datafiles, "incident_statins.parquet"))
riskscoring_pop <- merge(riskscoring_pop,incident_statins[, c("patid", "first_statin_dt")], by="patid", all.x=TRUE)

rm(incident_statins)
#create matching prop dates column
riskscoring_pop[, prop_dates := ceiling_date(risk_date, "month") -1]

setorder(riskscoring_pop, patid, prop_dates)

#create data table with month and years of study period using patients eligible for primary prevention
month_year_dt <- riskscoring_pop[,c("patid", "start_fup", "enddate")]
month_year_dt <- unique(month_year_dt) #delete repeating patids
uniqueN(month_year_dt, by = "patid") #8,088,498

month_year_dt <- month_year_dt[rep(seq_len(nrow(month_year_dt)), each = 153), ]
prop_dates <- seq(as.Date("2009/05/01"), as.Date("2022/01/01"), "month", format="%d/%m/%Y") -1
month_year_dt <- cbind(month_year_dt, prop_dates)
setorder(month_year_dt, patid, prop_dates)
#remove irrelevant rows
#remove where start_fup is after prop_date
month_year_dt <- month_year_dt[start_fup<prop_dates]
#remove rows in which enddate is before prop_date
month_year_dt <- month_year_dt[enddate>prop_dates]
month_year_dt[, .N] #71,321,592


#append proportion dates and primary prevention dataset
cvdrisk_pop_1<- rbind(riskscoring_pop,month_year_dt, fill=T)
uniqueN(cvdrisk_pop_1, by = "patid") #909,491
rm(riskscoring_pop)
constant_vars <- c("risk_tool", "risk_cat",  "risk_score", "first_statin_dt",
                   "risk_date", "excl_date1", "primary_prev")

#rm(month_year_dt)
#fill in constant variables
setorder(cvdrisk_pop_1, patid, prop_dates)
id_change = cvdrisk_pop_1[, c(TRUE, patid[-1] != patid[-.N])]
cvdrisk_pop_1[, (constant_vars) := lapply(.SD, function(x) x[cummax(((!is.na(x)) | id_change) * .I)]), .SDcols=constant_vars]

#cvdrisk_pop_1[, month_year := format(as.Date(prop_dates), "%Y-%m")]
rm(id_change)
#if more than one row in per person in same month keep latest date 
cvdrisk_pop_1[, year := lubridate::year(prop_dates)]
cvdrisk_pop_1[, month := lubridate::month(prop_dates)]

###record eligible if patient had a recorded risk score within 5 years of prop date
#create grace period of within 5 years of risk date 
cvdrisk_pop_1[, grace_period := risk_date+1826.25] 
cvdrisk_pop_1[prop_dates>grace_period, risk_record := 0]
cvdrisk_pop_1[prop_dates<=grace_period, risk_record := 1]
#if risk_date is missing code risk record as 0 - dates before first risk dates are missing risk date
cvdrisk_pop_1[is.na(risk_date), risk_record := 0]

#exclude those with exclusion code on the date of, or before month/year proportion date
cvdrisk_pop_1[!is.na(excl_date1) & excl_date1<=prop_dates, exclude := 1]
#exclude those with statin prescription before prop date
#cvdrisk_pop_1[!is.na(first_statin_dt) & first_statin_dt<=prop_dates, exclude := 1]
#rm(patients_dt)
cvdrisk_pop_1[, c("risk_score", "excl_flag", "primary_prev", "year", "month", "risk_record", "grace_period")]
cvdrisk_pop_1 <- cvdrisk_pop_1[exclude==0|is.na(exclude)]

#merge with study population dataset for grouped variables
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))

cvdrisk_pop_1 <- merge(cvdrisk_pop_1, studypop_all[,c("patid", "gender", "dob", "ethnicity", "deprivation")], by="patid", all.x=TRUE)

uniqueN(cvdrisk_pop_1, by = "patid")  

#remove where start_fup is after prop_date

#####
cvdrisk_pop_1[, start := pmax(start_fup,as.Date("2009-04-30"), na.rm=T)]
cvdrisk_pop_1[, exit := pmin(enddate,excl_date1,first_statin_dt, na.rm=T)]


cvdrisk_pop_1[, .N]
cvdrisk_pop_1 <- cvdrisk_pop_1[start_fup<prop_dates]
#remove rows in which enddate is before prop_date
cvdrisk_pop_1 <- cvdrisk_pop_1[enddate>prop_dates]

cvdrisk_pop_1[, total_pop := 1]
#generate age 
cvdrisk_pop_1[, age := as.numeric(prop_dates - dob) / 365.25]

##keep only those aged 84 years and younger 25-84
summary(cvdrisk_pop_1$age)
cvdrisk_pop_1 <- cvdrisk_pop_1[age<85]

cvdrisk_pop_1[, age_group := cut(age, breaks = c(25, 40, 75, 85),
                                     labels = c("25-39", "40-74", "75-84"))]

cvdrisk_pop_1[, .N, by=risk_record]

cvdrisk_pop_1 <- cvdrisk_pop_1[, c("patid", "risk_record", "total_pop", "prop_dates", "risk_tool", "risk_cat", "ethnicity", "deprivation", "gender", "age_group")]

###Save parquet
write_parquet(cvdrisk_pop_1, paste0(datafiles, "cvdrisk_pop_1.parquet"))
rm(cvdrisk_pop_1)

#############PART 2
cvdrisk_2 <- merge(cvdrisk_2,exclusions_all, by="patid", all.x=TRUE)
uniqueN(cvdrisk_2, by = "patid") #5,000,000
cvdrisk_2 <- cvdrisk_2[,c("patid", "risk_tool", "risk_cat", "risk_score", "risk_date", "excl_date1", "excl_flag", "start_fup", "enddate")]

#######risk scoring eligibility
#exclude those with cvd event before start of follow up
cvdrisk_2[excl_date1 <= start_fup & excl_flag==1, exclude := 1]
cvdrisk_2[is.na(excl_date1), exclude := 0]
cvdrisk_2[exclude==0|is.na(exclude), primary_prev := 1]
riskscoring_pop <- cvdrisk_2[exclude==0|is.na(exclude)]
rm(cvdrisk_2, cvdrisk_all)

uniqueN(riskscoring_pop, by = "patid") #4,583,824

##merge those with an ever statin prescription
incident_statins <- read_parquet(paste0(datafiles, "incident_statins.parquet"))
riskscoring_pop <- merge(riskscoring_pop,incident_statins[, c("patid", "first_statin_dt")], by="patid", all.x=TRUE)

rm(incident_statins)
#create matching prop dates column
riskscoring_pop[, prop_dates := ceiling_date(risk_date, "month") -1]

setorder(riskscoring_pop, patid, prop_dates)

#create data table with month and years of study period using patients eligible for primary prevention
month_year_dt <- riskscoring_pop[,c("patid", "start_fup", "enddate")]
month_year_dt <- unique(month_year_dt) #delete repeating patids
uniqueN(month_year_dt, by = "patid") #8,088,498

month_year_dt <- month_year_dt[rep(seq_len(nrow(month_year_dt)), each = 153), ]
prop_dates <- seq(as.Date("2009/05/01"), as.Date("2022/01/01"), "month", format="%d/%m/%Y") -1
month_year_dt <- cbind(month_year_dt, prop_dates)
setorder(month_year_dt, patid, prop_dates)
#remove irrelevant rows
#remove where start_fup is after prop_date
month_year_dt <- month_year_dt[start_fup<prop_dates]
#remove rows in which enddate is before prop_date
month_year_dt <- month_year_dt[enddate>prop_dates]
month_year_dt[, .N] #71,321,592


#append proportion dates and primary prevention dataset
cvdrisk_pop_2<- rbind(riskscoring_pop,month_year_dt, fill=T)
uniqueN(cvdrisk_pop_2, by = "patid") #909,491
rm(riskscoring_pop)
constant_vars <- c("risk_tool", "risk_cat",  "risk_score", "first_statin_dt",
                   "risk_date", "excl_date1", "primary_prev")

rm(month_year_dt)
#fill in constant variables
setorder(cvdrisk_pop_2, patid, prop_dates)
id_change = cvdrisk_pop_2[, c(TRUE, patid[-1] != patid[-.N])]
cvdrisk_pop_2[, (constant_vars) := lapply(.SD, function(x) x[cummax(((!is.na(x)) | id_change) * .I)]), .SDcols=constant_vars]

#cvdrisk_pop_2[, month_year := format(as.Date(prop_dates), "%Y-%m")]
rm(id_change)
#if more than one row in per person in same month keep latest date 
cvdrisk_pop_2[, year := lubridate::year(prop_dates)]
cvdrisk_pop_2[, month := lubridate::month(prop_dates)]

###record eligible if patient had a recorded risk score within 5 years of prop date
#create grace period of within 5 years of risk date 
cvdrisk_pop_2[, grace_period := risk_date+1826.25] 
cvdrisk_pop_2[prop_dates>grace_period, risk_record := 0]
cvdrisk_pop_2[prop_dates<=grace_period, risk_record := 1]
#if risk_date is missing code risk record as 0 - dates before first risk dates are missing risk date
cvdrisk_pop_2[is.na(risk_date), risk_record := 0]

#exclude those with exclusion code on the date of, or before month/year proportion date
cvdrisk_pop_2[!is.na(excl_date1) & excl_date1<=prop_dates, exclude := 1]
#exclude those with statin prescription before prop date
#cvdrisk_pop_2[!is.na(first_statin_dt) & first_statin_dt<=prop_dates, exclude := 1]
rm(split_parts, cvdrisk_pop_1)
cvdrisk_pop_2[, c("risk_score", "excl_flag", "primary_prev", "year", "month", "risk_record", "grace_period")]
cvdrisk_pop_2 <- cvdrisk_pop_2[exclude==0|is.na(exclude)]

#merge with study population dataset for grouped variables
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))

cvdrisk_pop_2 <- merge(cvdrisk_pop_2, studypop_all[,c("patid", "gender", "dob", "ethnicity", "deprivation")], by="patid", all.x=TRUE)


uniqueN(cvdrisk_pop_2, by = "patid")  

#####
cvdrisk_pop_2[, start := pmax(start_fup,as.Date("2009-04-30"), na.rm=T)]
cvdrisk_pop_2[, exit := pmin(enddate,excl_date1,first_statin_dt, na.rm=T)]


cvdrisk_pop_2[, .N]
cvdrisk_pop_2 <- cvdrisk_pop_2[start_fup<prop_dates]
#remove rows in which enddate is before prop_date
cvdrisk_pop_2 <- cvdrisk_pop_2[enddate>prop_dates]

cvdrisk_pop_2[, total_pop := 1]
#generate age 
cvdrisk_pop_2[, age := as.numeric(prop_dates - dob) / 365.25]

##keep only those aged 84 years and younger 25-84
summary(cvdrisk_pop_2$age)
cvdrisk_pop_2 <- cvdrisk_pop_2[age<85]

cvdrisk_pop_2[, age_group := cut(age, breaks = c(25, 40, 75, 85),
                                 labels = c("25-39", "40-74", "75-84"))]

cvdrisk_pop_2[, .N, by=risk_record]

cvdrisk_pop_2 <- cvdrisk_pop_2[, c("patid", "risk_record", "total_pop", "prop_dates", "risk_tool", "risk_cat", "ethnicity", "deprivation", "gender", "age_group")]

###Save parquet
write_parquet(cvdrisk_pop_2, paste0(datafiles, "cvdrisk_pop_2.parquet"))
rm(cvdrisk_pop_2)

#####PART 3

cvdrisk_3 <- merge(cvdrisk_3,exclusions_all, by="patid", all.x=TRUE)
uniqueN(cvdrisk_3, by = "patid") #5,000,000
cvdrisk_3 <- cvdrisk_3[,c("patid", "risk_tool", "risk_cat", "risk_score", "risk_date", "excl_date1", "excl_flag", "start_fup", "enddate")]

#######risk scoring eligibility
#exclude those with cvd event before start of follow up
cvdrisk_3[excl_date1 <= start_fup & excl_flag==1, exclude := 1]
cvdrisk_3[is.na(excl_date1), exclude := 0]
cvdrisk_3[exclude==0|is.na(exclude), primary_prev := 1]
riskscoring_pop <- cvdrisk_3[exclude==0|is.na(exclude)]
rm(cvdrisk_3, cvdrisk_all)

uniqueN(riskscoring_pop, by = "patid") #4,583,824

##merge those with an ever statin prescription
incident_statins <- read_parquet(paste0(datafiles, "incident_statins.parquet"))
riskscoring_pop <- merge(riskscoring_pop,incident_statins[, c("patid", "first_statin_dt")], by="patid", all.x=TRUE)

rm(incident_statins)
#create matching prop dates column
riskscoring_pop[, prop_dates := ceiling_date(risk_date, "month") -1]

setorder(riskscoring_pop, patid, prop_dates)

#create data table with month and years of study period using patients eligible for primary prevention
month_year_dt <- riskscoring_pop[,c("patid", "start_fup", "enddate")]
month_year_dt <- unique(month_year_dt) #delete repeating patids
uniqueN(month_year_dt, by = "patid") #8,088,498

month_year_dt <- month_year_dt[rep(seq_len(nrow(month_year_dt)), each = 153), ]
prop_dates <- seq(as.Date("2009/05/01"), as.Date("2022/01/01"), "month", format="%d/%m/%Y") -1
month_year_dt <- cbind(month_year_dt, prop_dates)
setorder(month_year_dt, patid, prop_dates)
#remove irrelevant rows
#remove where start_fup is after prop_date
month_year_dt <- month_year_dt[start_fup<prop_dates]
#remove rows in which enddate is before prop_date
month_year_dt <- month_year_dt[enddate>prop_dates]
month_year_dt[, .N] #71,321,592


#append proportion dates and primary prevention dataset
cvdrisk_pop_3<- rbind(riskscoring_pop,month_year_dt, fill=T)
uniqueN(cvdrisk_pop_3, by = "patid") #909,491
rm(riskscoring_pop)
constant_vars <- c("risk_tool", "risk_cat",  "risk_score", "first_statin_dt",
                   "risk_date", "excl_date1", "primary_prev")

rm(month_year_dt)
#fill in constant variables
setorder(cvdrisk_pop_3, patid, prop_dates)
id_change = cvdrisk_pop_3[, c(TRUE, patid[-1] != patid[-.N])]
cvdrisk_pop_3[, (constant_vars) := lapply(.SD, function(x) x[cummax(((!is.na(x)) | id_change) * .I)]), .SDcols=constant_vars]

#cvdrisk_pop_3[, month_year := format(as.Date(prop_dates), "%Y-%m")]
rm(id_change)
#if more than one row in per person in same month keep latest date 
cvdrisk_pop_3[, year := lubridate::year(prop_dates)]
cvdrisk_pop_3[, month := lubridate::month(prop_dates)]

# setorder(cvdrisk_pop_3, patid, -prop_dates, -risk_date)
# cvdrisk_pop_3[, dupl := seq_len(.N), by = list(patid, prop_dates)]
# cvdrisk_pop_3 <- cvdrisk_pop_3[dupl==1]
# cvdrisk_pop_3 <- cvdrisk_pop_3[, .SD[dupl == 1], by = list(patid, prop_dates)]
# unique_dt <- unique(dt, by = c("patid", "prop_dates"))

###record eligible if patient had a recorded risk score within 5 years of prop date
#create grace period of within 5 years of risk date 
cvdrisk_pop_3[, grace_period := risk_date+1826.25] 
cvdrisk_pop_3[prop_dates>grace_period, risk_record := 0]
cvdrisk_pop_3[prop_dates<=grace_period, risk_record := 1]
#if risk_date is missing code risk record as 0 - dates before first risk dates are missing risk date
cvdrisk_pop_3[is.na(risk_date), risk_record := 0]

#exclude those with exclusion code on the date of, or before month/year proportion date
cvdrisk_pop_3[!is.na(excl_date1) & excl_date1<=prop_dates, exclude := 1]
#exclude those with statin prescription before prop date
#cvdrisk_pop_3[!is.na(first_statin_dt) & first_statin_dt<=prop_dates, exclude := 1]
#rm(patients_dt)
cvdrisk_pop_3[, c("risk_score", "excl_flag", "primary_prev", "year", "month", "risk_record", "grace_period")]
cvdrisk_pop_3 <- cvdrisk_pop_3[exclude==0|is.na(exclude)]

#merge with study population dataset for grouped variables
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))

cvdrisk_pop_3 <- merge(cvdrisk_pop_3, studypop_all[,c("patid", "gender", "dob", "ethnicity", "deprivation")], by="patid", all.x=TRUE)

# cvdrisk_pop_3 <- cvdrisk_pop_3[prop_dates>="2009-05-31"] 
# cvdrisk_pop_3 <- cvdrisk_pop_3[prop_dates<="2021-10-31"] 

uniqueN(cvdrisk_pop_3, by = "patid")  

#remove where start_fup is after prop_date

#####
cvdrisk_pop_3[, start := pmax(start_fup,as.Date("2009-04-30"), na.rm=T)]
cvdrisk_pop_3[, exit := pmin(enddate,excl_date1,first_statin_dt, na.rm=T)]


cvdrisk_pop_3[, .N]
cvdrisk_pop_3 <- cvdrisk_pop_3[start_fup<prop_dates]
#remove rows in which enddate is before prop_date
cvdrisk_pop_3 <- cvdrisk_pop_3[enddate>prop_dates]

cvdrisk_pop_3[, total_pop := 1]
#generate age 
cvdrisk_pop_3[, age := as.numeric(prop_dates - dob) / 365.25]

##keep only those aged 84 years and younger 25-84
summary(cvdrisk_pop_3$age)
cvdrisk_pop_3 <- cvdrisk_pop_3[age<85]

cvdrisk_pop_3[, age_group := cut(age, breaks = c(25, 40, 75, 85),
                                 labels = c("25-39", "40-74", "75-84"))]

cvdrisk_pop_3[, .N, by=risk_record]

cvdrisk_pop_3 <- cvdrisk_pop_3[, c("patid", "risk_record", "total_pop", "prop_dates", "risk_tool", "risk_cat", "ethnicity", "deprivation", "gender", "age_group")]

###Save parquet
write_parquet(cvdrisk_pop_3, paste0(datafiles, "cvdrisk_pop_3.parquet"))
rm(cvdrisk_pop_3)

#####PART 4
cvdrisk_4 <- merge(cvdrisk_4,exclusions_all, by="patid", all.x=TRUE)
uniqueN(cvdrisk_4, by = "patid") #5,000,000
cvdrisk_4 <- cvdrisk_4[,c("patid", "risk_tool", "risk_cat", "risk_score", "risk_date", "excl_date1", "excl_flag", "start_fup", "enddate")]

#######risk scoring eligibility
#exclude those with cvd event before start of follow up
cvdrisk_4[excl_date1 <= start_fup & excl_flag==1, exclude := 1]
cvdrisk_4[is.na(excl_date1), exclude := 0]
cvdrisk_4[exclude==0|is.na(exclude), primary_prev := 1]
riskscoring_pop <- cvdrisk_4[exclude==0|is.na(exclude)]
rm(cvdrisk_4, cvdrisk_all)

uniqueN(riskscoring_pop, by = "patid") #4,583,824

##merge those with an ever statin prescription
incident_statins <- read_parquet(paste0(datafiles, "incident_statins.parquet"))
riskscoring_pop <- merge(riskscoring_pop,incident_statins[, c("patid", "first_statin_dt")], by="patid", all.x=TRUE)

rm(incident_statins)
#create matching prop dates column
riskscoring_pop[, prop_dates := ceiling_date(risk_date, "month") -1]

setorder(riskscoring_pop, patid, prop_dates)

#create data table with month and years of study period using patients eligible for primary prevention
month_year_dt <- riskscoring_pop[,c("patid", "start_fup", "enddate")]
month_year_dt <- unique(month_year_dt) #delete repeating patids
uniqueN(month_year_dt, by = "patid") #8,088,498

month_year_dt <- month_year_dt[rep(seq_len(nrow(month_year_dt)), each = 153), ]
prop_dates <- seq(as.Date("2009/05/01"), as.Date("2022/01/01"), "month", format="%d/%m/%Y") -1
month_year_dt <- cbind(month_year_dt, prop_dates)
setorder(month_year_dt, patid, prop_dates)
#remove irrelevant rows
#remove where start_fup is after prop_date
month_year_dt <- month_year_dt[start_fup<prop_dates]
#remove rows in which enddate is before prop_date
month_year_dt <- month_year_dt[enddate>prop_dates]
month_year_dt[, .N] #71,321,592


#append proportion dates and primary prevention dataset
cvdrisk_pop_4<- rbind(riskscoring_pop,month_year_dt, fill=T)
uniqueN(cvdrisk_pop_4, by = "patid") #909,491
rm(riskscoring_pop)
constant_vars <- c("risk_tool", "risk_cat",  "risk_score", "first_statin_dt",
                   "risk_date", "excl_date1", "primary_prev")

rm(month_year_dt)
#fill in constant variables
setorder(cvdrisk_pop_4, patid, prop_dates)
id_change = cvdrisk_pop_4[, c(TRUE, patid[-1] != patid[-.N])]
cvdrisk_pop_4[, (constant_vars) := lapply(.SD, function(x) x[cummax(((!is.na(x)) | id_change) * .I)]), .SDcols=constant_vars]

#cvdrisk_pop_4[, month_year := format(as.Date(prop_dates), "%Y-%m")]
rm(id_change)
#if more than one row in per person in same month keep latest date 
cvdrisk_pop_4[, year := lubridate::year(prop_dates)]
cvdrisk_pop_4[, month := lubridate::month(prop_dates)]

###record eligible if patient had a recorded risk score within 5 years of prop date
#create grace period of within 5 years of risk date 
cvdrisk_pop_4[, grace_period := risk_date+1826.25] 
cvdrisk_pop_4[prop_dates>grace_period, risk_record := 0]
cvdrisk_pop_4[prop_dates<=grace_period, risk_record := 1]
#if risk_date is missing code risk record as 0 - dates before first risk dates are missing risk date
cvdrisk_pop_4[is.na(risk_date), risk_record := 0]

#exclude those with exclusion code on the date of, or before month/year proportion date
cvdrisk_pop_4[!is.na(excl_date1) & excl_date1<=prop_dates, exclude := 1]
#exclude those with statin prescription before prop date
#cvdrisk_pop_4[!is.na(first_statin_dt) & first_statin_dt<=prop_dates, exclude := 1]
#rm(patients_dt)
cvdrisk_pop_4[, c("risk_score", "excl_flag", "primary_prev", "year", "month", "risk_record", "grace_period")]
cvdrisk_pop_4 <- cvdrisk_pop_4[exclude==0|is.na(exclude)]

#merge with study population dataset for grouped variables
studypop_all <- read_parquet(paste0(datafiles, "studypop_all.parquet"))

cvdrisk_pop_4 <- merge(cvdrisk_pop_4, studypop_all[,c("patid", "gender", "dob", "ethnicity", "deprivation")], by="patid", all.x=TRUE)

uniqueN(cvdrisk_pop_4, by = "patid")  


#remove where start_fup is after prop_date

#####
cvdrisk_pop_4[, start := pmax(start_fup,as.Date("2009-04-30"), na.rm=T)]
cvdrisk_pop_4[, exit := pmin(enddate,excl_date1,first_statin_dt, na.rm=T)]


cvdrisk_pop_4[, .N]
cvdrisk_pop_4 <- cvdrisk_pop_4[start_fup<prop_dates]
#remove rows in which enddate is before prop_date
cvdrisk_pop_4 <- cvdrisk_pop_4[enddate>prop_dates]

cvdrisk_pop_4[, total_pop := 1]
#generate age 
cvdrisk_pop_4[, age := as.numeric(prop_dates - dob) / 365.25]

##keep only those aged 84 years and younger 25-84
summary(cvdrisk_pop_4$age)
cvdrisk_pop_4 <- cvdrisk_pop_4[age<85]

cvdrisk_pop_4[, age_group := cut(age, breaks = c(25, 40, 75, 85),
                                 labels = c("25-39", "40-74", "75-84"))]

cvdrisk_pop_4[, .N, by=risk_record]

cvdrisk_pop_4 <- cvdrisk_pop_4[, c("patid", "risk_record", "total_pop", "prop_dates", "risk_tool", "risk_cat", "ethnicity", "deprivation", "gender", "age_group")]

###Save parquet
write_parquet(cvdrisk_pop_4, paste0(datafiles, "cvdrisk_pop_4.parquet"))




