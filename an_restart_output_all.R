##########################################################################################################
###################FIGURE x: FACTORS ASSOCIATED WITH restart
#############################PRIMARY AND SECONDARY PREVENTION##################################
#############################################################################################
rm(list=ls())

library(broom)
library(openxlsx)
library(cowplot)

lapply(c("ggpubr", "grid", "gridExtra", "forcats"), require, character.only=T)

##open parquet
prim_restart <- read_parquet(paste0(datafiles, "prim_restart_fact.parquet"))

##################################################################################################
###########################################PRIMARY PREVENTION#####################################
##################################################################################################

prim_restart <- prim_restart[ethnicity!="missing" & deprivation!="missing"]

prim_restart[, c("age_cat", "ethnicity", "gender", "deprivation") := lapply(.SD, as.factor), .SDcols = c("age_cat", "ethnicity", "gender", "deprivation")]
prim_restart[ethnicity=="missing", ethnicity := NA]
prim_restart[deprivation=="missing", deprivation := NA]

prim_restart[, age_cat := factor(age_cat, levels = c("25-39", "40-49", "50-59", "60-69", "70+"))]
prim_restart[, gender := factor(gender, levels = c("Male", "Female"))]
prim_restart[, ethnicity := factor(ethnicity, levels = c("White", "Black", "South Asian", "Mixed", "Other"))]
prim_restart[, deprivation := factor(deprivation, levels = c("1 - least deprived", "2", "3", "4", "5 - most deprived"))]

prim_restart[, .N, by="age_cat"]
prim_restart[, .N, by="ethnicity"]
prim_restart[, .N, by="gender"]
prim_restart[, .N, by="deprivation"]

cox_restart_young<- prim_restart[age_cat=="25-39"]
cox_restart_white<- prim_restart[ethnicity=="White"]
cox_restart_male<- prim_restart[gender=="Male"]
cox_restart_depri<- prim_restart[deprivation=="1 - least deprived"]

N_young_adj <- nrow(cox_restart_young)
N_white_adj <- nrow(cox_restart_white)
N_male_adj <- nrow(cox_restart_male)
N_depri_adj <- nrow(cox_restart_depri)

#Adjusted Cox model
cox_model1 <- coxph(Surv(time, restart) ~ age_cat + gender + ethnicity + deprivation,
                    data=prim_restart)
summary(cox_model1)

#model.matrix(cox_model1)

adjusted_df <- tidy(cox_model1, conf.int=T, exponentiate=T)

#add rows for reference groups
adjusted_df <- adjusted_df %>% add_row(term = "25-39", estimate = 1, conf.low=1,conf.high=1)
adjusted_df <- adjusted_df %>% add_row(term = "White", estimate = 1, conf.low=1,conf.high=1)
adjusted_df <- adjusted_df %>% add_row(term = "1 - least deprived", estimate = 1, conf.low=1,conf.high=1)
adjusted_df <- adjusted_df %>% add_row(term = "Men", estimate = 1, conf.low=1,conf.high=1)

N <- nobs(cox_model1)

adj_n <-  cox_model1 %>% model.matrix() %>% 
  colSums()%>%                           # show the total number in each level
  as.data.frame() %>% t %>% as.data.frame() %>% 
  dplyr::select(c("age_cat40-49", "age_cat50-59", "age_cat60-69", "age_cat70+", "genderFemale", "ethnicityBlack", "ethnicitySouth Asian",
                  "ethnicityMixed", "ethnicityOther", "deprivation2", "deprivation3",
                  "deprivation4", "deprivation5 - most deprived"))

adj_n <- pivot_longer(adj_n, cols = everything(), names_to = "term", values_to = "N")

adjusted_df <- merge(adjusted_df, adj_n, by="term", all = T)
##add N of reference groups
adjusted_df <- as.data.table(adjusted_df)

adjusted_df[term=="25-39" & is.na(N), N := N_young_adj]
adjusted_df[term=="White" & is.na(N), N := N_white_adj]
adjusted_df[term=="1 - least deprived" & is.na(N), N := N_depri_adj]
adjusted_df[term=="Men" & is.na(N), N := N_male_adj]

cols_round <- c("estimate", "conf.low", "conf.high")

for (col_name in cols_round) {
  adjusted_df[[col_name]] <- round(adjusted_df[[col_name]], 2)
}

adjusted_df <- adjusted_df %>%
  mutate(HR = paste(adjusted_df$estimate," (",adjusted_df$conf.low," - ",
                    adjusted_df$conf.high,")",sep = ""))

#reference categories HR
adjusted_df[term=="25-39", HR := "1 (reference)"]
adjusted_df[term=="White", HR := "1 (reference)"]
adjusted_df[term=="1 - least deprived", HR := "1 (reference)"]
adjusted_df[term=="Men", HR := "1 (reference)"]


#####person years and events
pyears_age <- pyears(Surv(time, restart) ~ age_cat, data = prim_restart, model = T)
pyears_age <- tidy(pyears_age)
colnames(pyears_age)[colnames(pyears_age) == "n"] <- "N"

#pyears_age$term <- c('25-39', 'age_cat40-49', 'age_cat50-59', 'age_cat60-69', 'age_cat70+')

pyears_sex <- pyears(Surv(time, restart) ~ gender, data = prim_restart, model = T)
pyears_sex <- tidy(pyears_sex)
colnames(pyears_sex)[colnames(pyears_sex) == "n"] <- "N"

pyears_eth <- pyears(Surv(time, restart) ~ ethnicity, data = prim_restart, model = T)
pyears_eth <- tidy(pyears_eth)
colnames(pyears_eth)[colnames(pyears_eth) == "n"] <- "N"

pyears_depr <- pyears(Surv(time, restart) ~ deprivation, data = prim_restart, model = T)
pyears_depr <- tidy(pyears_depr)
colnames(pyears_depr)[colnames(pyears_depr) == "n"] <- "N"
prim_restart[, .N, by="deprivation"]
pyears_depr$term <- c('1 - least deprived', 'deprivation2', 'deprivation3', 'deprivation4', 'deprivation5 - most deprived')
adjusted_df <- as.data.table(adjusted_df)

adjusted_df <- merge(adjusted_df, pyears_age, by="N", all.x=T)
adjusted_df <- merge(adjusted_df, pyears_sex, by="N", all.x=T)
adjusted_df <- merge(adjusted_df, pyears_eth, by="N", all.x=T)
adjusted_df[, event := coalesce(event.x, event.y, event)]
adjusted_df[, c("event.x", "event.y", "pyears", "pyears.x", "pyears.y") := NULL]
adjusted_df <- merge(adjusted_df, pyears_depr, by=c("N", "term"), all.x=T)
adjusted_df[, event := coalesce(event.x, event.y)]

adjusted_df <- adjusted_df[, c("event.x", "event.y", "pyears", "N") := NULL]
setnames(adjusted_df, "event", "N")

##########FOREST PLOT########
adjusted_df[term=="age_cat40-49", term := "40-49"]
adjusted_df[term=="age_cat50-59", term := "50-59"]
adjusted_df[term=="age_cat60-69", term := "60-69"]
adjusted_df[term=="age_cat70+", term := "70+"]
adjusted_df[term=="genderFemale", term := "Women"]
adjusted_df[term=="ethnicitySouth Asian", term := "South Asian"]
adjusted_df[term=="ethnicityBlack", term := "Black"]
adjusted_df[term=="ethnicityMixed", term := "Mixed"]
adjusted_df[term=="ethnicityOther", term := "Other"]
adjusted_df[term=="ethnicityWhite", term := "White"]
adjusted_df[term=="deprivation2", term := "2"]
adjusted_df[term=="deprivation3", term := "3"]
adjusted_df[term=="deprivation4", term := "4"]
adjusted_df[term=="deprivation5 - most deprived", term := "5 - most"]

adjusted_df[term=="1 - least deprived", term := "1 - least"]
adjusted_df[term=="1 - least deprived"|term=="2"|term=="3"|term=="4"|term=="5 - most",
            deprivation := term]
adjusted_df[term=="White"|term=="Black"|term=="Mixed"|term=="Other"|term=="South Asian",
            ethnicity := term]
adjusted_df[term=="Men"|term=="Women", gender := term]
adjusted_df[term=="25-39"|term=="40-49"|term=="50-59"|term=="60-69"|term=="70+",
            age_cat := term]


adjusted_df[!is.na(age_cat), Factor := "Age Group"]
adjusted_df[!is.na(gender), Factor := "Gender"]
adjusted_df[!is.na(deprivation), Factor := "Deprivation"]
adjusted_df[!is.na(ethnicity), Factor := "Ethnicity"]

adjusted_df[term=="1 - least", Factor := "Deprivation"]
adjusted_df[term=="5 - most", Factor := "Deprivation"]

adjusted_df[!is.na(age_cat), factor_order := 1]
adjusted_df[!is.na(gender), factor_order := 2]
adjusted_df[!is.na(ethnicity), factor_order := 3]
adjusted_df[!is.na(deprivation), factor_order := 4]

plot_data <- adjusted_df %>%
  mutate(factor_order=case_when(
    Factor == "Age Group" ~ 1,
    Factor == "Gender" ~ 2,
    Factor == "Ethnicity" ~ 3,
    Factor == "Deprivation" ~ 4,                    # for facet order
    TRUE ~ NA_integer_  # Handle other cases or missing values
  ))

plot_data <- plot_data %>%
  dplyr::select(factor_order, Factor, term, N, estimate, conf.low, conf.high, HR)

plot_data[,1:4] <- map(plot_data[,1:4], as.factor)                      # factor coercion
plot_data$Factor %>% factor(
  labels = c("Deprivation", "Ethnicity", "Gender", "Age Group"))    # labelled for facet
plot_data

###order correctly
plot_data[HR=="1 (reference)", ref := 1] 
plot_data[HR!="1 (reference)", ref := 0]

setorder(plot_data, factor_order, -ref, term)

#depr_order <- c("1 - least deprived", "2", "3", "4", "5 - most deprived")
#plot_data[, term := factor(term, levels=depr_order)]

plot_data$order <- 17:1
plot_data$order %>% as.factor
plot_data$order %>% levels()

#################################################################################
##############################PLOT###################################

prim_data <- plot_data
##ALL
prim_main_plot <- ggplot(prim_data, aes(y=order, x=estimate, xmin=conf.low, xmax=conf.high, color=term)) +
  geom_point() +
  geom_errorbarh(height=.3) +
  geom_vline(xintercept=1, lty=2) +
  ggtitle(" ") +
  theme_classic2( ) +
  theme(strip.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  xlab("Hazard ratio")

###Plot text
prim_data[term!="25-39"&term!="White"&term!="1 - least"&term!="Men", Factor := " "] 

prim_left_plot <- ggplot(data = prim_data, aes(y = order)) + xlim(0,0.8) +
  geom_text(aes(x = 0, label = Factor ),lineheight = 0.001, hjust = 0, size = 4, colour = "black") + 
  geom_text(aes(x = 0.35, label = term),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") +
  geom_text(aes(x = 0.68, label = N),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") + 
  theme_classic2() + 
  ggtitle("   Category                                                                                  Events")+  
  theme(plot.title = element_text(size = 10, lineheight=.001, face="bold", vjust=-8),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        axis.title.x = element_text(colour = "white", size = 8), # make the text invisible
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "white"),               # make the lines invisible
        axis.text.x = element_text(colour = "white", size = 8),  # make the text invisible
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + xlab(" ") 

###Plot HR and 95% CI

prim_right_plot <- ggplot(data = prim_data,  aes(y = order)) + xlim(0,0.5) +
  geom_text(aes(x = 0, label = HR),lineheight = 0.01, hjust = 0, size = 4, colour = "black") + 
  theme_classic2( ) +
  ggtitle("        HR (95%CI)")+
  theme(plot.title = element_text(size = 10 ,lineheight=.001, face="bold", vjust=-8),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        axis.title.x = element_text(colour = "white", size = 12),# make the text invisible
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "white", size = 12), # make the text invisible
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "white"),              # make the lines invisible
        panel.background=element_blank(),
        plot.background=element_blank()) + xlab(" ") 

#grid.arrange(prim_left_plot, prim_main_plot, prim_right_plot,  ncol=3)

prim_restart_plot <- plot_grid(prim_left_plot, prim_main_plot, prim_right_plot, labels=c('A. Primary Prevention'), align="center", label_size = 12, label_y=0.99, ncol=3)
#ggsave("prim_restart_prim.png", plot = last_plot(), width = 12, height = 8, units = "in", dpi = 300, path = graphs)


##################################################################################################
#########################################SECONDARY PREVENTION#####################################
##################################################################################################
##open parquet
sec_restart <- read_parquet(paste0(datafiles, "sec_restart_fact.parquet"))

uniqueN(sec_restart, "patid")
sec_restart[, .N, by="CVD"]

sec_restart <- sec_restart[ethnicity!="missing" & deprivation!="missing"]

sec_restart[, c("ethnicity", "gender", "deprivation", "age_cat") := lapply(.SD, as.factor), .SDcols = c("ethnicity", "gender", "deprivation", "age_cat")]
sec_restart[ethnicity=="missing", ethnicity := NA]
sec_restart[deprivation=="missing", deprivation := NA]

sec_restart[, age_cat := factor(age_cat, levels = c("25-39", "40-49", "50-59", "60-69", "70+"))]
sec_restart[, gender := factor(gender, levels = c("Male", "Female"))]
sec_restart[, ethnicity := factor(ethnicity, levels = c("White", "Black", "South Asian", "Mixed", "Other"))]
sec_restart[, deprivation := factor(deprivation, levels = c("1 - least deprived", "2", "3", "4", "5 - most deprived"))]

#sec_restart <- sec_restart[CVD!="Acute Coronary Syndrome (non-specific)" & CVD!="Coronary Heart Disease (non-specific)"]
sec_restart[, .N, by="age_cat"]
sec_restart[, .N, by="ethnicity"]
sec_restart[, .N, by="gender"]
sec_restart[, .N, by="deprivation"]
sec_restart[, .N, by="CVD"]
cox_restart_young<- sec_restart[age_cat=="25-39"]
cox_restart_white<- sec_restart[ethnicity=="White"]
cox_restart_male<- sec_restart[gender=="Male"]
cox_restart_depri<- sec_restart[deprivation=="1 - least deprived"]

N_young_adj <- nrow(cox_restart_young)
N_white_adj <- nrow(cox_restart_white)
N_male_adj <- nrow(cox_restart_male)
N_depri_adj <- nrow(cox_restart_depri)

#Adjusted Cox model
cox_model1 <- coxph(Surv(time, restart) ~ age_cat + gender + ethnicity + deprivation,
                    data=sec_restart)
summary(cox_model1)

#model.matrix(cox_model1)

adjusted_df <- tidy(cox_model1, conf.int=T, exponentiate=T)

#add rows for reference groups
adjusted_df <- adjusted_df %>% add_row(term = "25-39", estimate = 1, conf.low=1,conf.high=1)
adjusted_df <- adjusted_df %>% add_row(term = "White", estimate = 1, conf.low=1,conf.high=1)
adjusted_df <- adjusted_df %>% add_row(term = "1 - least deprived", estimate = 1, conf.low=1,conf.high=1)
adjusted_df <- adjusted_df %>% add_row(term = "Men", estimate = 1, conf.low=1,conf.high=1)

N <- nobs(cox_model1)

adj_n <-  cox_model1 %>% model.matrix() %>% 
  colSums()%>%                           # show the total number in each level
  as.data.frame() %>% t %>% as.data.frame() %>% 
  dplyr::select(c("age_cat40-49", "age_cat50-59", "age_cat60-69", "age_cat70+", "genderFemale", "ethnicityBlack", "ethnicitySouth Asian",
                  "ethnicityMixed", "ethnicityOther", "deprivation2", "deprivation3",
                  "deprivation4", "deprivation5 - most deprived"))

adj_n <- pivot_longer(adj_n, cols = everything(), names_to = "term", values_to = "N")

adjusted_df <- merge(adjusted_df, adj_n, by="term", all = T)
##add N of reference groups
adjusted_df <- as.data.table(adjusted_df)

adjusted_df[term=="25-39" & is.na(N), N := N_young_adj]
adjusted_df[term=="White" & is.na(N), N := N_white_adj]
adjusted_df[term=="1 - least deprived" & is.na(N), N := N_depri_adj]
adjusted_df[term=="Men" & is.na(N), N := N_male_adj]

cols_round <- c("estimate", "conf.low", "conf.high")

for (col_name in cols_round) {
  adjusted_df[[col_name]] <- round(adjusted_df[[col_name]], 2)
}

adjusted_df <- adjusted_df %>%
  mutate(HR = paste(adjusted_df$estimate," (",adjusted_df$conf.low," - ",
                    adjusted_df$conf.high,")",sep = ""))

#reference categories HR
adjusted_df[term=="25-39", HR := "1 (reference)"]
adjusted_df[term=="White", HR := "1 (reference)"]
adjusted_df[term=="1 - least deprived", HR := "1 (reference)"]
adjusted_df[term=="Men", HR := "1 (reference)"]


#####person years and events
pyears_age <- pyears(Surv(time, restart) ~ age_cat, data = sec_restart, model = T)
pyears_age <- tidy(pyears_age)
colnames(pyears_age)[colnames(pyears_age) == "n"] <- "N"

#pyears_age$term <- c('25-39', 'age_cat40-49', 'age_cat50-59', 'age_cat60-69', 'age_cat70+')

pyears_sex <- pyears(Surv(time, restart) ~ gender, data = sec_restart, model = T)
pyears_sex <- tidy(pyears_sex)
colnames(pyears_sex)[colnames(pyears_sex) == "n"] <- "N"

pyears_eth <- pyears(Surv(time, restart) ~ ethnicity, data = sec_restart, model = T)
pyears_eth <- tidy(pyears_eth)
colnames(pyears_eth)[colnames(pyears_eth) == "n"] <- "N"

pyears_depr <- pyears(Surv(time, restart) ~ deprivation, data = sec_restart, model = T)
pyears_depr <- tidy(pyears_depr)
colnames(pyears_depr)[colnames(pyears_depr) == "n"] <- "N"

#pyears_depr$term <- c('1 - least deprived', 'deprivation2', 'deprivation3', 'deprivation4', 'deprivation5 - most deprived')
adjusted_df <- as.data.table(adjusted_df)

adjusted_df <- merge(adjusted_df, pyears_age, by="N", all.x=T)
adjusted_df <- merge(adjusted_df, pyears_sex, by="N", all.x=T)
adjusted_df <- merge(adjusted_df, pyears_eth, by="N", all.x=T)
adjusted_df[, event := coalesce(event.x, event.y, event)]
adjusted_df[, c("event.x", "event.y", "pyears", "pyears.x", "pyears.y") := NULL]
adjusted_df <- merge(adjusted_df, pyears_depr, by="N", all.x=T)
adjusted_df[, event := coalesce(event.x, event.y)]

adjusted_df <- adjusted_df[, c("event.x", "event.y", "pyears", "N") := NULL]
setnames(adjusted_df, "event", "N")

##########FOREST PLOT########
adjusted_df[term=="age_cat40-49", term := "40-49"]
adjusted_df[term=="age_cat50-59", term := "50-59"]
adjusted_df[term=="age_cat60-69", term := "60-69"]
adjusted_df[term=="age_cat70+", term := "70+"]
adjusted_df[term=="genderFemale", term := "Women"]
adjusted_df[term=="ethnicitySouth Asian", term := "South Asian"]
adjusted_df[term=="ethnicityBlack", term := "Black"]
adjusted_df[term=="ethnicityMixed", term := "Mixed"]
adjusted_df[term=="ethnicityOther", term := "Other"]
adjusted_df[term=="ethnicityWhite", term := "White"]
adjusted_df[term=="deprivation2", term := "2"]
adjusted_df[term=="deprivation3", term := "3"]
adjusted_df[term=="deprivation4", term := "4"]
adjusted_df[term=="deprivation5 - most deprived", term := "5 - most"]

adjusted_df[term=="1 - least deprived", term := "1 - least"]
adjusted_df[term=="1 - least deprived"|term=="2"|term=="3"|term=="4"|term=="5 - most",
            deprivation := term]
adjusted_df[term=="White"|term=="Black"|term=="Mixed"|term=="Other"|term=="South Asian",
            ethnicity := term]
adjusted_df[term=="Men"|term=="Women", gender := term]
adjusted_df[term=="25-39"|term=="40-49"|term=="50-59"|term=="60-69"|term=="70+",
            age_cat := term]


adjusted_df[!is.na(age_cat), Factor := "Age Group"]
adjusted_df[!is.na(gender), Factor := "Gender"]
adjusted_df[!is.na(deprivation), Factor := "Deprivation"]
adjusted_df[!is.na(ethnicity), Factor := "Ethnicity"]

adjusted_df[term=="1 - least", Factor := "Deprivation"]
adjusted_df[term=="5 - most", Factor := "Deprivation"]

adjusted_df[!is.na(age_cat), factor_order := 1]
adjusted_df[!is.na(gender), factor_order := 2]
adjusted_df[!is.na(ethnicity), factor_order := 3]
adjusted_df[!is.na(deprivation), factor_order := 4]

plot_data <- adjusted_df %>%
  mutate(factor_order=case_when(
    Factor == "Age Group" ~ 1,
    Factor == "Gender" ~ 2,
    Factor == "Ethnicity" ~ 3,
    Factor == "Deprivation" ~ 4,                    # for facet order
    TRUE ~ NA_integer_  # Handle other cases or missing values
  ))

plot_data <- plot_data %>%
  dplyr::select(factor_order, Factor, term, N, estimate, conf.low, conf.high, HR)

plot_data[,1:4] <- map(plot_data[,1:4], as.factor)                      # factor coercion
plot_data$Factor %>% factor(
  labels = c("Deprivation", "Ethnicity", "Gender", "Age Group"))    # labelled for facet
plot_data

###order correctly
plot_data[HR=="1 (reference)", ref := 1] 
plot_data[HR!="1 (reference)", ref := 0]

setorder(plot_data, factor_order, -ref, term)

#depr_order <- c("1 - least deprived", "2", "3", "4", "5 - most deprived")
#plot_data[, term := factor(term, levels=depr_order)]

plot_data$order <- 17:1
plot_data$order %>% as.factor
plot_data$order %>% levels()

#################################################################################
##############################PLOT###################################

sec_data <- plot_data
##ALL
main_plot <- ggplot(sec_data, aes(y=order, x=estimate, xmin=conf.low, xmax=conf.high, color=term)) +
  geom_point() +
  geom_errorbarh(height=.3) +
  geom_vline(xintercept=1, lty=2) +
  ggtitle(" ") +
  theme_classic2( ) +
  theme(strip.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  xlab("Hazard ratio")

###Plot text
sec_data[term!="25-39"&term!="White"&term!="1 - least"&term!="Men", Factor := " "] 

left_plot <- ggplot(data = sec_data, aes(y = order)) + xlim(0,0.8) +
  geom_text(aes(x = 0, label = Factor ),lineheight = 0.001, hjust = 0, size = 4, colour = "black") + 
  geom_text(aes(x = 0.35, label = term),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") +
  geom_text(aes(x = 0.68, label = N),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") + 
  theme_classic2() + 
  ggtitle("   Category                                                                                  Events")+  
  theme(plot.title = element_text(size = 10, lineheight=.001, face="bold", vjust=-8),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        axis.title.x = element_text(colour = "white", size = 8), # make the text invisible
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "white"),               # make the lines invisible
        axis.text.x = element_text(colour = "white", size = 8),  # make the text invisible
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + xlab(" ") 

###Plot HR and 95% CI

right_plot <- ggplot(data = sec_data,  aes(y = order)) + xlim(0,0.5) +
  geom_text(aes(x = 0, label = HR),lineheight = 0.01, hjust = 0, size = 4, colour = "black") + 
  theme_classic2( ) +
  ggtitle("        HR (95%CI)")+
  theme(plot.title = element_text(size = 10 ,lineheight=.001, face="bold", vjust=-8),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        axis.title.x = element_text(colour = "white", size = 12),# make the text invisible
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "white", size = 12), # make the text invisible
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "white"),              # make the lines invisible
        panel.background=element_blank(),
        plot.background=element_blank()) + xlab(" ") 


sec_restart_plot <- plot_grid(left_plot, main_plot, right_plot, labels=c('B. Secondary Prevention'), align="center", label_size = 12, label_y=0.99, ncol=3)

#########################################################################################
####COMBINE WITH PRIMARY PREVENTION PLOT
total_comb <- plot_grid(prim_restart_plot, sec_restart_plot, ncol=1)
print(total_comb)
ggsave("prim_sec_restart_plot.pdf", plot = last_plot(), width=15, height=15, units = "in", dpi = 400, path = final_graphs)


###############################COX ASSUMPTION###############################
#haven::write_dta(prim_restart, paste0(datafiles, "prim_restart.dta"))
# 
# 
# 
# library("survival")
# #Adjusted Cox model
# prim_model <- coxph(Surv(time, restart) ~ gender + ethnicity + deprivation,
#                     data=prim_restart)
# prim_model
# 
# prim.ph <- cox.zph(prim_model)
# prim.ph
# 
# # chisq df       p
# # gender       1.04  1    0.31
# # ethnicity   67.72  4 6.9e-14
# # deprivation  2.62  4    0.62
# # GLOBAL      80.80  9 1.1e-13
# 
# ggcoxzph(prim.ph)
# 
# 
# #survival plots
# restart.cox <- coxph(Surv(time, restart)~as.factor(ethnicity),data=prim_restart)
# summary(restart.cox)
# restart.survfit=survfit(restart.cox,newdata=data.frame(ethnicity=c("Black","South Asian", "Mixed", "Other")))
# summary(restart.survfit)
# plot(restart.survfit,mark.time=F,col=c("black","grey"),xlab="Time",ylab="Estimated survivor function")
# 
# library(survival)
# #log log plots
# prim_cox <- coxph(Surv(time, restart) ~ ethnicity,
#                   data=prim_restart)
# prim_cox
# # Plot the log-log survival curve
# # Log-Log Survival Curve with Kaplan-Meier
# 
# km=survfit(Surv(time,restart)~ethnicity,data=prim_restart)
# 
# prim.km <- survfit(Surv(time,restart)~as.factor(ethnicity),data=prim_restart)
# ggsurvplot(prim.km, data = prim_restart,conf.int = T,censor=F,legend.title="",
#            legend.labs = c("White","Black", "South Asian", "Mixed", "Other"))
# 
# plot(prim.km,fun="cloglog",xlab="time (log scale)",ylab="log(-log S(t))",
#      col=c("blue","red", "green", "orange","yellow"))
# legend(0.02,0,c("White","Black", "South Asian", "Mixed", "Other"),col=c("blue","red", "green", "orange", "yellow"),lty=1,cex=0.5)
# 
# ggsurvplot(prim.km, data = prim_restart,conf.int = T,fun="cloglog",censor=F,legend.title="",
#            legend.labs = c("White","Black", "South Asian", "Mixed", "Other"))
# 
# ggsurvplot(prim.km, data = prim_restart,conf.int = F,fun="cloglog",loglog =T, #xlim = c(1, 12),
#            legend.title="",
#            legend.labs = c("White","Black", "South Asian", "Mixed", "Other"))
# 
# 
# #Interaction term between ethnicity and time
# prim.split=survSplit(Surv(time,restart)~., data=prim_restart, cut=6, end="time", event="restart", start="time0", episode="time_period")
# 
# prim.cox.t2new=coxph(Surv(time0,time,restart)~as.factor(ethnicity)*as.factor(time_period),data=prim.split)
# summary(prim.cox.t2new)
# 
# 
# 
# ########################SECONDARY PREVENTION
# sec_model <- coxph(Surv(time, restart) ~ gender + ethnicity + deprivation,
#                    data=sec_restart)
# sec_model
# 
# sec_model <- coxph(Surv(time, restart) ~ gender,
#                    data=sec_restart)
# 
# 
# sec.ph <- cox.zph(sec_model)
# sec.ph
# 
# # chisq df       p
# # gender       90.67  1 < 2e-16
# # ethnicity    19.42  4 0.00065
# # deprivation   1.26  4 0.86793
# # GLOBAL      109.35  9 < 2e-16
# 
# #gender
# km=survfit(Surv(time,restart)~gender,data=sec_restart)
# 
# sec.km <- survfit(Surv(time, restart)~as.factor(gender),data=sec_restart)
# ggsurvplot(sec.km, data = sec_restart,conf.int = T,censor=F,legend.title="",
#            legend.labs = c("Female","Male"))
# 
# plot(sec.km,fun="cloglog",xlab="time (log scale)",ylab="log(-log S(t))",
#      col=c("red","blue"))
# legend(0.02,0,c("Female","Male"),col=c("red","blue"),lty=1,cex=0.5)
# 
# ggsurvplot(sec.km, data = sec_restart,conf.int = T,fun="cloglog",censor=F,legend.title="",
#            legend.labs = c("Female","Male"))
# 
# 
# 
# km=survfit(Surv(time,restart)~ethnicity,data=sec_restart)
# 
# sec.km <- survfit(Surv(time,restart)~as.factor(ethnicity),data=sec_restart)
# ggsurvplot(sec.km, data = sec_restart,conf.int = T,censor=F,legend.title="",
#            legend.labs = c("White","Black", "South Asian", "Mixed", "Other"))
# 
# plot(sec.km,fun="cloglog",xlab="time (log scale)",ylab="log(-log S(t))",
#      col=c("blue","red", "green", "orange","yellow"))
# legend(0.02,0,c("White","Black", "South Asian", "Mixed", "Other"),col=c("blue","red", "green", "orange", "yellow"),lty=1,cex=0.5)
# 
# ggsurvplot(sec.km, data = sec_restart,conf.int = T,fun="cloglog",censor=F,legend.title="",
#            legend.labs = c("White","Black", "South Asian", "Mixed", "Other"))
# 
# 
