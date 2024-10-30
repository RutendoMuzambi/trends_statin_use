##########################################################################################################
###################FIGURE x: FACTORS ASSOCIATED WITH DISCONTINUATION
#############################PRIMARY AND SECONDARY PREVENTION##################################
#############################################################################################
rm(list=ls())

library(broom)
library(openxlsx)
library(cowplot)

lapply(c("ggpubr", "grid", "gridExtra", "forcats"), require, character.only=T)

##open parquet
prim_disc <- read_parquet(paste0(datafiles, "prim_discont_fact.parquet"))

##################################################################################################
###########################################PRIMARY PREVENTION#####################################
##################################################################################################

prim_disc <- prim_disc[ethnicity!="missing" & deprivation!="missing"]

prim_disc[, c("age_cat", "ethnicity", "gender", "deprivation") := lapply(.SD, as.factor), .SDcols = c("age_cat", "ethnicity", "gender", "deprivation")]
prim_disc[ethnicity=="missing", ethnicity := NA]
prim_disc[deprivation=="missing", deprivation := NA]

prim_disc[, age_cat := factor(age_cat, levels = c("25-39", "40-49", "50-59", "60-69", "70+"))]
prim_disc[, gender := factor(gender, levels = c("Male", "Female"))]
prim_disc[, ethnicity := factor(ethnicity, levels = c("White", "Black", "South Asian", "Mixed", "Other"))]
prim_disc[, deprivation := factor(deprivation, levels = c("1 - least deprived", "2", "3", "4", "5 - most deprived"))]

prim_disc[, .N, by="age_cat"]
prim_disc[, .N, by="ethnicity"]
prim_disc[, .N, by="gender"]
prim_disc[, .N, by="deprivation"]

cox_disc_young<- prim_disc[age_cat=="25-39"]
cox_disc_white<- prim_disc[ethnicity=="White"]
cox_disc_male<- prim_disc[gender=="Male"]
cox_disc_depri<- prim_disc[deprivation=="1 - least deprived"]

N_young_adj <- nrow(cox_disc_young)
N_white_adj <- nrow(cox_disc_white)
N_male_adj <- nrow(cox_disc_male)
N_depri_adj <- nrow(cox_disc_depri)

#Adjusted Cox model
cox_model1 <- coxph(Surv(time, discontinuation) ~ age_cat + gender + ethnicity + deprivation,
                    data=prim_disc)
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
pyears_age <- pyears(Surv(time, discontinuation) ~ age_cat, data = prim_disc, model = T)
pyears_age <- tidy(pyears_age)
colnames(pyears_age)[colnames(pyears_age) == "n"] <- "N"

#pyears_age$term <- c('25-39', 'age_cat40-49', 'age_cat50-59', 'age_cat60-69', 'age_cat70+')

pyears_sex <- pyears(Surv(time, discontinuation) ~ gender, data = prim_disc, model = T)
pyears_sex <- tidy(pyears_sex)
colnames(pyears_sex)[colnames(pyears_sex) == "n"] <- "N"

pyears_eth <- pyears(Surv(time, discontinuation) ~ ethnicity, data = prim_disc, model = T)
pyears_eth <- tidy(pyears_eth)
colnames(pyears_eth)[colnames(pyears_eth) == "n"] <- "N"

pyears_depr <- pyears(Surv(time, discontinuation) ~ deprivation, data = prim_disc, model = T)
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

prim_disc_plot <- plot_grid(prim_left_plot, prim_main_plot, prim_right_plot, labels=c('A. Primary Prevention'), align="center", label_size = 12, label_y=0.99, ncol=3)
#ggsave("prim_discont_prim.png", plot = last_plot(), width = 12, height = 8, units = "in", dpi = 300, path = graphs)


##################################################################################################
#########################################SECONDARY PREVENTION#####################################
##################################################################################################
##open parquet
sec_disc <- read_parquet(paste0(datafiles, "sec_discont_fact.parquet"))

uniqueN(sec_disc, "patid")
sec_disc[, .N, by="CVD"]

sec_disc <- sec_disc[ethnicity!="missing" & deprivation!="missing"]

sec_disc[, c("ethnicity", "gender", "deprivation", "age_cat") := lapply(.SD, as.factor), .SDcols = c("ethnicity", "gender", "deprivation", "age_cat")]
sec_disc[ethnicity=="missing", ethnicity := NA]
sec_disc[deprivation=="missing", deprivation := NA]

sec_disc[, age_cat := factor(age_cat, levels = c("25-39", "40-49", "50-59", "60-69", "70+"))]
sec_disc[, gender := factor(gender, levels = c("Male", "Female"))]
sec_disc[, ethnicity := factor(ethnicity, levels = c("White", "Black", "South Asian", "Mixed", "Other"))]
sec_disc[, deprivation := factor(deprivation, levels = c("1 - least deprived", "2", "3", "4", "5 - most deprived"))]

#sec_disc <- sec_disc[CVD!="Acute Coronary Syndrome (non-specific)" & CVD!="Coronary Heart Disease (non-specific)"]
sec_disc[, .N, by="age_cat"]
sec_disc[, .N, by="ethnicity"]
sec_disc[, .N, by="gender"]
sec_disc[, .N, by="deprivation"]
sec_disc[, .N, by="CVD"]
cox_disc_young<- sec_disc[age_cat=="25-39"]
cox_disc_white<- sec_disc[ethnicity=="White"]
cox_disc_male<- sec_disc[gender=="Male"]
cox_disc_depri<- sec_disc[deprivation=="1 - least deprived"]

N_young_adj <- nrow(cox_disc_young)
N_white_adj <- nrow(cox_disc_white)
N_male_adj <- nrow(cox_disc_male)
N_depri_adj <- nrow(cox_disc_depri)

#Adjusted Cox model
cox_model1 <- coxph(Surv(time, discontinuation) ~ age_cat + gender + ethnicity + deprivation,
                    data=sec_disc)
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
pyears_age <- pyears(Surv(time, discontinuation) ~ age_cat, data = sec_disc, model = T)
pyears_age <- tidy(pyears_age)
colnames(pyears_age)[colnames(pyears_age) == "n"] <- "N"

#pyears_age$term <- c('25-39', 'age_cat40-49', 'age_cat50-59', 'age_cat60-69', 'age_cat70+')

pyears_sex <- pyears(Surv(time, discontinuation) ~ gender, data = sec_disc, model = T)
pyears_sex <- tidy(pyears_sex)
colnames(pyears_sex)[colnames(pyears_sex) == "n"] <- "N"

pyears_eth <- pyears(Surv(time, discontinuation) ~ ethnicity, data = sec_disc, model = T)
pyears_eth <- tidy(pyears_eth)
colnames(pyears_eth)[colnames(pyears_eth) == "n"] <- "N"

pyears_depr <- pyears(Surv(time, discontinuation) ~ deprivation, data = sec_disc, model = T)
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


sec_disc_plot <- plot_grid(left_plot, main_plot, right_plot, labels=c('B. Secondary Prevention'), align="center", label_size = 12, label_y=0.99, ncol=3)

#########################################################################################
####COMBINE WITH PRIMARY PREVENTION PLOT
total_comb <- plot_grid(prim_disc_plot, sec_disc_plot, ncol=1)
print(total_comb)
ggsave("prim_sec_disc_plot.pdf", plot = last_plot(), width=15, height=15, units = "in", dpi = 400, path = final_graphs)

