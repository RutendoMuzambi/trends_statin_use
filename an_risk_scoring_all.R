
##########################################################################################################
###################FIGURE x: FACTORS ASSOCIATED WITH RISK SCORING##################################
#############################################################################################
rm(list=ls())

library("broom")
library("openxlsx")
library("cowplot")

lapply(c("ggpubr", "grid", "gridExtra", "forcats"), require, character.only=T)

##open parqeut
cox_risk_dt <- read_parquet(paste0(datafiles, "risk_scoring_fact.parquet"))

##################################################################################################
##################################################################################################
cox_risk_dt <- cox_risk_dt[ethnicity!="missing" & deprivation!="missing"]

cox_risk_dt[, .N, by="age_cat"]
cox_risk_dt[, .N, by="ethnicity"]
cox_risk_dt[, .N, by="gender"]
cox_risk_dt[, .N, by="deprivation"]
cox_risk_forty<- cox_risk_dt[age_cat=="40-49"]
cox_risk_white<- cox_risk_dt[ethnicity=="White"]
cox_risk_male<- cox_risk_dt[gender=="Male"]
cox_risk_depri<- cox_risk_dt[deprivation=="1 - least deprived"]

N_forty_adj <- nrow(cox_risk_forty)
N_white_adj <- nrow(cox_risk_white)
N_male_adj <- nrow(cox_risk_male)
N_depri_adj <- nrow(cox_risk_depri)

#Adjusted Cox model
cox_model1 <- coxph(Surv(time, risk_scoring) ~ age_cat + gender + ethnicity + deprivation,
                    data=cox_risk_dt)
summary(cox_model1)

model.matrix(cox_model1)

adjusted_df <- tidy(cox_model1, conf.int=T, exponentiate=T)

#add rows for reference groups
adjusted_df <- adjusted_df %>% add_row(term = "40-49", estimate = 1, conf.low=1,conf.high=1)
adjusted_df <- adjusted_df %>% add_row(term = "White", estimate = 1, conf.low=1,conf.high=1)
adjusted_df <- adjusted_df %>% add_row(term = "1 - least deprived", estimate = 1, conf.low=1,conf.high=1)
adjusted_df <- adjusted_df %>% add_row(term = "Men", estimate = 1, conf.low=1,conf.high=1)

N <- nobs(cox_model1)

adj_n <-  cox_model1 %>% model.matrix() %>% 
  colSums()%>%                           # show the total number in each level
  as.data.frame() %>% t %>% as.data.frame() %>% 
  dplyr::select(c("age_cat50-59", "age_cat60-69", "age_cat70-74", "genderFemale", "ethnicityBlack", "ethnicitySouth Asian",
                  "ethnicityMixed", "ethnicityOther", "deprivation2", "deprivation3",
                  "deprivation4", "deprivation5 - most deprived"))

adjusted_df <- adjusted_df %>% mutate(N=case_when(
  term == "age_cat50-59" ~ adj_n$age_cat50-59,
  term == "age_cat60-69" ~ adj_n$age_cat60-69,
  term == "age_cat70-74" ~ adj_n$age_cat70-74,
  term == "genderFemale" ~ adj_n$genderFemale,
  term == "ethnicityBlack" ~ adj_n$ethnicityBlack,
  term == "ethnicitySouth Asian" ~ adj_n$"ethnicitySouth Asian",
  term == "ethnicityMixed" ~ adj_n$ethnicityMixed,
  term == "ethnicityOther" ~ adj_n$ethnicityOther,
  term == "deprivation2" ~ adj_n$deprivation2,
  term == "deprivation3" ~ adj_n$deprivation3,
  term == "deprivation4" ~ adj_n$deprivation4,
  term == "deprivation5 - most deprived" ~ adj_n$"deprivation5 - most deprived"
))

rm(adj_n)

##add N of reference groups
adjusted_df <- as.data.table(adjusted_df)

adjusted_df[term=="40-49" & is.na(N), N := N_forty_adj]
adjusted_df[term=="White" & is.na(N), N := N_white_adj]
adjusted_df[term=="1 - least deprived" & is.na(N), N := N_depri_adj]
adjusted_df[term=="Men" & is.na(N), N := N_male_adj]

cols_round <- c("estimate", "conf.low", "conf.high")

for (col_name in cols_round) {
  adjusted_df[[col_name]] <- round(adjusted_df[[col_name]], 2)
}

#adjusted_df$p.value <- round(adjusted_df$p.value, 3)
#adjusted_df <- adjusted_df[1:9, ]

adjusted_df <- adjusted_df %>%
  mutate(HR = paste(adjusted_df$estimate," (",adjusted_df$conf.low," - ",
                    adjusted_df$conf.high,")",sep = ""))

#reference categories HR
adjusted_df[term=="40-49", HR := "1 (reference)"]
adjusted_df[term=="White", HR := "1 (reference)"]
adjusted_df[term=="1 - least deprived", HR := "1 (reference)"]
adjusted_df[term=="Men", HR := "1 (reference)"]

##########FOREST PLOT########
adjusted_df[term=="age_cat50-59", term := "50-59"]
adjusted_df[term=="age_cat60-69", term := "60-69"]
adjusted_df[term=="age_cat70-74", term := "70-74"]
adjusted_df[term=="genderFemale", term := "Women"]
adjusted_df[term=="ethnicitySouth Asian", term := "South Asian"]
adjusted_df[term=="ethnicityBlack", term := "Black"]
adjusted_df[term=="ethnicityMixed", term := "Mixed"]
adjusted_df[term=="ethnicityOther", term := "Other"]
adjusted_df[term=="ethnicityWhite", term := "White"]
adjusted_df[term=="deprivation2", term := "2"]
adjusted_df[term=="deprivation3", term := "3"]
adjusted_df[term=="deprivation4", term := "4"]
adjusted_df[term=="deprivation5 - most deprived", term := "5 - most deprived"]

adjusted_df[term=="1 - least deprived"|term=="2"|term=="3"|term=="4"|term=="5 - most deprived",
            deprivation := term]
adjusted_df[term=="White"|term=="Black"|term=="Mixed"|term=="Other"|term=="South Asian",
            ethnicity := term]
adjusted_df[term=="Men"|term=="Women", gender := term]
adjusted_df[term=="40-49"|term=="50-59"|term=="60-69"|term=="70-74",
            age_cat := term]

adjusted_df[!is.na(age_cat), Factor := "Age Group"]
adjusted_df[!is.na(gender), Factor := "Gender"]
adjusted_df[!is.na(deprivation), Factor := "Deprivation"]
adjusted_df[!is.na(ethnicity), Factor := "Ethnicity"]

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
plot_data[term=="40-49" |term=="White" | term=="1 - least deprived" | term=="Men", ref := 1] 
plot_data[term!="40-49" & term!="White"&term!="1 - least deprived"&term!="Men", ref := 0] 

setorder(plot_data, factor_order, -ref)

plot_data$order <- 16:1
plot_data$order %>% as.factor
plot_data$order %>% levels()

#################################################################################
##############################PLOT###################################



main_plot <- ggplot(plot_data, aes(y=order, x=estimate, xmin=conf.low, xmax=conf.high, color=term)) +
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
plot_data[term!="White"&term!="1 - least deprived"&term!="Men"&term!="40-49", Factor := " "] 
plot_data[term == "1 - least deprived", term := "1 - least"]
plot_data[term == "5 - most deprived", term := "5 - most"]

left_plot <- ggplot(data = plot_data, aes(y = order)) + xlim(0,0.8) +
  geom_text(aes(x = 0, label = Factor ),lineheight = 0.001, hjust = 0, size = 4, colour = "black") + 
  geom_text(aes(x = 0.3, label = term),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") +
  geom_text(aes(x = 0.6, label = N),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") + 
  theme_classic2() + 
  ggtitle(" Category                                           Events")+   
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

right_plot <- ggplot(data = plot_data,  aes(y = order)) + xlim(0,0.5) +
  geom_text(aes(x = 0, label = HR),lineheight = 0.01, hjust = 0, size = 4, colour = "black") + 
  theme_classic2( ) +
  ggtitle("    HR (95%CI)")+
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


risk_plot <- plot_grid(left_plot, main_plot, right_plot, ncol=3)
ggsave("risk_scoring_plot_all.pdf", plot=last_plot(), width=10, height=8, units = "in", dpi = 400, path = final_graphs)


