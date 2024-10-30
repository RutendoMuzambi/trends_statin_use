##########################################################################################################
###################TABLE AND FIGURES: FACTORS ASSOCIATED WITH INITIATION###################
#####################################SECONDARY PREVENTION BY CVD SUBTYPE##################################
#############################################################################################
#rm(list=ls())
library(epiDisplay)
library(epitools)
library(finalfit)
library(openxlsx)
library(broom)
library(stringr)
library(zoo)

lapply(c("ggpubr", "grid", "gridExtra", "forcats"), require, character.only=T)

##open parquet
sec_init_dt <- read_parquet(paste0(datafiles, "sec_initiators.parquet"))

##format bmi
sec_init_dt[bmi_cat=="Normal weight (18.5-24.9)", bmi_cat := "Normal weight"]
sec_init_dt[bmi_cat=="Overweight (25.0-29.9)", bmi_cat := "Overweight"]
sec_init_dt[bmi_cat=="Obese (30.0>)", bmi_cat := "Obese"]
sec_init_dt[bmi_cat=="Underweight (<18.5)", bmi_cat := "Underweight"]

sec_init_dt[, bmi_cat := relevel(factor(bmi_cat), ref="Normal weight")]

sec_init_dt[, .N, by="bmi_cat"]
##########################MYOCARDIAL INFARCTION######################################
init_mi_dt <- sec_init_dt[CVD=="MI"]
uniqueN(init_mi_dt , by = "patid") #77,039
init_mi_dt[, .N, by=initiation]

mi_lg <- glm(initiation ~ age_cat + gender + ethnicity + deprivation + bmi_cat + smokingstatus + 
               treatedhyp + af_cat + ckd_cat + ra_cat + t2dm_cat, data = init_mi_dt, family = binomial())
logistic.display(mi_lg)

explanatory= c("age_cat", "gender", "ethnicity",  "deprivation", "bmi_cat", "smokingstatus",
               "treatedhyp","af_cat", "ckd_cat", "ra_cat", "t2dm_cat")
dependent="initiation"
init_mi_dt %>%
  finalfit(dependent, explanatory, confint_sep=" - ", confint_level=0.95) %>%
  ff_remove_p() -> mi_lg

########modify for table and plot
setnames(mi_lg, "Dependent: initiation", "Factor")
setnames(mi_lg, "yes", "events")
setnames(mi_lg, "OR (univariable)", "crude_OR")
setnames(mi_lg, "OR (multivariable)", "adjusted_OR")
colnames(mi_lg)[2] <- "term"

mi_lg <- as.data.table(mi_lg)

###extract values to different column
#crude
mi_lg[, crd_OR := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
mi_lg[, crd_lci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
mi_lg[, crd_uci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]


mi_lg[, adj_OR := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
mi_lg[, adj_lci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
mi_lg[, adj_uci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]
mi_lg[, events := str_extract(events, "(\\d+) \\((.*)\\)", group=1)]

#rename reference groups
mi_lg[crude_OR=="-", c("crd_OR", "crd_lci", "crd_uci") := list(1, 1, 1)]
mi_lg[crude_OR=="-", crude_OR := "1 (reference)"]

mi_lg[adjusted_OR=="-", c("adj_OR", "adj_lci", "adj_uci") := list(1, 1, 1)]
mi_lg[adjusted_OR=="-", adjusted_OR := "1 (reference)"]

# Remove 'missing' rows
mi_lg <- mi_lg[term != "missing"]
mi_lg <- mi_lg[term !="Normal weight (18.5-24.9)"& term !="Overweight (25.0-29.9)"&
                 term!="Obese (30.0>)"&term !="Underweight (<18.5)"]
mi_lg[term=="Normal weight", Factor := "BMI"]

################################################################################
###################################PLOT#########################################
################################################################################

mi_plot <- mi_lg
mi_plot[, factor2 := Factor]
mi_plot[factor2 == "", factor2 := NA]
mi_plot[, factor2 := na.locf(factor2)]
mi_plot[, factor_order := rleid(factor2)]

###order correctly
mi_plot[adjusted_OR=="1 (reference)", ref := 1] 
mi_plot[adjusted_OR!="1 (reference)", ref := 0]

setorder(mi_plot, factor_order, -ref)
mi_plot[term==1, term := "Yes"]
mi_plot[term==0, term := "No"]

mi_plot <- mi_plot %>%
  dplyr::select(factor_order, Factor, term, events, adj_OR, adj_lci, adj_uci, adjusted_OR, crd_OR, crd_lci, crd_uci, crude_OR )
mi_plot[,1:3] <- map(mi_plot[,1:3], as.factor)    
mi_plot[, order := .N:1]
str(mi_plot)

mi_plot[, adj_OR := as.numeric(adj_OR)]
mi_plot[, adj_lci := as.numeric(adj_lci)]
mi_plot[, adj_uci:= as.numeric(adj_uci)]

mi_plot[, crd_OR := as.numeric(crd_OR)]
mi_plot[, crd_lci := as.numeric(crd_lci)]
mi_plot[, crd_uci:= as.numeric(crd_uci)]

mi_long_plot <- mi_plot %>%
  pivot_longer(cols=c(crude_OR, crd_OR, crd_lci, crd_uci,
                      adjusted_OR, adj_OR, adj_lci, adj_uci),
               names_to=c(".value", "Model"),
               names_pattern="(\\w+)(\\w+)$")
crd_long_plot <- mi_plot[, c("Factor", "term", "crd_OR", "crude_OR", "crd_OR", "crd_lci", "crd_uci", "factor_order", "events", "order")]
crd_long_plot[, model := 1]

adj_long_plot <- mi_plot[, c("Factor", "term", "adjusted_OR", "adj_OR", "adj_lci", "adj_uci", "factor_order", "events", "order")]
adj_long_plot[, model := 2]

final_mi_plot <- rbind(crd_long_plot, adj_long_plot, fill=T)

final_mi_plot[model==1, OR := crd_OR]
final_mi_plot[model==2, OR := adj_OR]
final_mi_plot[model==1, lci := crd_lci]
final_mi_plot[model==2, lci := adj_lci]
final_mi_plot[model==1, uci := crd_uci]
final_mi_plot[model==2, uci := adj_uci]
final_mi_plot[model==1, OR_full := crude_OR]
final_mi_plot[model==2, OR_full := adjusted_OR]

final_mi_plot <- final_mi_plot[, c("Factor", "term", "OR", "lci", "uci", "OR_full", "events", "factor_order", "order", "model")]

setorder(final_mi_plot, -order, model)

final_mi_plot <- final_mi_plot[(OR==1.00 & model==1) | OR!=1.00]
final_mi_plot[model==2, term :=" "] 
final_mi_plot[model==2, events :=" "] 
final_mi_plot[, order_2 := .N:1]

final_mi_plot[model==1, Model := "Crude"]
final_mi_plot[model==2, Model := "Adjusted"]

main_mi_plot <- ggplot(final_mi_plot, aes(y=order_2, x=OR, xmin=lci, xmax=uci, color=Model)) +
  geom_point() +
  geom_errorbarh(height=.3) +
  geom_vline(xintercept=1, lty=2) +
  ggtitle(" ") +
  theme_classic2( ) +
  theme(strip.text.y = element_blank(),
        #legend.position="bottom",
        #legend.justification = "right",
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  #guides(color=guide_legend(reverse=TRUE)) +
  scale_color_manual(breaks = c("Crude", "Adjusted"),
                     values=c("skyblue", "navy"),
                     labels=c("0"="crude", "1"="adjusted")) +
  xlab("Odds ratio") 

#rename for plot
final_mi_plot[Factor == "age_cat", Factor := "Age group"]
final_mi_plot[Factor == "ethnicity", Factor := "Ethnicity"]
final_mi_plot[Factor == "gender", Factor := "Gender"]
final_mi_plot[Factor == "deprivation", Factor := "Deprivation"]
final_mi_plot[Factor == "smokingstatus", Factor := "Smoking status"]
final_mi_plot[Factor == "treatedhyp", Factor := "Treated hypertension"]
final_mi_plot[Factor == "af_cat", Factor := "Atrial Fibrillation"]
final_mi_plot[Factor == "ckd_cat", Factor := "CKD"]
final_mi_plot[Factor == "ra_cat", Factor := "RA"]
final_mi_plot[Factor == "t2dm_cat", Factor := "Type II diabetes"]
final_mi_plot[Factor == "BMI" & term==" ", Factor := " "]
final_mi_plot[term == "Male", term := "Men"]
final_mi_plot[term == "Female", term := "Women"]

left_mi_plot <- ggplot(data = final_mi_plot, aes(y = order_2)) + xlim(0,0.8) +
  geom_text(aes(x = 0, label = Factor ),lineheight = 0.001, hjust = 0, size = 4, colour = "black") + 
  geom_text(aes(x = 0.37, label = term),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") +
  geom_text(aes(x = 0.70, label = events),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") + 
  theme_classic2() + 
  ggtitle("     Factor                                                      Events")+   
  theme(plot.title = element_text(size = 12, lineheight=.001, face="bold", vjust=-13),
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

###Plot OR and 95% CI
right_mi_plot <- ggplot(data = final_mi_plot,  aes(y = order_2)) + xlim(0,0.5) +
  geom_text(aes(x = 0, label = OR_full),lineheight = 0.01, hjust = 0, size = 4, colour = "black") + 
  theme_classic2( ) +
  ggtitle("  OR (95% CI)")+
  theme(plot.title = element_text(size = 12 ,lineheight=.001, face="bold", vjust=-13),
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

library(cowplot)

init_mi_plot <- plot_grid(left_mi_plot, main_mi_plot, right_mi_plot, labels=c('MI'), label_size =12,   ncol=3)
ggsave("sec_init_mi_plot.png", plot=init_mi_plot, width=10, height=12, path = graphs)

##########################STROKE######################################
init_stroke_dt <- sec_init_dt[CVD=="Stroke"]
uniqueN(init_stroke_dt , by = "patid") #77,039
init_stroke_dt[, .N, by=initiation]

stroke_lg <- glm(initiation ~ age_cat + gender + ethnicity + deprivation + bmi_cat + smokingstatus + 
                   treatedhyp + af_cat + ckd_cat + ra_cat + t2dm_cat, data = init_stroke_dt, family = binomial())
logistic.display(stroke_lg)

explanatory= c("age_cat", "gender", "ethnicity",  "deprivation", "bmi_cat", "smokingstatus",
               "treatedhyp","af_cat", "ckd_cat", "ra_cat", "t2dm_cat")
dependent="initiation"
init_stroke_dt %>%
  finalfit(dependent, explanatory, confint_sep=" - ", confint_level=0.95) %>%
  ff_remove_p() -> stroke_lg

########modify for table and plot
setnames(stroke_lg, "Dependent: initiation", "Factor")
setnames(stroke_lg, "yes", "events")
setnames(stroke_lg, "OR (univariable)", "crude_OR")
setnames(stroke_lg, "OR (multivariable)", "adjusted_OR")
colnames(stroke_lg)[2] <- "term"

stroke_lg <- as.data.table(stroke_lg)

###extract values to different column
#crude
stroke_lg[, crd_OR := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
stroke_lg[, crd_lci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
stroke_lg[, crd_uci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]


stroke_lg[, adj_OR := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
stroke_lg[, adj_lci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
stroke_lg[, adj_uci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]
stroke_lg[, events := str_extract(events, "(\\d+) \\((.*)\\)", group=1)]

#rename reference groups
stroke_lg[crude_OR=="-", c("crd_OR", "crd_lci", "crd_uci") := list(1, 1, 1)]
stroke_lg[crude_OR=="-", crude_OR := "1 (reference)"]

stroke_lg[adjusted_OR=="-", c("adj_OR", "adj_lci", "adj_uci") := list(1, 1, 1)]
stroke_lg[adjusted_OR=="-", adjusted_OR := "1 (reference)"]

# Remove 'missing' rows
stroke_lg <- stroke_lg[term != "missing"]
stroke_lg <- stroke_lg[term !="Normal weight (18.5-24.9)"& term !="Overweight (25.0-29.9)"&
                         term!="Obese (30.0>)"&term !="Underweight (<18.5)"]
stroke_lg[term=="Normal weight", Factor := "BMI"]

################################################################################
###################################PLOT#########################################
################################################################################

stroke_plot <- stroke_lg
stroke_plot[, factor2 := Factor]
stroke_plot[factor2 == "", factor2 := NA]
stroke_plot[, factor2 := na.locf(factor2)]
stroke_plot[, factor_order := rleid(factor2)]

###order correctly
stroke_plot[adjusted_OR=="1 (reference)", ref := 1] 
stroke_plot[adjusted_OR!="1 (reference)", ref := 0]

setorder(stroke_plot, factor_order, -ref)
stroke_plot[term==1, term := "Yes"]
stroke_plot[term==0, term := "No"]

stroke_plot <- stroke_plot %>%
  dplyr::select(factor_order, Factor, term, events, adj_OR, adj_lci, adj_uci, adjusted_OR, crd_OR, crd_lci, crd_uci, crude_OR )
stroke_plot[,1:3] <- map(stroke_plot[,1:3], as.factor)    
stroke_plot[, order := .N:1]
str(stroke_plot)

stroke_plot[, adj_OR := as.numeric(adj_OR)]
stroke_plot[, adj_lci := as.numeric(adj_lci)]
stroke_plot[, adj_uci:= as.numeric(adj_uci)]

stroke_plot[, crd_OR := as.numeric(crd_OR)]
stroke_plot[, crd_lci := as.numeric(crd_lci)]
stroke_plot[, crd_uci:= as.numeric(crd_uci)]

stroke_long_plot <- stroke_plot %>%
  pivot_longer(cols=c(crude_OR, crd_OR, crd_lci, crd_uci,
                      adjusted_OR, adj_OR, adj_lci, adj_uci),
               names_to=c(".value", "Model"),
               names_pattern="(\\w+)(\\w+)$")
crd_long_plot <- stroke_plot[, c("Factor", "term", "crd_OR", "crude_OR", "crd_OR", "crd_lci", "crd_uci", "factor_order", "events", "order")]
crd_long_plot[, model := 1]

adj_long_plot <- stroke_plot[, c("Factor", "term", "adjusted_OR", "adj_OR", "adj_lci", "adj_uci", "factor_order", "events", "order")]
adj_long_plot[, model := 2]

final_stroke_plot <- rbind(crd_long_plot, adj_long_plot, fill=T)

final_stroke_plot[model==1, OR := crd_OR]
final_stroke_plot[model==2, OR := adj_OR]
final_stroke_plot[model==1, lci := crd_lci]
final_stroke_plot[model==2, lci := adj_lci]
final_stroke_plot[model==1, uci := crd_uci]
final_stroke_plot[model==2, uci := adj_uci]
final_stroke_plot[model==1, OR_full := crude_OR]
final_stroke_plot[model==2, OR_full := adjusted_OR]

final_stroke_plot <- final_stroke_plot[, c("Factor", "term", "OR", "lci", "uci", "OR_full", "events", "factor_order", "order", "model")]

setorder(final_stroke_plot, -order, model)

final_stroke_plot <- final_stroke_plot[(OR==1.00 & model==1) | OR!=1.00]
final_stroke_plot[model==2, term :=" "] 
final_stroke_plot[model==2, events :=" "] 
final_stroke_plot[, order_2 := .N:1]

final_stroke_plot[model==1, Model := "Crude"]
final_stroke_plot[model==2, Model := "Adjusted"]

main_stroke_plot <- ggplot(final_stroke_plot, aes(y=order_2, x=OR, xmin=lci, xmax=uci, color=Model)) +
  geom_point() +
  geom_errorbarh(height=.3) +
  geom_vline(xintercept=1, lty=2) +
  ggtitle(" ") +
  theme_classic2( ) +
  theme(strip.text.y = element_blank(),
        #legend.position="bottom",
        #legend.justification = "right",
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  #guides(color=guide_legend(reverse=TRUE)) +
  scale_color_manual(breaks = c("Crude", "Adjusted"),
                     values=c("skyblue", "navy"),
                     labels=c("0"="crude", "1"="adjusted")) +
  xlab("Odds ratio") 

#rename for plot
final_stroke_plot[Factor == "age_cat", Factor := "Age group"]
final_stroke_plot[Factor == "ethnicity", Factor := "Ethnicity"]
final_stroke_plot[Factor == "gender", Factor := "Gender"]
final_stroke_plot[Factor == "deprivation", Factor := "Deprivation"]
final_stroke_plot[Factor == "smokingstatus", Factor := "Smoking status"]
final_stroke_plot[Factor == "treatedhyp", Factor := "Treated hypertension"]
final_stroke_plot[Factor == "af_cat", Factor := "Atrial Fibrillation"]
final_stroke_plot[Factor == "ckd_cat", Factor := "CKD"]
final_stroke_plot[Factor == "ra_cat", Factor := "RA"]
final_stroke_plot[Factor == "t2dm_cat", Factor := "Type II diabetes"]
final_stroke_plot[Factor == "BMI" & term==" ", Factor := " "]
final_stroke_plot[term == "Male", term := "Men"]
final_stroke_plot[term == "Female", term := "Women"]

left_stroke_plot <- ggplot(data = final_stroke_plot, aes(y = order_2)) + xlim(0,0.8) +
  geom_text(aes(x = 0, label = Factor ),lineheight = 0.001, hjust = 0, size = 4, colour = "black") + 
  geom_text(aes(x = 0.37, label = term),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") +
  geom_text(aes(x = 0.70, label = events),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") + 
  theme_classic2() + 
  ggtitle("     Factor                                      Events")+   
  theme(plot.title = element_text(size = 12, lineheight=.001, face="bold", vjust=-13),
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

###Plot OR and 95% CI
right_stroke_plot <- ggplot(data = final_stroke_plot,  aes(y = order_2)) + xlim(0,0.5) +
  geom_text(aes(x = 0, label = OR_full),lineheight = 0.01, hjust = 0, size = 4, colour = "black") + 
  theme_classic2( ) +
  ggtitle("  OR (95% CI)")+
  theme(plot.title = element_text(size = 12 ,lineheight=.001, face="bold", vjust=-13),
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

library(cowplot)

init_stroke_plot <- plot_grid(left_stroke_plot, main_stroke_plot, right_stroke_plot, labels=c('Stroke'), label_size =12,   ncol=3)
ggsave("sec_init_stroke_plot.png", plot=init_stroke_plot, width=10, height=12, path = graphs)


##########################ANGINA######################################
init_angina_dt <- sec_init_dt[CVD=="Angina"]
uniqueN(init_angina_dt , by = "patid") #77,039
init_angina_dt[, .N, by=initiation]

angina_lg <- glm(initiation ~ age_cat + gender + ethnicity + deprivation + bmi_cat + smokingstatus + 
                   treatedhyp + af_cat + ckd_cat + ra_cat + t2dm_cat, data = init_angina_dt, family = binomial())
logistic.display(angina_lg)

explanatory= c("age_cat", "gender", "ethnicity",  "deprivation", "bmi_cat", "smokingstatus",
               "treatedhyp","af_cat", "ckd_cat", "ra_cat", "t2dm_cat")
dependent="initiation"
init_angina_dt %>%
  finalfit(dependent, explanatory, confint_sep=" - ", confint_level=0.95) %>%
  ff_remove_p() -> angina_lg

########modify for table and plot
setnames(angina_lg, "Dependent: initiation", "Factor")
setnames(angina_lg, "yes", "events")
setnames(angina_lg, "OR (univariable)", "crude_OR")
setnames(angina_lg, "OR (multivariable)", "adjusted_OR")
colnames(angina_lg)[2] <- "term"

angina_lg <- as.data.table(angina_lg)

###extract values to different column
#crude
angina_lg[, crd_OR := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
angina_lg[, crd_lci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
angina_lg[, crd_uci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]


angina_lg[, adj_OR := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
angina_lg[, adj_lci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
angina_lg[, adj_uci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]
angina_lg[, events := str_extract(events, "(\\d+) \\((.*)\\)", group=1)]

#rename reference groups
angina_lg[crude_OR=="-", c("crd_OR", "crd_lci", "crd_uci") := list(1, 1, 1)]
angina_lg[crude_OR=="-", crude_OR := "1 (reference)"]

angina_lg[adjusted_OR=="-", c("adj_OR", "adj_lci", "adj_uci") := list(1, 1, 1)]
angina_lg[adjusted_OR=="-", adjusted_OR := "1 (reference)"]

# Remove 'missing' rows
angina_lg <- angina_lg[term != "missing"]
angina_lg <- angina_lg[term !="Normal weight (18.5-24.9)"& term !="Overweight (25.0-29.9)"&
                         term!="Obese (30.0>)"&term !="Underweight (<18.5)"]
angina_lg[term=="Normal weight", Factor := "BMI"]

################################################################################
###################################PLOT#########################################
################################################################################

angina_plot <- angina_lg
angina_plot[, factor2 := Factor]
angina_plot[factor2 == "", factor2 := NA]
angina_plot[, factor2 := na.locf(factor2)]
angina_plot[, factor_order := rleid(factor2)]

###order correctly
angina_plot[adjusted_OR=="1 (reference)", ref := 1] 
angina_plot[adjusted_OR!="1 (reference)", ref := 0]

setorder(angina_plot, factor_order, -ref)
angina_plot[term==1, term := "Yes"]
angina_plot[term==0, term := "No"]

angina_plot <- angina_plot %>%
  dplyr::select(factor_order, Factor, term, events, adj_OR, adj_lci, adj_uci, adjusted_OR, crd_OR, crd_lci, crd_uci, crude_OR )
angina_plot[,1:3] <- map(angina_plot[,1:3], as.factor)    
angina_plot[, order := .N:1]
str(angina_plot)

angina_plot[, adj_OR := as.numeric(adj_OR)]
angina_plot[, adj_lci := as.numeric(adj_lci)]
angina_plot[, adj_uci:= as.numeric(adj_uci)]

angina_plot[, crd_OR := as.numeric(crd_OR)]
angina_plot[, crd_lci := as.numeric(crd_lci)]
angina_plot[, crd_uci:= as.numeric(crd_uci)]

angina_long_plot <- angina_plot %>%
  pivot_longer(cols=c(crude_OR, crd_OR, crd_lci, crd_uci,
                      adjusted_OR, adj_OR, adj_lci, adj_uci),
               names_to=c(".value", "Model"),
               names_pattern="(\\w+)(\\w+)$")
crd_long_plot <- angina_plot[, c("Factor", "term", "crd_OR", "crude_OR", "crd_OR", "crd_lci", "crd_uci", "factor_order", "events", "order")]
crd_long_plot[, model := 1]

adj_long_plot <- angina_plot[, c("Factor", "term", "adjusted_OR", "adj_OR", "adj_lci", "adj_uci", "factor_order", "events", "order")]
adj_long_plot[, model := 2]

final_angina_plot <- rbind(crd_long_plot, adj_long_plot, fill=T)

final_angina_plot[model==1, OR := crd_OR]
final_angina_plot[model==2, OR := adj_OR]
final_angina_plot[model==1, lci := crd_lci]
final_angina_plot[model==2, lci := adj_lci]
final_angina_plot[model==1, uci := crd_uci]
final_angina_plot[model==2, uci := adj_uci]
final_angina_plot[model==1, OR_full := crude_OR]
final_angina_plot[model==2, OR_full := adjusted_OR]

final_angina_plot <- final_angina_plot[, c("Factor", "term", "OR", "lci", "uci", "OR_full", "events", "factor_order", "order", "model")]

setorder(final_angina_plot, -order, model)

final_angina_plot <- final_angina_plot[(OR==1.00 & model==1) | OR!=1.00]
final_angina_plot[model==2, term :=" "] 
final_angina_plot[model==2, events :=" "] 
final_angina_plot[, order_2 := .N:1]

final_angina_plot[model==1, Model := "Crude"]
final_angina_plot[model==2, Model := "Adjusted"]

main_angina_plot <- ggplot(final_angina_plot, aes(y=order_2, x=OR, xmin=lci, xmax=uci, color=Model)) +
  geom_point() +
  geom_errorbarh(height=.3) +
  geom_vline(xintercept=1, lty=2) +
  ggtitle(" ") +
  theme_classic2( ) +
  theme(strip.text.y = element_blank(),
        #legend.position="bottom",
        #legend.justification = "right",
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  #guides(color=guide_legend(reverse=TRUE)) +
  scale_color_manual(breaks = c("Crude", "Adjusted"),
                     values=c("skyblue", "navy"),
                     labels=c("0"="crude", "1"="adjusted")) +
  xlab("Odds ratio") 

#rename for plot
final_angina_plot[Factor == "age_cat", Factor := "Age group"]
final_angina_plot[Factor == "ethnicity", Factor := "Ethnicity"]
final_angina_plot[Factor == "gender", Factor := "Gender"]
final_angina_plot[Factor == "deprivation", Factor := "Deprivation"]
final_angina_plot[Factor == "smokingstatus", Factor := "Smoking status"]
final_angina_plot[Factor == "treatedhyp", Factor := "Treated hypertension"]
final_angina_plot[Factor == "af_cat", Factor := "Atrial Fibrillation"]
final_angina_plot[Factor == "ckd_cat", Factor := "CKD"]
final_angina_plot[Factor == "ra_cat", Factor := "RA"]
final_angina_plot[Factor == "t2dm_cat", Factor := "Type II diabetes"]
final_angina_plot[Factor == "BMI" & term==" ", Factor := " "]
final_angina_plot[term == "Male", term := "Men"]
final_angina_plot[term == "Female", term := "Women"]

left_angina_plot <- ggplot(data = final_angina_plot, aes(y = order_2)) + xlim(0,0.8) +
  geom_text(aes(x = 0, label = Factor ),lineheight = 0.001, hjust = 0, size = 4, colour = "black") + 
  geom_text(aes(x = 0.37, label = term),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") +
  geom_text(aes(x = 0.70, label = events),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") + 
  theme_classic2() + 
  ggtitle("     Factor                                                      Events")+   
  theme(plot.title = element_text(size = 12, lineheight=.001, face="bold", vjust=-13),
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

###Plot OR and 95% CI
right_angina_plot <- ggplot(data = final_angina_plot,  aes(y = order_2)) + xlim(0,0.5) +
  geom_text(aes(x = 0, label = OR_full),lineheight = 0.01, hjust = 0, size = 4, colour = "black") + 
  theme_classic2( ) +
  ggtitle("  OR (95% CI)")+
  theme(plot.title = element_text(size = 12 ,lineheight=.001, face="bold", vjust=-13),
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

library(cowplot)

init_angina_plot <- plot_grid(left_angina_plot, main_angina_plot, right_angina_plot, labels=c('Angina'), label_size =12,   ncol=3)
ggsave("sec_init_angina_plot.png", plot=init_angina_plot, width=10, height=12, path = graphs)


##########################TRANSIENT ISCHAEMIC ATTACK######################################
init_tia_dt <- sec_init_dt[CVD=="TIA"]
uniqueN(init_tia_dt , by = "patid") #77,039
init_tia_dt[, .N, by=initiation]

tia_lg <- glm(initiation ~ age_cat + gender + ethnicity + deprivation + bmi_cat + smokingstatus + 
                treatedhyp + af_cat + ckd_cat + ra_cat + t2dm_cat, data = init_tia_dt, family = binomial())
logistic.display(tia_lg)

explanatory= c("age_cat", "gender", "ethnicity",  "deprivation", "bmi_cat", "smokingstatus",
               "treatedhyp","af_cat", "ckd_cat", "ra_cat", "t2dm_cat")
dependent="initiation"
init_tia_dt %>%
  finalfit(dependent, explanatory, confint_sep=" - ", confint_level=0.95) %>%
  ff_remove_p() -> tia_lg

########modify for table and plot
setnames(tia_lg, "Dependent: initiation", "Factor")
setnames(tia_lg, "yes", "events")
setnames(tia_lg, "OR (univariable)", "crude_OR")
setnames(tia_lg, "OR (multivariable)", "adjusted_OR")
colnames(tia_lg)[2] <- "term"

tia_lg <- as.data.table(tia_lg)

###extract values to different column
#crude
tia_lg[, crd_OR := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
tia_lg[, crd_lci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
tia_lg[, crd_uci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]


tia_lg[, adj_OR := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
tia_lg[, adj_lci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
tia_lg[, adj_uci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]
tia_lg[, events := str_extract(events, "(\\d+) \\((.*)\\)", group=1)]

#rename reference groups
tia_lg[crude_OR=="-", c("crd_OR", "crd_lci", "crd_uci") := list(1, 1, 1)]
tia_lg[crude_OR=="-", crude_OR := "1 (reference)"]

tia_lg[adjusted_OR=="-", c("adj_OR", "adj_lci", "adj_uci") := list(1, 1, 1)]
tia_lg[adjusted_OR=="-", adjusted_OR := "1 (reference)"]

# Remove 'missing' rows
tia_lg <- tia_lg[term != "missing"]
tia_lg <- tia_lg[term !="Normal weight (18.5-24.9)"& term !="Overweight (25.0-29.9)"&
                   term!="Obese (30.0>)"&term !="Underweight (<18.5)"]
tia_lg[term=="Normal weight", Factor := "BMI"]

################################################################################
###################################PLOT#########################################
################################################################################

tia_plot <- tia_lg
tia_plot[, factor2 := Factor]
tia_plot[factor2 == "", factor2 := NA]
tia_plot[, factor2 := na.locf(factor2)]
tia_plot[, factor_order := rleid(factor2)]

###order correctly
tia_plot[adjusted_OR=="1 (reference)", ref := 1] 
tia_plot[adjusted_OR!="1 (reference)", ref := 0]

setorder(tia_plot, factor_order, -ref)
tia_plot[term==1, term := "Yes"]
tia_plot[term==0, term := "No"]

tia_plot <- tia_plot %>%
  dplyr::select(factor_order, Factor, term, events, adj_OR, adj_lci, adj_uci, adjusted_OR, crd_OR, crd_lci, crd_uci, crude_OR )
tia_plot[,1:3] <- map(tia_plot[,1:3], as.factor)    
tia_plot[, order := .N:1]
str(tia_plot)

tia_plot[, adj_OR := as.numeric(adj_OR)]
tia_plot[, adj_lci := as.numeric(adj_lci)]
tia_plot[, adj_uci:= as.numeric(adj_uci)]

tia_plot[, crd_OR := as.numeric(crd_OR)]
tia_plot[, crd_lci := as.numeric(crd_lci)]
tia_plot[, crd_uci:= as.numeric(crd_uci)]

tia_long_plot <- tia_plot %>%
  pivot_longer(cols=c(crude_OR, crd_OR, crd_lci, crd_uci,
                      adjusted_OR, adj_OR, adj_lci, adj_uci),
               names_to=c(".value", "Model"),
               names_pattern="(\\w+)(\\w+)$")
crd_long_plot <- tia_plot[, c("Factor", "term", "crd_OR", "crude_OR", "crd_OR", "crd_lci", "crd_uci", "factor_order", "events", "order")]
crd_long_plot[, model := 1]

adj_long_plot <- tia_plot[, c("Factor", "term", "adjusted_OR", "adj_OR", "adj_lci", "adj_uci", "factor_order", "events", "order")]
adj_long_plot[, model := 2]

final_tia_plot <- rbind(crd_long_plot, adj_long_plot, fill=T)

final_tia_plot[model==1, OR := crd_OR]
final_tia_plot[model==2, OR := adj_OR]
final_tia_plot[model==1, lci := crd_lci]
final_tia_plot[model==2, lci := adj_lci]
final_tia_plot[model==1, uci := crd_uci]
final_tia_plot[model==2, uci := adj_uci]
final_tia_plot[model==1, OR_full := crude_OR]
final_tia_plot[model==2, OR_full := adjusted_OR]

final_tia_plot <- final_tia_plot[, c("Factor", "term", "OR", "lci", "uci", "OR_full", "events", "factor_order", "order", "model")]

setorder(final_tia_plot, -order, model)

final_tia_plot <- final_tia_plot[(OR==1.00 & model==1) | OR!=1.00]
final_tia_plot[model==2, term :=" "] 
final_tia_plot[model==2, events :=" "] 
final_tia_plot[, order_2 := .N:1]

final_tia_plot[model==1, Model := "Crude"]
final_tia_plot[model==2, Model := "Adjusted"]

main_tia_plot <- ggplot(final_tia_plot, aes(y=order_2, x=OR, xmin=lci, xmax=uci, color=Model)) +
  geom_point() +
  geom_errorbarh(height=.3) +
  geom_vline(xintercept=1, lty=2) +
  ggtitle(" ") +
  theme_classic2( ) +
  theme(strip.text.y = element_blank(),
        #legend.position="bottom",
        #legend.justification = "right",
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  #guides(color=guide_legend(reverse=TRUE)) +
  scale_color_manual(breaks = c("Crude", "Adjusted"),
                     values=c("skyblue", "navy"),
                     labels=c("0"="crude", "1"="adjusted")) +
  xlab("Odds ratio") 

#rename for plot
final_tia_plot[Factor == "age_cat", Factor := "Age group"]
final_tia_plot[Factor == "ethnicity", Factor := "Ethnicity"]
final_tia_plot[Factor == "gender", Factor := "Gender"]
final_tia_plot[Factor == "deprivation", Factor := "Deprivation"]
final_tia_plot[Factor == "smokingstatus", Factor := "Smoking status"]
final_tia_plot[Factor == "treatedhyp", Factor := "Treated hypertension"]
final_tia_plot[Factor == "af_cat", Factor := "Atrial Fibrillation"]
final_tia_plot[Factor == "ckd_cat", Factor := "CKD"]
final_tia_plot[Factor == "ra_cat", Factor := "RA"]
final_tia_plot[Factor == "t2dm_cat", Factor := "Type II diabetes"]
final_tia_plot[Factor == "BMI" & term==" ", Factor := " "]
final_tia_plot[term == "Male", term := "Men"]
final_tia_plot[term == "Female", term := "Women"]

left_tia_plot <- ggplot(data = final_tia_plot, aes(y = order_2)) + xlim(0,0.8) +
  geom_text(aes(x = 0, label = Factor ),lineheight = 0.001, hjust = 0, size = 4, colour = "black") + 
  geom_text(aes(x = 0.35, label = term),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") +
  geom_text(aes(x = 0.70, label = events),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") + 
  theme_classic2() + 
  ggtitle("     Factor                                                      Events")+   
  theme(plot.title = element_text(size = 12, lineheight=.001, face="bold",  vjust=-13),
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

###Plot OR and 95% CI
right_tia_plot <- ggplot(data = final_tia_plot,  aes(y = order_2)) + xlim(0,0.5) +
  geom_text(aes(x = 0, label = OR_full),lineheight = 0.01, hjust = 0, size = 4, colour = "black") + 
  theme_classic2( ) +
  ggtitle("  OR (95% CI)")+
  theme(plot.title = element_text(size = 12 ,lineheight=.001, face="bold", vjust=-13),
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

library(cowplot)

init_tia_plot <- plot_grid(left_tia_plot, main_tia_plot, right_tia_plot, labels=c('TIA'), label_size =12,   ncol=3)
ggsave("sec_init_tia_plot.png", plot=init_tia_plot, width=10, height=12, path = graphs)


##########################PERIPHERAL ARTERIAL DISEASE######################################
init_pad_dt <- sec_init_dt[CVD=="PAD"]
uniqueN(init_pad_dt , by = "patid") #77,039
init_pad_dt[, .N, by=initiation]

pad_lg <- glm(initiation ~ age_cat + gender + ethnicity + deprivation + bmi_cat + smokingstatus + 
                treatedhyp + af_cat + ckd_cat + ra_cat + t2dm_cat, data = init_pad_dt, family = binomial())
logistic.display(pad_lg)

explanatory= c("age_cat", "gender", "ethnicity",  "deprivation", "bmi_cat", "smokingstatus",
               "treatedhyp","af_cat", "ckd_cat", "ra_cat", "t2dm_cat")
dependent="initiation"
init_pad_dt %>%
  finalfit(dependent, explanatory, confint_sep=" - ", confint_level=0.95) %>%
  ff_remove_p() -> pad_lg

########modify for table and plot
setnames(pad_lg, "Dependent: initiation", "Factor")
setnames(pad_lg, "yes", "events")
setnames(pad_lg, "OR (univariable)", "crude_OR")
setnames(pad_lg, "OR (multivariable)", "adjusted_OR")
colnames(pad_lg)[2] <- "term"

pad_lg <- as.data.table(pad_lg)

###extract values to different column
#crude
pad_lg[, crd_OR := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
pad_lg[, crd_lci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
pad_lg[, crd_uci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]


pad_lg[, adj_OR := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
pad_lg[, adj_lci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
pad_lg[, adj_uci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]
pad_lg[, events := str_extract(events, "(\\d+) \\((.*)\\)", group=1)]

#rename reference groups
pad_lg[crude_OR=="-", c("crd_OR", "crd_lci", "crd_uci") := list(1, 1, 1)]
pad_lg[crude_OR=="-", crude_OR := "1 (reference)"]

pad_lg[adjusted_OR=="-", c("adj_OR", "adj_lci", "adj_uci") := list(1, 1, 1)]
pad_lg[adjusted_OR=="-", adjusted_OR := "1 (reference)"]

# Remove 'missing' rows
pad_lg <- pad_lg[term != "missing"]
pad_lg <- pad_lg[term !="Normal weight (18.5-24.9)"& term !="Overweight (25.0-29.9)"&
                   term!="Obese (30.0>)"&term !="Underweight (<18.5)"]
pad_lg[term=="Normal weight", Factor := "BMI"]

################################################################################
###################################PLOT#########################################
################################################################################

pad_plot <- pad_lg
pad_plot[, factor2 := Factor]
pad_plot[factor2 == "", factor2 := NA]
pad_plot[, factor2 := na.locf(factor2)]
pad_plot[, factor_order := rleid(factor2)]

###order correctly
pad_plot[adjusted_OR=="1 (reference)", ref := 1] 
pad_plot[adjusted_OR!="1 (reference)", ref := 0]

setorder(pad_plot, factor_order, -ref)
pad_plot[term==1, term := "Yes"]
pad_plot[term==0, term := "No"]

pad_plot <- pad_plot %>%
  dplyr::select(factor_order, Factor, term, events, adj_OR, adj_lci, adj_uci, adjusted_OR, crd_OR, crd_lci, crd_uci, crude_OR )
pad_plot[,1:3] <- map(pad_plot[,1:3], as.factor)    
pad_plot[, order := .N:1]
str(pad_plot)

pad_plot[, adj_OR := as.numeric(adj_OR)]
pad_plot[, adj_lci := as.numeric(adj_lci)]
pad_plot[, adj_uci:= as.numeric(adj_uci)]

pad_plot[, crd_OR := as.numeric(crd_OR)]
pad_plot[, crd_lci := as.numeric(crd_lci)]
pad_plot[, crd_uci:= as.numeric(crd_uci)]

pad_long_plot <- pad_plot %>%
  pivot_longer(cols=c(crude_OR, crd_OR, crd_lci, crd_uci,
                      adjusted_OR, adj_OR, adj_lci, adj_uci),
               names_to=c(".value", "Model"),
               names_pattern="(\\w+)(\\w+)$")
crd_long_plot <- pad_plot[, c("Factor", "term", "crd_OR", "crude_OR", "crd_OR", "crd_lci", "crd_uci", "factor_order", "events", "order")]
crd_long_plot[, model := 1]

adj_long_plot <- pad_plot[, c("Factor", "term", "adjusted_OR", "adj_OR", "adj_lci", "adj_uci", "factor_order", "events", "order")]
adj_long_plot[, model := 2]

final_pad_plot <- rbind(crd_long_plot, adj_long_plot, fill=T)

final_pad_plot[model==1, OR := crd_OR]
final_pad_plot[model==2, OR := adj_OR]
final_pad_plot[model==1, lci := crd_lci]
final_pad_plot[model==2, lci := adj_lci]
final_pad_plot[model==1, uci := crd_uci]
final_pad_plot[model==2, uci := adj_uci]
final_pad_plot[model==1, OR_full := crude_OR]
final_pad_plot[model==2, OR_full := adjusted_OR]

final_pad_plot <- final_pad_plot[, c("Factor", "term", "OR", "lci", "uci", "OR_full", "events", "factor_order", "order", "model")]

setorder(final_pad_plot, -order, model)

final_pad_plot <- final_pad_plot[(OR==1.00 & model==1) | OR!=1.00]
final_pad_plot[model==2, term :=" "] 
final_pad_plot[model==2, events :=" "] 
final_pad_plot[, order_2 := .N:1]

final_pad_plot[model==1, Model := "Crude"]
final_pad_plot[model==2, Model := "Adjusted"]

main_pad_plot <- ggplot(final_pad_plot, aes(y=order_2, x=OR, xmin=lci, xmax=uci, color=Model)) +
  geom_point() +
  geom_errorbarh(height=.3) +
  geom_vline(xintercept=1, lty=2) +
  ggtitle(" ") +
  theme_classic2( ) +
  theme(strip.text.y = element_blank(),
        #legend.position="bottom",
        #legend.justification = "right",
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  #guides(color=guide_legend(reverse=TRUE)) +
  scale_color_manual(breaks = c("Crude", "Adjusted"),
                     values=c("skyblue", "navy"),
                     labels=c("0"="crude", "1"="adjusted")) +
  xlab("Odds ratio") 

#rename for plot
final_pad_plot[Factor == "age_cat", Factor := "Age group"]
final_pad_plot[Factor == "ethnicity", Factor := "Ethnicity"]
final_pad_plot[Factor == "gender", Factor := "Gender"]
final_pad_plot[Factor == "deprivation", Factor := "Deprivation"]
final_pad_plot[Factor == "smokingstatus", Factor := "Smoking status"]
final_pad_plot[Factor == "treatedhyp", Factor := "Treated hypertension"]
final_pad_plot[Factor == "af_cat", Factor := "Atrial Fibrillation"]
final_pad_plot[Factor == "ckd_cat", Factor := "CKD"]
final_pad_plot[Factor == "ra_cat", Factor := "RA"]
final_pad_plot[Factor == "t2dm_cat", Factor := "Type II diabetes"]
final_pad_plot[Factor == "BMI" & term==" ", Factor := " "]
final_pad_plot[term == "Male", term := "Men"]
final_pad_plot[term == "Female", term := "Women"]
final_pad_plot[term == "1 - least deprived", term := "1 - least"]
final_pad_plot[term == "5 - most deprived", term := "5 - most"]

left_pad_plot <- ggplot(data = final_pad_plot, aes(y = order_2)) + xlim(0,0.8) +
  geom_text(aes(x = 0, label = Factor ),lineheight = 0.001, hjust = 0, size = 4, colour = "black") + 
  geom_text(aes(x = 0.35, label = term),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") +
  geom_text(aes(x = 0.70, label = events),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") + 
  theme_classic2() + 
  ggtitle("     Factor                                                      Events")+   
  theme(plot.title = element_text(size = 12, lineheight=.001, face="bold", vjust=-13),
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

###Plot OR and 95% CI
right_pad_plot <- ggplot(data = final_pad_plot,  aes(y = order_2)) + xlim(0,0.5) +
  geom_text(aes(x = 0, label = OR_full),lineheight = 0.01, hjust = 0, size = 4, colour = "black") + 
  theme_classic2( ) +
  ggtitle("    OR (95% CI)")+
  theme(plot.title = element_text(size = 12 ,lineheight=.001, face="bold", vjust=-13),
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


init_pad_plot <- plot_grid(left_pad_plot, main_pad_plot, right_pad_plot, labels=c('PAD'), label_size =12,   ncol=3)


####################################################ALL COMBINED#####################################
#sec_comb_plot <- plot_grid(init_prepan_plot, init_pan_plot, init_mi_plot, init_stroke_plot, init_angina_plot, init_tia_plot, init_pad_plot, labels=c('Pre-pandemic', 'Pandemic'), align="center", ncol=2)
sec_comb_plot <- plot_grid(init_mi_plot, init_stroke_plot, init_angina_plot, init_tia_plot, init_pad_plot, align="center", ncol=2)
sec_comb_plot <- plot_grid(legend, sec_comb_plot, ncol = 1, rel_heights = c(1, 25), labels=c('Secondary Prevention'), vjust = 8, label_size =16 )

ggsave("sec_init_subtype.pdf", plot = last_plot(), width=25, height=45, units = "in", dpi = 300, path = final_graphs)


