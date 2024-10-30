##########################################################################################################
###################TABLE AND FIGURES: FACTORS ASSOCIATED WITH INITIATION###################
#####################################PRIMARY AND SECONDARY PREVENTION ##################################
#############################################################################################
rm(list=ls())

library(epiDisplay)
library(epitools)
library(finalfit)
library(openxlsx)
library(broom)
library(stringr)
library(zoo)
library(cowplot)

lapply(c("ggpubr", "grid", "gridExtra", "forcats"), require, character.only=T)

##open parquet
analysis_prim <- read_parquet(paste0(datafiles, "prim_init_rec_risk.parquet"))

test <- analysis_prim[is.na(mean_sbp)]

#explore
library(crosstable)
crosstable(analysis_prim, c(age_cat), by=initiation) %>%
  as_flextable(keep_id=TRUE)
#40-49 initiation 1104 (16.93%), no 5417 (83.07%)
#50-59 initiation 3236 (13.15%), no 21377 (86.85%)

table(analysis_prim$initiation)
prop.table(table(analysis_prim$initiation))

##format bmi
analysis_prim[bmi_cat=="Normal weight (18.5-24.9)", bmi_cat := "Normal weight"]
analysis_prim[bmi_cat=="Overweight (25.0-29.9)", bmi_cat := "Overweight"]
analysis_prim[bmi_cat=="Obese (30.0>)", bmi_cat := "Obese"]
analysis_prim[bmi_cat=="Underweight (<18.5)", bmi_cat := "Underweight"]

analysis_prim[, bmi_cat := relevel(factor(bmi_cat), ref="Normal weight")]

analysis_prim[, .N, by="bmi_cat"]

####################Primary prevention###################################
uniqueN(analysis_prim , by = "patid") #77,039
analysis_prim[, .N, by=initiation]

#cholesterol variables have high missingness so don't include. sbp and dbp have less missingness
#than the cholesterol variables
prim_lg <- glm(initiation ~ age_cat + gender + ethnicity + deprivation + bmi_cat + smokingstatus + 
                treatedhyp + af_cat + ra_cat + t2dm_cat, data = analysis_prim, family = binomial())
logistic.display(prim_lg)

explanatory= c("age_cat", "gender", "ethnicity",  "deprivation", "bmi_cat", "smokingstatus",
               "treatedhyp","af_cat", "ra_cat", "t2dm_cat")
dependent="initiation"
analysis_prim %>%
  finalfit(dependent, explanatory, confint_sep=" - ", confint_level=0.95) %>%
  ff_remove_p() -> prim_lg

########modify for table and plot
setnames(prim_lg, "Dependent: initiation", "Factor")
setnames(prim_lg, "yes", "events")
setnames(prim_lg, "OR (univariable)", "crude_OR")
setnames(prim_lg, "OR (multivariable)", "adjusted_OR")
colnames(prim_lg)[2] <- "term"

prim_lg <- as.data.table(prim_lg)

###extract values to different column
#crude
prim_lg[, crd_OR := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
prim_lg[, crd_lci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
prim_lg[, crd_uci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]


prim_lg[, adj_OR := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
prim_lg[, adj_lci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
prim_lg[, adj_uci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]
prim_lg[, events := str_extract(events, "(\\d+) \\((.*)\\)", group=1)]

#rename reference groups
prim_lg[crude_OR=="-", c("crd_OR", "crd_lci", "crd_uci") := list(1, 1, 1)]
prim_lg[crude_OR=="-", crude_OR := "1 (reference)"]

prim_lg[adjusted_OR=="-", c("adj_OR", "adj_lci", "adj_uci") := list(1, 1, 1)]
prim_lg[adjusted_OR=="-", adjusted_OR := "1 (reference)"]

# Remove 'missing' rows
prim_lg <- prim_lg[term != "missing"]
prim_lg <- prim_lg[term !="Normal weight (18.5-24.9)"& term !="Overweight (25.0-29.9)"&
                         term!="Obese (30.0>)"&term !="Underweight (<18.5)"]
prim_lg[term=="Normal weight", Factor := "BMI"]

################################################################################
###################################PLOT#########################################
################################################################################

prim_plot <- prim_lg
prim_plot[, factor2 := Factor]
prim_plot[factor2 == "", factor2 := NA]
prim_plot[, factor2 := na.locf(factor2)]
prim_plot[, factor_order := rleid(factor2)]

###order correctly
prim_plot[adjusted_OR=="1 (reference)", ref := 1] 
prim_plot[adjusted_OR!="1 (reference)", ref := 0]

setorder(prim_plot, factor_order, -ref)
prim_plot[term==1, term := "Yes"]
prim_plot[term==0, term := "No"]

prim_plot <- prim_plot %>%
  dplyr::select(factor_order, Factor, term, events, adj_OR, adj_lci, adj_uci, adjusted_OR, crd_OR, crd_lci, crd_uci, crude_OR )
prim_plot[,1:3] <- map(prim_plot[,1:3], as.factor)    
prim_plot[, order := .N:1]
str(prim_plot)

prim_plot[, adj_OR := as.numeric(adj_OR)]
prim_plot[, adj_lci := as.numeric(adj_lci)]
prim_plot[, adj_uci:= as.numeric(adj_uci)]

prim_plot[, crd_OR := as.numeric(crd_OR)]
prim_plot[, crd_lci := as.numeric(crd_lci)]
prim_plot[, crd_uci:= as.numeric(crd_uci)]

prepan_long_plot <- prim_plot %>%
  pivot_longer(cols=c(crude_OR, crd_OR, crd_lci, crd_uci,
                      adjusted_OR, adj_OR, adj_lci, adj_uci),
               names_to=c(".value", "Model"),
               names_pattern="(\\w+)(\\w+)$")
crd_long_plot <- prim_plot[, c("Factor", "term", "crd_OR", "crude_OR", "crd_OR", "crd_lci", "crd_uci", "factor_order", "events", "order")]
crd_long_plot[, model := 1]

adj_long_plot <- prim_plot[, c("Factor", "term", "adjusted_OR", "adj_OR", "adj_lci", "adj_uci", "factor_order", "events", "order")]
adj_long_plot[, model := 2]

final_prim_plot <- rbind(crd_long_plot, adj_long_plot, fill=T)

final_prim_plot[model==1, OR := crd_OR]
final_prim_plot[model==2, OR := adj_OR]
final_prim_plot[model==1, lci := crd_lci]
final_prim_plot[model==2, lci := adj_lci]
final_prim_plot[model==1, uci := crd_uci]
final_prim_plot[model==2, uci := adj_uci]
final_prim_plot[model==1, OR_full := crude_OR]
final_prim_plot[model==2, OR_full := adjusted_OR]

final_prim_plot <- final_prim_plot[, c("Factor", "term", "OR", "lci", "uci", "OR_full", "events", "factor_order", "order", "model")]

setorder(final_prim_plot, -order, model)

final_prim_plot <- final_prim_plot[(OR==1.00 & model==1) | OR!=1.00]
final_prim_plot[model==2, term :=" "] 
final_prim_plot[model==2, events :=" "] 
final_prim_plot[, order_2 := .N:1]

final_prim_plot[model==1, Model := "Crude"]
final_prim_plot[model==2, Model := "Adjusted"]

main_prim_plot <- ggplot(final_prim_plot, aes(y=order_2, x=OR, xmin=lci, xmax=uci, color=Model)) +
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
final_prim_plot[Factor == "age_cat", Factor := "Age group"]
final_prim_plot[Factor == "ethnicity", Factor := "Ethnicity"]
final_prim_plot[Factor == "gender", Factor := "Gender"]
final_prim_plot[Factor == "deprivation", Factor := "Deprivation"]
final_prim_plot[Factor == "smokingstatus", Factor := "Smoking status"]
final_prim_plot[Factor == "treatedhyp", Factor := "Treated hypertension"]
final_prim_plot[Factor == "af_cat", Factor := "Atrial Fibrillation"]
final_prim_plot[Factor == "ra_cat", Factor := "RA"]
final_prim_plot[Factor == "t2dm_cat", Factor := "Type II diabetes"]
final_prim_plot[term == "1 - least deprived", term := "1 - least"]
final_prim_plot[term == "5 - most deprived", term := "5 - most"]
final_prim_plot[term == "Male", term := "Men"]
final_prim_plot[term == "Female", term := "Women"]

left_prim_plot <- ggplot(data = final_prim_plot, aes(y = order_2)) + xlim(0,0.8) +
  geom_text(aes(x = 0, label = Factor ),lineheight = 0.001, hjust = 0, size = 4, colour = "black") + 
  geom_text(aes(x = 0.37, label = term),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") +
  geom_text(aes(x = 0.70, label = events),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") + 
  theme_classic2() + 
  ggtitle("     Factor                                                         Events")+         
  theme(plot.title = element_text(size = 11, lineheight=.001,  face="bold", vjust=-10),
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
right_prim_plot <- ggplot(data = final_prim_plot,  aes(y = order_2)) + xlim(0,0.5) +
  geom_text(aes(x = 0, label = OR_full),lineheight = 0.01, hjust = 0, size = 4, colour = "black") + 
  theme_classic2( ) +
  ggtitle("    OR (95% CI)")+
  theme(plot.title = element_text(size = 11,lineheight=.001,  face="bold", vjust=-10),
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

init_prim_plot <- plot_grid(left_prim_plot, main_prim_plot, right_prim_plot, labels=c('A. Primary Prevention'),  label_size =12, label_y=0.98,  ncol=3) 
ggsave("initiation_prim_plot.png", plot=init_prim_plot, width=10, height=12, path = graphs)

##create legend
legend <- get_legend(
  main_prim_plot + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top",
          legend.text = element_text(size=12)))

prim_comb_plot <- plot_grid(legend, init_prim_plot, ncol = 1, rel_heights = c(1, 45))
##################################################################################################################
#####################################SECONDARY PREVENTION##################################
###################################################################################################################################

##open parquet
sec_init_dt <- read_parquet(paste0(datafiles, "sec_initiators.parquet"))

##format bmi
sec_init_dt[bmi_cat=="Normal weight (18.5-24.9)", bmi_cat := "Normal weight"]
sec_init_dt[bmi_cat=="Overweight (25.0-29.9)", bmi_cat := "Overweight"]
sec_init_dt[bmi_cat=="Obese (30.0>)", bmi_cat := "Obese"]
sec_init_dt[bmi_cat=="Underweight (<18.5)", bmi_cat := "Underweight"]

sec_init_dt[, bmi_cat := relevel(factor(bmi_cat), ref="Normal weight")]

####################Pre-pandemic###################################
uniqueN(sec_init_dt , by = "patid") #77,039
sec_init_dt[, .N, by=initiation]
#initiation     N
#1:        yes 10196
#2:         no 72702


####PRE-PANDEMIC###
#cholesterol variables have high missingness so don't include. sbp and dbp have less missingness
#than the cholesterol variables
sec_lg <- glm(initiation ~ age_cat + gender + ethnicity + deprivation + bmi_cat + smokingstatus + 
                treatedhyp + af_cat + ckd_cat + ra_cat + t2dm_cat, data = sec_init_dt, family = binomial())
logistic.display(sec_lg)

explanatory= c("age_cat", "gender", "ethnicity",  "deprivation", "bmi_cat", "smokingstatus",
               "treatedhyp","af_cat", "ckd_cat", "ra_cat", "t2dm_cat")
dependent="initiation"
sec_init_dt %>%
  finalfit(dependent, explanatory, confint_sep=" - ", confint_level=0.95) %>%
  ff_remove_p() -> sec_lg

########modify for table and plot
setnames(sec_lg, "Dependent: initiation", "Factor")
setnames(sec_lg, "yes", "events")
setnames(sec_lg, "OR (univariable)", "crude_OR")
setnames(sec_lg, "OR (multivariable)", "adjusted_OR")
colnames(sec_lg)[2] <- "term"

sec_lg <- as.data.table(sec_lg)

###extract values to different column
#crude
sec_lg[, crd_OR := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
sec_lg[, crd_lci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
sec_lg[, crd_uci := str_extract(crude_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]


sec_lg[, adj_OR := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=1)]
sec_lg[, adj_lci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=2)]
sec_lg[, adj_uci := str_extract(adjusted_OR, "(\\d+\\.?\\d*) \\((\\d+\\.?\\d*) - (\\d+\\.?\\d*)\\)", group=3)]
sec_lg[, events := str_extract(events, "(\\d+) \\((.*)\\)", group=1)]

#rename reference groups
sec_lg[crude_OR=="-", c("crd_OR", "crd_lci", "crd_uci") := list(1, 1, 1)]
sec_lg[crude_OR=="-", crude_OR := "1 (reference)"]

sec_lg[adjusted_OR=="-", c("adj_OR", "adj_lci", "adj_uci") := list(1, 1, 1)]
sec_lg[adjusted_OR=="-", adjusted_OR := "1 (reference)"]

# Remove 'missing' rows
sec_lg <- sec_lg[term != "missing"]
sec_lg <- sec_lg[term !="Normal weight (18.5-24.9)"& term !="Overweight (25.0-29.9)"&
                         term!="Obese (30.0>)"&term !="Underweight (<18.5)"]
sec_lg[term=="Normal weight", Factor := "BMI"]

################################################################################
###################################PLOT#########################################
################################################################################

sec_plot <- sec_lg
sec_plot[, factor2 := Factor]
sec_plot[factor2 == "", factor2 := NA]
sec_plot[, factor2 := na.locf(factor2)]
sec_plot[, factor_order := rleid(factor2)]

###order correctly
sec_plot[adjusted_OR=="1 (reference)", ref := 1] 
sec_plot[adjusted_OR!="1 (reference)", ref := 0]

setorder(sec_plot, factor_order, -ref)
sec_plot[term==1, term := "Yes"]
sec_plot[term==0, term := "No"]

sec_plot <- sec_plot %>%
  dplyr::select(factor_order, Factor, term, events, adj_OR, adj_lci, adj_uci, adjusted_OR, crd_OR, crd_lci, crd_uci, crude_OR )
sec_plot[,1:3] <- map(sec_plot[,1:3], as.factor)    
sec_plot[, order := .N:1]
str(sec_plot)

sec_plot[, adj_OR := as.numeric(adj_OR)]
sec_plot[, adj_lci := as.numeric(adj_lci)]
sec_plot[, adj_uci:= as.numeric(adj_uci)]

sec_plot[, crd_OR := as.numeric(crd_OR)]
sec_plot[, crd_lci := as.numeric(crd_lci)]
sec_plot[, crd_uci:= as.numeric(crd_uci)]

prepan_long_plot <- sec_plot %>%
  pivot_longer(cols=c(crude_OR, crd_OR, crd_lci, crd_uci,
                      adjusted_OR, adj_OR, adj_lci, adj_uci),
               names_to=c(".value", "Model"),
               names_pattern="(\\w+)(\\w+)$")
crd_long_plot <- sec_plot[, c("Factor", "term", "crd_OR", "crude_OR", "crd_OR", "crd_lci", "crd_uci", "factor_order", "events", "order")]
crd_long_plot[, model := 1]

adj_long_plot <- sec_plot[, c("Factor", "term", "adjusted_OR", "adj_OR", "adj_lci", "adj_uci", "factor_order", "events", "order")]
adj_long_plot[, model := 2]

final_sec_plot <- rbind(crd_long_plot, adj_long_plot, fill=T)

final_sec_plot[model==1, OR := crd_OR]
final_sec_plot[model==2, OR := adj_OR]
final_sec_plot[model==1, lci := crd_lci]
final_sec_plot[model==2, lci := adj_lci]
final_sec_plot[model==1, uci := crd_uci]
final_sec_plot[model==2, uci := adj_uci]
final_sec_plot[model==1, OR_full := crude_OR]
final_sec_plot[model==2, OR_full := adjusted_OR]

final_sec_plot <- final_sec_plot[, c("Factor", "term", "OR", "lci", "uci", "OR_full", "events", "factor_order", "order", "model")]

setorder(final_sec_plot, -order, model)

final_sec_plot <- final_sec_plot[(OR==1.00 & model==1) | OR!=1.00]
final_sec_plot[model==2, term :=" "] 
final_sec_plot[model==2, events :=" "] 
final_sec_plot[, order_2 := .N:1]

final_sec_plot[model==1, Model := "Crude"]
final_sec_plot[model==2, Model := "Adjusted"]

main_sec_plot <- ggplot(final_sec_plot, aes(y=order_2, x=OR, xmin=lci, xmax=uci, color=Model)) +
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
final_sec_plot[Factor == "age_cat", Factor := "Age group"]
final_sec_plot[Factor == "ethnicity", Factor := "Ethnicity"]
final_sec_plot[Factor == "gender", Factor := "Gender"]
final_sec_plot[Factor == "deprivation", Factor := "Deprivation"]
final_sec_plot[Factor == "smokingstatus", Factor := "Smoking status"]
final_sec_plot[Factor == "treatedhyp", Factor := "Treated hypertension"]
final_sec_plot[Factor == "af_cat", Factor := "Atrial Fibrillation"]
final_sec_plot[Factor == "ckd_cat", Factor := "CKD"]
final_sec_plot[Factor == "ra_cat", Factor := "RA"]
final_sec_plot[Factor == "t2dm_cat", Factor := "Type II diabetes"]
final_sec_plot[Factor == "BMI" & term==" ", Factor := " "]
final_sec_plot[term == "Male", term := "Men"]
final_sec_plot[term == "Female", term := "Women"]
final_sec_plot[term == "1 - least deprived", term := "1 - least"]
final_sec_plot[term == "5 - most deprived", term := "5 - most"]

left_sec_plot <- ggplot(data = final_sec_plot, aes(y = order_2)) + xlim(0,0.8) +
  geom_text(aes(x = 0, label = Factor ),lineheight = 0.001, hjust = 0, size = 4, colour = "black") + 
  geom_text(aes(x = 0.37, label = term),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") +
  geom_text(aes(x = 0.70, label = events),lineheight = 0.001, hjust = 0 ,size = 4, colour = "black") + 
  theme_classic2() + 
  ggtitle("     Factor                                                         Events")+     
  theme(plot.title = element_text(size = 11, lineheight=.001, face="bold", vjust=-10),
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
right_sec_plot <- ggplot(data = final_sec_plot,  aes(y = order_2)) + xlim(0,0.5) +
  geom_text(aes(x = 0, label = OR_full),lineheight = 0.01, hjust = 0, size = 4, colour = "black") + 
  theme_classic2( ) +
  ggtitle("   OR (95% CI)")+
  theme(plot.title = element_text(size = 11 ,lineheight=.001, face="bold", vjust=-10),
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

sec_comb_plot<- plot_grid(left_sec_plot, main_sec_plot, right_sec_plot, labels=c('B. Secondary Prevention'),  label_size =12, label_y=0.98, hjust=-0.8, ncol=3)
ggsave("sec_init_plot.pdf", plot=sec_comb_plot, width=10, height=12, path = final_graphs)

#########################################################################################
####COMBINE WITH PRIMARY PREVENTION PLOT
total_comb <- plot_grid(prim_comb_plot, sec_comb_plot, ncol=1)
print(total_comb)
ggsave("prim_sec_init_plot.pdf", plot = last_plot(), width=12, height=20, units = "in", dpi = 400, path = final_graphs)

