
#================
# Load libraries
#================
library(survminer)
library(survival)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(dplyr)

#===========================================================================
# 1、Test the predictive value of iLSPS for survival in ICI-treated cohorts
#===========================================================================

# Load the iLSPS values of seven independent ICI-treated cohorts
ICI_cohort <- read.csv("step4.iLSPS_ICIcohorts.csv")

# Obtain the optimal cut-point for each cohort from Figure S5D
cut_off <- data.frame(
  study = c("Snyder 2014", "Van Allen 2015", "Hugo 2016", "Riaz 2017",
            "Mariathasan 2018", "Liu 2019", "Wang 2022"),
  cut.off = c(43.5, 38.5, 49.5, 49, 51.5, 38, 44)) 

# Survival analysis Use a for-loop  
splots <- list() 

for (i in c("Snyder 2014", "Mariathasan 2018", "Van Allen 2015", 
            "Liu 2019", "Hugo 2016", "Wang 2022", "Riaz 2017")) {
  # Transform the survival data set           
  surv <- ICI_cohort %>% filter(study == i) %>% 
    mutate(iLSPS = ifelse(iLSPS > cut_off[cut_off$study == i, "cut.off"], "High","Low")) %>%
    mutate(iLSPS = factor(iLSPS, levels = c("Low", "High")))
  
  # Fit a Cox proportional hazards regression model for the survival data set
  hr <- coxph(Surv(OS, censor == 1) ~ iLSPS, surv) %>% summary() 
  
  # Create Kaplan–Meier survival curves of iLSPS-high and iLSPS-low arms
  fit <- survfit(Surv(OS, censor == 1) ~ iLSPS, surv)   
  
  # compare the Kaplan–Meier survival curves using the log-rank test   
  p <- surv_pvalue(fit, surv) 
  
  # ggsurvplot
  psurv <- ggsurvplot(fit, surv, fun="pct", surv.median.line = "hv", censor.shape = (124),
                      legend.labs = c("iLSPS-Low", "iLSPS-High"), 
                      xlab = "Time in month", ylab ='OS (%)',
                      title = paste0(i," (N = ",nrow(surv),")"),)
  psurv$plot <- psurv$plot + annotate(
    "text", x=Inf, y=Inf, hjust = 1, vjust = 1, 
    label = paste0("HR=", round(hr$conf.int[1], 2), "(", 
                   round(hr$conf.int[3], 2), "-", round(hr$conf.int[4], 2), 
                   ")", "\n p = ", sprintf("%0.4f", p$pval))) +
    theme(plot.title = element_text(size = 10))
  splots[[i]] <- psurv
}

#-----------------------
#   Figure 4E
#-----------------------  
# Arrange multiple ggsurvplots
psurv7 <- arrange_ggsurvplots(splots, print = F, ncol = 4, nrow = 2)
psurv7
