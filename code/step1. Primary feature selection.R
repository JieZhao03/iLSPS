
#================
# Load libraries
#================
library(stats)


#=============================================================
# 1、Primary feature selection--Spearman correlation analysis
#=============================================================

# Load the panel of 107 biomarkers and 1-year OS rates per cancer type 
feature107 <- read.csv("step1.input_feature107.csv", row.names = 1)

# Conduct Spearman correlation analysis between biomarkers and 1-year OS rates following ICI
y <- as.numeric(feature107$ICI.OS_1y_rate)
var <- colnames(feature107)[-1:-2]
spearman <- data.frame(var)
for (i in 1:length(var)){
  test <- cor.test(as.numeric(feature107[,i+2]),y,method = "spearman")
  spearman[i,2] <- test$estimate
  spearman[i,3] <- test$p.value
}
spearman[,4] <- p.adjust(spearman[,3], method = "BH")
names(spearman)[c(2:4)] <- c("correlation with ICI OS","pvalue (ICI OS)","p.adjust")

# Conduct Spearman correlation analysis between biomarkers and TCGA prognostic 1-year OS rates
y <- as.numeric(feature107$TCGA.OS_1y_rate)
for (i in 1:length(var)){
  test <- cor.test(as.numeric(feature107[,i+2]),y,method = "spearman")
  spearman[i,5] <- test$estimate
  spearman[i,6] <- test$p.value
}
names(spearman)[c(5:6)] <- c("correlation with TCGA OS","pvalue (TCGA OS)")

#-----------------------
#   Table S5
#-----------------------
spearman <- spearman[order(spearman$`pvalue (ICI OS)`, spearman$p.adjust),]
write.csv(spearman, file = "step1.output_feature_1st.csv", row.names = F)
