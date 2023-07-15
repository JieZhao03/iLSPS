
#================
# Load libraries
#================
library(randomForest)
library(mlr)
library(DMwR)
library(ggplot2)
library(dplyr)

# Load top25 features data from primary feature selection via Spearman's correlation
top25 <- read.csv("top25features.csv", row.names = 1)


#====================================================
# 1、Secondary feature selection--random forest (RF)
#====================================================

# Setting a seed to create reproducible results
set.seed(209) 

# Impute missing values using proximity from randomForest
top25_rf <- rfImpute(OS_1y_rate ~ ., data = top25, iter = 8, ntree = 390)

# Utilized an RF of 200 trees to identify the importance of the features
rf <- randomForest(OS_1y_rate ~ ., data = top25_rf, ntree= 200, importance = T) 

# Assess the feature importance using the measure of %IncMSE
rfvip <- data.frame(vip = rf$importance[,"%IncMSE"])

# Transform data for plotting
rfvip <- rfvip %>% mutate(variable = rownames(rfvip)) %>% 
  arrange(desc(vip)) %>% 
  mutate(variable = factor(variable, levels = variable))

# Obtain the top 10 important features in randomForest
rf10_feature <- head(rfvip$variable, 10)

#-----------------------
#   Figure 3C
#-----------------------
ggplot(rfvip, aes(x = variable, y = vip, fill = variable)) + 
  geom_bar(stat = "identity", width=0.5) + 
  ylab("%IncMSE (feature importance)") + 
  theme_bw() + theme(panel.grid=element_blank(), legend.position = "",
                     axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))


#===========================================================================
# 2、Secondary feature selection--GPR feature selection with wrapper method
#===========================================================================

# Setting a seed to create reproducible results
set.seed(4032)

# Impute missing values using the KNN algorithm with k=10 and normalize features
top25_gpr <- knnImputation(top25, k = 10) %>%
  normalizeFeatures(target = "OS_1y_rate")

# Select a feature subset leading to the best learner performance using the wrapper method
select_task = makeRegrTask(id = "gausspr", top25_gpr, target = "OS_1y_rate")
lrn = makeFeatSelWrapper("regr.gausspr",
                         resampling = makeResampleDesc("CV", iters = 10),
                         control = makeFeatSelControlSequential(method = "sbs"),
                         show.info = FALSE)
select = train(lrn, task = select_task)

# Obtain the feature subset from GPR wrapper method
gpr_feature <- getFeatSelResult(select)$x

#-----------------------
#   Figure 3D
#-----------------------
final_feature <- intersect(rf10_feature, gpr_feature)
final_feature 
# plot Figure 3D using adobe illustrator







