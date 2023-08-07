
#================
# Load libraries
#================
library(mlr)
library(DMwR)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(stringr)
library(dplyr)


#======================================================
# 1、Train the gaussian process regression (GPR) model
#======================================================

# Setting a seed to create reproducible results
set.seed(3011)

# Load train dataset
# Impute missing values using the KNN with k=10 and normalize features
trainset <- read.csv("trainset.csv", row.names = 1) %>% 
  knnImputation() %>% 
  normalizeFeatures(target = "OS_1y_rate")

# Create a gaussian process regression train task    
train.task = makeRegrTask(id = "gausspr", trainset, target = "OS_1y_rate")

# Create a GPR learner 
# Tune the parameter "sigma" for the "mse" using a random search in bootstrap resampling strategy
lrn = makeTuneWrapper(makeLearner("regr.gausspr"), 
                      resampling = makeResampleDesc("Bootstrap"), 
                      measures = mse, 
                      par.set = makeParamSet(makeNumericParam("sigma", lower = 0, upper = 40)), 
                      control = makeTuneControlRandom(maxit = 300L))

# Train the GPR model with the tuned bagged learner                 
model_gpr = train(lrn, train.task)

# Predict the survival for the train set and obtain the performance with "mse" and "rmse"
pred_train <- predict(model_gpr, train.task)
performance(pred_train, list("MSE" = mse, "RMSE" = rmse))

# Obtain the predicted iLSPS outcomes for the train set
train_result <- pred_train$data 
train_result <- train_result %>% mutate(cancer = rownames(train_result)) 

# Measure the performance of GPR model with Pearson R between the iLSPS and 1y OS rates 
cor.test(train_result$truth, train_result$response, method = "pearson")

#-----------------------
#   Figure 4A
#-----------------------    
ggplot(train_result, aes(x = response*100, y = truth*100, color = cancer)) + 
  geom_point() + 
  geom_smooth(method = lm, se = T, color = "black", lwd = 0.5, lty = 2) +
  geom_text_repel(aes(response*100, truth*100, label = cancer)) + 
  labs(x = "iLSPS (predict)", y = "1-year OS rate (%) (truth)") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "")



#=====================================================
# 2、Test the gaussian process regression (GPR) model
#=====================================================

# Load independent test dataset 
# Impute missing values using the KNN algorithm with k=10 and normalize features
testset <- read.csv("testset.csv", row.names = 1) %>% 
  knnImputation() %>% 
  normalizeFeatures(target = "OS_1y_rate")

# Create a gaussian process regression test task    
test.task = makeRegrTask(id = "gausspr", testset, target = "OS_1y_rate")

# Predict the survival for the test set and obtain the performance with "mse" and "rmse"
pred_test <- predict(model_gpr, test.task)
performance(pred_test, list("MES" = mse, "RMSE" = rmse))

# Obtain the predicted iLSPS outcomes for the test set
test_result <- pred_test$data 
test_result <- test_result %>% mutate(cancer = rownames(test_result)) %>%
  mutate(cancer = str_split(cancer,"_", simplify = T)[,1])

# Measure the performance of GPR model with Pearson R between the iLSPS and 1y OS rates     
cor.test(test_result$truth, test_result$response, method = "pearson")

#-----------------------
#   Figure S5C
#-----------------------  
ggplot(test_result, aes(x = response*100, y = truth*100, color = cancer)) + 
  geom_point() + 
  geom_smooth(method = lm, se = T, color = "black", lwd = 0.5, lty = 2) +
  labs(x = "iLSPS (predict)", y = "1-year OS rate (%) (truth)") +
  scale_color_ucscgb(guide = guide_legend(byrow = F, ncol = 2))


#===============================================
# 3、Calculate the iLSPS of ICI-treated cohorts
#===============================================

# Load seven independent ICI-treated cohorts
ICI_cohort <- read.csv("ICI-treated cohorts.csv", row.names = 1)

# Obtain the predicted iLSPS outcomes for each ICI-treated cohort using a for-loop
ICI_cohort_iLSPS <- list()
for (i in c("Snyder 2014", "Van Allen 2015", "Hugo 2016", "Riaz 2017",
            "Mariathasan 2018", "Liu 2019", "Wang 2022")) {
  testset <- ICI_cohort[ICI_cohort$study == i, c(-1,-7)] %>% 
    rename("OS_1y_rate" = "OS") %>%
    knnImputation() %>% 
    normalizeFeatures(target = "OS_1y_rate")
  test.task = makeRegrTask(id = "gausspr", testset, target = "OS_1y_rate")
  pred_test <- predict(model_gpr, test.task)
  test_result <- pred_test$data %>% mutate(response = response*100)
  test_result$Patient <- rownames(test_result)
  ICI_cohort_iLSPS[[i]] <- test_result
}  
ICI_cohort_iLSPS <- do.call(rbind, ICI_cohort_iLSPS) %>% 
  rename("iLSPS" = "response")
ICI_cohort$Patient <- rownames(ICI_cohort)
ICI_cohort <- merge(ICI_cohort, ICI_cohort_iLSPS[,c("Patient","iLSPS")])
write.csv(ICI_cohort, file = "step4.iLSPS_ICIcohorts.csv", row.names = F)
