# iLSPS
The following experiments were performed with R (version 4.0.2)

- **Step0**,  install the packages used for the experiments using the _"step0. list of R packages.R"_ script in R (version 4.0.2).
- **Step1**,  primary feature selection through **Spearman correlation** should be run using the _"step1. Primary feature selection.R"_ script, with "step1.input_feature107.csv" as input data, which will output the "step1.output_feature_1st.csv" file (Table S5 of the paper). The top 20 variables according to the "p.adjust" and top 5 variables of mutation-based biomarkers were selected for the following secondary feature selection.
- **Step2**,  secondary feature selection through **Random forest (RF)** algorithm and **Gaussian process regression (GPR)** algorithm-based wrapper method should be run using the _"step2. Secondary feature selection.R"_ script, with "step2.input_top25features.csv" file.\
Intersections of the top 10 important features in RF and the selected feature subset from GPR were obtained for subsequent model construction.
- **Step3**,  the train and test of the **GPR model** with TCGA multi-omics datasets and ICI-treated clinical trials of 25 cancer types should be run using the _"step3. Model_construction.R"_ script [Point 1 and 2], with "step3.input_trainset.csv" and "step3.input_testset.csv" files. Then, the iLSPS of six external and one internal ICI-treated cohorts should be calculated using the _"step3. Model_construction.R"_ script [Point 3], with "step3.input_ICI-treated cohorts.csv" file, which will output the "step4.iLSPS_ICIcohorts.csv" file.
- **Step4**,  the validation of the predictive value of iLSPS should be run using the _"step4. Model_utility.R"_ script with the "step4.iLSPS_ICIcohorts.csv" file.
