# iLSPS
The following experiments were performed with R (version 4.0.2)

- Feature selection through **Random forest (RF)** algorithm and **Gaussian process regression (GPR)** algorithm-based wrapper method using _feature_selection.R_\
 Intersections of the top 10 important features in RF and the selected feature subset from GPR were obtained for subsequent model establishment.
- Train and test the **GPR** model with TCGA multi-omics datasets and ICI-treated clinical trials of 25 cancer types to calculate iLSPS, validate the predictive value of iLSPS in six external and one internal ICI-treated cohorts using _model_construction_and_utility.R_
