Multi-kingdom microbiota analyses identify bacterial-fungal interactions and biomarkers of colorectal cancer across cohorts

###Data Preprocess
##major step could be found in 1_raw sequence process.sh, mainly include:
##quality control with kneaddata, taxonomic profile with Kraken2 and Braken, functional profile by assembly (megahit) and gene prediction with Prodigal. Then gene abundance is estimated with Coverm

###Profile_preprocess, code in 2_Profile_preprocess.R
##Profiles were converted into relative abundance and only that with prevalence above 20% samples were filtered for further analysis 

###Differential species/functions, code in 3_MMUPHin_differential_analysis.R
##Since the microbial profiles are compositional and sparse and heterogeneity exist among different cohorts, Meta-analysis Methods with a Uniform Pipeline for Heterogeneity (MMUPHin)76 was performed to identify CRC-related differential microbial species, which enables the normalization and combination of multiple microbial community studies. In the MMUPHin analysis, microbial community batch effect among cohorts was corrected via a Combat-like extended method. Microbial profile was arcsine square root (AST) transformed and age, gender and BMI of subjects were treated as covariates.

###The random forest-based model construction
##Feature selection, code in 4_feature_selection.R
##feature selection via ‘Boruta’ package in R with default parameters (pValue=0.05, mcAdj=T, maxRuns=1000) which iteratively removes the features that are proved by a statistical test to be less relevant than random probes. Correlations between “Confirmed features” identified by Boruta were then calculated and only features with correlation less than 0.7 were selected to further model construction to avoid co-linearity issue

##model construction, code in 5_classification_model_cv_repeat.R
##to construct predictive models, we tuned hyperparameters (e.g. mtry, ntree, nodesize, maxnodes) using ‘caret’ packages. Finally, with the best combination of hyperparameters, five-fold cross-validation model was constructed to avoid overfitting issues, which was constructed with each cohort and repeated 20 times.
##Model’s significance was accessed with 1000 permutations via “A3” package. code in 6_classification_model_significance.R

##Evaluation generalization of microbial markers, code in 7_classification_model_studytostudy_loco.R
##To further test the generalization of CRC microbial markers across technical and geographic differences in multiple populations, we extensively validated the diagnostic models with cohort-to-cohort transfer validation and leave-one-cohort-out (LOCO) validation as described previously

##Independent validation for robustness of multi-kingdom microbial markers with additional datasets,code in 8_classification_model_independent_validation.R

##specificity of features in non-CRC diseases, code in 9_classification_model_nonCRC_model_specificity.R
##To avoid false positives in clinical diagnoses, we estimated the specificity of microbial biomarkers for CRC via testing the AUROC values of models constructed with the best panel of features aganist non-CRC disease.
