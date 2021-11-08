library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(vegan)
setwd('/Users/profile/')

### meta_analysis
#load data
meta.all <- read.csv(file = '../metadata/CRC_metadata_crcdis_CHI_JAP_com.csv',stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                     check.name = FALSE)
rownames(meta.all) <- meta.all$Run
meta.all$StudyID <- factor(meta.all$country)

feat.abu <- read.csv("Profile/CRC_Discovery_Bacteria_Species_table_ra.csv",stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                     check.name = FALSE)
feat.abu <- feat.abu[,rownames(meta.all)]

feat.abu[is.na(feat.abu)] <- 0
feat.abu <- feat.abu/100

fit_meta <- lm_meta(feature_abd = feat.abu,
                    exposure = "Group",
                    batch = "StudyID",
                    covariates = c("Gender", "Age","BMI"),
                    control = list(rma_method="HS",transform="AST"),
                    data = meta.all)
meta_results <- fit_meta$meta_fits
maaslin_results <- fit_meta$maaslin_fits


write.csv(meta_results, 'Results/CRC_Bacteria_Species_ra_unadj_MMUPHin_result.csv')
write.csv(maaslin_results, 'Results/CRC_Bacteria_Species_ra_unadj_MMUPHin_maaslin_result.csv')

meta_results %>% 
  filter(pval < 0.05) %>% 
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature)) %>% 
  ggplot(aes(y = coef, x = feature)) +
  geom_bar(stat = "identity") +
  coord_flip()



