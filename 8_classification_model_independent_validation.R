#!/usr/bin/env Rscript
# Copyright 2021-2026 Na Jiao
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-W", "--workplace"), type="character", default="~/Project/CRC_kingdom/",
                help="Input workplace [default %default]"),
    make_option(c("-m", "--metadata"), type="character", default="CRC_metadata_crcdis_CHI_JAP_com.csv",
                help="metadata file [default %default]"),
    make_option(c("-f", "--feature_profile"), type="character", default="Results/CRC_Discovery_Bacteria_Species_ra_unadj_abunt_rename.csv",
                help="profile of features [default %default]"),
    make_option(c("-s", "--select_feature"), type="character", default="Results/CRC_Bacteria_feature.csv",
                help="MMUPHin file [default %default]"),
    make_option(c("-v", "--validation_profile"), type="character", default="Results/CRC_Validation_Bacteria_Species_table_ra_rename.csv",
                help="validation profile file [default %default]"),
    make_option(c("-p", "--model_para"), type="character", default="Results/CRC_Discovery_Bacteria_all_selfcv_bestmodel.csv",
                help="model parameter file [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix, include pdf and txt [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  prefix = gsub(".txt$", "", opts$input, perl = T)
  if (opts$output==""){
    opts$output=paste(prefix, "_", opts$type, sep = "")}

  # 显示输入输出确认是否正确
  print(paste("The workplace is ", opts$workplace,  sep = ""))
  print(paste("The metadata file is ", opts$metadata,  sep = ""))
  print(paste("The profile of features is ", opts$feature_profile,  sep = ""))
  print(paste("The selected feature file ", opts$select_feature,  sep = ""))
  print(paste("The model parameter file ", opts$model_para,  sep = ""))
  print(paste("The profile of validation is ", opts$validation_profile,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}

package_list = c("dplyr","ggplot2","randomForest","caret","A3","Boruta","ROCR","ImageGP","corrplot","vegan","Hmisc")

###CHECK PACKAGES
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

###load file
model_rf <- read.csv(file = paste(opts$workplace,opts$model_para,sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
meta.all <- read.csv(file = paste(opts$workplace,opts$metadata,sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
feat.all <- read.csv(file = paste(opts$workplace,opts$feature_profile,sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
feature <- read.csv(file = paste(opts$workplace,opts$select_feature,sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
feature <- rownames(feature)
feat.all[is.na(feat.all)] <- 0
rownames(meta.all) <- meta.all$Run
meta.all$StudyID <- factor(meta.all$country)
feat.sig <- t(feat.all[feature,meta.all$Run])
data_use <- cbind(feat.sig,data.frame(meta.all[rownames(feat.sig),'Group']))
colnames(data_use)[ncol(data_use)] <- 'Group'
Group <- factor(data_use$Group)
data_val <- read.csv(file = paste(opts$workplace,opts$validation_profile,sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
meta_val <- read.csv("~/Project/CRC_kingdom/CRC_count_metadata_CHI_JAP_validation.csv",stringsAsFactors = FALSE, header = TRUE, row.names = 1,check.name = FALSE)

meta_val <- filter(meta_val,Group=='CTR'|Group=='CRC')
meta_val$Run <- rownames(meta_val)
data_val[is.na(data_val)] <- 0
feat_val <- data.frame(t(data_val[feature,]))
data_val_use <- cbind(feat_val,meta_val[rownames(feat_val),'Group'])
colnames(data_val_use)[ncol(data_val_use)] <- "Group"

# define customRF
customRF <- list(type = "Classification", 
                 library = "randomForest", 
                 loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), 
                                  class = rep("numeric", 2), 
                                  label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

### cv fold
set.seed(123)
data.SH <- data_val_use[filter(meta_val,Study=="CHISH")$Run,]
Group <- data.SH$Group
print(Group)
SH_rf <- tune_model(data.SH,Group)
SH_AUC <- get_self_cv(data.SH,Group,best_rf=SH_rf)
mean(SH_AUC$metrics[,'AUC'])
write.csv(SH_AUC$metrics,paste(opts$workplace,"CRC_Validation_",opts$output,"_SH_selfcv_metric.csv",sep=''))

data.USA <- data_val_use[filter(meta_val,Study=="USA")$Run,]
Group <- data.USA$Group
USA_rf <- tune_model(data.USA,Group)
USA_AUC <- get_self_cv(data.USA,Group,best_rf=USA_rf)
mean(USA_AUC$metrics[,'AUC'])
write.csv(USA_AUC$metrics,paste(opts$workplace,"CRC_Validation_",opts$output,"_USA_selfcv_metric.csv",sep=''))


###Study-to-study validation

Valid_STS <- get_study_to_study_validation(data_use,data_val_use,model_rf)
write.csv(Valid_STS,paste(opts$workplace,"Result_new/CRC_Validation_",opts$output,"_studytostudy_metric.csv",sep=''))
print('Validation done')

