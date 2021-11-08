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
    make_option(c("-W", "--workplace"), type="character", default="~/CRC_kingdom/",
                help="Input workplace [default %default]"),
    make_option(c("-m", "--metadata"), type="character", default="CRC_metadata_crcdis_CHI_JAP_com.csv",
                help="metadata file [default %default]"),
    make_option(c("-f", "--feature_profile"), type="character", default="Results/CRC_Discovery_Bacteria_Species_ra_unadj_abunt_rename.csv",
                help="profile of features [default %default]"),
    make_option(c("-s", "--select_feature"), type="character", default="Results/CRC_Bacteria_feature.csv",
                help="MMUPHin file [default %default]"),
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
  print(paste("The output file prefix is ", opts$output, sep = ""))
}

#package_list = c("dplyr","ggplot2","randomForest","caret","A3","Boruta","ROCR","corrplot","vegan","Hmisc")
package_list = c("ROSE","dplyr","ggplot2","randomForest","caret","A3","Boruta","ROCR","corrplot","vegan","Hmisc")
###CHECK PACKAGES
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

###load file
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

# define customRF
customRF <- list(type = "Classification",
                 library = "randomForest",
                 loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree","nodesize","maxnodes"),
                                  class = rep("numeric", 4),
                                  label = c("mtry", "ntree","nodesize","maxnodes"))
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


###study-to-study
metrics <- matrix(NA,20,12)
colnames(metrics) <- c("Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision",
                       "Recall","F1","Prevalence","Detection Rate","Detection Prevalence",
                       "Balanced Accuracy","AUC")

rownames(metrics) <- c("AUS_CHI","AUS_FRA","AUS_GER","AUS_JAP","CHI_AUS","CHI_FRA","CHI_GER","CHI_JAP",
                       "FRA_AUS","FRA_CHI","FRA_GER","FRA_JAP","GER_AUS","GER_CHI","GER_FRA","GER_JAP",
                       "JAP_AUS","JAP_CHI","JAP_FRA","JAP_GER")
sites <- c('AUS',"CHI","FRA","GER","JAP")
for (i in 1:length(sites)){
  print(sites[i])
  meta_train <- filter(meta.all,country==sites[i])
  print(dim(meta_train))
  fold_train <- data_use[meta_train$Run,]
  #print(dim(fold_train_raw))
  #fold_train <- ovun.sample(Group ~ ., data = fold_train_raw, method = "under", seed=2021)$data
  Group.train <- fold_train$Group
  sites_test <- sites[-i]
  print(sites_test)
  for (j in 1:length(sites_test)){
    print(sites_test[j])
    meta_test <- filter(meta.all,country==sites_test[j])
    print(dim(meta_test))
    fold_test <- data_use[meta_test$Run,]
    #fold_test <- ovun.sample(Group ~ ., data = fold_test_raw, method = "under", seed=2021)$data
    set.seed(0606)
    fold_fit <- randomForest(as.factor(Group.train)~., data=fold_train,mtry=model_rf[sites[i],'mtry'],
			     ntree=model_rf[sites[i],'ntree'],nodesize=model_rf[sites[i],'nodesize'],
			     maxnodes=model_rf[sites[i],'maxnodes'],importance=TRUE)
    fold_pred <- predict(fold_fit,fold_test)
    result <- confusionMatrix(factor(as.vector(fold_pred)),as.factor(fold_test$Group),mode = "everything")
    metrics[paste(sites[i],sites_test[j],sep = '_'),1:11] <- as.numeric(result$byClass)
    predicted_probs <- predict(fold_fit, fold_test, type = 'prob')
    pred <- prediction(predicted_probs[, 'CTR'], fold_test$Group)
    auc <- performance(pred, 'auc')
    metrics[paste(sites[i],sites_test[j],sep = '_'),12] <- auc@y.values[[1]]
  }
}

write.csv(metrics,paste(opts$workplace,"CRC_Discovery_",opts$output,"_all_studytostudy_metric_single_com_final.csv",sep=''))

### leave-one-cohort-out
set.seed(123)

###tune loco model
data_LOCO_AUS <- data_use[filter(meta.all,country!='AUS')$Run,]
Group <- data_LOCO_AUS$Group
AUS_loco_rf <- tune_loco_model(data_LOCO_AUS,Group)

data_LOCO_CHI <- data_use[filter(meta.all,country!='CHI')$Run,]
Group <- data_LOCO_CHI$Group
CHI_LOCO_rf <- tune_LOCO_model(data_LOCO_CHI,Group)

data_LOCO_FRA <- data_use[filter(meta.all,country!='FRA')$Run,]
Group <- data_LOCO_FRA$Group
FRA_LOCO_rf <- tune_LOCO_model(data_LOCO_FRA,Group)

data_LOCO_GER <- data_use[filter(meta.all,country!='GER')$Run,]
Group <- data_LOCO_GER$Group
GER_LOCO_rf <- tune_LOCO_model(data_LOCO_GER,Group)

data_LOCO_JAP <- data_use[filter(meta.all,country!='JAP')$Run,]
Group <- data_LOCO_JAP$Group
JAP_LOCO_rf <- tune_LOCO_model(data_LOCO_JAP,Group)

model_LOCO_rf <- rbind(AUS_LOCO_rf,CHI_LOCO_rf,FRA_LOCO_rf,GER_LOCO_rf,JAP_LOCO_rf)
rownames(model_LOCO_rf) <- c("AUS","CHI","FRA","GER","JAP")

Discov_LOCO <- get_LOCO_model(data_use,model_LOCO_rf)
write.csv(Discov_LOCO,paste(opts$workplace,"CRC_Discovery_",opts$output,"_all_LOCO_metric_single_com_final.csv",sep=''))

print("model construction done")


