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
    make_option(c("-s", "--select_feature"), type="character", default="Results/CRC_Bacteria_feature_com.csv",
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

package_list = c("dplyr","ggplot2","randomForest","caret","A3","Boruta","ROCR","corrplot","vegan","Hmisc")

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
print(feature)
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

### cv fold
set.seed(123)
tune_model <- function(data,Group,num){
  print("tune model")
  set.seed(2021*num)
  control <- trainControl(method="repeatedcv", number = 5,classProbs=T,summaryFunction=twoClassSummary)
  tunegrid <- expand.grid(.mtry=c(1:ceiling(sqrt(length(feature)))),
                          .ntree=c(301,501,801,1001),
                          .nodesize=c(50,100,150),
                          .maxnodes=c(5,10,15,20))
  rf_default <- train(as.factor(Group)~., data=data,
                      method=customRF,
                      metric="ROC",
                      preProcess = c("center", "scale"),
                      tuneGrid=tunegrid,
                      trControl=control)
  best_rf <- rf_default$bestTune
  print(rf_default)
  #print("get_model metrix") 
  return(best_rf)
}

get_self_cv <- function(data,Group,best_rf,num){
  set.seed(2021*num)
  fold <- createFolds(y = Group, k=5)
  #return(fold)
  metrics <- matrix(NA,5,12)
  colnames(metrics) <- c("Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision",
                         "Recall","F1","Prevalence","Detection Rate","Detection Prevalence",
                         "Balanced Accuracy","AUC")
  for(j in 1:5){
    #print(j)
    fold_test <- data[fold[[j]],]
    fold_train <- data[-fold[[j]],]
    print(table(fold_train$Group))
    Group <- fold_train$Group
    set.seed(0518)
    fold_fit <- randomForest(as.factor(Group)~., data=fold_train,mtry=best_rf$mtry,
                             ntree=best_rf$ntree,importance=TRUE)
    fold_pred <- predict(fold_fit,fold_test)
    #print(fold_pred)
    result <- confusionMatrix(factor(as.vector(fold_pred)),as.factor(fold_test$Group),mode = "everything",positive='CRC')
    metrics[j,1:11] <- as.numeric(result$byClass)
    predicted_probs <- predict(fold_fit, fold_test, type = 'prob')
    pred <- prediction(predicted_probs[, 'CTR'], fold_test$Group)
    auc <- performance(pred, 'auc')
    #print(auc@y.values[[1]])
    metrics[j,12] <- auc@y.values[[1]]

  }
  best_index <- which(metrics[,'AUC']==max(metrics[,'AUC']))[1]
  fold_test_best <- data[fold[[best_index]],]
  fold_train_best <- data[-fold[[best_index]],]
  Group <- fold_train_best$Group
  best_model <- randomForest(as.factor(Group)~., data=fold_train_best,mtry=best_rf$mtry,
                             ntree=best_rf$ntree,importance=TRUE)
  result_list <- list("model" = best_model,"metrics" = metrics,"best_fold"=best_index)
  return(result_list)
}

sites <- c('AUS','FRA','GER','CHI','JAP')
for (s in sites){
  print(s)
  meta.s <- filter(meta.all, country==s)
  data.s <- data_use[meta.s$Run,]
  Group.s <- factor(data_use[meta.s$Run,]$Group)
  model_s <- matrix(NA,20,4)
  colnames(model_s) <- c('mtry','ntree','nodesize','maxnodes')
  s_matrix <- matrix(NA,100,12)
  colnames(s_matrix) <- c("Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision",
                          "Recall","F1","Prevalence","Detection Rate","Detection Prevalence",
                          "Balanced Accuracy","AUC")

  s_matrix_max <- matrix(NA,20,12)
  colnames(s_matrix_max) <- c("Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision",
                          "Recall","F1","Prevalence","Detection Rate","Detection Prevalence",
                          "Balanced Accuracy","AUC")
  for (i in (1:20)){
    print(i)
    s_rf <- tune_model(data.s,Group.s,i)
    model_s[i,] <- c(s_rf$mtry,s_rf$ntree,s_rf$nodesize,s_rf$maxnodes)
    s_AUC <- get_self_cv(data.s,Group.s,best_rf=s_rf,i)
    s_matrix[(i*5-4):(i*5),] <- s_AUC$metrics
    s_matrix_max[i,] <- s_AUC$metrics[s_AUC$best_fold,]
    mean(s_matrix_max[,'AUC'])
    write.csv(s_matrix,paste(opts$workplace,"CRC_Discovery_",opts$output,"_selfcv_repeat_metric_",s,".csv",sep=''))
    write.csv(s_matrix_max,paste(opts$workplace,"CRC_Discovery_",opts$output,"_selfcv_repeat_max_",s,".csv",sep=''))
    write.csv(model_s,paste(opts$workplace,"CRC_Discovery_",opts$output,"_selfcv_repeat_max_para_",s,".csv",sep=''))
  }
}

print('Self repeat done')
