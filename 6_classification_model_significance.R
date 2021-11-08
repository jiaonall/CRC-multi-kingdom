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

get_para <- function(i){
  model_rf <- data.frame()

  para_AUS <- read.csv(file = paste(opts$workplace,'CRC_Discovery_',opts$output,'_AUS_selfcv_repeat_max_para_com.csv',sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
  model_rf <- rbind(model_rf,para_AUS[i,])
  para_FRA <- read.csv(file = paste(opts$workplace,'CRC_Discovery_',opts$output,'_FRA_selfcv_repeat_max_para_com.csv',sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
  model_rf <- rbind(model_rf,para_FRA[i,])
  para_GER <- read.csv(file = paste(opts$workplace,'CRC_Discovery_',opts$output,'_GER_selfcv_repeat_max_para_com.csv',sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
  model_rf <- rbind(model_rf,para_GER[i,])
  para_CHI <- read.csv(file = paste(opts$workplace,'CRC_Discovery_',opts$output,'_CHI_selfcv_repeat_max_para_com.csv',sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
  model_rf <- rbind(model_rf,para_CHI[i,])
  para_JAP <- read.csv(file = paste(opts$workplace,'CRC_Discovery_',opts$output,'_JAP_selfcv_repeat_max_para_com.csv',sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
  model_rf <- rbind(model_rf,para_JAP[i,])
  rownames(model_rf) <- c('AUS','FRA','GER','CHI','JAP')
  return(model_rf)
}

sites <- c('AUS','FRA','GER','CHI','JAP')
for (s in sites){
  print(s)
  meta.s <- filter(meta.all, country==s)
  data.s <- data_use[meta.s$Run,]
  Group.s <- factor(data_use[meta.s$Run,]$Group)
  model_s <- get_para(1)
  timestart<-Sys.time()
  set.seed(2021)
  model.pval <- a3(Group.s ~., data = data.s,model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree=model_s[s,'ntree'],mtry=model_s[s,'mtry'],nodesize=model_s[s,'nodesize'],maxnodes=model_s[s,'maxnodes']))
  write.csv(model.pval$table,file = paste(opts$workplace,'CRC_Discovery_',opts$output,'_model_A3_sig_',s,'.csv',sep=''))
  timeend<-Sys.time()
  runningtime <- timeend-timestart
  print(runningtime)
}
print('All done')
		   


