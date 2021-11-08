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
    make_option(c("-p", "--pval_profile"), type="character", default="Results/CRC_Bacteria_feature.csv",
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
  print(paste("The MMUPHin results file ", opts$pval_profile,  sep = ""))
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

###
meta.all <- read.csv(file = paste(opts$workplace,opts$metadata,sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
feat.all <- read.csv(file = paste(opts$workplace,opts$feature_profile,sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
pval.all <- read.csv(file = paste(opts$workplace,opts$pval_profile,sep=''),stringsAsFactors = FALSE, header = TRUE, row.names =1, check.name = FALSE)
#pval.sig <- subset(pval.all,pval<0.05)
diff <- rownames(pval.all)
feat.all[is.na(feat.all)] <- 0
rownames(meta.all) <- meta.all$Run
meta.all$StudyID <- factor(meta.all$country)
feat.sig <- t(feat.all[diff,meta.all$Run])
data.sig <- cbind(feat.sig,data.frame(meta.all[rownames(feat.sig),'Group']))
colnames(data.sig)[ncol(data.sig)] <- 'Group'
Group <- factor(data.sig$Group)

###feature_selection
timestart<-Sys.time()
set.seed(1)
boruta <- Boruta(x=feat.sig, y=Group, pValue=0.05, mcAdj=T,
                 maxRuns=1000)
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)
boruta
print(table(boruta$finalDecision))
#eatract feature
boruta.finalVarsWithTentative <- data.frame(Item=getSelectedAttributes(boruta, withTentative = F), Type="Boruta_with_tentative")
feature <- boruta.finalVarsWithTentative$Item
boruta.variable.imp.use <- boruta$ImpHistory[,feature]
feature_importance <- apply(boruta.variable.imp.use,2,mean)
feature_importance <- data.frame(sort(feature_importance,decreasing = TRUE))
feature <- rownames(feature_importance)
write.csv(feature_importance,paste(opts$workplace,'Results/CRC_',opts$output,'_feature_com.csv',sep=''))

##importance
boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]

  variableGrp <- data.frame(Variable=names(x$finalDecision),
                            finalDecision=x$finalDecision)

  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))

  variableGrp <- rbind(variableGrp, showGrp)

  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)

  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>%
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)


  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)

  invisible(boruta.variable.imp)
}
boruta.variable.imp <- boruta.imp(boruta)
#head(boruta.variable.imp)
feature_impor_plot <- sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
           legend_variable = "finalDecision", legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
           xtics_angle = 90, coordinate_flip =T)

pdf(paste(opts$workplace,'Plot/CRC_Discovery_',opts$output,'_feature_selection_importance_com.pdf',sep=''), useDingbats = FALSE,width = 8, height = 6)
plot(feature_impor_plot)
dev.off()
##correlation
up_CorMatrix <- function(cor,p) {ut <- upper.tri(cor) 
data.frame(row = rownames(cor)[row(cor)[ut]] ,
           column = rownames(cor)[col(cor)[ut]], 
           cor =(cor)[ut] ) }

corrplot(res$r, type="upper",  p.mat = res$P, sig.level = 0.01, 
         insig = "blank",tl.cex=0.4,tl.col="black",addCoef.col = "grey",addCoefasPercent = TRUE,number.cex=0.4)


print('feature selection done')

