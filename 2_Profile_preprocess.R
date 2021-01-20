###Packages
library("tidyverse")
library("coin")
library("pROC")

###workplace
setwd('xxxx')

###load data
file_count <- read.csv('CRC_merged_species.csv',stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                         check.name = FALSE)
###convert relative abundance
All_counts <- data.frame(apply(file_count,2, sum))
colnames(All_counts) <- 'Sum'

file_rela <- sweep(file_count,2,All_counts[match(colnames(file_count),rownames(All_counts)),"Sum"], "/")
write.csv(file_rela,'CRC_merged_rela.csv')

### remove rare species

# read in raw data

feat.all.rela <- read.csv("CRC_Bacteria_merged_species_rela.csv",stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                          check.name = FALSE)

feat.all.rela = as.matrix(feat.all.rela)
cat('Successfully loaded data...\n')
cat('Feature dimensions:', dim(feat.all.rela), '\n')

# metadata

meta.all <- read.csv(file = '../metadata/CRC_metadata_crcdis_CHI_JAP_com_delect.csv',stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                     check.name = FALSE)


meta.all <- filter(meta.all,Group != 'adenoma' &Group != 'HS')
feat.all.rela <- feat.all.rela[,meta.all$Run]
dim(feat.all.rela)

filter.f <- function(dat, Num){
  SD <- apply(dat,1,sd)
  num_0 <- apply(dat, 1, function(x) length(which(x == 0)))
  ave_abun <- apply(dat,1,mean)
  tmp <- cbind(dat,data.frame(SD),data.frame(num_0),data.frame(ave_abun))
  colnames(tmp)[(ncol(dat)+1):(ncol(dat)+3)] <- c("sd","count0",'aveabun')
  #dat_filter <- tmp %>% filter(count0 <= as.numeric(Num*0.9) & sd >0) 
  dat_filter <- tmp[(tmp$count0 <= as.numeric(Num*0.8)) & (tmp$sd >0) & (tmp$aveabun > 0.0001),]
  #dat_filter <- tmp[(tmp$count0 <= as.numeric(Num*0.8)) & (tmp$sd >0),]
  return(dat_filter[,1:ncol(dat)])
}


temp.mean = function(data,meta){t(sapply(row.names(data),
                                         FUN=function(marker){sapply(unique(meta$country),
                                                                     FUN=function(site, marker){
                                                                       mean.ab = mean(data[marker, which(meta$country == site)])
                                                                     },
                                                                     marker=marker)}))
}

temp.mean.crc <- temp.mean(feat.all.rela,meta.all)

f.idx.crc = rowSums(temp.mean.crc >= 0.00001) >= 3 &
  row.names(feat.all.rela) != '-1'

feat.crc_filter = feat.all.rela[f.idx.crc,]
cat('Retaining', sum(f.idx.crc), 'features after low-abundance filtering...\n')

feat.filter <- filter.f(feat.crc_filter,ncol(feat.crc_filter))
dim(feat.filter)

# save cleaned data
fn.filter <- '../Result_new/CRC_Bacteria_rela_select_CHI_JAP_abunt.csv'
write.csv(feat.filter, file=fn.filter)
