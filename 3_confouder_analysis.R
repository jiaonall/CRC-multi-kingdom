# Packages
library("tidyverse")
library("coin")

###workplace
setwd('/Users/CRC_fungi/NewProfile/')
###load data
meta.all <- read.csv('../metadata/CRC_metadata_crcdis.csv',  stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                     check.name = FALSE)

feat.all <- read.csv('CRC_4domain_merged_species_rela.csv',stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                     check.name = FALSE)
feat.all <- t(feat.all)
# block for colonoscopy and study as well
meta.all <- meta.all %>%
  filter(!is.na(Sampling_rel_to_colonoscopy)) %>%
  mutate(block=ifelse(country!='CHI_CRC', country, 
                      paste0(country, '_', Sampling_rel_to_colonoscopy)))

feat.all <- feat.all[,meta.all$Run]
###
ss.disease <- apply(feat.all, 1, FUN=function(x, label){
  rank.x <- rank(x)/length(x)
  ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
  ss.o.i <- sum(vapply(unique(label), function(l){
    sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
  }, FUN.VALUE = double(1)))/length(rank.x)
  return(1-ss.o.i/ss.tot)
}, label=meta.all %>% pull(Group))


t.mean <- apply(feat.all, 1, mean, trim=0.05)

df.plot.all <- tibble(
  species=rownames(feat.all),
  disease=ss.disease,
  t.mean=t.mean)


# ##############################################################################
# Test all possible confounder variables
df.list <- list()
for (meta.var in c('age_factor','bmi_factor', 'Gender', 'country','Group','size_factor','Sampling_rel_to_colonoscopy')){
  
  cat('###############################\n', meta.var, '\n')
  metadata.c <- meta.all %>%
    filter(!is.na(eval(parse(text=meta.var))))
  
  cat('After filtering, the distribution of variables is:\n')
  print(table(metadata.c$Group, metadata.c %>% pull(meta.var)))
  print(table(metadata.c$Study))
  feat.red <- feat.all[,metadata.c$Run]
  
  cat('Calculating variance explained by meta-variable...\n')
  ss.var <- apply(feat.red, 1, FUN=function(x, label){
    rank.x <- rank(x)/length(x)
    ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
    ss.o.i <- sum(vapply(unique(label), function(l){
      sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
    }, FUN.VALUE = double(1)))/length(rank.x)
    return(1 - ss.o.i/ss.tot)
  }, label=metadata.c %>% pull(meta.var))
  df.plot.all[[meta.var]] <- ss.var
  
  cat('Calculating association with the meta-variable...\n')
  if (metadata.c %>% pull(meta.var) %>% unique %>% length > 2){
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      kruskal.test(x~as.factor(var))$p.value
    }, var=metadata.c %>% pull(meta.var))
  } else {
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      wilcox.test(x~as.factor(var))$p.value
    }, var=metadata.c %>% pull(meta.var))
  }
  meta.significance <- p.adjust(meta.significance, method='fdr')
  df.plot.all[[paste0(meta.var, '.significance')]] <- meta.significance
  cat('\n')
}

write.csv(df.plot.all,"../Results/CRC_confounder_4domain_all.csv")
