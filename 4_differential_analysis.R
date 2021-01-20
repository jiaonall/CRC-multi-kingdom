library("tidyverse")
library("coin")
library("pROC")

setwd('CRC_fungi/NewProfile/')

meta.all <- read.csv(file = '../metadata/CRC_metadata_crcdis.csv',stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                     check.name = FALSE)
feat.all <- read.csv("CRC_Fungi_rela_select_CHI_JAP_abunt.csv",stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                     check.name = FALSE)

meta.all <- filter(meta.all,Group != 'adenoma' & Group != 'HS' )


# block for colonoscopy and study as well
meta.all <- meta.all %>%
  filter(!is.na(Sampling_rel_to_colonoscopy)) %>%
  mutate(block=ifelse(country!='CHI', country, 
                      paste0(country, '_', Sampling_rel_to_colonoscopy)))


feat.all <- feat.all[,meta.all$Run]

#feat.trans <- asin(sqrt(feat.all))
feat.trans <- feat.all
stopifnot(all(meta.all$Run %in% colnames(feat.all)))

sites <- meta.all %>% pull(country) %>% unique

# calculate pval, gFC, AUROC
log.n0 <- 1e-05
p.val <- matrix(NA, nrow=nrow(feat.trans), ncol=length(sites)+1, 
                dimnames=list(row.names(feat.trans), c(sites, 'all')))
fc <- p.val
fc2 <- fc
aucs.mat <- p.val
p.adj <- p.val
aucs.all  <- vector('list', nrow(feat.trans))

cat("Calculating effect size for every feature...\n")
pb <- txtProgressBar(max=nrow(feat.trans), style=3)

# caluclate wilcoxon test and effect size for each feature and study
for (f in row.names(feat.trans)) {
  # for each sites
  for (s in sites) {
    
    x <- as.numeric(feat.trans[f, meta.all %>% filter(country==s) %>% 
                                 filter(Group=='CRC') %>% pull(Run)])
    y <- as.numeric(feat.trans[f, meta.all %>% filter(country==s) %>% 
                                 filter(Group=='CTR') %>% pull(Run)])
    
    # Wilcoxon
    p.val[f,s] <- wilcox.test(x, y, exact=FALSE)$p.value
    #p.adj[f,s] <- p.adjust(wilcox.test(x, y, exact=FALSE)$p.value, method='fdr')
    # AUC
    aucs.all[[f]][[s]] <- c(roc(controls=y, cases=x, 
                                direction='<', ci=TRUE, auc=TRUE)$ci)
    aucs.mat[f,s] <- c(roc(controls=y, cases=x, 
                           direction='<', ci=TRUE, auc=TRUE)$ci)[2]
    
    # FC
    q.p <- quantile(log10(x+log.n0), probs=seq(.1, .9, .05))
    q.n <- quantile(log10(y+log.n0), probs=seq(.1, .9, .05))
    fc[f, s] <- sum(q.p - q.n)/length(q.p)
    ## FC_new
    fc2[f,s] <- log2((median(x)+log.n0)/(median(y)+log.n0))
  }
  
  # calculate effect size for all studies combined
  # Wilcoxon + blocking factor
  d <- data.frame(y=as.numeric(feat.trans[f,]), 
                  x=meta.all$Group, block=as.factor(meta.all$block))
  p.val[f,'all'] <- pvalue(wilcox_test(y ~ x | block, data=d))
  p.adj[f,'all'] <- p.adjust(as.numeric(p.val[f,'all']), method='fdr')
  # other metrics
  x <- as.numeric(feat.trans[f, meta.all %>% filter(Group=='CRC') %>% pull(Run)])
  y <- as.numeric(feat.trans[f, meta.all %>% filter(Group=='CTR') %>% pull(Run)])
  # FC
  fc[f, 'all'] <- mean(fc[f, sites])
  fc2[f,'all'] <- log2((median(x)+log.n0)/(median(y)+log.n0))
  # AUC
  aucs.mat[f,'all'] <- c(roc(controls=y, cases=x, 
                             direction='<', ci=TRUE, auc=TRUE)$ci)[2]
  
  # progressbar
  setTxtProgressBar(pb, (pb$getVal()+1))
}
cat('\n')

## statistic
feat.crc <- feat.all[,meta.all %>% filter(Group=='CRC') %>% pull(Run)]
feat.ctr <- feat.all[,meta.all %>% filter(Group=='CTR') %>% pull(Run)]
stat <- cbind(apply(feat.crc,1,mean),apply(feat.ctr,1,mean),
              apply(feat.crc,1,sd),apply(feat.ctr,1,sd))
colnames(stat) <- c('CRC_mean','CTR_mean','CRC_sd','CTR_sd')
rownames(stat) <- rownames(feat.trans)
 multiple hypothesis correction
p.adj <- data.frame(apply(p.val, MARGIN=2, FUN=p.adjust, method='fdr'),
                    check.names = FALSE)


stat_p <- cbind(p.adj,stat,fc)


# ##############################################################################
# save results
write.csv(fc2, file='../Results/CRC_Fungi_block_logfc_abunt.csv')
write.csv(stat_p, file='../Results/CRC_Fungi_block_stats_abunt.csv')
write.csv(fc, file='../Results/CRC_Fungi_block_fc_abunt.csv')
write.csv(p.adj, file = '../Results/CRC_Fungi_padj_block_abunt.csv')
write.csv(p.val, '../Results/CRC_Fungi_pvalue_block_abunt.csv')
write.csv(aucs.mat, '../Results/CRC_Fungi_aucs_block_abunt.csv')
