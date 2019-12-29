
```R
#######################################################

eQTL mapping in AA and EA samples in the GENOA study

#######################################################

> dim(AA_table)
[1] 14511338       12
> dim(EA_table)
[1] 8521801      12
> length(unique(as.character(AA_table$rs)))
[1] 6432684
> length(unique(as.character(EA_table$rs)))
[1] 3818520

> dim(gene_AA_anno_order)
[1] 17616    38
> dim(gene_EA_anno_order)
[1] 17360    38
> dim(gene_AA_anno_use)
[1] 17572    42
> dim(gene_EA_anno_use)
[1] 17343    42

> summary(gene_AA_anno_use$SNPnum)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    1.0   502.0   722.0   825.8   984.0 19808.0 
> sd(gene_AA_anno_use$SNPnum)
[1] 649.7872
>  summary(gene_EA_anno_use$SNPnum)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    1.0   282.5   418.0   491.4   589.0 13379.0 
> sd(gene_EA_anno_use$SNPnum)
[1] 410.7874


dim(gene_AA_anno_eGene)
dim(gene_EA_anno_eGene)

> dim(gene_AA_anno_eGene)
[1] 5475   42
> dim(gene_EA_anno_eGene)
[1] 4402   42

length(intersect(gene_AA_anno_eGene$GENE, gene_EA_anno_eGene$GENE))
> length(intersect(gene_AA_anno_eGene$GENE, gene_EA_anno_eGene$GENE))
[1] 3048

3048/(5475+4402-3048)
> 3048/(5475+4402-3048)
[1] 0.4463318

dim(AA_table_subset)
dim(EA_table_subset)

> dim(AA_table_subset)
[1] 354931     12
> dim(EA_table_subset)
[1] 371309     12

common_eSNP = merge(AA_table_subset, EA_table_subset, by=c("GENE","rs"))
dim(common_eSNP)
> dim(common_eSNP)
[1] 112316     22
112316/(354931+371309-112316)
> 112316/(354931+371309-112316)
[1] 0.1829477


   1    2    3    4    5    6    7    8    9 
3725 1203  368  104   49   14    7    4    1 
   1    2    3    4    5    6    7 
3577  690  108   18    7    1    1 


#######################################################

Characteristics of eGenes and eQTLs 

#######################################################


cor.test(gene_AA_anno_use$topSNP_af, abs(gene_AA_anno_use$topSNP_beta))$estimate
cor.test(gene_AA_anno_use$topSNP_af, abs(gene_AA_anno_use$topSNP_beta))$p.value


> cor.test(gene_AA_anno_use$topSNP_af, abs(gene_AA_anno_use$topSNP_beta))$estimate
       cor 
-0.4698748 
> cor.test(gene_AA_anno_use$topSNP_af, abs(gene_AA_anno_use$topSNP_beta))$p.value
[1] 0


cor.test(gene_EA_anno_use$topSNP_af, abs(gene_EA_anno_use$topSNP_beta))$estimate
cor.test(gene_EA_anno_use$topSNP_af, abs(gene_EA_anno_use$topSNP_beta))$p.value


> cor.test(gene_EA_anno_use$topSNP_af, abs(gene_EA_anno_use$topSNP_beta))$estimate
       cor 
-0.4595592 
> cor.test(gene_EA_anno_use$topSNP_af, abs(gene_EA_anno_use$topSNP_beta))$p.value
[1] 0

############
cor.test(gene_AA_anno_use$topSNP_af, -log10(gene_AA_anno_use$topSNP_pval))$estimate
cor.test(gene_AA_anno_use$topSNP_af, -log10(gene_AA_anno_use$topSNP_pval))$p.value

cor.test(gene_EA_anno_use$topSNP_af, -log10(gene_EA_anno_use$topSNP_pval))$estimate
cor.test(gene_EA_anno_use$topSNP_af, -log10(gene_EA_anno_use$topSNP_pval))$p.value


> cor.test(gene_AA_anno_use$topSNP_af, -log10(gene_AA_anno_use$topSNP_pval))$estimate
      cor 
0.2483079 
> cor.test(gene_AA_anno_use$topSNP_af, -log10(gene_AA_anno_use$topSNP_pval))$p.value
[1] 3.817783e-245
> 
> cor.test(gene_EA_anno_use$topSNP_af, -log10(gene_EA_anno_use$topSNP_pval))$estimate
      cor 
0.2045359 
> cor.test(gene_EA_anno_use$topSNP_af, -log10(gene_EA_anno_use$topSNP_pval))$p.value
[1] 3.543658e-163


common_eSNP$allele1.x = as.character(common_eSNP$allele1.x)
common_eSNP$allele0.x = as.character(common_eSNP$allele0.x)
common_eSNP$allele1.y = as.character(common_eSNP$allele1.y)
common_eSNP$allele0.y = as.character(common_eSNP$allele0.y)
ind1 = which(common_eSNP$allele1.x != common_eSNP$allele1.y)
common_eSNP[ind1,]$beta.x = -common_eSNP[ind1,]$beta.x 
cor.test(common_eSNP$beta.x, common_eSNP$beta.y)$estimate
cor.test(common_eSNP$beta.x, common_eSNP$beta.y)$p.value
> cor.test(common_eSNP$beta.x, common_eSNP$beta.y)$estimate
      cor 
0.9028881 
> cor.test(common_eSNP$beta.x, common_eSNP$beta.y)$p.value
[1] 0

sum(common_eSNP$beta.x*common_eSNP$beta.y>0)/dim(common_eSNP)[1]

> sum(common_eSNP$beta.x*common_eSNP$beta.y>0)/dim(common_eSNP)[1]
[1] 0.9734499

cor.test(-log10(common_eSNP$p_wald.x), -log10(common_eSNP$p_wald.y))$estimate
cor.test(-log10(common_eSNP$p_wald.x), -log10(common_eSNP$p_wald.y))$p.value
> cor.test(-log10(common_eSNP$p_wald.x), -log10(common_eSNP$p_wald.y))$estimate
      cor 
0.5871749 
> cor.test(-log10(common_eSNP$p_wald.x), -log10(common_eSNP$p_wald.y))$p.value
[1] 0


common_eSNP_maf = common_eSNP[-ind,]
cor.test(common_eSNP_maf$af.x, common_eSNP_maf$af.y)$estimate
cor.test(common_eSNP_maf$af.x, common_eSNP_maf$af.y)$p.value

> cor.test(common_eSNP_maf$af.x, common_eSNP_maf$af.y)$estimate
      cor 
0.3362597 
> cor.test(common_eSNP_maf$af.x, common_eSNP_maf$af.y)$p.value
[1] 0


library(dplyr)

AA_table_eQTL = AA_table_subset %>%  group_by(GENE) %>%
summarise(eSNPnum = length(p_wald))
AA_table_eQTL_num = merge(gene_AA_anno_eGene,AA_table_eQTL , by = "GENE")
cor.test(AA_table_eQTL_num$distance, AA_table_eQTL_num$eSNPnum)$estimate
cor.test(AA_table_eQTL_num$distance, AA_table_eQTL_num$eSNPnum)$p.value

> cor.test(AA_table_eQTL_num$distance, AA_table_eQTL_num$eSNPnum)$estimate
       cor 
0.04534255 
> cor.test(AA_table_eQTL_num$distance, AA_table_eQTL_num$eSNPnum)$p.value
[1] 0.0007907798

EA_table_eQTL = EA_table_subset %>%  group_by(GENE) %>%
summarise(eSNPnum = length(p_wald))
EA_table_eQTL_num = merge(gene_EA_anno_eGene,EA_table_eQTL , by = "GENE")
cor.test(EA_table_eQTL_num$distance, EA_table_eQTL_num$eSNPnum)$estimate
cor.test(EA_table_eQTL_num$distance, EA_table_eQTL_num$eSNPnum)$p.value

> cor.test(EA_table_eQTL_num$distance, EA_table_eQTL_num$eSNPnum)$estimate
       cor 
0.05194163 
> cor.test(EA_table_eQTL_num$distance, EA_table_eQTL_num$eSNPnum)$p.value
[1] 0.0005656972


cor.test(AA_table_eQTL_num$distance, AA_table_eQTL_num$eSNPnum/AA_table_eQTL_num$SNPnum)$estimate
cor.test(AA_table_eQTL_num$distance, AA_table_eQTL_num$eSNPnum/AA_table_eQTL_num$SNPnum)$p.value
cor.test(EA_table_eQTL_num$distance, EA_table_eQTL_num$eSNPnum/EA_table_eQTL_num$SNPnum)$estimate
cor.test(EA_table_eQTL_num$distance, EA_table_eQTL_num$eSNPnum/EA_table_eQTL_num$SNPnum)$p.value

> cor.test(AA_table_eQTL_num$distance, AA_table_eQTL_num$eSNPnum/AA_table_eQTL_num$SNPnum)$estimate
      cor 
-0.113886 
> cor.test(AA_table_eQTL_num$distance, AA_table_eQTL_num$eSNPnum/AA_table_eQTL_num$SNPnum)$p.value
[1] 2.854189e-17
> cor.test(EA_table_eQTL_num$distance, EA_table_eQTL_num$eSNPnum/EA_table_eQTL_num$SNPnum)$estimate
       cor 
-0.1196792 
> cor.test(EA_table_eQTL_num$distance, EA_table_eQTL_num$eSNPnum/EA_table_eQTL_num$SNPnum)$p.value
[1] 1.6281e-15


cor.test(AA_table_eQTL_num$distance, AA_table_eQTL_num$eSNPnum/AA_table_eQTL_num$distance)$estimate
cor.test(AA_table_eQTL_num$distance, AA_table_eQTL_num$eSNPnum/AA_table_eQTL_num$distance)$p.value
cor.test(EA_table_eQTL_num$distance, EA_table_eQTL_num$eSNPnum/EA_table_eQTL_num$distance)$estimate
cor.test(EA_table_eQTL_num$distance, EA_table_eQTL_num$eSNPnum/EA_table_eQTL_num$distance)$p.value


#########

conservation

#########

cor.test(gene_AA_anno_eGene$AA_phylop, gene_AA_anno_eGene$eSNPnum)$estimate
cor.test(gene_AA_anno_eGene$distance, gene_AA_anno_eGene$eSNPnum)$p.value

common_eGene = intersect(gene_AA_anno_eGene$GENE, gene_EA_anno_eGene$GENE)
gene_AA_anno_use$TYPE = "noneGene"
gene_AA_anno_use$TYPE[gene_AA_anno_use$eGene==1] = "AAspecific"
gene_AA_anno_use$TYPE[which(gene_AA_anno_use$GENE %in% common_eGene)] = "common"
gene_EA_anno_use$TYPE = "noneGene"
gene_EA_anno_use$TYPE[gene_EA_anno_use$eGene==1] = "EAspecific"
gene_EA_anno_use$TYPE[which(gene_EA_anno_use$GENE %in% common_eGene)] = "common"

x1=gene_AA_anno_use[,c(1,14,15,16)]
colnames(x1)=c("GENE", "dNdS","phylop","phastcon")
x2=gene_EA_anno_use[,c(1,14,15,16)]
colnames(x2)=c("GENE", "dNdS","phylop","phastcon")
gene_all = rbind(x1,x2)
gene_all = unique(gene_all)
allnoneGenes = intersect(gene_AA_anno_use$GENE[gene_AA_anno_use$eGene==0],gene_EA_anno_use$GENE[gene_EA_anno_use$eGene==0])
noneGenes = gene_all[which(gene_all$GENE %in% allnoneGenes),]

library(clinfun)

trendtest = list()
trendtest$Background = noneGenes$phylop
trendtest$EA_unique = gene_EA_anno_use$EA_phylop[gene_EA_anno_use$TYPE=="EAspecific"]
trendtest$AA_unique = gene_AA_anno_use$AA_phylop[gene_AA_anno_use$TYPE=="AAspecific"]
trendtest$Common = gene_AA_anno_use$AA_phylop[gene_AA_anno_use$TYPE=="common"]
pieces<-trendtest
n<-c(length(trendtest$Background),length(trendtest$EA_unique),length(trendtest$AA_unique),length(trendtest$Common))
grp<-as.ordered(factor(rep(1:length(n),n))) 
res = jonckheere.test(unlist(pieces),grp,alternative="decreasing")
# phylop
>  res$p.value
1.86292e-07


mean(trendtest$Background)
mean(trendtest$EA_unique)
mean(trendtest$AA_unique)
mean(trendtest$Common)

> mean(trendtest$Background)
[1] 0.255344
> mean(trendtest$EA_unique)
[1] 0.2065692
> mean(trendtest$AA_unique)
[1] 0.2278803
> mean(trendtest$Common)
[1] 0.1656216


trendtest = list()
trendtest$Background = noneGenes$phastcon
trendtest$EA_unique = gene_EA_anno_use$EA_phastcon[gene_EA_anno_use$TYPE=="EAspecific"]
trendtest$AA_unique = gene_AA_anno_use$AA_phastcon[gene_AA_anno_use$TYPE=="AAspecific"]
trendtest$Common = gene_AA_anno_use$AA_phastcon[gene_AA_anno_use$TYPE=="common"]
pieces<-trendtest
n<-c(length(trendtest$Background),length(trendtest$EA_unique),length(trendtest$AA_unique),length(trendtest$Common))
grp<-as.ordered(factor(rep(1:length(n),n))) 
res = jonckheere.test(unlist(pieces),grp,alternative="decreasing")
# phastcon
res$p.value
# phastcon
> res$p.value
[1] 1.599104e-35
mean(trendtest$Background)
mean(trendtest$EA_unique)
mean(trendtest$AA_unique)
mean(trendtest$Common)

> mean(trendtest$Background)
[1] 0.1508371
> mean(trendtest$EA_unique)
[1] 0.1345587
> mean(trendtest$AA_unique)
[1] 0.1334839
> mean(trendtest$Common)
[1] 0.1170731



trendtest = list()
trendtest$Background = noneGenes$dNdS
trendtest$EA_unique = gene_EA_anno_use$dNdS_ratio[gene_EA_anno_use$TYPE=="EAspecific"]
trendtest$AA_unique = gene_AA_anno_use$dNdS_ratio[gene_AA_anno_use$TYPE=="AAspecific"]
trendtest$Common = gene_AA_anno_use$dNdS_ratio[gene_AA_anno_use$TYPE=="common"]
pieces<-trendtest
n<-c(length(trendtest$Background),length(trendtest$EA_unique),length(trendtest$AA_unique),length(trendtest$Common))
grp<-as.ordered(factor(rep(1:length(n),n))) 
res = jonckheere.test(unlist(pieces),grp,alternative="increasing")
# dNdS
res$p.value

> res$p.value
[1] 2.936078e-10

mean(na.omit(trendtest$Background))
mean(na.omit(trendtest$EA_unique))
mean(na.omit(trendtest$AA_unique))
mean(na.omit(trendtest$Common))

> mean(na.omit(trendtest$Background))
[1] 0.1347604
> mean(na.omit(trendtest$EA_unique))
[1] 0.1274798
> mean(na.omit(trendtest$AA_unique))
[1] 0.1312477
> mean(na.omit(trendtest$Common))
[1] 0.1439989
> 




gene_AA_anno_use_tmp = gene_AA_anno_use[gene_AA_anno_use$PVEsd_bslmm_combined != 0,]
gene_AA_anno_use_tmp$PVE_bslmm_combined_z = gene_AA_anno_use_tmp$PVE_bslmm_combined/gene_AA_anno_use_tmp$PVEsd_bslmm_combined
gene_AA_anno_use_tmp$PVE_bslmm_combined_p  = 2*pnorm(abs(gene_AA_anno_use_tmp$PVE_bslmm_combined_z), lower.tail = F)
porder = gene_AA_anno_use_tmp$PVE_bslmm_combined_p[order(gene_AA_anno_use_tmp$PVE_bslmm_combined_p)]
PK = c(1:length(porder))*0.1/length(porder)
tail(porder[which(porder<=PK)])
combined_pve_FDR <- p.adjust(gene_AA_anno_use_tmp$PVE_bslmm_combined_p, method="BH")
gene_AA_anno_use_tmp$combined_pve_FDR = combined_pve_FDR


summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR<0.05])
summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR>=0.05])

> summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR<0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.08485 0.20108 0.24596 0.26676 0.30759 0.79457 
> summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR>=0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01418 0.04209 0.06257 0.07616 0.09977 0.27725 

> sd(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR<0.05])
[1] 0.1049609

> sum(gene_AA_anno_use_tmp$combined_pve_FDR<0.05)
[1] 1986

sum(gene_AA_anno_use_tmp$combined_pve_FDR<0.05)/dim(gene_AA_anno_use_tmp)[1]
[1] 0.1130722


gene_EA_anno_use_tmp = gene_EA_anno_use[gene_EA_anno_use$PVEsd_bslmm_combined != 0,]
gene_EA_anno_use_tmp$PVE_bslmm_combined_z = gene_EA_anno_use_tmp$PVE_bslmm_combined/gene_EA_anno_use_tmp$PVEsd_bslmm_combined
gene_EA_anno_use_tmp$PVE_bslmm_combined_p  = 2*pnorm(abs(gene_EA_anno_use_tmp$PVE_bslmm_combined_z), lower.tail = F)
porder = gene_EA_anno_use_tmp$PVE_bslmm_combined_p[order(gene_EA_anno_use_tmp$PVE_bslmm_combined_p)]
PK = c(1:length(porder))*0.1/length(porder)
tail(porder[which(porder<=PK)])
combined_pve_FDR <- p.adjust(gene_EA_anno_use_tmp$PVE_bslmm_combined_p, method="BH")
gene_EA_anno_use_tmp$combined_pve_FDR = combined_pve_FDR
sum(gene_EA_anno_use_tmp$combined_pve_FDR<0.05)/dim(gene_EA_anno_use_tmp)[1]
[1] 0.08305456


> summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR<0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.09813 0.22104 0.26780 0.28246 0.32222 0.78694 
> summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR>=0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01433 0.04669 0.06896 0.08395 0.10919 0.28467 


> sum(gene_EA_anno_use_tmp$combined_pve_FDR<0.05)
[1] 1440
sd(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR<0.05])
> sd(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR<0.05])
[1] 0.09881507


gene_AAEA = merge(gene_AA_anno_use, gene_EA_anno_use,by="GENE")
> dim(gene_AAEA)
[1] 17220    85

cor.test(gene_AAEA$PVE_bslmm_combined.x, gene_AAEA$PVE_bslmm_combined.y)$estimate
cor.test(gene_AAEA$PVE_bslmm_combined.x, gene_AAEA$PVE_bslmm_combined.y)$p.value
> cor.test(gene_AAEA$PVE_bslmm_combined.x, gene_AAEA$PVE_bslmm_combined.y)$estimate
      cor 
0.5666613 
> cor.test(gene_AAEA$PVE_bslmm_combined.x, gene_AAEA$PVE_bslmm_combined.y)$p.value
[1] 0

gene_AAEA_eGene  = merge(gene_AA_anno_eGene, gene_EA_anno_eGene,by="GENE")
cor.test(gene_AAEA_eGene$PVE_bslmm_combined.x, gene_AAEA_eGene$PVE_bslmm_combined.y)$estimate
cor.test(gene_AAEA_eGene$PVE_bslmm_combined.x, gene_AAEA_eGene$PVE_bslmm_combined.y)$p.value
> cor.test(gene_AAEA_eGene$PVE_bslmm_combined.x, gene_AAEA_eGene$PVE_bslmm_combined.y)$estimate
      cor 
0.5019758 
> cor.test(gene_AAEA_eGene$PVE_bslmm_combined.x, gene_AAEA_eGene$PVE_bslmm_combined.y)$p.value
[1] 2.68525e-194


x=unlist(gene_AA_anno_use$PVE_bslmm_combined[gene_AA_anno_use$eGene==1])
y=unlist(gene_AA_anno_use$PVE_bslmm_combined[gene_AA_anno_use$eGene==0])
wilcox.test(x, y, alternative = "greater")$p.value
summary(x)
summary(y)
> wilcox.test(x, y, alternative = "greater")$p.value
[1] 0
> summary(x)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.02427 0.09893 0.14194 0.16804 0.21080 0.79457 
> summary(y)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000089 0.0382187 0.0528382 0.0659095 0.0781831 0.7829927 

x=unlist(gene_EA_anno_use$PVE_bslmm_combined[gene_EA_anno_use$eGene==1])
y=unlist(gene_EA_anno_use$PVE_bslmm_combined[gene_EA_anno_use$eGene==0])
wilcox.test(x, y, alternative = "greater")$p.value
summary(x)
summary(y)

> wilcox.test(x, y, alternative = "greater")$p.value
[1] 0
> summary(x)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.03231 0.10269 0.14697 0.17132 0.21380 0.78694 
> summary(y)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01273 0.04350 0.06030 0.07631 0.09169 0.47233 

#save(gene_EA_anno_use, file = "gene_EA_anno_use.RData")
#save(gene_AA_anno_use, file = "gene_AA_anno_use.RData")


summary(gene_AA_anno_use$por_cis_combined)
sd(gene_AA_anno_use$por_cis_combined)

> summary(gene_AA_anno_use$por_cis_combined)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.000000 0.005726 0.010288 0.030894 0.032113 0.725383 
> sd(gene_AA_anno_use$por_cis_combined)
[1] 0.0544459


summary(gene_EA_anno_use$por_cis_combined)
sd(gene_EA_anno_use$por_cis_combined)

> summary(gene_EA_anno_use$por_cis_combined)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000003 0.0060379 0.0101317 0.0263195 0.0243969 0.7054134 
> sd(gene_EA_anno_use$por_cis_combined)
[1] 0.04845148


summary(gene_AA_anno_eGene$por_cis_combined)
sd(gene_AA_anno_eGene$por_cis_combined)


> summary(gene_AA_anno_eGene$por_cis_combined)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00385 0.03418 0.05364 0.07929 0.09372 0.72538 
> sd(gene_AA_anno_eGene$por_cis_combined)
[1] 0.07761329


summary(gene_EA_anno_eGene$por_cis_combined)
sd(gene_EA_anno_eGene$por_cis_combined)


> summary(gene_EA_anno_eGene$por_cis_combined)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.004951 0.033180 0.050897 0.076115 0.086598 0.705413 
> sd(gene_EA_anno_eGene$por_cis_combined)
[1] 0.07633261



cor.test(gene_AAEA_eGene$por_cis_combined.x, gene_AAEA_eGene$por_cis_combined.y)$estimate
cor.test(gene_AAEA_eGene$por_cis_combined.x, gene_AAEA_eGene$por_cis_combined.y)$p.value

> cor.test(gene_AAEA_eGene$por_cis_combined.x, gene_AAEA_eGene$por_cis_combined.y)$estimate
      cor 
0.6227689 
> cor.test(gene_AAEA_eGene$por_cis_combined.x, gene_AAEA_eGene$por_cis_combined.y)$p.value
[1] 0


cor.test(gene_AAEA_eGene$por_trans_combined.x, gene_AAEA_eGene$por_trans_combined.y)$estimate
cor.test(gene_AAEA_eGene$por_trans_combined.x, gene_AAEA_eGene$por_trans_combined.y)$p.value

> cor.test(gene_AAEA_eGene$por_trans_combined.x, gene_AAEA_eGene$por_trans_combined.y)$estimate
      cor 
0.1758956 
> cor.test(gene_AAEA_eGene$por_trans_combined.x, gene_AAEA_eGene$por_trans_combined.y)$p.value
[1] 1.327372e-22


> load("/net/mulan/home/shanglu/GENOA/analysis/conditional/gene_anno_AA_eGene.RData")
> load("/net/mulan/home/shanglu/GENOA/analysis/conditional/gene_anno_EA_eGene.RData")

summary(gene_anno_AA_eGene$pve_primary/gene_anno_AA_eGene$por_cis_combined)

> summary(gene_anno_AA_eGene$pve_primary/gene_anno_AA_eGene$por_cis_combined)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0189  0.4527  0.6452  0.6309  0.7926  9.2368 


summary(gene_anno_EA_eGene$pve_primary/gene_anno_EA_eGene$por_cis_combined)
> summary(gene_anno_EA_eGene$pve_primary/gene_anno_EA_eGene$por_cis_combined)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  0.01615   0.65905   0.80433   0.82259   0.91107 146.39379 


summary(gene_anno_AA_eGene$pve_indep/gene_anno_AA_eGene$por_cis_combined)
summary(gene_anno_EA_eGene$pve_indep/gene_anno_EA_eGene$por_cis_combined)
> summary(gene_anno_AA_eGene$pve_indep/gene_anno_AA_eGene$por_cis_combined)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.2436  0.6485  0.7783  0.7787  0.8907  4.2782 
> summary(gene_anno_EA_eGene$pve_indep/gene_anno_EA_eGene$por_cis_combined)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.3025  0.7548  0.8628  0.8596  0.9484  3.8640 


> sum(gene_anno_AA_eGene$indep_snp)
[1] 8070
> sum(gene_anno_EA_eGene$indep_snp)
[1] 5401


sum(gene_anno_AA_eGene$indep_snp)-dim(gene_anno_AA_eGene)[1]
> sum(gene_anno_AA_eGene$indep_snp)-dim(gene_anno_AA_eGene)[1]
[1] 2595

sum(gene_anno_EA_eGene$indep_snp)-dim(gene_anno_EA_eGene)[1]
> sum(gene_anno_EA_eGene$indep_snp)-dim(gene_anno_EA_eGene)[1]
[1] 999

> sum(gene_anno_AA_eGene$indep_snp==1)/dim(gene_anno_AA_eGene)[1]
[1] 0.6803653
sum(gene_anno_EA_eGene$indep_snp==1)/dim(gene_anno_EA_eGene)[1]
> sum(gene_anno_EA_eGene$indep_snp==1)/dim(gene_anno_EA_eGene)[1]
[1] 0.8125852


fisher.test(rbind(c(2595, 999),c((dim(gene_anno_AA_eGene)[1]-2595),(dim(gene_anno_EA_eGene)[1]-999))), alternative="greater")$p.value
> fisher.test(rbind(c(2595, 999),c((dim(gene_anno_AA_eGene)[1]-2595),(dim(gene_anno_EA_eGene)[1]-999))), alternative="greater")$p.value
[1] 8.530872e-146


sum(gene_anno_AA_eGene$indep_snp==2)/dim(gene_anno_AA_eGene)[1]
sum(gene_anno_EA_eGene$indep_snp==2)/dim(gene_anno_EA_eGene)[1]

> sum(gene_anno_AA_eGene$indep_snp==2)/dim(gene_anno_AA_eGene)[1]
[1] 0.219726
> sum(gene_anno_EA_eGene$indep_snp==2)/dim(gene_anno_EA_eGene)[1]
[1] 0.1567469

fisher.test(rbind(c(sum(gene_anno_AA_eGene$indep_snp==2), sum(gene_anno_EA_eGene$indep_snp==2)),c((dim(gene_anno_AA_eGene)[1]-sum(gene_anno_AA_eGene$indep_snp==2)),(dim(gene_anno_EA_eGene)[1]-sum(gene_anno_EA_eGene$indep_snp==2)))), alternative="greater")$p.value
[1] 1.108535e-15

sum(gene_anno_AA_eGene$indep_snp>=3)/dim(gene_anno_AA_eGene)[1]
sum(gene_anno_EA_eGene$indep_snp>=3)/dim(gene_anno_EA_eGene)[1]
> sum(gene_anno_AA_eGene$indep_snp>=3)/dim(gene_anno_AA_eGene)[1]
[1] 0.09990868
> sum(gene_anno_EA_eGene$indep_snp>=3)/dim(gene_anno_EA_eGene)[1]
[1] 0.03066788

fisher.test(rbind(c(sum(gene_anno_AA_eGene$indep_snp>=3), sum(gene_anno_EA_eGene$indep_snp>=3)),c((dim(gene_anno_AA_eGene)[1]-sum(gene_anno_AA_eGene$indep_snp>=3)),(dim(gene_anno_EA_eGene)[1]-sum(gene_anno_EA_eGene$indep_snp>=3)))), alternative="greater")$p.value
[1] 6.942577e-45

sum(gene_anno_AA_eGene$indep_snp)/dim(gene_anno_AA_eGene)[1]
sum(gene_anno_EA_eGene$indep_snp)/dim(gene_anno_EA_eGene)[1]
> sum(gene_anno_AA_eGene$indep_snp)/dim(gene_anno_AA_eGene)[1]
[1] 1.473973
> sum(gene_anno_EA_eGene$indep_snp)/dim(gene_anno_EA_eGene)[1]
[1] 1.226942


wilcox.test(gene_anno_AA_eGene$indep_snp, gene_anno_EA_eGene$indep_snp, alternative = "greater")$p.value
[1] 2.044394e-56



AAeGene = merge( AA_table_eQTL, gene_anno_AA_eGene, by="GENE")
EAeGene = merge( AA_table_eQTL, gene_anno_EA_eGene, by="GENE")

cor.test(AAeGene$indep_snp, AAeGene$eSNPnum)$estimate
cor.test(AAeGene$indep_snp, AAeGene$eSNPnum)$p.value
> cor.test(AAeGene$indep_snp, AAeGene$eSNPnum)$estimate
      cor 
0.4075766 
> cor.test(AAeGene$indep_snp, AAeGene$eSNPnum)$p.value
[1] 3.342018e-218


cor.test(EAeGene$indep_snp, EAeGene$eSNPnum)$estimate
cor.test(EAeGene$indep_snp, EAeGene$eSNPnum)$p.value

> cor.test(EAeGene$indep_snp, EAeGene$eSNPnum)$estimate
      cor 
0.1752742 
> cor.test(EAeGene$indep_snp, EAeGene$eSNPnum)$p.value
[1] 1.876865e-22



cor.test(log10(AAeGene$distance), AAeGene$indep_snp)$estimate
cor.test(log10(AAeGene$distance), AAeGene$indep_snp)$p.value
cor.test(log10(EAeGene$distance), EAeGene$indep_snp)$estimate
cor.test(log10(EAeGene$distance), EAeGene$indep_snp)$p.value

> cor.test(log10(AAeGene$distance), AAeGene$indep_snp)$estimate
        cor 
0.007617861 
> cor.test(log10(AAeGene$distance), AAeGene$indep_snp)$p.value
[1] 0.5730606
> cor.test(log10(EAeGene$distance), EAeGene$indep_snp)$estimate
       cor 
0.03410713 
> cor.test(log10(EAeGene$distance), EAeGene$indep_snp)$p.value
[1] 0.05972941


cor.test(AAeGene$PVE_bslmm_combined, AAeGene$indep_snp)$estimate
cor.test(AAeGene$PVE_bslmm_combined, AAeGene$indep_snp)$p.value
cor.test(EAeGene$PVE_bslmm_combined, EAeGene$indep_snp)$estimate
cor.test(EAeGene$PVE_bslmm_combined, EAeGene$indep_snp)$p.value
> cor.test(AAeGene$PVE_bslmm_combined, AAeGene$indep_snp)$estimate
      cor 
0.5683793 
> cor.test(AAeGene$PVE_bslmm_combined, AAeGene$indep_snp)$p.value
[1] 0
> cor.test(EAeGene$PVE_bslmm_combined, EAeGene$indep_snp)$estimate
      cor 
0.4095875 
> cor.test(EAeGene$PVE_bslmm_combined, EAeGene$indep_snp)$p.value
[1] 1.213508e-123


cor.test(AAeGene$dNdS_ratio, AAeGene$indep_snp)$estimate
cor.test(AAeGene$dNdS_ratio, AAeGene$indep_snp)$p.value
cor.test(EAeGene$dNdS_ratio, EAeGene$indep_snp)$estimate
cor.test(EAeGene$dNdS_ratio, EAeGene$indep_snp)$p.value


> cor.test(AAeGene$dNdS_ratio, AAeGene$indep_snp)$estimate
       cor 
0.04744324 
> cor.test(AAeGene$dNdS_ratio, AAeGene$indep_snp)$p.value
[1] 0.0005000723
> cor.test(EAeGene$dNdS_ratio, EAeGene$indep_snp)$estimate
        cor 
0.001289464 
> cor.test(EAeGene$dNdS_ratio, EAeGene$indep_snp)$p.value
[1] 0.9437463


cor.test(AAeGene$AA_phylop, AAeGene$indep_snp)$estimate
cor.test(AAeGene$AA_phylop, AAeGene$indep_snp)$p.value
cor.test(EAeGene$EA_phylop, EAeGene$indep_snp)$estimate
cor.test(EAeGene$EA_phylop, EAeGene$indep_snp)$p.value

> cor.test(AAeGene$AA_phylop, AAeGene$indep_snp)$estimate
        cor 
-0.09151548 
> cor.test(AAeGene$AA_phylop, AAeGene$indep_snp)$p.value
[1] 1.16715e-11
> cor.test(EAeGene$EA_phylop, EAeGene$indep_snp)$estimate
       cor 
-0.1059256 
> cor.test(EAeGene$EA_phylop, EAeGene$indep_snp)$p.value
[1] 4.570181e-09

cor.test(AAeGene$AA_phastcon, AAeGene$indep_snp)$estimate
cor.test(AAeGene$AA_phastcon, AAeGene$indep_snp)$p.value
cor.test(EAeGene$EA_phastcon, EAeGene$indep_snp)$estimate
cor.test(EAeGene$EA_phastcon, EAeGene$indep_snp)$p.value
> cor.test(AAeGene$AA_phastcon, AAeGene$indep_snp)$estimate
        cor 
-0.07797363 
> cor.test(AAeGene$AA_phastcon, AAeGene$indep_snp)$p.value
[1] 7.604382e-09
> cor.test(EAeGene$EA_phastcon, EAeGene$indep_snp)$estimate
        cor 
-0.08586032 
> cor.test(EAeGene$EA_phastcon, EAeGene$indep_snp)$p.value
[1] 2.063952e-06








```
