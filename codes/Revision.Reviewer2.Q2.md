

# Page 9 eQTL detection: 
The authors detect an overlap of ~30% of eSNPs between populations. 
I am not sure that I would refer to that overlap as small as the authors did. 
The eQTL detection is based on a P-value (or better FDR) cutoff, 
which is highly dependent (as also the authors acknowledge) on the MAF. 
I am wondering how much of that overlap - or lack thereof - is simply caused by allele-frequency differences 
between population, 
i.e. a significant eQTL in population A is not found in population B because of a much smaller MAF for the SNP in B 
(The authors describe that even on page 10 showing that eSNPs have larger Fsts). 

First, we displayed -log10(p-value) difference between the two populations vs MAF difference between the two populations. The results show that…. 
```R

setwd("/net/mulan/home/shanglu/GENOA/analysis/figure")

pathAA="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
pathEA="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"

load(paste0(pathAA,"/AA_table.RData"))
load(paste0(pathEA,"/EA_table.RData"))

dim(AA_table)
dim(EA_table)


AA_thr = 6.245907e-05
EA_thr = 0.0001385504

AA_table$rs = as.character(AA_table$rs)
AA_table$GENE = as.character(AA_table$GENE)
EA_table$rs = as.character(EA_table$rs)
EA_table$GENE = as.character(EA_table$GENE)

AA_table_subset = AA_table[AA_table$p_wald <= AA_thr,]
EA_table_subset = EA_table[EA_table$p_wald <= EA_thr,]

common_table_all = merge(AA_table, EA_table, by = c("GENE","rs"))

common_table_all$allele1.x = as.character(common_table_all$allele1.x)
common_table_all$allele0.x = as.character(common_table_all$allele0.x)
common_table_all$allele1.y = as.character(common_table_all$allele1.y)
common_table_all$allele0.y = as.character(common_table_all$allele0.y)
ind1 = which(common_table_all$allele1.x == common_table_all$allele1.y)
ind2 = which(common_table_all$allele0.x == common_table_all$allele0.y)
ind = intersect(ind1, ind2)
common_table_use_all = common_table_all[ind,]

dim(common_table_all)
dim(common_table_use_all)
#> dim(common_table_all)
#[1] 6060046      22
#> dim(common_table_use_all)
#[1] 5171951      22

common_table_use_all$log10p.x = -log10(common_table_use_all$p_wald.x)
common_table_use_all$log10p.y = -log10(common_table_use_all$p_wald.y)

log10p_diff = common_table_use_all$log10p.x - common_table_use_all$log10p.y
maf_diff = common_table_use_all$af.x - common_table_use_all$af.y

/*
> summary(log10p_diff)
      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
-136.19928   -0.46620   -0.01764   -0.10295    0.38975  177.91264 
> summary(maf_diff)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-0.43300 -0.09600 -0.01700 -0.01053  0.07000  0.48600 

cor.test(log10p_diff, maf_diff, method = "spearman")$estimate
cor.test(log10p_diff, maf_diff, method = "spearman")$p.value
       rho 
0.08418018 
> cor.test(log10p_diff, maf_diff, method = "spearman")$p.value
[1] 0

*/

#----------------
# scatter plot
#----------------
library(ggplot2)

# pdf too large, need to use tiff to tune the resolution
path_R2Q2 = "/net/mulan/home/shanglu/GENOA/analysis/revision/R2Q2"

mydata = data.frame(log10p_diff, maf_diff)
tiff(paste0(path_R2Q2,"/diff_log10p_maf.tiff"), units="in", width=8, height=8, res=150)
ggplot(mydata, aes(x=maf_diff, y=log10p_diff)) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=0.5)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  geom_abline(intercept = 0)+
  #ylim(0,1)+
  #xlim(0,1)+
  labs(title=paste0(""),
       x="MAF diff between AA and EA", y = "-log10(p-value) diff between AA and EA")
dev.off()



```


Second, we draw boxplots or violin plots displaying the distribution of MAF difference between the two populations for common eSNPs, population unique eSNPs, and non-eSNPs. The results show that…  

```R
common_table_use_all$maf_diff = common_table_use_all$af.x - common_table_use_all$af.y

rsAA_significant = unique(common_table_use_all$rs[common_table_use_all$p_wald.x <= AA_thr]) # using updated FDR
rsEA_significant = unique(common_table_use_all$rs[common_table_use_all$p_wald.y <= EA_thr])  # using updated FDR

common_table_use_all$signif.x = 0
common_table_use_all$signif.x[which(common_table_use_all$p_wald.x <= AA_thr)]=1
common_table_use_all$signif.y = 0
common_table_use_all$signif.y[which(common_table_use_all$p_wald.y <= EA_thr)]=1
common_table_use_all$signif_x_y = common_table_use_all$signif.x + common_table_use_all$signif.y
common_table_use_all$spec_AA = common_table_use_all$signif.x * common_table_use_all$signif_x_y
common_table_use_all$spec_AA = common_table_use_all$signif.x * common_table_use_all$signif_x_y
common_table_use_all$spec_AA[common_table_use_all$spec_AA==2] = 0
common_table_use_all$spec_EA = common_table_use_all$signif.y * common_table_use_all$signif_x_y
common_table_use_all$spec_EA[common_table_use_all$spec_EA==2] = 0


> table(common_table_use_all$signif_x_y)
      0       1       2 
4832908  251680   87363 

> sum(common_table_use_all$spec_AA)
[1] 92725
#> sum(common_table_use_all$spec_EA)
# [1] 158955


common_table_use_all$SNP_type = ""
common_table_use_all$SNP_type[which(common_table_use_all$signif_x_y==0)] = "non-eSNPs"
common_table_use_all$SNP_type[which(common_table_use_all$signif_x_y==2)] = "common eSNPs"
common_table_use_all$SNP_type[which(common_table_use_all$spec_AA==1)] = "AA specific"
common_table_use_all$SNP_type[which(common_table_use_all$spec_EA==1)] = "EA specific"

#> table(common_table_use_all$SNP_type)
# AA specific common eSNPs  EA specific    non-eSNPs 
#       92725        87363       158955      4832908 
      
common_table_use_all$maf_diff = common_table_use_all$af.x -  common_table_use_all$af.y

MAF_diff = common_table_use_all$maf_diff
SNP_type = common_table_use_all$SNP_type
dat = data.frame(MAF_diff, SNP_type)
Class = factor(SNP_type, levels = c("non-eSNPs","AA specific","EA specific","common eSNPs"),order=T)
dat = data.frame(MAF_diff,SNP_type, Class)


tiff(paste0(path_R2Q2,"/maf_diff_4groups_SNPs.tiff"), units="in", width=8, height=8, res=150)
dp <- ggplot(dat, aes(x=Class, y=MAF_diff)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="MAF difference between AA and EA",x="", y = "Difference in MAF")+
  scale_fill_brewer(palette="RdBu") + 
  theme_minimal(base_size = 22)
dp
dev.off()

  library(dplyr)

maf_diff_4groups = common_table_use_all %>%  group_by(SNP_type) %>%
summarise(average_maf_diff = mean(maf_diff), 
          largest_maf_diff = max(maf_diff),
	  smallest_maf_diff = min(maf_diff))
	  
maf_diff_4groups	  

#> maf_diff_4groups  
## A tibble: 4 x 4
#  SNP_type     average_maf_diff largest_maf_diff smallest_maf_diff
#  <chr>                   <dbl>            <dbl>             <dbl>
#1 AA specific           0.0580             0.482            -0.429
#2 common eSNPs         -0.00287            0.467            -0.42 
#3 EA specific          -0.0554             0.468            -0.426
#4 non-eSNPs            -0.0105             0.486            -0.433


wilcox.test(common_table_use_all$maf_diff[common_table_use_all$SNP_type=="non-eSNPs"],  alternative="two.sided")$p.value
wilcox.test(common_table_use_all$maf_diff[common_table_use_all$SNP_type=="AA specific"],  alternative="two.sided")$p.value
wilcox.test(common_table_use_all$maf_diff[common_table_use_all$SNP_type=="EA specific"],  alternative="two.sided")$p.value
wilcox.test(common_table_use_all$maf_diff[common_table_use_all$SNP_type=="common eSNPs"], alternative="two.sided")$p.value

/*
> wilcox.test(common_table_use_all$maf_diff[common_table_use_all$SNP_type=="non-eSNPs"],  alternative="two.sided")$p.value

[1] 0
> wilcox.test(common_table_use_all$maf_diff[common_table_use_all$SNP_type=="AA specific"],  alternative="two.sided")$p.value
[1] 0
> wilcox.test(common_table_use_all$maf_diff[common_table_use_all$SNP_type=="EA specific"],  alternative="two.sided")$p.value
[1] 0
> wilcox.test(common_table_use_all$maf_diff[common_table_use_all$SNP_type=="common eSNPs"], alternative="two.sided")$p.value
[1] 9.152245e-27
*/


```

Third, we display MAF in AA vs MAF in EA for all four categories of eSNPs. 

```R
library(patchwork)
mydata = common_table_use_all[common_table_use_all$SNP_type == "non-eSNPs",]
p1 = ggplot(mydata, aes(x=af.x, y=af.y)) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=0.5)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  geom_abline(intercept = 0)+
  ylim(0,0.5)+
  xlim(0,0.5)+
  labs(title=paste0("non-eSNPs"),
       x="MAF in AA", y = "MAF in EA")

mydata = common_table_use_all[common_table_use_all$SNP_type == "AA specific",]
p2 = ggplot(mydata, aes(x=af.x, y=af.y)) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=0.5)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  geom_abline(intercept = 0)+
  ylim(0,0.5)+
  xlim(0,0.5)+
  labs(title=paste0("AA specific"),
       x="MAF in AA", y = "MAF in EA")
       
mydata = common_table_use_all[common_table_use_all$SNP_type == "EA specific",]
p3 = ggplot(mydata, aes(x=af.x, y=af.y)) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=0.5)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  geom_abline(intercept = 0)+
  ylim(0,0.5)+
  xlim(0,0.5)+
  labs(title=paste0("EA specific"),
       x="MAF in AA", y = "MAF in EA")
       
mydata = common_table_use_all[common_table_use_all$SNP_type == "common eSNPs",]
p4 = ggplot(mydata, aes(x=af.x, y=af.y)) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=0.5)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  geom_abline(intercept = 0)+
  ylim(0,0.5)+
  xlim(0,0.5)+
  labs(title=paste0("common eSNPs"),
       x="MAF in AA", y = "MAF in EA")
    
tiff(paste0(path_R2Q2,"/maf_in_4groups_SNPs_patch.tiff"), units="in", width=8, height=8, res=150)    
(p1 | p2 ) /
      (p3 | p4)
dev.off()

# show 2d density plot should be better

tiff(paste0(path_R2Q2,"/MAF_in_noneSNPs.tiff"), units="in", width=7, height=6, res=150) 
mydata = common_table_use_all[common_table_use_all$SNP_type == "non-eSNPs",]      
p5 = ggplot(mydata, aes(x=af.x, y=af.y)) +    
  geom_point(size=0.1)+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=1) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs(title=paste0("non-eSNPs"),
       x="MAF in AA", y = "MAF in EA")+
  theme_bw(base_size=22)+  
  xlim(0,0.5)+
  ylim(0,0.5) 
p5
dev.off()


mydata = common_table_use_all[common_table_use_all$SNP_type == "AA specific",]
tiff(paste0(path_R2Q2,"/MAF_in_AAspecificSNP.tiff"), units="in", width=7, height=6, res=150) 
p6 = ggplot(mydata, aes(x=af.x, y=af.y)) +    
  geom_point(size=0.1)+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=1) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs(title=paste0("AA specific"),
       x="MAF in AA", y = "MAF in EA")+
  theme_bw(base_size=22)+       
  xlim(0,0.5)+
  ylim(0,0.5) 
p6
dev.off()

 
mydata = common_table_use_all[common_table_use_all$SNP_type == "EA specific",]
tiff(paste0(path_R2Q2,"/MAF_in_EAspecificSNP.tiff"), units="in", width=7, height=6, res=150) 
p7 = ggplot(mydata, aes(x=af.x, y=af.y)) +    
  geom_point(size=0.1)+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=1) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs(title=paste0("EA specific"),
       x="MAF in AA", y = "MAF in EA")+
  theme_bw(base_size=22)+      
  xlim(0,0.5)+
  ylim(0,0.5) 
p7
dev.off()
  
mydata = common_table_use_all[common_table_use_all$SNP_type == "common eSNPs",]
tiff(paste0(path_R2Q2,"/MAF_in_commonSNP.tiff"), units="in", width=7, height=6, res=150) 
p8 = ggplot(mydata, aes(x=af.x, y=af.y)) +    
  geom_point(size=0.1)+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=1) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs(title=paste0("common eSNPs"),
       x="MAF in AA", y = "MAF in EA")+
  theme_bw(base_size=22)+   
  xlim(0,0.5)+
  ylim(0,0.5) 
p8
dev.off()

```

Fourth, we compare the effect sizes between AA and EA for common eSNPs, population specific eSNPs, as well as non-eSNPs. We did observe that EA unique eSNPs tends to have large effect size in EA as compared to the effect size in AA; and vice versa.
```R
> summary(common_table_use_all$beta.x)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
-1.946808 -0.060525 -0.001166 -0.001511  0.058112  1.716955 
> summary(common_table_use_all$beta.y)
      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
-1.8799640 -0.0677804 -0.0009995 -0.0011039  0.0655924  1.9299760 


mydata = common_table_use_all[common_table_use_all$SNP_type == "non-eSNPs",] 
p9 = ggplot(mydata, aes(x=beta.x, y=beta.y)) +    
  geom_point(shape=19, fill="#999999", color="#999999", size=0.5)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  labs(title=paste0("non-eSNPs"),
       x="Effect size in AA", y = "Effect size in EA")+
  xlim(-2,2)+
  ylim(-2,2)  


mydata = common_table_use_all[common_table_use_all$SNP_type == "AA specific",]
p10 = ggplot(mydata, aes(x=beta.x, y=beta.y)) +    
  geom_point(shape=19, fill="tan1", color="tan1", size=0.5)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  labs(title=paste0("AA specific"),
       x="Effect size in AA", y = "Effect size in EA")+
  xlim(-2,2)+
  ylim(-2,2)


mydata = common_table_use_all[common_table_use_all$SNP_type == "EA specific",]
p11 = ggplot(mydata, aes(x=beta.x, y=beta.y)) +    
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=0.5)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  labs(title=paste0("EA specific"),
       x="Effect size in AA", y = "Effect size in EA")+
  xlim(-2,2)+
  ylim(-2,2)  

  
mydata = common_table_use_all[common_table_use_all$SNP_type == "common eSNPs",]
p12 = ggplot(mydata, aes(x=beta.x, y=beta.y)) +    
  geom_point(shape=19, fill="violet", color="violet", size=0.5)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  labs(title=paste0("common eSNPs"),
       x="Effect size in AA", y = "Effect size in EA")+
  xlim(-2,2)+
  ylim(-2,2) 

tiff(paste0(path_R2Q2,"/effectsize_in_4groups_SNPs.tiff"), units="in", width=10, height=10, res=150)    
(p9 | p10 ) /
      (p11 | p12)
dev.off()



```



Fifth, we also calculated the number of “primary eQTLs (based on P value) are in fact also having the largest effect sizes across all SNPs close to a given gene”. 

```R
load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order_correct.RData")
load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order_correct.RData")


gene_AA_eGene = gene_AA_anno_order[which(gene_AA_anno_order$eGene==1),]
gene_EA_eGene = gene_EA_anno_order[which(gene_EA_anno_order$eGene==1),]

sum(gene_AA_eGene$topSNP_rs == gene_AA_eGene$maxbeta_rs)
sum(gene_EA_eGene$topSNP_rs == gene_EA_eGene$maxbeta_rs)

sum(gene_AA_eGene$topSNP_rs == gene_AA_eGene$maxbeta_rs)/dim(gene_AA_eGene)[1]
sum(gene_EA_eGene$topSNP_rs == gene_EA_eGene$maxbeta_rs)/dim(gene_EA_eGene)[1]

> sum(gene_AA_eGene$topSNP_rs == gene_AA_eGene$maxbeta_rs)
[1] 428
> sum(gene_EA_eGene$topSNP_rs == gene_EA_eGene$maxbeta_rs)
[1] 488
> 
> sum(gene_AA_eGene$topSNP_rs == gene_AA_eGene$maxbeta_rs)/dim(gene_AA_eGene)[1]
[1] 0.07817352
> sum(gene_EA_eGene$topSNP_rs == gene_EA_eGene$maxbeta_rs)/dim(gene_EA_eGene)[1]
[1] 0.1108587

```
    
In AA, xxx (xxx%) out of the xxx primary eQTLs have the largest effect sizes across all SNPs mapped to a given gene (by chance alone you would expect xxx%). In EA, xxx (xxx%) out of xxx primary eQTLs also have the largest effect sizes across all SNPs mapped to a given gene (by chance alone you would expect xxx%). 

```R

gene_EA_eGene$one_over_numSNP = 1/gene_EA_eGene$SNPnum
mean(gene_EA_eGene$one_over_numSNP)
summary(gene_EA_eGene$one_over_numSNP)


gene_AA_eGene$one_over_numSNP = 1/gene_AA_eGene$SNPnum
mean(gene_AA_eGene$one_over_numSNP)
summary(gene_AA_eGene$one_over_numSNP)

> mean(gene_AA_eGene$one_over_numSNP)
[1] 0.001511355
> summary(gene_AA_eGene$one_over_numSNP)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
5.048e-05 9.533e-04 1.250e-03 1.511e-03 1.681e-03 4.167e-02 

```

In addition, we calculated the number of primary eQTLs (based on P value) having the largest MAF across all SNPs close to a given gene. We found that… Therefore, all these new evidence suggest that …..

```R
library(dplyr)

AA_table_eQTL = AA_table_subset %>%  group_by(GENE) %>%
summarise(topSNP = which.min(p_wald), 
          largest_maf = which.max(af))

sum(AA_table_eQTL$topSNP == AA_table_eQTL$largest_maf)
> sum(AA_table_eQTL$topSNP == AA_table_eQTL$largest_maf)
[1] 1105

sum(AA_table_eQTL$topSNP == AA_table_eQTL$largest_maf)/dim(gene_AA_eGene)[1]
[1] 0.2018265


EA_table_eQTL = EA_table_subset %>%  group_by(GENE) %>%
summarise(topSNP = which.min(p_wald), 
          largest_maf = which.max(af))

sum(EA_table_eQTL$topSNP == EA_table_eQTL$largest_maf)
> sum(EA_table_eQTL$topSNP == EA_table_eQTL$largest_maf)
[1] 715
715/dim(gene_EA_eGene)[1]
[1] 0.1624262



```
From Jonathan prichard Genetics paper, 2019:
Variable eQTL-detection power by allele frequency, however, did impact our calls. One might expect eQTL-based analyses to be limited to relatively common regulatory variants. Indeed, our eQTLs were, on average, more common than candidate SNPs. 

I added another figure similar to Figure 4B in Jonathan paper.
```

common_table_use_all$AA_maf_class = "<0.1"
common_table_use_all$AA_maf_class[which(common_table_use_all$af.x <=0.1)] = "<0.1"
common_table_use_all$AA_maf_class[which(common_table_use_all$af.x >0.1 & common_table_use_all$af.x <=0.15)] = "<0.15"
common_table_use_all$AA_maf_class[which(common_table_use_all$af.x >0.15 & common_table_use_all$af.x <=0.2)] = "<0.2"
common_table_use_all$AA_maf_class[which(common_table_use_all$af.x >0.2 & common_table_use_all$af.x <=0.25)] = "<0.25"
common_table_use_all$AA_maf_class[which(common_table_use_all$af.x >0.25 & common_table_use_all$af.x <=0.3)] = "<0.3"
common_table_use_all$AA_maf_class[which(common_table_use_all$af.x >0.3 & common_table_use_all$af.x <=0.35)] = "<0.35"
common_table_use_all$AA_maf_class[which(common_table_use_all$af.x >0.35 & common_table_use_all$af.x <=0.4)] = "<0.4"
common_table_use_all$AA_maf_class[which(common_table_use_all$af.x >0.4 & common_table_use_all$af.x <=0.45)] = "<0.45"
common_table_use_all$AA_maf_class[which(common_table_use_all$af.x >0.45 & common_table_use_all$af.x <=0.5)] = "<0.5"

> sum(common_table_use_all$signif.x)
[1] 180088

AA_bar = common_table_use_all %>%  group_by(AA_maf_class) %>%
summarise(proportion = sum(signif.x)/180088)

> table(common_table_use_all$AA_maf_class)/dim(common_table_use_all)[1]

      <0.1      <0.15       <0.2      <0.25       <0.3      <0.35       <0.4 
0.37572688 0.13693711 0.11332629 0.09211978 0.07615250 0.06431577 0.05480388 
     <0.45       <0.5 
0.04623072 0.04038708 

> AA_bar
# A tibble: 9 x 2
  AA_maf_class proportion
  <chr>             <dbl>
1 <0.1             0.168 
2 <0.15            0.122 
3 <0.2             0.132 
4 <0.25            0.116 
5 <0.3             0.113 
6 <0.35            0.107 
7 <0.4             0.0944
8 <0.45            0.0756
9 <0.5             0.0712






common_table_use_all$EA_maf_class = "<0.1"
common_table_use_all$EA_maf_class[which(common_table_use_all$af.y <=0.1)] = "<0.1"
common_table_use_all$EA_maf_class[which(common_table_use_all$af.y >0.1 & common_table_use_all$af.y <=0.15)] = "<0.15"
common_table_use_all$EA_maf_class[which(common_table_use_all$af.y >0.15 & common_table_use_all$af.y <=0.2)] = "<0.2"
common_table_use_all$EA_maf_class[which(common_table_use_all$af.y >0.2 & common_table_use_all$af.y <=0.25)] = "<0.25"
common_table_use_all$EA_maf_class[which(common_table_use_all$af.y >0.25 & common_table_use_all$af.y <=0.3)] = "<0.3"
common_table_use_all$EA_maf_class[which(common_table_use_all$af.y >0.3 & common_table_use_all$af.y <=0.35)] = "<0.35"
common_table_use_all$EA_maf_class[which(common_table_use_all$af.y >0.35 & common_table_use_all$af.y <=0.4)] = "<0.4"
common_table_use_all$EA_maf_class[which(common_table_use_all$af.y >0.4 & common_table_use_all$af.y <=0.45)] = "<0.45"
common_table_use_all$EA_maf_class[which(common_table_use_all$af.y >0.45 & common_table_use_all$af.y <=0.5)] = "<0.5"

> sum(common_table_use_all$signif.y)
[1] 246318

EA_bar = common_table_use_all %>%  group_by(EA_maf_class) %>%
summarise(proportion = sum(signif.y)/246318)

> table(common_table_use_all$EA_maf_class)/dim(common_table_use_all)[1]

      <0.1      <0.15       <0.2      <0.25       <0.3      <0.35       <0.4 
0.34047248 0.14504004 0.11468071 0.09579016 0.07980567 0.06957491 0.05961309 
     <0.45       <0.5 
0.05183228 0.04319066 

> EA_bar
# A tibble: 9 x 2
  EA_maf_class proportion
  <chr>             <dbl>
1 <0.1             0.169 
2 <0.15            0.132 
3 <0.2             0.127 
4 <0.25            0.119 
5 <0.3             0.110 
6 <0.35            0.103 
7 <0.4             0.0893
8 <0.45            0.0824
9 <0.5             0.0669

all_snps_propor = c(table(common_table_use_all$AA_maf_class)/dim(common_table_use_all)[1])
esnps_propor = AA_bar$proportion
propor = c(all_snps_propor, esnps_propor)
class = c(rep("All SNPs",length(all_snps_propor)),rep( "eSNPs",length(esnps_propor)))
maf = rep(c( "<0.1", "<0.15" , "<0.2" , "<0.25" , "<0.3", "<0.35"  ,"<0.4", "<0.45",  "<0.5"),2)
MAF= factor(maf, levels = c( "<0.1", "<0.15" , "<0.2" , "<0.25" , "<0.3", "<0.35"  ,"<0.4", "<0.45",  "<0.5"),order=T)
dat = data.frame(propor, class,maf)


p1 = ggplot(data=dat, aes(x=maf, y=propor, fill=class)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  labs(title="African American",x="Allele frequency", y = "Proportion")+
  theme_minimal(base_size = 22)+
  theme(legend.position="bottom") 

all_snps_propor = c(table(common_table_use_all$EA_maf_class)/dim(common_table_use_all)[1])
esnps_propor = EA_bar$proportion
propor = c(all_snps_propor, esnps_propor)
class = c(rep("All SNPs",length(all_snps_propor)),rep( "eSNPs",length(esnps_propor)))
maf = rep(c( "<0.1", "<0.15" , "<0.2" , "<0.25" , "<0.3", "<0.35"  ,"<0.4", "<0.45",  "<0.5"),2)
MAF= factor(maf, levels = c( "<0.1", "<0.15" , "<0.2" , "<0.25" , "<0.3", "<0.35"  ,"<0.4", "<0.45",  "<0.5"),order=T)
dat = data.frame(propor, class,maf)

p2 = ggplot(data=dat, aes(x=maf, y=propor, fill=class)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  labs(title="European American",x="Allele frequency", y = "Proportion")+
  theme_minimal(base_size = 22)+
  theme(legend.position="bottom") 



pdf(paste0(path_R2Q2,"/Bar_proportion_MAF.pdf"),width=14, height=5)
(p1 | p2 )
dev.off()

```

Another analysis that is affected by this issue is the assignment of a primary eQTL. 
It would be helpful to understand the magnitude of this effect by, for example, 
testing how many primary eQTLs (based on P value) are in fact also having the 
largest effect sizes across all SNPs close to a given gene.

```
# CHECK ea and aa, how many primary eQTL have largest abs effect size (two groups )
# compare between EA and AA, compare 
#    - largest abs effect size overlap 
#    - non-largest abs effect size overlap 

# AA
path_AA="/net/mulan/home/shanglu/GENOA/data/AA"
path_AA_result="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping/result"
path_AA_output="/net/mulan/home/shanglu/GENOA/data/AA/output/"
load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order.RData")
SNPnum = c()
topSNP_ps = c()
topSNP_allele1 = c()
topSNP_allele0 = c()
topSNP_af = c()
topSNP_rs = c()
topSNP_beta = c()
topSNP_pval = c()
maxbeta = c()
maxbeta_snp = c()
maxbeta_pval = c()
count = 0
for(i in 1:22){
  print(i)
  load(paste0(path_AA_result,"/AA_eqtl_res_chr_",i,"pheno2.RData"))
  for(j in 1:length(res)){
    count = count + 1
    if(sum(dim(res[[j]])[1])!=0){
    top_index = which.min(res[[j]]$p_wald)
    maxbeta_index = which.max(abs(res[[j]]$beta))
    maxbeta[count] = res[[j]]$beta[maxbeta_index]
    maxbeta_snp[count] = as.character(res[[j]]$rs[maxbeta_index])
    maxbeta_pval[count] = res[[j]]$p_wald[maxbeta_index]
    SNPnum[count] = dim(res[[j]])[1]
    topSNP_ps[count] = res[[j]]$ps[top_index]
    topSNP_allele1[count] = as.character(res[[j]]$allele1[top_index])
    topSNP_allele0[count] = as.character(res[[j]]$allele0[top_index])
    topSNP_af[count] = res[[j]]$af[top_index]
    topSNP_rs[count] = as.character(res[[j]]$rs[top_index])
    topSNP_pval[count] = res[[j]]$p_wald[top_index]
    topSNP_beta[count] = res[[j]]$beta[top_index]
    }
   }
}

gene_AA_anno_order$maxbeta_snp = maxbeta_snp
gene_AA_anno_order$maxbeta_pval = maxbeta_pval
gene_AA_anno_order$maxbeta = maxbeta
gene_AA_anno_order$topSNP_beta = topSNP_beta


# EA
path_EA="/net/mulan/home/shanglu/GENOA/data/EA"
path_EA_result="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping/result"
path_EA_output="/net/mulan/home/shanglu/GENOA/data/EA/output/"
load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order.RData")
SNPnum = c()
topSNP_ps = c()
topSNP_allele1 = c()
topSNP_allele0 = c()
topSNP_af = c()
topSNP_rs = c()
topSNP_beta = c()
topSNP_pval = c()
maxbeta = c()
maxbeta_snp = c()
maxbeta_pval = c()
count = 0
for(i in 1:22){
  print(i)
  load(paste0(path_EA_result,"/EA_eqtl_res_chr_",i,"pheno2.RData"))
  for(j in 1:length(res)){
    count = count + 1
    if(sum(dim(res[[j]])[1])!=0){
    top_index = which.min(res[[j]]$p_wald)
    maxbeta_index = which.max(abs(res[[j]]$beta))
    maxbeta[count] = res[[j]]$beta[maxbeta_index]
    maxbeta_snp[count] = as.character(res[[j]]$rs[maxbeta_index])
    maxbeta_pval[count] = res[[j]]$p_wald[maxbeta_index]
    SNPnum[count] = dim(res[[j]])[1]
    topSNP_ps[count] = res[[j]]$ps[top_index]
    topSNP_allele1[count] = as.character(res[[j]]$allele1[top_index])
    topSNP_allele0[count] = as.character(res[[j]]$allele0[top_index])
    topSNP_af[count] = res[[j]]$af[top_index]
    topSNP_rs[count] = as.character(res[[j]]$rs[top_index])
    topSNP_pval[count] = res[[j]]$p_wald[top_index]
    topSNP_beta[count] = res[[j]]$beta[top_index]
    }
   }
}

gene_EA_anno_order$maxbeta_snp = maxbeta_snp
gene_EA_anno_order$maxbeta_pval = maxbeta_pval
gene_EA_anno_order$maxbeta = maxbeta
gene_EA_anno_order$topSNP_beta = topSNP_beta


sum(as.character(gene_anno_AA_eGene$maxbeta_rs)==as.character(gene_anno_AA_eGene$topSNP_rs))
dim(gene_anno_AA_eGene)[1]
sum(as.character(gene_anno_AA_eGene$maxbeta_rs)==as.character(gene_anno_AA_eGene$topSNP_rs))/dim(gene_anno_AA_eGene)[1]
mean(1/gene_anno_AA_eGene$SNPnum)


> sum(as.character(gene_anno_AA_eGene$maxbeta_rs)==as.character(gene_anno_AA_eGene$topSNP_rs))
[1] 428
> dim(gene_anno_AA_eGene)[1]
[1] 5475
> sum(as.character(gene_anno_AA_eGene$maxbeta_rs)==as.character(gene_anno_AA_eGene$topSNP_rs))/dim(gene_anno_AA_eGene)[1]
[1] 0.07817352
> mean(1/gene_anno_AA_eGene$SNPnum)
[1] 0.001511355


sum(as.character(gene_anno_EA_eGene$maxbeta_rs)==as.character(gene_anno_EA_eGene$topSNP_rs))
dim(gene_anno_EA_eGene)[1]
sum(as.character(gene_anno_EA_eGene$maxbeta_rs)==as.character(gene_anno_EA_eGene$topSNP_rs))/dim(gene_anno_EA_eGene)[1]
mean(1/gene_anno_EA_eGene$SNPnum)

> sum(as.character(gene_anno_EA_eGene$maxbeta_rs)==as.character(gene_anno_EA_eGene$topSNP_rs))
[1] 488
> dim(gene_anno_EA_eGene)[1]
[1] 4402
> sum(as.character(gene_anno_EA_eGene$maxbeta_rs)==as.character(gene_anno_EA_eGene$topSNP_rs))/dim(gene_anno_EA_eGene)[1]
[1] 0.1108587
> mean(1/gene_anno_EA_eGene$SNPnum)
[1] 0.002679471






```









