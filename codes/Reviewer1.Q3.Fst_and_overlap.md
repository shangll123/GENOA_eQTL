
In Table 1, you show that only a small percentage of eSNPs overlap between the two populations (28-30%), 
but a larger percentage of eGenes overlap (53-63%). 
Do the eSNPs that overlap have a lower Fst than those that don't overlap? 

This would indicate the differences in eSNPs are largely due to allele frequency differences between populations. 
You've done this analysis at the eGene level, but not the eSNP level. 

Also, I'm assuming the overlap increases if you relax your significance threshold, 
but it would be helpful to show this and point out that some of the lack of overlap is a power issue. 



## 1. overlap at different FDR
## 2. Fst compare eSNPs vs noneSNPs, eGenes vs non-eGenes, at different FDR
## use FDR: 0.01, 0.05, 0.1, 0.2


```R
cd /net/mulan/home/shanglu/GENOA/analysis/compare

# make AA full table

pathAA="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
load(paste0(pathAA,"/AA_full_table_chr_1.RData"))
AA_table = AA_full_table
for(i in 2:22){
print(i)
load(paste0(pathAA,"/AA_full_table_chr_",i,".RData"))
AA_table = rbind(AA_table, AA_full_table)
}
save(AA_table, file = paste0(pathAA,"/AA_table.RData"))

# make EA full table

pathEA="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"
load(paste0(pathEA,"/EA_full_table_chr_1.RData"))
EA_table = EA_full_table
for(i in 2:22){
print(i)
load(paste0(pathEA,"/EA_full_table_chr_",i,".RData"))
EA_table = rbind(EA_table, EA_full_table)
}
save(EA_table, file = paste0(pathEA,"/EA_table.RData"))

#-----------------------
> dim(AA_table)
[1] 14511338       12
> 
> dim(EA_table)
[1] 8521801      12


#-----------------------
```

```R
> summary(AA_table$af)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0100  0.0240  0.0575  0.1169  0.1690  0.5000 
> summary(EA_table$af)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0100  0.0600  0.1540  0.1899  0.3060  0.5000  
 ```

# FDR 0.01
```R

AA_thr =  9.479415e-06
EA_thr = 2.069628e-05

AA_table_subset = AA_table[AA_table$p_wald <= AA_thr,]
EA_table_subset = EA_table[EA_table$p_wald <= EA_thr,]

#-------- eGenes
# AA
eGene_AA_num = length(unique(AA_table_subset$GENE))
gene_AA_num = length(unique(AA_table$GENE))
gene_percent_AA = eGene_AA_num/gene_AA_num
# EA
eGene_EA_num = length(unique(EA_table_subset$GENE))
gene_EA_num = length(unique(EA_table$GENE))
gene_percent_EA = eGene_EA_num/gene_EA_num
# overlap
eGene_AAEA_num = length(intersect(unique(AA_table_subset$GENE), unique(EA_table_subset$GENE)))
gene_percent_AA_eGene = eGene_AAEA_num/eGene_AA_num
gene_percent_EA_eGene = eGene_AAEA_num/eGene_EA_num

#-------- eSNPs
# AA
eSNP_AA_num = dim(AA_table_subset)[1]
SNP_AA_num = dim(AA_table)[1]
snp_percent_AA = eSNP_AA_num/SNP_AA_num
# EA
eSNP_EA_num = dim(EA_table_subset)[1]
SNP_EA_num = dim(EA_table)[1]
snp_percent_EA = eSNP_EA_num/SNP_EA_num
# overlap
common_table = merge(AA_table_subset, EA_table_subset, by = c("GENE","rs"))
eSNP_overlap_num = dim(common_table)[1]
snp_percent_AA_eSNP = eSNP_overlap_num/eSNP_AA_num
snp_percent_EA_eSNP = eSNP_overlap_num/eSNP_EA_num


eGene_AA_num
gene_AA_num
gene_percent_AA
eGene_EA_num
gene_EA_num
gene_percent_EA
eGene_AAEA_num
gene_percent_AA_eGene
gene_percent_EA_eGene
eSNP_AA_num
SNP_AA_num
snp_percent_AA
eSNP_EA_num
SNP_EA_num
snp_percent_EA
eSNP_overlap_num
snp_percent_AA_eSNP
snp_percent_EA_eSNP

/*
> eGene_AA_num
[1] 4433
> gene_AA_num
[1] 17572
> gene_percent_AA
[1] 0.2522763
> eGene_EA_num
[1] 3517
> gene_EA_num
[1] 17343
> gene_percent_EA
[1] 0.2027908
> eGene_AAEA_num
[1] 2350
> gene_percent_AA_eGene
[1] 0.530115
> gene_percent_EA_eGene
[1] 0.6681831
> eSNP_AA_num
[1] 278473
> SNP_AA_num
[1] 14511338
> snp_percent_AA
[1] 0.01919003
> eSNP_EA_num
[1] 294878
> SNP_EA_num
[1] 8521801
> snp_percent_EA
[1] 0.03460278
> eSNP_overlap_num
[1] 88971
> snp_percent_AA_eSNP
[1] 0.319496
> snp_percent_EA_eSNP
[1] 0.3017214
*/
```

# FDR 0.05
```R
AA_thr = 6.245907e-05
EA_thr = 0.0001385504
AA_table_subset = AA_table[AA_table$p_wald <= AA_thr,]
EA_table_subset = EA_table[EA_table$p_wald <= EA_thr,]

#-------- eGenes
# AA
eGene_AA_num = length(unique(AA_table_subset$GENE))
gene_AA_num = length(unique(AA_table$GENE))
gene_percent_AA = eGene_AA_num/gene_AA_num
# EA
eGene_EA_num = length(unique(EA_table_subset$GENE))
gene_EA_num = length(unique(EA_table$GENE))
gene_percent_EA = eGene_EA_num/gene_EA_num
# overlap
eGene_AAEA_num = length(intersect(unique(AA_table_subset$GENE), unique(EA_table_subset$GENE)))
gene_percent_AA_eGene = eGene_AAEA_num/eGene_AA_num
gene_percent_EA_eGene = eGene_AAEA_num/eGene_EA_num

#-------- eSNPs
# AA
eSNP_AA_num = dim(AA_table_subset)[1]
SNP_AA_num = dim(AA_table)[1]
snp_percent_AA = eSNP_AA_num/SNP_AA_num
# EA
eSNP_EA_num = dim(EA_table_subset)[1]
SNP_EA_num = dim(EA_table)[1]
snp_percent_EA = eSNP_EA_num/SNP_EA_num
# overlap
common_table = merge(AA_table_subset, EA_table_subset, by = c("GENE","rs"))
eSNP_overlap_num = dim(common_table)[1]
snp_percent_AA_eSNP = eSNP_overlap_num/eSNP_AA_num
snp_percent_EA_eSNP = eSNP_overlap_num/eSNP_EA_num


eGene_AA_num
gene_AA_num
gene_percent_AA
eGene_EA_num
gene_EA_num
gene_percent_EA
eGene_AAEA_num
gene_percent_AA_eGene
gene_percent_EA_eGene
eSNP_AA_num
SNP_AA_num
snp_percent_AA
eSNP_EA_num
SNP_EA_num
snp_percent_EA
eSNP_overlap_num
snp_percent_AA_eSNP
snp_percent_EA_eSNP

/*
> eGene_AA_num
[1] 5475
> gene_AA_num
[1] 17572
> gene_percent_AA
[1] 0.3115752
> eGene_EA_num
[1] 4402
> gene_EA_num
[1] 17343
> gene_percent_EA
[1] 0.25382
> eGene_AAEA_num
[1] 3048
> gene_percent_AA_eGene
[1] 0.5567123
> gene_percent_EA_eGene
[1] 0.6924125
> eSNP_AA_num
[1] 354931
> SNP_AA_num
[1] 14511338
> snp_percent_AA
[1] 0.02445887
> eSNP_EA_num
[1] 371309
> SNP_EA_num
[1] 8521801
> snp_percent_EA
[1] 0.04357166
> eSNP_overlap_num
[1] 112316
> snp_percent_AA_eSNP
[1] 0.3164446
> snp_percent_EA_eSNP
[1] 0.3024866
*/
```

# FDR 0.1
```R
AA_thr = 0.0001526481
EA_thr = 0.0003483179

AA_table_subset = AA_table[AA_table$p_wald <= AA_thr,]
EA_table_subset = EA_table[EA_table$p_wald <= EA_thr,]

#-------- eGenes
# AA
eGene_AA_num = length(unique(AA_table_subset$GENE))
gene_AA_num = length(unique(AA_table$GENE))
gene_percent_AA = eGene_AA_num/gene_AA_num
# EA
eGene_EA_num = length(unique(EA_table_subset$GENE))
gene_EA_num = length(unique(EA_table$GENE))
gene_percent_EA = eGene_EA_num/gene_EA_num
# overlap
eGene_AAEA_num = length(intersect(unique(AA_table_subset$GENE), unique(EA_table_subset$GENE)))
gene_percent_AA_eGene = eGene_AAEA_num/eGene_AA_num
gene_percent_EA_eGene = eGene_AAEA_num/eGene_EA_num

#-------- eSNPs
# AA
eSNP_AA_num = dim(AA_table_subset)[1]
SNP_AA_num = dim(AA_table)[1]
snp_percent_AA = eSNP_AA_num/SNP_AA_num
# EA
eSNP_EA_num = dim(EA_table_subset)[1]
SNP_EA_num = dim(EA_table)[1]
snp_percent_EA = eSNP_EA_num/SNP_EA_num
# overlap
common_table = merge(AA_table_subset, EA_table_subset, by = c("GENE","rs"))
eSNP_overlap_num = dim(common_table)[1]
snp_percent_AA_eSNP = eSNP_overlap_num/eSNP_AA_num
snp_percent_EA_eSNP = eSNP_overlap_num/eSNP_EA_num


eGene_AA_num
gene_AA_num
gene_percent_AA
eGene_EA_num
gene_EA_num
gene_percent_EA
eGene_AAEA_num
gene_percent_AA_eGene
gene_percent_EA_eGene
eSNP_AA_num
SNP_AA_num
snp_percent_AA
eSNP_EA_num
SNP_EA_num
snp_percent_EA
eSNP_overlap_num
snp_percent_AA_eSNP
snp_percent_EA_eSNP

/*
> 
> eGene_AA_num
[1] 6172
> gene_AA_num
[1] 17572
> gene_percent_AA
[1] 0.3512406
> eGene_EA_num
[1] 5121
> gene_EA_num
[1] 17343
> gene_percent_EA
[1] 0.2952776
> eGene_AAEA_num
[1] 3540
> gene_percent_AA_eGene
[1] 0.573558
> gene_percent_EA_eGene
[1] 0.6912712
> eSNP_AA_num
[1] 404971
> SNP_AA_num
[1] 14511338
> snp_percent_AA
[1] 0.02790721
> eSNP_EA_num
[1] 421583
> SNP_EA_num
[1] 8521801
> snp_percent_EA
[1] 0.04947112
> eSNP_overlap_num
[1] 126564
> snp_percent_AA_eSNP
[1] 0.3125261
> snp_percent_EA_eSNP
[1] 0.3002113
*/

```



# FDR 0.2
```R
AA_thr = 0.0003909549
EA_thr =0.0009466665

AA_table_subset = AA_table[AA_table$p_wald <= AA_thr,]
EA_table_subset = EA_table[EA_table$p_wald <= EA_thr,]

#-------- eGenes
# AA
eGene_AA_num = length(unique(AA_table_subset$GENE))
gene_AA_num = length(unique(AA_table$GENE))
gene_percent_AA = eGene_AA_num/gene_AA_num
# EA
eGene_EA_num = length(unique(EA_table_subset$GENE))
gene_EA_num = length(unique(EA_table$GENE))
gene_percent_EA = eGene_EA_num/gene_EA_num
# overlap
eGene_AAEA_num = length(intersect(unique(AA_table_subset$GENE), unique(EA_table_subset$GENE)))
gene_percent_AA_eGene = eGene_AAEA_num/eGene_AA_num
gene_percent_EA_eGene = eGene_AAEA_num/eGene_EA_num

#-------- eSNPs
# AA
eSNP_AA_num = dim(AA_table_subset)[1]
SNP_AA_num = dim(AA_table)[1]
snp_percent_AA = eSNP_AA_num/SNP_AA_num
# EA
eSNP_EA_num = dim(EA_table_subset)[1]
SNP_EA_num = dim(EA_table)[1]
snp_percent_EA = eSNP_EA_num/SNP_EA_num
# overlap
common_table = merge(AA_table_subset, EA_table_subset, by = c("GENE","rs"))
eSNP_overlap_num = dim(common_table)[1]
snp_percent_AA_eSNP = eSNP_overlap_num/eSNP_AA_num
snp_percent_EA_eSNP = eSNP_overlap_num/eSNP_EA_num


eGene_AA_num
gene_AA_num
gene_percent_AA
eGene_EA_num
gene_EA_num
gene_percent_EA
eGene_AAEA_num
gene_percent_AA_eGene
gene_percent_EA_eGene
eSNP_AA_num
SNP_AA_num
snp_percent_AA
eSNP_EA_num
SNP_EA_num
snp_percent_EA
eSNP_overlap_num
snp_percent_AA_eSNP
snp_percent_EA_eSNP

/*

> eGene_AA_num
[1] 7256
> gene_AA_num
[1] 17572
> gene_percent_AA
[1] 0.4129297
> eGene_EA_num
[1] 6216
> gene_EA_num
[1] 17343
> gene_percent_EA
[1] 0.3584155
> eGene_AAEA_num
[1] 4292
> gene_percent_AA_eGene
[1] 0.5915105
> gene_percent_EA_eGene
[1] 0.6904762
> eSNP_AA_num
[1] 473459
> SNP_AA_num
[1] 14511338
> snp_percent_AA
[1] 0.03262683
> eSNP_EA_num
[1] 494518
> SNP_EA_num
[1] 8521801
> snp_percent_EA
[1] 0.05802975
> eSNP_overlap_num
[1] 146684
> snp_percent_AA_eSNP
[1] 0.3098135
> snp_percent_EA_eSNP
[1] 0.2966201
*/



```


# use FDR 0.05

Between populations, the effect sizes of the identified eSNPs shared between EAs and AAs are highly correlated with 
each other (Pearson coefficient =xxx; p-value=xxx; Figure 1C), with xxx% of gene-SNP pairs sharing the 
same effect sign between AA and EA.

```R
# check effect size correlation
common_table$allele1.x = as.character(common_table$allele1.x)
common_table$allele1.y = as.character(common_table$allele1.y)
common_table$allele0.x = as.character(common_table$allele0.x)
common_table$allele0.y = as.character(common_table$allele0.y)
sum(common_table$allele1.x == common_table$allele1.y)
sum(common_table$allele0.x == common_table$allele0.y)
[1] 87363

ind = which(common_table$allele1.x != common_table$allele1.y)
common_table[ind,]$beta.x = -common_table[ind,]$beta.x
cor.test(common_table$beta.x, common_table$beta.y)
/*
> cor.test(common_table$beta.x, common_table$beta.y)$estimate
      cor 
0.9028881 
> cor.test(common_table$beta.x, common_table$beta.y)$p.value
[1] 0
*/

sum(common_table$beta.x*common_table$beta.y>0)
/*
> sum(common_table$beta.x*common_table$beta.y>0)
[1] 109334
> sum(common_table$beta.x*common_table$beta.y>0)
[1] 109334
> length(common_table$beta.x)
[1] 112316
> 109334/112316
[1] 0.9734499
*/



```

# Fst analysis
```R


# pathAA="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
# load(paste0(pathAA,"/AA_table.RData"))
# pathEA="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"
# load(paste0(pathEA,"/EA_table.RData"))
# load("/net/mulan/home/shanglu/GENOA/analysis/compare/Fst_result.RData")

AA_thr = 6.245907e-05
EA_thr = 0.0001385504

AA_table_subset = AA_table[AA_table$p_wald <= AA_thr,]
EA_table_subset = EA_table[EA_table$p_wald <= EA_thr,]

#> dim(AA_table_subset)
#[1] 354931     12
#> dim(EA_table_subset)
#[1] 371309     12


common_table_all = merge(AA_table, EA_table, by = c("GENE","rs"))
dim(common_table_all)
#> dim(common_table_all)
#[1] 6060046      22

# only use snps where AA and EA have same ref allele and alt allele:

common_table_all$allele1.x = as.character(common_table_all$allele1.x)
common_table_all$allele0.x = as.character(common_table_all$allele0.x)
common_table_all$allele1.y = as.character(common_table_all$allele1.y)
common_table_all$allele0.y = as.character(common_table_all$allele0.y)
common_table_all$SNP = common_table_all$rs

ind = which(common_table_all$allele1.x == common_table_all$allele1.y)

common_table_all_sameref = common_table_all[ind,]
#dim(common_table_all_sameref)
#[1] 5171951      23

common_table_all_Fst = merge(common_table_all_sameref,Fst_result, by="SNP")

dim(common_table_all_Fst)
[1] 5171951      29

common_table_all_Fst$signif.x = 0
common_table_all_Fst$signif.x[which(common_table_all_Fst$p_wald.x <= AA_thr)]=1
common_table_all_Fst$signif.y = 0
common_table_all_Fst$signif.y[which(common_table_all_Fst$p_wald.y <= EA_thr)]=1
common_table_all_Fst$signif_x_y = common_table_all_Fst$signif.x + common_table_all_Fst$signif.y
common_table_all_Fst$spec_AA = common_table_all_Fst$signif.x * common_table_all_Fst$signif_x_y
common_table_all_Fst$spec_AA = common_table_all_Fst$signif.x * common_table_all_Fst$signif_x_y
common_table_all_Fst$spec_AA[common_table_all_Fst$spec_AA==2] = 0
common_table_all_Fst$spec_EA = common_table_all_Fst$signif.y * common_table_all_Fst$signif_x_y
common_table_all_Fst$spec_EA[common_table_all_Fst$spec_EA==2] = 0

########################################################
# first view each eSNP in the context of gene-SNP pairs
########################################################

common_table_all_Fst$SNP_type = "non-eSNPs"
common_table_all_Fst$SNP_type[common_table_all_Fst$spec_AA==1] = "AA specific"
common_table_all_Fst$SNP_type[common_table_all_Fst$spec_EA==1] = "EA specific"
common_table_all_Fst$SNP_type[common_table_all_Fst$signif_x_y==2] = "common eSNPs"

#> table(common_table_all_Fst$SNP_type )
# AA specific common eSNPs  EA specific    non-eSNPs 
#       92725        87363       158955      4832908 
       
library(qvalue)
library(data.table)
library(dplyr)

common_table_all_Fst_result = common_table_all_Fst %>%  group_by(SNP_type) %>%
summarise(meanFst = mean(Fst), 
          medianFst = median(Fst),
          minFst = min(Fst),
          maxFst = max(Fst))

common_table_all_Fst_result
> common_table_all_Fst_result
# A tibble: 4 x 5
  SNP_type     meanFst medianFst minFst maxFst
  <chr>          <dbl>     <dbl>  <dbl>  <dbl>
1 AA specific   0.0622    0.0324      0  0.442
2 common eSNPs  0.0576    0.0311      0  0.413
3 EA specific   0.0608    0.0367      0  0.403
4 non-eSNPs     0.0563    0.0314      0  0.443




########################################################
# then view use unique SNPs to represent eSNPs 
########################################################

AA_specific = as.character(unique(common_table_all_Fst[which(common_table_all_Fst$SNP_type=="AA specific"),]$rs))
EA_specific = as.character(unique(common_table_all_Fst[which(common_table_all_Fst$SNP_type=="EA specific"),]$rs))
noeSNP = as.character(unique(common_table_all_Fst[which(common_table_all_Fst$SNP_type=="non-eSNPs"),]$rs))
common_signif_rs = as.character(unique(common_table_all_Fst[which(common_table_all_Fst$SNP_type=="common eSNPs"),]$rs))

Fst_result$SNP = as.character(Fst_result$SNP)
Fst_result_tmp = Fst_result[which(Fst_result$SNP %in% c(AA_specific,EA_specific,noeSNP,common_signif_rs)),]
#> dim(Fst_result_tmp)
#[1] 2367612       8

Fst_result_tmp$TYPE[which(Fst_result_tmp$SNP  %in% AA_specific)] = "AA specific"
Fst_result_tmp$TYPE[which(Fst_result_tmp$SNP  %in% AA_specific)] = "AA specific"
Fst_result_tmp$TYPE[which(Fst_result_tmp$SNP  %in% EA_specific)] = "EA specific"
Fst_result_tmp$TYPE[which(Fst_result_tmp$SNP  %in% common_signif_rs)] = "common eSNPs"

#table(Fst_result_tmp$TYPE)
#> table(Fst_result_tmp$TYPE)
# AA specific common eSNPs  EA specific    non eSNPs 
#       68008        78567       132088      2088949 
       
Fst_result_tmp_result = Fst_result_tmp %>%  group_by(TYPE) %>%
summarise(meanFst = mean(Fst), 
          medianFst = median(Fst),
          minFst = min(Fst),
          maxFst = max(Fst))

> Fst_result_tmp_result
# A tibble: 4 x 5
  TYPE         meanFst medianFst minFst maxFst
  <chr>          <dbl>     <dbl>  <dbl>  <dbl>
1 AA specific   0.0649    0.0347      0  0.442
2 common eSNPs  0.0579    0.0318      0  0.413
3 EA specific   0.0603    0.0366      0  0.403
4 non eSNPs     0.0560    0.0311      0  0.443


```


# if use all gene-SNP pairs rather than common gene-SNP pairs
```R

AA_table$signif = (AA_table$p_wald<=AA_thr)*1
EA_table$signif = (EA_table$p_wald<=EA_thr)*1
common_signif = merge(AA_table_subset,EA_table_subset, by=c("GENE","rs") )
common_signif_rs = common_signif$rs


> length(common_signif_rs)
[1] 112316


AA_specific = unique(AA_table_subset[-which(AA_table_subset$rs %in% common_signif_rs),]$rs)
EA_specific = unique(EA_table_subset[-which(EA_table_subset$rs %in% common_signif_rs),]$rs)
noeSNP = intersect(AA_table[AA_table$signif==0,]$rs,EA_table[EA_table$signif==0,]$rs)

Fst_result$TYPE = "non eSNPs"
Fst_result$TYPE[which(Fst_result$SNP  %in% AA_specific)] = "AA specific"
Fst_result$TYPE[which(Fst_result$SNP  %in% EA_specific)] = "EA specific"
Fst_result$TYPE[which(Fst_result$SNP  %in% common_signif_rs)] = "common eSNPs"

Fst_result_result = Fst_result %>%  group_by(TYPE) %>%
summarise(meanFst = mean(Fst), 
          medianFst = median(Fst),
          minFst = min(Fst),
          maxFst = max(Fst))
          
> table(Fst_result$TYPE)

 AA specific common eSNPs  EA specific    non eSNPs 
      197188       101285       216653     24391818 
      
> Fst_result_result
# A tibble: 4 x 5
  TYPE         meanFst medianFst minFst maxFst
  <chr>          <dbl>     <dbl>  <dbl>  <dbl>
1 AA specific   0.0895   0.0559       0  0.750
2 common eSNPs  0.0888   0.0472       0  0.654
3 EA specific   0.0790   0.0440       0  0.694
4 non eSNPs     0.0338   0.00902      0  0.807


```


```R
# figures

Fst = Fst_result$Fst
Class = factor(Fst_result$TYPE,
levels = c("non eSNPs","AA specific","EA specific","common eSNPs"),order=T)

dat = data.frame(Fst,Class)

setwd("/net/mulan/home/shanglu/GENOA/analysis/compare")
library(ggplot2)
 
tiff(paste0("Fst_boxplot.tiff"), units="in", width=6, height=5, res=150)
ggplot(dat, aes(x=Class, y=Fst)) +
    geom_boxplot(alpha=0.5,fill="white") +
    ylim(0,1)+
    stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme(legend.position="none") +theme_bw(base_size = 15)+
    scale_fill_brewer(palette="Set3")+labs(title="",x="", y = "Fst")
dev.off()


dat = data.frame(Fst,Class)
dat$logFst = log10(dat$Fst + 1)
tiff(paste0("Fst_boxplot_log10plus1.tiff"), units="in", width=6, height=5, res=150)
ggplot(dat, aes(x=Class, y=logFst)) +
    geom_boxplot(alpha=0.5,fill="white") +
   # ylim(0,1)+
    stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme(legend.position="none") +theme_bw(base_size = 15)+
    scale_fill_brewer(palette="Set3")+labs(title="",x="", y = "log10(Fst+1)")
dev.off()



```


```R
# other analysis
AA_spec_fst = Fst_result$Fst[Fst_result$TYPE == "AA specific"]
EA_spec_fst = Fst_result$Fst[Fst_result$TYPE == "EA specific"]
common_fst = Fst_result$Fst[Fst_result$TYPE == "common eSNPs"]
noneSNP_fst = Fst_result$Fst[Fst_result$TYPE == "non eSNPs"]


wilcox.test(common_fst , AA_spec_fst, alternative = "less")$p.value
wilcox.test(common_fst , EA_spec_fst, alternative = "greater")$p.value
wilcox.test(common_fst , noneSNP_fst, alternative = "greater")$p.value

> wilcox.test(common_fst , AA_spec_fst, alternative = "less")$p.value
[1] 1.250164e-272
> wilcox.test(common_fst , EA_spec_fst, alternative = "greater")$p.value
[1] 1.604164e-40
> wilcox.test(common_fst , noneSNP_fst, alternative = "greater")$p.value
[1] 0


library(clinfun)
# order: non-eSNPs,  EA specific, AA specific, common 
Fst_trendtest = list()
Fst_trendtest$Background = noneSNP_fst
Fst_trendtest$AA_unique = AA_spec_fst
Fst_trendtest$EA_unique = EA_spec_fst
Fst_trendtest$Common = common_fst
pieces<-Fst_trendtest
n<-c(length(Fst_trendtest$Background),length(Fst_trendtest$AA_unique),length(Fst_trendtest$EA_unique),length(Fst_trendtest$Common))
grp<-as.ordered(factor(rep(1:length(n),n))) 
res = jonckheere.test(unlist(pieces),grp,alternative="increasing")

#data:  
#JT = 9.0865e+12, p-value < 2.2e-16
#> res$p.value
#[1] 0

# order: non-eSNPs, EA specific, common, AA specific
Fst_trendtest = list()
Fst_trendtest$Background = noneSNP_fst
Fst_trendtest$EA_unique = EA_spec_fst
Fst_trendtest$Common = common_fst
Fst_trendtest$AA_unique = AA_spec_fst
pieces<-Fst_trendtest
n<-c(length(Fst_trendtest$Background),length(Fst_trendtest$EA_unique),length(Fst_trendtest$Common),length(Fst_trendtest$AA_unique))
grp<-as.ordered(factor(rep(1:length(n),n))) 
res = jonckheere.test(unlist(pieces),grp,alternative="increasing")

#data
#JT = 9.0929e+12, p-value < 2.2e-16
#> res$p.value
#[1] 0

```
