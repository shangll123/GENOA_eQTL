
```
# from MESA paper:

poplist2 <- c('AFA','HIS','CAU','MEX','YRI','GEU','FHS')   #j, testin, pop2
poplist <- c('AFA','HIS','CAU','AFHI','ALL')   #i, refin, pop1

if(poplist[i] == poplist2[j]){
          pi1 <- 1
          }else{
          pop1 <- refinfile 
          pop2 <- testinfile
          pop1fdr05 <- dplyr::filter(pop1, FDR < 0.05)
          fdr<-dim(pop1fdr05)
          pop2tested <- inner_join(pop1fdr05,pop2,by=c("snps","gene"))
          over<-dim(pop2tested)
          pop2pvals <- pop2tested$pvalue.y
          qobjpop2 <- qvalue(p = pop2pvals)
          pi1 <- 1 - qobjpop2$pi0
}

```

```

library(qvalue)
library(data.table)
library(dplyr)

###########################
# AA - EA
###########################

pathAA="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
pathEA="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"

load(paste0(pathAA,"/AA_table.RData"))
load(paste0(pathEA,"/EA_table.RData"))
AA_table$snps = as.character(AA_table$rs)
AA_table$gene = as.character(AA_table$GENE)
EA_table$snps = as.character(EA_table$rs)
EA_table$gene = as.character(EA_table$GENE)


AA_thr = 6.245907e-05
EA_thr = 0.0001385504


pop1 <- AA_table # refinfile
pop2 <- EA_table # testinfile
pop1fdr05 <- dplyr::filter(pop1, p_wald <= AA_thr) 
fdr<-dim(pop1fdr05)
pop2tested <- inner_join(pop1fdr05,pop2,by=c("snps","gene"))
over<-dim(pop2tested)
pop2pvals <- pop2tested$p_wald.y
qobjpop2 <- qvalue(p = pop2pvals)
pi1 <- 1 - qobjpop2$pi0
pi1
pop2pvals_2 = pop2pvals
qobjpop2=qvalue_truncp(pop2pvals_2)
pi1 <- 1 - qobjpop2$pi0
pi1

#---------------
# qvalue version
AA-EA
> pi1
[1] 0.8340495
# qvalue_truncp version
> pi1
[1] 0.8340495
#---------------

###########################
# EA - AA
###########################
pop1 <- EA_table # refinfile
pop2 <- AA_table # testinfile
pop1fdr05 <- dplyr::filter(pop1, p_wald <= EA_thr) 
fdr<-dim(pop1fdr05)
pop2tested <- inner_join(pop1fdr05,pop2,by=c("snps","gene"))
over<-dim(pop2tested)
pop2pvals <- pop2tested$p_wald.y
qobjpop2 <- qvalue(p = pop2pvals)
pi1 <- 1 - qobjpop2$pi0
pi1
pop2pvals_2 = pop2pvals
qobjpop2=qvalue_truncp(pop2pvals_2)
pi1 <- 1 - qobjpop2$pi0
pi1

#---------------
# qvalue version
EA-AA
> pi1
[1] 0.7525856
#---------------
# qvalue_truncp version
EA-AA
> pi1
[1] 0.7524728
#---------------

#################################################################
# try to download YRI and EUR results from Geuvadis website:

YRI = read.table("/net/mulan/home/shanglu/GENOA/analysis/compare/YRI89.gene.cis.FDR5.all.rs137.txt",header=T)
YRI$gene = as.character(YRI$GENE_ID)
gene_name = as.character(YRI$gene)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
YRI$gene=unlist(gene_name02)
YRI$rs = as.character(YRI$SNP_ID)

EUR = read.table("/net/mulan/home/shanglu/GENOA/analysis/compare/EUR373.gene.cis.FDR5.all.rs137.txt",header=T)
EUR$gene = as.character(EUR$GENE_ID)
gene_name = as.character(EUR$gene)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
EUR$gene=unlist(gene_name02)
EUR$rs = as.character(EUR$SNP_ID)



AA_table$rs = as.character(AA_table$rs)
AA_table$gene = as.character(AA_table$gene)

pop1 <- AA_table
pop2 <- YRI
pop1fdr05_yri <- dplyr::filter(pop1, p_wald <= AA_thr) 
pop2tested <- inner_join(pop1fdr05_yri,pop2,by=c("rs","gene"))
pop2pvals_1 <- pop2tested$pvalue  
qobjpop2=qvalue_truncp(pop2pvals_1)
pi1 <- 1 - qobjpop2$pi0
pi1

> pi1
[1] 0.8033737
> length(pop2pvals_1)
[1] 6590



EA_table$rs = as.character(EA_table$rs)
EA_table$gene = as.character(EA_table$gene)

pop1 <- EA_table
pop2 <- YRI
pop1fdr05_eur <- dplyr::filter(pop1, p_wald <= EA_thr) 
pop2tested <- inner_join(pop1fdr05_eur,pop2,by=c("rs","gene"))
pop2pvals_1 <- pop2tested$pvalue  
qobjpop2=qvalue_truncp(pop2pvals_1)
pi1 <- 1 - qobjpop2$pi0
pi1

> pi1
[1] 0.8130049

> length(pop2pvals_1)
[1] 4901

# AA EUR
pop1 <- AA_table
pop2 <- EUR
pop1fdr05_yri <- dplyr::filter(pop1, p_wald <= AA_thr) 
pop2tested <- inner_join(pop1fdr05_yri,pop2,by=c("rs","gene"))
pop2pvals_1 <- pop2tested$pvalue  
qobjpop2=qvalue_truncp(pop2pvals_1)
pi1 <- 1 - qobjpop2$pi0
pi1

[1] 0.9310375
> length(pop2pvals_1)
[1] 63630

# EA EUR
pop1 <- EA_table
pop2 <- EUR
pop1fdr05_yri <- dplyr::filter(pop1, p_wald <= EA_thr) 
pop2tested <- inner_join(pop1fdr05_yri,pop2,by=c("rs","gene"))
pop2pvals_1 <- pop2tested$pvalue  
qobjpop2=qvalue_truncp(pop2pvals_1)
pi1 <- 1 - qobjpop2$pi0
pi1

[1] 0.9248046
> length(pop2pvals_1)
[1] 101022

#################################################################
# AA vs YRI & EUR

pop2pvals_1_yri = c()
pop2pvals_2_eur = c()
for(chr in 1:22){
print(chr)
path_AA="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
load(paste0(path_AA,"/AA_full_table_chr_",chr,".RData"))
AA_full_table$snps = as.character(AA_full_table$rs)
AA_full_table$gene = as.character(AA_full_table$GENE)
path_Geuvadis = "/net/mulan/home/shanglu/GENOA/analysis/compare_study"
load(paste0(path_Geuvadis,"/geu_yri_chr_",chr,".RData"))
load(paste0(path_Geuvadis,"/geu_eur_chr_",chr,".RData"))

geu_yri_chr$gene = as.character(geu_yri_chr$GENE_ID)
gene_name = as.character(geu_yri_chr$GENE_ID)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
geu_yri_chr$gene=unlist(gene_name02)
geu_yri_chr$ps=geu_yri_chr$SNPpos

geu_eur_chr$gene = as.character(geu_eur_chr$GENE_ID)
gene_name = as.character(geu_eur_chr$GENE_ID)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
geu_eur_chr$gene=unlist(gene_name02)
geu_eur_chr$ps=geu_eur_chr$SNPpos

pop1 <- AA_full_table
pop2 <- geu_yri_chr
pop1fdr05_yri <- dplyr::filter(pop1, p_wald <= AA_thr) 
pop2tested_yri <- inner_join(pop1fdr05_yri,pop2,by=c("ps","gene"))
pop2pvals_1_yri <- c(pop2pvals_1_yri, pop2tested_yri$pvalue)

pop1 <- AA_full_table
pop2 <- geu_eur_chr
pop1fdr05_eur <- dplyr::filter(pop1, p_wald <= AA_thr) 
pop2tested_eur <- inner_join(pop1fdr05_eur,pop2,by=c("ps","gene"))
pop2pvals_2_eur <- c(pop2pvals_2_eur, pop2tested_eur$pvalue)

}


#qobjpop2 <- qvalue(p = pop2pvals_1)
# Error in pi0est(p, ...) : 
#  ERROR: maximum p-value is smaller than lambda range. Change the range of lambda or use qvalue_truncp() for truncated p values.
  
# AA - YRI
qobjpop2=qvalue_truncp(pop2pvals_1_yri)
pi1 <- 1 - qobjpop2$pi0
pi1
> pi1
[1] 0.6928576

# AA - EUR
qobjpop2=qvalue_truncp(pop2pvals_2_eur)
pi1 <- 1 - qobjpop2$pi0
pi1

> pi1
[1] 0.9190936

# the difference might due to the number of matched gene-SNP pairs in YRI and EUR
> length(pop2pvals_1_yri)
[1] 29051
> length(pop2pvals_2_eur)
[1] 94038




#################################################################
# EA vs YRI & EUR

pop2pvals_1 = c()
pop2pvals_2 = c()
for(chr in 1:22){
print(i)
path_EA="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"
load(paste0(path_EA,"/EA_full_table_chr_",chr,".RData"))
EA_full_table$snps = as.character(EA_full_table$rs)
EA_full_table$gene = as.character(EA_full_table$GENE)
path_Geuvadis = "/net/mulan/home/shanglu/GENOA/analysis/compare_study"
load(paste0(path_Geuvadis,"/geu_yri_chr_",chr,".RData"))
load(paste0(path_Geuvadis,"/geu_eur_chr_",chr,".RData"))

geu_yri_chr$gene = as.character(geu_yri_chr$GENE_ID)
gene_name = as.character(geu_yri_chr$GENE_ID)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
geu_yri_chr$gene=unlist(gene_name02)
geu_yri_chr$ps=geu_yri_chr$SNPpos

geu_eur_chr$gene = as.character(geu_eur_chr$GENE_ID)
gene_name = as.character(geu_eur_chr$GENE_ID)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
geu_eur_chr$gene=unlist(gene_name02)
geu_eur_chr$ps=geu_eur_chr$SNPpos

pop1 <- EA_full_table
pop2 <- geu_yri_chr
pop1fdr05 <- dplyr::filter(pop1, p_wald <= 0.0001367967) 
pop2tested <- inner_join(pop1fdr05,pop2,by=c("ps","gene"))
pop2pvals_1 <- c(pop2pvals_1, pop2tested$pvalue)

pop1 <- EA_full_table
pop2 <- geu_eur_chr
pop1fdr05 <- dplyr::filter(pop1, p_wald <= 0.0001367967) 
pop2tested <- inner_join(pop1fdr05,pop2,by=c("ps","gene"))
pop2pvals_2 <- c(pop2pvals_2, pop2tested$pvalue)


}

# EA - YRI
 #qobjpop2 <- qvalue(p = pop2pvals_1)
qobjpop2=qvalue_truncp(pop2pvals_1)
pi1 <- 1 - qobjpop2$pi0
pi1
> pi1
[1] 0.6989223

# EA - EUR
#qobjpop2 <- qvalue(p = pop2pvals_2)
qobjpop2=qvalue_truncp(pop2pvals_2)
pi1 <- 1 - qobjpop2$pi0
pi1
> pi1
[1] 0.9092186

#################################################################
# AA vs AFA, CAU, HIS
pop2pvals_1 = c()
pop2pvals_2 = c()
pop2pvals_3 = c()

for(chr in 1:22){
print(chr)
path_AA="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
load(paste0(path_AA,"/AA_full_table_chr_",chr,".RData"))
AA_full_table$snps = as.character(AA_full_table$rs)
AA_full_table$gene = as.character(AA_full_table$GENE)

path_AFA="/net/mulan/home/shanglu/GENOA/analysis/prediction/AFA"
load(paste0(path_AFA,"/AFA_to_use_chr_",chr,".RData"))
path_CAU="/net/mulan/home/shanglu/GENOA/analysis/prediction/CAU"
load(paste0(path_CAU,"/CAU_to_use_chr_",chr,".RData"))
path_HIS="/net/mulan/home/shanglu/GENOA/analysis/prediction/HIS"
load(paste0(path_HIS,"/HIS_to_use_chr_",chr,".RData"))

AFA_to_use$snps = as.character(AFA_to_use$snps)
AFA_to_use$gene = as.character(AFA_to_use$gene)
CAU_to_use$snps = as.character(CAU_to_use$snps)
CAU_to_use$gene = as.character(CAU_to_use$gene)
HIS_to_use$snps = as.character(HIS_to_use$snps)
HIS_to_use$gene = as.character(HIS_to_use$gene)

pop1 <- AA_full_table
pop1fdr05 <- dplyr::filter(pop1, p_wald <= AA_thr) 

pop2 <- AFA_to_use
pop2tested <- inner_join(pop1fdr05,pop2,by=c("snps","gene"))
pop2pvals_1 <- c(pop2pvals_1, pop2tested$pvalue)


pop2 <- CAU_to_use
pop2tested <- inner_join(pop1fdr05,pop2,by=c("snps","gene"))
pop2pvals_2 <- c(pop2pvals_2, pop2tested$pvalue)


pop2 <- HIS_to_use
pop2tested <- inner_join(pop1fdr05,pop2,by=c("snps","gene"))
pop2pvals_3 <- c(pop2pvals_3, pop2tested$pvalue)


}

# AA - AFA
# qobjpop2 <- qvalue(p = pop2pvals_1)
qobjpop2=qvalue_truncp(pop2pvals_1)
pi1 <- 1 - qobjpop2$pi0
pi1
> pi1
[1] 0.5930053


# AA - CAU
# qobjpop2 <- qvalue(p = pop2pvals_2)
qobjpop2=qvalue_truncp(pop2pvals_2)
pi1 <- 1 - qobjpop2$pi0
pi1
> pi1
[1] 0.6896177

# AA - HIS
# qobjpop2 <- qvalue(p = pop2pvals_3)
qobjpop2=qvalue_truncp(pop2pvals_3)
pi1 <- 1 - qobjpop2$pi0
pi1

> pi1
[1] 0.632631

length(pop2pvals_1)
length(pop2pvals_2)
length(pop2pvals_3)

> length(pop2pvals_1)
[1] 164471
> length(pop2pvals_2)
[1] 128755
> length(pop2pvals_3)
[1] 160899

#################################################################
# EA vs AFA, CAU, HIS
pop2pvals_1 = c()
pop2pvals_2 = c()
pop2pvals_3 = c()

for(chr in 1:22){
print(chr)
path_EA="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"
load(paste0(path_EA,"/EA_full_table_chr_",chr,".RData"))
EA_full_table$snps = as.character(EA_full_table$rs)
EA_full_table$gene = as.character(EA_full_table$GENE)

path_AFA="/net/mulan/home/shanglu/GENOA/analysis/prediction/AFA"
load(paste0(path_AFA,"/AFA_to_use_chr_",chr,".RData"))
path_CAU="/net/mulan/home/shanglu/GENOA/analysis/prediction/CAU"
load(paste0(path_CAU,"/CAU_to_use_chr_",chr,".RData"))
path_HIS="/net/mulan/home/shanglu/GENOA/analysis/prediction/HIS"
load(paste0(path_HIS,"/HIS_to_use_chr_",chr,".RData"))

AFA_to_use$snps = as.character(AFA_to_use$snps)
AFA_to_use$gene = as.character(AFA_to_use$gene)
CAU_to_use$snps = as.character(CAU_to_use$snps)
CAU_to_use$gene = as.character(CAU_to_use$gene)
HIS_to_use$snps = as.character(HIS_to_use$snps)
HIS_to_use$gene = as.character(HIS_to_use$gene)


pop1 <- EA_full_table
pop1fdr05 <- dplyr::filter(pop1, p_wald <= EA_thr) 

pop2 <- AFA_to_use

pop2tested <- inner_join(pop1fdr05,pop2,by=c("snps","gene"))
pop2pvals_1 <- c(pop2pvals_1, pop2tested$pvalue)


pop2 <- CAU_to_use
pop2tested <- inner_join(pop1fdr05,pop2,by=c("snps","gene"))
pop2pvals_2 <- c(pop2pvals_2, pop2tested$pvalue)


pop2 <- HIS_to_use
pop2tested <- inner_join(pop1fdr05,pop2,by=c("snps","gene"))
pop2pvals_3 <- c(pop2pvals_3, pop2tested$pvalue)


}

length(pop2pvals_1)
length(pop2pvals_2)
length(pop2pvals_3)
> length(pop2pvals_1)
[1] 185339
> length(pop2pvals_2)
[1] 192689
> length(pop2pvals_3)
[1] 194417


# EA - AFA
# qobjpop2 <- qvalue(p = pop2pvals_1)
qobjpop2=qvalue_truncp(pop2pvals_1)
pi1 <- 1 - qobjpop2$pi0
pi1
[1] 0.513567


# EA - CAU
#qobjpop2 <- qvalue(p = pop2pvals_2)
qobjpop2=qvalue_truncp(pop2pvals_2)
pi1 <- 1 - qobjpop2$pi0
pi1
[1] 0.7267128


# EA - HIS
 # qobjpop2 <- qvalue(p = pop2pvals_3)
qobjpop2=qvalue_truncp(pop2pvals_3)
pi1 <- 1 - qobjpop2$pi0
pi1

[1] 0.7051797


```


# try AFA, CAU, HIS vs YRI and EUR
```R
# AFA vs EUR, YRI

pop2pvals_1 = c()
pop2pvals_2 = c()
for(chr in 1:22){
print(chr)
path_Geuvadis = "/net/mulan/home/shanglu/GENOA/analysis/compare_study"
load(paste0(path_Geuvadis,"/geu_yri_chr_",chr,".RData"))
load(paste0(path_Geuvadis,"/geu_eur_chr_",chr,".RData"))

geu_yri_chr$gene = as.character(geu_yri_chr$GENE_ID)
gene_name = as.character(geu_yri_chr$GENE_ID)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
geu_yri_chr$gene=unlist(gene_name02)
geu_yri_chr$ps=geu_yri_chr$SNPpos

geu_eur_chr$gene = as.character(geu_eur_chr$GENE_ID)
gene_name = as.character(geu_eur_chr$GENE_ID)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
geu_eur_chr$gene=unlist(gene_name02)
geu_eur_chr$ps=geu_eur_chr$SNPpos

path_AFA="/net/mulan/home/shanglu/GENOA/analysis/prediction/AFA"
load(paste0(path_AFA,"/AFA_to_use_chr_",chr,".RData"))
AFA_to_use$snps = as.character(AFA_to_use$snps)
AFA_to_use$gene = as.character(AFA_to_use$gene)
AFA_to_use$ps = AFA_to_use$pos_snps

pop1 <- AFA_to_use
pop2 <- geu_yri_chr
pop1fdr05 <- dplyr::filter(pop1, FDR<=0.05) 
pop2tested <- inner_join(pop1fdr05,pop2,by=c("ps","gene"))
pop2pvals_1 <- c(pop2pvals_1, pop2tested$pvalue.y)

pop1 <- AFA_to_use
pop2 <- geu_eur_chr
pop1fdr05 <- dplyr::filter(pop1, FDR<=0.05) 
pop2tested <- inner_join(pop1fdr05,pop2,by=c("ps","gene"))
pop2pvals_2 <- c(pop2pvals_2, pop2tested$pvalue.y)


}

# AFA - YRI
qobjpop2=qvalue_truncp(pop2pvals_1)
pi1 <- 1 - qobjpop2$pi0
pi1
[1] 0.7494572

# AFA - EUR
qobjpop2=qvalue_truncp(pop2pvals_2)
pi1 <- 1 - qobjpop2$pi0
pi1
[1] 0.8909875


#-------------
# CAU vs YRI and EUR
#-------------

pop2pvals_1 = c()
pop2pvals_2 = c()
for(chr in 1:22){
print(chr)
path_Geuvadis = "/net/mulan/home/shanglu/GENOA/analysis/compare_study"
load(paste0(path_Geuvadis,"/geu_yri_chr_",chr,".RData"))
load(paste0(path_Geuvadis,"/geu_eur_chr_",chr,".RData"))

geu_yri_chr$gene = as.character(geu_yri_chr$GENE_ID)
gene_name = as.character(geu_yri_chr$GENE_ID)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
geu_yri_chr$gene=unlist(gene_name02)
geu_yri_chr$ps=geu_yri_chr$SNPpos

geu_eur_chr$gene = as.character(geu_eur_chr$GENE_ID)
gene_name = as.character(geu_eur_chr$GENE_ID)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
geu_eur_chr$gene=unlist(gene_name02)
geu_eur_chr$ps=geu_eur_chr$SNPpos

path_CAU="/net/mulan/home/shanglu/GENOA/analysis/prediction/CAU"
load(paste0(path_CAU,"/CAU_to_use_chr_",chr,".RData"))
CAU_to_use$snps = as.character(CAU_to_use$snps)
CAU_to_use$gene = as.character(CAU_to_use$gene)
CAU_to_use$ps = CAU_to_use$pos_snps

pop1 <- CAU_to_use
pop2 <- geu_yri_chr
pop1fdr05 <- dplyr::filter(pop1, FDR<=0.05) 
pop2tested <- inner_join(pop1fdr05,pop2,by=c("ps","gene"))
pop2pvals_1 <- c(pop2pvals_1, pop2tested$pvalue.y)

pop1 <- CAU_to_use
pop2 <- geu_eur_chr
pop1fdr05 <- dplyr::filter(pop1, FDR<=0.05) 
pop2tested <- inner_join(pop1fdr05,pop2,by=c("ps","gene"))
pop2pvals_2 <- c(pop2pvals_2, pop2tested$pvalue.y)


}
# CAU - YRI
qobjpop2=qvalue_truncp(pop2pvals_1)
pi1 <- 1 - qobjpop2$pi0
pi1
[1] 0.6800789

# CAU - EUR
qobjpop2=qvalue_truncp(pop2pvals_2)
pi1 <- 1 - qobjpop2$pi0
pi1
[1] 0.8386803

#-------------
# HIS vs YRI and EUR
#-------------
pop2pvals_1 = c()
pop2pvals_2 = c()
for(chr in 1:22){
print(chr)
path_Geuvadis = "/net/mulan/home/shanglu/GENOA/analysis/compare_study"
load(paste0(path_Geuvadis,"/geu_yri_chr_",chr,".RData"))
load(paste0(path_Geuvadis,"/geu_eur_chr_",chr,".RData"))

geu_yri_chr$gene = as.character(geu_yri_chr$GENE_ID)
gene_name = as.character(geu_yri_chr$GENE_ID)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
geu_yri_chr$gene=unlist(gene_name02)
geu_yri_chr$ps=geu_yri_chr$SNPpos

geu_eur_chr$gene = as.character(geu_eur_chr$GENE_ID)
gene_name = as.character(geu_eur_chr$GENE_ID)
gene_name01 = strsplit(gene_name,"[.]")
gene_name02 = c()
for(i in 1:length(gene_name01)){gene_name02[i]=paste0(gene_name01[[i]][1])}
geu_eur_chr$gene=unlist(gene_name02)
geu_eur_chr$ps=geu_eur_chr$SNPpos

path_HIS="/net/mulan/home/shanglu/GENOA/analysis/prediction/HIS"
load(paste0(path_HIS,"/HIS_to_use_chr_",chr,".RData"))
HIS_to_use$snps = as.character(HIS_to_use$snps)
HIS_to_use$gene = as.character(HIS_to_use$gene)
HIS_to_use$ps = HIS_to_use$pos_snps

pop1 <- HIS_to_use
pop2 <- geu_yri_chr
pop1fdr05 <- dplyr::filter(pop1, FDR<=0.05) 
pop2tested <- inner_join(pop1fdr05,pop2,by=c("ps","gene"))
pop2pvals_1 <- c(pop2pvals_1, pop2tested$pvalue.y)

pop1 <- HIS_to_use
pop2 <- geu_eur_chr
pop1fdr05 <- dplyr::filter(pop1, FDR<=0.05) 
pop2tested <- inner_join(pop1fdr05,pop2,by=c("ps","gene"))
pop2pvals_2 <- c(pop2pvals_2, pop2tested$pvalue.y)

}



# HIS - YRI
qobjpop2=qvalue_truncp(pop2pvals_1)
pi1 <- 1 - qobjpop2$pi0
pi1
[1] 0.7327489

# HIS - EUR
qobjpop2=qvalue_truncp(pop2pvals_2)
pi1 <- 1 - qobjpop2$pi0
pi1
[1] 0.8520439
```
