

# Gene length as a potential source of bias: 
The authors conduct a couple of analyses that use statistics that are dependent on the lengths of genes 
and/or the number of SNPs linked to a given gene. 
First, they define eGenes as genes that harbor at least one significant eQTL. 

Genes with a larger number of SNPs in close proximity (possibly these genes are also on average longer) 
are under neutrality having lower "lowest p-values" (assuming a uniform distribution of P values). 
The authors even report this phenomena in paragraph 2 page 10. 

I find this observation troubling since it looks like a technical, rather than a biological cause for it. 

In order to possibly investigate whether the source of that observation is technical or biological 
one could for example 

(a) use extreme value statistics and model the data accordingly or 

(b) use a sub-sampling approach (of for example equal number of SNPs per gene) to 
correct this bias and check whether the lowest P values these sub-set samples are still correlated with gene length 
and SNP density. While I am not certain that there's a general solution for that issue, 
I - at a minimum - would expect that this problem is acknowledged in the text. 
Also, all subsequent analyses, regarding conservation, etc. will of course carry along that initial potential bias, 
which would possibly question the robustness of these results as well. 

# version 1: use same number of SNPs
# AA
```R
AA_thr = 6.245907e-05
EA_thr = 0.0001385504

path_AA_result="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping/result"
path_EA_result="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping/result"


# get genes 
load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order_correct.RData")
load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order_correct.RData")


summary(gene_AA_anno_order$SNPnum)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    1.0   502.0   722.0   825.8   984.0 19808.0      44 
    

gene_AA_anno_subset = gene_AA_anno_order[which(gene_AA_anno_order$SNPnum >=median(na.omit(gene_AA_anno_order$SNPnum))),]

>  dim(gene_AA_anno_subset)
[1] 8802   38


FDR_AA_mat = matrix(0,  dim(gene_AA_anno_subset)[1], 10)

for(k in 1:10){

print(k)

pheno=2
res_sub = list()
sample_index = list()
count = 0
for(i in 1:22){
print(i)
load(paste0(path_AA_result,"/AA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_AA_anno_subset[which(gene_AA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
sample_index[[count]] = sample(1:dim(res[[j]])[1], 722)
res_sub[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}
}

res_permu = list()
pheno=3
print(pheno)
count = 0
for(i in 1:22){
print(i)
load(paste0(path_AA_result,"/AA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_AA_anno_subset[which(gene_AA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}
}
for(pheno in 4:12){
count = 0
print(pheno)
for(i in 1:22){
print(i)
load(paste0(path_AA_result,"/AA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_AA_anno_subset[which(gene_AA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = c(res_permu[[count]], min(res[[j]]$p_wald[sample_index[[count]]]))
}
}
}
min_p_AA = unlist(res_sub)
order_p_AA = min_p_AA[order(min_p_AA)]
min_p_AA_permu_mat = matrix(unlist(res_permu), ncol = length(order_p_AA), byrow = TRUE) 
count = 0 
num_real = NULL 
FDR_AA = NULL 
for(threshold in min_p_AA){ 
count = count + 1 
num_real[count] = sum(min_p_AA<=threshold) 
FDR_AA[count] = mean(apply(min_p_AA_permu_mat, 1, function(x) sum(x<=threshold)/num_real[count])) }

FDR_AA_mat[,k] = FDR_AA

}

eGene_sub_AA_mat = (FDR_AA_mat<0.05)*1
gene_AA_anno_subset = cbind(gene_AA_anno_subset,eGene_sub_AA_mat)
summary(colSums(eGene_sub_AA_mat))
table(rowMeans(eGene_sub_AA_mat))
sum(gene_AA_anno_subset$eGene)

gene_AA_anno_subset$res_sub_pval = unlist(res_sub)

# v1, gene length
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value
# v2, SNP density defined by SNPnum/distance
gene_AA_anno_subset$SNPdensity = gene_AA_anno_subset$SNPnum/gene_AA_anno_subset$distance
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
# v3, SNP density defined by SNPnum/distance
gene_AA_anno_subset$SNPdensity = gene_AA_anno_subset$SNPnum/gene_AA_anno_subset$distance
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value

cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value


save(FDR_AA_mat, file = "FDR_AA_mat_v1.RData")
save(sample_index, file = "sample_index_v1.RData")
save(gene_AA_anno_subset, file = "gene_AA_anno_subset_v1.RData")




```
# EA
```R
pathAAeqtl = "/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
pathEAeqtl = "/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"
pathcompare = "/net/mulan/home/shanglu/GENOA/analysis/compare"
path_EA_result="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping/result"
load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order_correct.RData")
> summary(gene_EA_anno_order$SNPnum)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    1.0   282.5   418.0   491.4   589.0 13379.0      17 
gene_EA_anno_subset = gene_EA_anno_order[which(gene_EA_anno_order$SNPnum >=median(na.omit(gene_EA_anno_order$SNPnum))),]
dim(gene_EA_anno_subset)
[1] 8694   38


FDR_EA_mat = matrix(0, dim(gene_EA_anno_subset)[1], 10)

#------------------------------------------------------------------------------------
# repeat 10 times to account for stochasticity in the down-sampling process
#------------------------------------------------------------------------------------
 
for(k in 1:10){

print(k)

pheno=2
res_sub = list()
sample_index = list()
count = 0
for(i in 1:22){
print(i)
load(paste0(path_EA_result,"/EA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_EA_anno_subset[which(gene_EA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
sample_index[[count]] = sample(1:dim(res[[j]])[1], median(na.omit(gene_EA_anno_order$SNPnum)))
res_sub[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}
}

res_permu = list()
pheno=3
print(pheno)
count = 0
for(i in 1:22){
print(i)
load(paste0(path_EA_result,"/EA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_EA_anno_subset[which(gene_EA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}
}
for(pheno in 4:12){
count = 0
print(pheno)
for(i in 1:22){
print(i)
load(paste0(path_EA_result,"/EA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_EA_anno_subset[which(gene_EA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = c(res_permu[[count]], min(res[[j]]$p_wald[sample_index[[count]]]))
}
}
}
min_p_EA = unlist(res_sub)
order_p_EA = min_p_EA[order(min_p_EA)]
min_p_EA_permu_mat = matrix(unlist(res_permu), ncol = length(order_p_EA), byrow = TRUE) 
count = 0 
num_real = NULL 
FDR_EA = NULL 
for(threshold in min_p_EA){ 
count = count + 1 
num_real[count] = sum(min_p_EA<=threshold) 
FDR_EA[count] = mean(apply(min_p_EA_permu_mat, 1, function(x) sum(x<=threshold)/num_real[count])) }

FDR_EA_mat[,k] = FDR_EA

}


eGene_sub_EA_mat = (FDR_EA_mat<0.05)*1
gene_EA_anno_subset = cbind(gene_EA_anno_subset,eGene_sub_EA_mat)
summary(colSums(eGene_sub_EA_mat))
table(rowMeans(eGene_sub_EA_mat))
sum(gene_EA_anno_subset$eGene)

gene_EA_anno_subset$res_sub_pval = unlist(res_sub)

# v1, gene length
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$estimate
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$p.value
# v2, SNP density defined by SNPnum/distance
gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
# v3, SNP density defined by SNPnum/distance
gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value

save(FDR_EA_mat, file = "FDR_EA_mat_v1.RData")
save(sample_index, file = "sample_index_v1.RData")
save(gene_EA_anno_subset, file = "gene_EA_anno_subset_v1.RData")

```


# version 2: use same density of SNPs

```R
gene_AA_anno_order$SNP_density = gene_AA_anno_order$SNPnum/gene_AA_anno_order$distance
> summary(gene_AA_anno_order$SNP_density )
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
0.00002 0.01017 0.02157 0.07757 0.05325 5.60156      44 

gene_AA_anno_subset = gene_AA_anno_order[which(gene_AA_anno_order$SNP_density >= 0.02157),]
dim(gene_AA_anno_subset)
[1] 8787   39
> summary(gene_AA_anno_subset$SNPnum)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   26.0   481.0   656.0   680.9   844.0  3851.0 
   

FDR_AA_mat = matrix(0, 8787, 10)

for(k in 1:10){

print(k)

pheno=2
res_sub = list()
sample_index = list()
count = 0
for(i in 1:22){
print(i)
load(paste0(path_AA_result,"/AA_eqtl_res_chr_",i,"pheno",pheno,".RData"))

tmp = gene_AA_anno_subset[which(gene_AA_anno_subset$chr == i),]
ind = tmp$index_j
count_j = 0
for(j in ind){
count_j = count_j + 1
count = count + 1
sample_index[[count]] = sample(1:dim(res[[j]])[1], floor(tmp[count_j,]$distance* 0.02157))
res_sub[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}

}


res_permu = list()
pheno=3
print(pheno)
count = 0
for(i in 1:22){
print(i)
load(paste0(path_AA_result,"/AA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_AA_anno_subset[which(gene_AA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}
}
for(pheno in 4:12){
count = 0
print(pheno)
for(i in 1:22){
print(i)
load(paste0(path_AA_result,"/AA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_AA_anno_subset[which(gene_AA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = c(res_permu[[count]], min(res[[j]]$p_wald[sample_index[[count]]]))
}
}
}
min_p_AA = unlist(res_sub)
order_p_AA = min_p_AA[order(min_p_AA)]
min_p_AA_permu_mat = matrix(unlist(res_permu), ncol = length(order_p_AA), byrow = TRUE) 
count = 0 
num_real = NULL 
FDR_AA = NULL 
for(threshold in min_p_AA){ 
count = count + 1 
num_real[count] = sum(min_p_AA<=threshold) 
FDR_AA[count] = mean(apply(min_p_AA_permu_mat, 1, function(x) sum(x<=threshold)/num_real[count])) }

FDR_AA_mat[,k] = FDR_AA

}

eGene_sub_AA_mat = (FDR_AA_mat<0.05)*1
gene_AA_anno_subset = cbind(gene_AA_anno_subset,eGene_sub_AA_mat)
summary(colSums(eGene_sub_AA_mat))
table(rowMeans(eGene_sub_AA_mat))
sum(gene_AA_anno_subset$eGene)

gene_AA_anno_subset$res_sub_pval = unlist(res_sub)

# v1, gene length
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value
# v2, SNP density defined by SNPnum/distance
gene_AA_anno_subset$SNPdensity = gene_AA_anno_subset$SNPnum/gene_AA_anno_subset$distance
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
# v3, SNP density defined by SNPnum/distance
gene_AA_anno_subset$SNPdensity = gene_AA_anno_subset$SNPnum/gene_AA_anno_subset$distance
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value

save(FDR_AA_mat, file = "FDR_AA_mat_v2.RData")
save(sample_index, file = "sample_index_v2.RData")
save(gene_AA_anno_subset, file = "gene_AA_anno_subset_v2.RData")


```

# EA
```R
gene_EA_anno_order$SNP_density = gene_EA_anno_order$SNPnum/gene_EA_anno_order$distance
> summary(gene_EA_anno_order$SNP_density )
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
0.000004 0.005860 0.012579 0.045805 0.031795 2.855072       17 

gene_EA_anno_subset = gene_EA_anno_order[which(gene_EA_anno_order$SNP_density >= median(na.omit(gene_EA_anno_order$SNP_density))),]
> dim(gene_EA_anno_subset)
[1] 8672   39

> summary(gene_EA_anno_subset$SNPnum)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   17.0   280.0   391.0   419.5   525.0  3595.0 
   

FDR_EA_mat = matrix(0, dim(gene_EA_anno_subset)[1], 10)

for(k in 1:10){

print(k)

pheno=2
res_sub = list()
sample_index = list()
count = 0
for(i in 1:22){
print(i)
load(paste0(path_EA_result,"/EA_eqtl_res_chr_",i,"pheno",pheno,".RData"))

tmp = gene_EA_anno_subset[which(gene_EA_anno_subset$chr == i),]
ind = tmp$index_j
count_j = 0
for(j in ind){
count_j = count_j + 1
count = count + 1
sample_index[[count]] = sample(1:dim(res[[j]])[1], floor(tmp[count_j,]$distance*median(na.omit(gene_EA_anno_order$SNP_density))))
res_sub[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}
}


res_permu = list()
pheno=3
print(pheno)
count = 0
for(i in 1:22){
print(i)
load(paste0(path_EA_result,"/EA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_EA_anno_subset[which(gene_EA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}
}
for(pheno in 4:12){
count = 0
print(pheno)
for(i in 1:22){
print(i)
load(paste0(path_EA_result,"/EA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_EA_anno_subset[which(gene_EA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = c(res_permu[[count]], min(res[[j]]$p_wald[sample_index[[count]]]))
}
}
}
min_p_EA = unlist(res_sub)
order_p_EA = min_p_EA[order(min_p_EA)]
min_p_EA_permu_mat = matrix(unlist(res_permu), ncol = length(order_p_EA), byrow = TRUE) 
count = 0 
num_real = NULL 
FDR_EA = NULL 
for(threshold in min_p_EA){ 
count = count + 1 
num_real[count] = sum(min_p_EA<=threshold) 
FDR_EA[count] = mean(apply(min_p_EA_permu_mat, 1, function(x) sum(x<=threshold)/num_real[count])) }

FDR_EA_mat[,k] = FDR_EA

}


eGene_sub_EA_mat = (FDR_EA_mat<0.05)*1
gene_EA_anno_subset = cbind(gene_EA_anno_subset,eGene_sub_EA_mat)
summary(colSums(eGene_sub_EA_mat))
table(rowMeans(eGene_sub_EA_mat))
sum(gene_EA_anno_subset$eGene)

gene_EA_anno_subset$res_sub_pval = unlist(res_sub)

# v1, gene length
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$estimate
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$p.value
# v2, SNP density defined by SNPnum/distance
gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
# v3, SNP density defined by SNPnum/distance
gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value

save(FDR_EA_mat, file = "FDR_EA_mat_v2.RData")
save(sample_index, file = "sample_index_v2.RData")
save(gene_EA_anno_subset, file = "gene_EA_anno_subset_v2.RData")
 
 

 
```

# version 3: use high density SNPs, but use same number of SNPs in each gene

```R

gene_AA_anno_order$SNP_density = gene_AA_anno_order$SNPnum/gene_AA_anno_order$distance
>  summary(gene_AA_anno_order$SNP_density )
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
0.00002 0.01017 0.02157 0.07757 0.05325 5.60156      44 

gene_AA_anno_subset = gene_AA_anno_order[which(gene_AA_anno_order$SNP_density >=0.02157),]
>  dim(gene_AA_anno_subset)
[1] 8787   39
> summary(gene_AA_anno_subset$SNPnum)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   26.0   481.0   656.0   680.9   844.0  3851.0 
   

FDR_AA_mat = matrix(0, dim(gene_AA_anno_subset)[1], 10)

for(k in 1:10){

print(k)

pheno=2
res_sub = list()
sample_index = list()
count = 0
for(i in 1:22){
print(i)
load(paste0(path_AA_result,"/AA_eqtl_res_chr_",i,"pheno",pheno,".RData"))

tmp = gene_AA_anno_subset[which(gene_AA_anno_subset$chr == i),]
ind = tmp$index_j
count_j = 0
for(j in ind){
count_j = count_j + 1
count = count + 1
sample_index[[count]] = sample(1:dim(res[[j]])[1], 26)
res_sub[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}

}


res_permu = list()
pheno=3
print(pheno)
count = 0
for(i in 1:22){
print(i)
load(paste0(path_AA_result,"/AA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_AA_anno_subset[which(gene_AA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}
}
for(pheno in 4:12){
count = 0
print(pheno)
for(i in 1:22){
print(i)
load(paste0(path_AA_result,"/AA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_AA_anno_subset[which(gene_AA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = c(res_permu[[count]], min(res[[j]]$p_wald[sample_index[[count]]]))
}
}
}
min_p_AA = unlist(res_sub)
order_p_AA = min_p_AA[order(min_p_AA)]
min_p_AA_permu_mat = matrix(unlist(res_permu), ncol = length(order_p_AA), byrow = TRUE) 
count = 0 
num_real = NULL 
FDR_AA = NULL 
for(threshold in min_p_AA){ 
count = count + 1 
num_real[count] = sum(min_p_AA<=threshold) 
FDR_AA[count] = mean(apply(min_p_AA_permu_mat, 1, function(x) sum(x<=threshold)/num_real[count])) }

FDR_AA_mat[,k] = FDR_AA

}


eGene_sub_AA_mat = (FDR_AA_mat<0.05)*1
gene_AA_anno_subset = cbind(gene_AA_anno_subset,eGene_sub_AA_mat)
summary(colSums(eGene_sub_AA_mat))
table(rowMeans(eGene_sub_AA_mat))
sum(gene_AA_anno_subset$eGene)

gene_AA_anno_subset$res_sub_pval = unlist(res_sub)

# v1, gene length
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value
# v2, SNP density defined by SNPnum/distance
gene_AA_anno_subset$SNPdensity = gene_AA_anno_subset$SNPnum/gene_AA_anno_subset$distance
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
# v3, SNP density defined by SNPnum/distance
gene_AA_anno_subset$SNPdensity = gene_AA_anno_subset$SNPnum/gene_AA_anno_subset$distance
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value

save(FDR_AA_mat, file = "FDR_AA_mat_v3.RData")
save(sample_index, file = "sample_index_v3.RData")
save(gene_AA_anno_subset, file = "gene_AA_anno_subset_v3.RData")


cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value
cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value

```

# EA
```R
gene_EA_anno_order$SNP_density = gene_EA_anno_order$SNPnum/gene_EA_anno_order$distance
>  summary(gene_EA_anno_order$SNP_density )
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
0.000004 0.005860 0.012579 0.045805 0.031795 2.855072       17

gene_EA_anno_subset = gene_EA_anno_order[which(gene_EA_anno_order$SNP_density >=0.012579),]
> dim(gene_EA_anno_subset)
[1] 8672   39
>  summary(gene_EA_anno_subset$SNPnum)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   17.0   280.0   391.0   419.5   525.0  3595.0 
   

FDR_EA_mat = matrix(0,  dim(gene_EA_anno_subset)[1], 10)

for(k in 1:10){

print(k)

pheno=2
res_sub = list()
sample_index = list()
count = 0
for(i in 1:22){
print(i)
load(paste0(path_EA_result,"/EA_eqtl_res_chr_",i,"pheno",pheno,".RData"))

tmp = gene_EA_anno_subset[which(gene_EA_anno_subset$chr == i),]
ind = tmp$index_j
count_j = 0
for(j in ind){
count_j = count_j + 1
count = count + 1
sample_index[[count]] = sample(1:dim(res[[j]])[1], min(gene_EA_anno_subset$SNPnum))
res_sub[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}
}



res_permu = list()
pheno=3
print(pheno)
count = 0
for(i in 1:22){
print(i)
load(paste0(path_EA_result,"/EA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_EA_anno_subset[which(gene_EA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = min(res[[j]]$p_wald[sample_index[[count]]])
}
}
for(pheno in 4:12){
count = 0
print(pheno)
for(i in 1:22){
print(i)
load(paste0(path_EA_result,"/EA_eqtl_res_chr_",i,"pheno",pheno,".RData"))
ind = gene_EA_anno_subset[which(gene_EA_anno_subset$chr == i),]$index_j
for(j in ind){
count = count + 1
res_permu[[count]] = c(res_permu[[count]], min(res[[j]]$p_wald[sample_index[[count]]]))
}
}
}
min_p_EA = unlist(res_sub)
order_p_EA = min_p_EA[order(min_p_EA)]
min_p_EA_permu_mat = matrix(unlist(res_permu), ncol = length(order_p_EA), byrow = TRUE) 
count = 0 
num_real = NULL 
FDR_EA = NULL 
for(threshold in min_p_EA){ 
count = count + 1 
num_real[count] = sum(min_p_EA<=threshold) 
FDR_EA[count] = mean(apply(min_p_EA_permu_mat, 1, function(x) sum(x<=threshold)/num_real[count])) }

FDR_EA_mat[,k] = FDR_EA

}


eGene_sub_EA_mat = (FDR_EA_mat<0.05)*1
gene_EA_anno_subset = cbind(gene_EA_anno_subset,eGene_sub_EA_mat)
summary(colSums(eGene_sub_EA_mat))
table(rowMeans(eGene_sub_EA_mat))
sum(gene_EA_anno_subset$eGene)

gene_EA_anno_subset$res_sub_pval = unlist(res_sub)

# v1, gene length
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$estimate
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$p.value
# v2, SNP density defined by SNPnum/distance
gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
# v3, SNP density defined by SNPnum/distance
gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value

save(FDR_EA_mat, file = "FDR_EA_mat_v3.RData")
save(sample_index, file = "sample_index_v3.RData")
save(gene_EA_anno_subset, file = "gene_EA_anno_subset_v3.RData")


```


AA v1
```R
>
> eGene_sub_AA_mat = (FDR_AA_mat<0.05)*1
> gene_AA_anno_subset = cbind(gene_AA_anno_subset,eGene_sub_AA_mat)
> summary(colSums(eGene_sub_AA_mat))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   3194    3208    3211    3213    3219    3231
> table(rowMeans(eGene_sub_AA_mat))

   0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9    1
5462    8   11   26   31   29   47   52   71  107 2958
> sum(gene_AA_anno_subset$eGene)
[1] 3263
>
> gene_AA_anno_subset$res_sub_pval = unlist(res_sub)
>
> # v1, gene length
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
       rho
-0.1072658
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value
[1] 6.038687e-24


> # v2, SNP density defined by SNPnum/distance
> gene_AA_anno_subset$SNPdensity = gene_AA_anno_subset$SNPnum/gene_AA_anno_subset$distance
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho
0.1033526
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 2.451263e-22


> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho
0.1352047
> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 3.475974e-37

> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
      rho
-0.146395
> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value
[1] 2.308678e-43

```



EA v1
```R
> eGene_sub_EA_mat = (FDR_EA_mat<0.05)*1
> gene_EA_anno_subset = cbind(gene_EA_anno_subset,eGene_sub_EA_mat)
> summary(colSums(eGene_sub_EA_mat))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   2631    2647    2650    2650    2656    2661
> table(rowMeans(eGene_sub_EA_mat))

   0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9    1
5957   10   10   15   19   21   20   45   47   66 2484
> sum(gene_EA_anno_subset$eGene)
[1] 2691
>
> gene_EA_anno_subset$res_sub_pval = unlist(res_sub)
>
> # v1, gene length
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$estimate
      rho
-0.109128
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$p.value
[1] 1.896946e-24
> # v2, SNP density defined by SNPnum/distance
> gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho
0.0937808
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 1.912509e-18
> # v3, SNP density defined by SNPnum/distance
> gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho
0.0937808
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 1.912509e-18

cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance, method = "spearman")$estimate
cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance, method = "spearman")$p.value
cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value

> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance, method = "spearman")$estimate
       rho 
-0.1413658 
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance,  :
  Cannot compute exact p-value with ties
> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance, method = "spearman")$p.value
[1] 4.768333e-40
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance,  :
  Cannot compute exact p-value with ties
> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho 
0.1176226 
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity,  :
  Cannot compute exact p-value with ties
> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 3.657028e-28
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity,  :
  Cannot compute exact p-value with ties
> 



  ```

AA v2
```R
> eGene_sub_AA_mat = (FDR_AA_mat<0.05)*1
> gene_AA_anno_subset = cbind(gene_AA_anno_subset,eGene_sub_AA_mat)
> summary(colSums(eGene_sub_AA_mat))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   2397    2409    2414    2411    2417    2420
> table(rowMeans(eGene_sub_AA_mat))

   0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9    1
6156   46   30   41   44   54   57   66   83  130 2080
> sum(gene_AA_anno_subset$eGene)
[1] 2446
>
> gene_AA_anno_subset$res_sub_pval = unlist(res_sub)
>
> # v1, gene length
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
       rho
-0.4581551
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value
[1] 0
> # v2, SNP density defined by SNPnum/distance
> gene_AA_anno_subset$SNPdensity = gene_AA_anno_subset$SNPnum/gene_AA_anno_subset$distance
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho
0.3728191
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 7.460189e-288
> # v3, SNP density defined by SNPnum/distance
> gene_AA_anno_subset$SNPdensity = gene_AA_anno_subset$SNPnum/gene_AA_anno_subset$distance
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho
0.3728191
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 7.460189e-288


cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value
cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value

> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
       rho 
-0.2668333 
Warning message:
In cor.test.default(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance,  :
  Cannot compute exact p-value with ties
> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value
[1] 3.994753e-143
Warning message:
In cor.test.default(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance,  :
  Cannot compute exact p-value with ties
> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho 
0.1576073 
Warning message:
In cor.test.default(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity,  :
  Cannot compute exact p-value with ties
> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 5.575266e-50
Warning message:
In cor.test.default(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity,  :
  Cannot compute exact p-value with ties
```

EA v2
```R
> eGene_sub_EA_mat = (FDR_EA_mat<0.05)*1
> gene_EA_anno_subset = cbind(gene_EA_anno_subset,eGene_sub_EA_mat)
> summary(colSums(eGene_sub_EA_mat))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   1987    2000    2004    2004    2008    2019
> table(rowMeans(eGene_sub_EA_mat))

   0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9    1
6508   32   37   27   31   24   47   51   46   83 1786
> sum(gene_EA_anno_subset$eGene)
[1] 2022
>
> gene_EA_anno_subset$res_sub_pval = unlist(res_sub)
>
> # v1, gene length
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$estimate
       rho
-0.4251828
Warning message:
In cor.test.default(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance,  :
  Cannot compute exact p-value with ties
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$p.value
[1] 0
Warning message:
In cor.test.default(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance,  :
  Cannot compute exact p-value with ties
> # v2, SNP density defined by SNPnum/distance
> gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
     rho
0.326034
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 6.908268e-214
> # v3, SNP density defined by SNPnum/distance
> gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
     rho
0.326034
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 6.908268e-214

> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance, method = "spearman")$estimate
       rho 
-0.2692169 
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance,  :
  Cannot compute exact p-value with ties
> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance, method = "spearman")$p.value
[1] 7.132667e-144
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance,  :
  Cannot compute exact p-value with ties
> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho 
0.1467193 
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity,  :
  Cannot compute exact p-value with ties
> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 6.227342e-43
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity,  :
  Cannot compute exact p-value with ties
  
  ```

AA v3
```R
> eGene_sub_AA_mat = (FDR_AA_mat<0.05)*1
> gene_AA_anno_subset = cbind(gene_AA_anno_subset,eGene_sub_AA_mat)
> summary(colSums(eGene_sub_AA_mat))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   1933    1974    1982    1985    1998    2039
> table(rowMeans(eGene_sub_AA_mat))

   0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9    1
5778  363  225  199  167  188  158  163  201  315 1030
> sum(gene_AA_anno_subset$eGene)
[1] 2446
>
> gene_AA_anno_subset$res_sub_pval = unlist(res_sub)
>
> # v1, gene length
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
       rho
-0.2143909
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value
[1] 6.872302e-92
> # v2, SNP density defined by SNPnum/distance
> gene_AA_anno_subset$SNPdensity = gene_AA_anno_subset$SNPnum/gene_AA_anno_subset$distance
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho
0.1616872
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 1.531689e-52
> # v3, SNP density defined by SNPnum/distance
> gene_AA_anno_subset$SNPdensity = gene_AA_anno_subset$SNPnum/gene_AA_anno_subset$distance
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho
0.1616872
> cor.test(gene_AA_anno_subset$res_sub_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 1.531689e-52

> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$estimate
       rho
-0.2668333
> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$distance, method = "spearman")$p.value
[1] 3.994753e-143

cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value

> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho 
0.1576073 
Warning message:
In cor.test.default(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity,  :
  Cannot compute exact p-value with ties
> cor.test(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 5.575266e-50
Warning message:
In cor.test.default(gene_AA_anno_subset$topSNP_pval, gene_AA_anno_subset$SNPdensity,  :
  Cannot compute exact p-value with ties
  
```

  EA v3
 ```R
  > eGene_sub_EA_mat = (FDR_EA_mat<0.05)*1
> gene_EA_anno_subset = cbind(gene_EA_anno_subset,eGene_sub_EA_mat)
> summary(colSums(eGene_sub_EA_mat))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   1813    1838    1844    1845    1856    1873
> table(rowMeans(eGene_sub_EA_mat))

   0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9    1
6153  215  160  120  117  106  135  117  162  241 1146
> sum(gene_EA_anno_subset$eGene)
[1] 2022
>
> gene_EA_anno_subset$res_sub_pval = unlist(res_sub)
>
> # v1, gene length
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$estimate
       rho
-0.2224948
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$distance, method = "spearman")$p.value
[1] 9.983728e-98
> # v2, SNP density defined by SNPnum/distance
> gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho
0.1557918
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 3.022623e-48
> # v3, SNP density defined by SNPnum/distance
> gene_EA_anno_subset$SNPdensity = gene_EA_anno_subset$SNPnum/gene_EA_anno_subset$distance
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho
0.1557918
> cor.test(gene_EA_anno_subset$res_sub_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 3.022623e-48


> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance, method = "spearman")$estimate
       rho 
-0.2692169 
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance,  :
  Cannot compute exact p-value with ties
> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance, method = "spearman")$p.value
[1] 7.132667e-144
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$distance,  :
  Cannot compute exact p-value with ties
> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$estimate
      rho 
0.1467193 
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity,  :
  Cannot compute exact p-value with ties
> cor.test(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity, method = "spearman")$p.value
[1] 6.227342e-43
Warning message:
In cor.test.default(gene_EA_anno_subset$topSNP_pval, gene_EA_anno_subset$SNPdensity,  :
  Cannot compute exact p-value with ties
  
```



