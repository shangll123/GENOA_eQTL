

```

pathAAeqtl = "/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
pathEAeqtl = "/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"
pathcompare = "/net/mulan/home/shanglu/GENOA/analysis/compare"

load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order_correct.RData")
load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order_correct.RData")
load(paste0(pathcompare, "/gene_anno_EU.RData"))
load(paste0(pathcompare, "/gene_anno_AA.RData"))

```

## AA gene annotation
```
path_AA="/net/mulan/home/shanglu/GENOA/data/AA"
path_AA_result="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping/result"
path_AA_output="/net/mulan/home/shanglu/GENOA/data/AA/output/"
load(paste0(path_AA,"/gene_AA_anno_protein.RData"))

#------------------------------------------------------------------------------
# add SNP number, top SNP's position, allele1 and allele0, MAF, pvalue to annotation
#------------------------------------------------------------------------------

load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order.RData")
SNPnum = c()
topSNP_ps = c()
topSNP_allele1 = c()
topSNP_allele0 = c()
topSNP_af = c()
topSNP_rs = c()
topSNP_pval = c()
topSNP_beta = c()
largest_maf = c()
largest_maf_rs = c()
largest_maf_pval = c()
maxbeta_index = c()
maxbeta_rs = c()
maxbeta_pval = c()
count = 0
for(i in 1:22){
  print(i)
  load(paste0(path_AA_result,"/AA_eqtl_res_chr_",i,"pheno2.RData"))
  for(j in 1:length(res)){
    count = count + 1
    if(sum(dim(res[[j]])[1])!=0){
    top_index = which.min(res[[j]]$p_wald)
    SNPnum[count] = dim(res[[j]])[1]
    topSNP_ps[count] = res[[j]]$ps[top_index]
    topSNP_allele1[count] = as.character(res[[j]]$allele1[top_index])
    topSNP_allele0[count] = as.character(res[[j]]$allele0[top_index])
    topSNP_af[count] = res[[j]]$af[top_index]
    topSNP_rs[count] = as.character(res[[j]]$rs[top_index])
    topSNP_pval[count] = res[[j]]$p_wald[top_index]
    topSNP_beta[count] = res[[j]]$beta[top_index]
    
    largest_maf_index = which.max(res[[j]]$af)
    largest_maf_rs[count] = as.character(res[[j]]$rs[largest_maf_index])
    largest_maf[count] = res[[j]]$af[largest_maf_index]
    largest_maf_pval[count] = res[[j]]$p_wald[largest_maf_index]
    
    maxbeta_index = which.max(abs(res[[j]]$beta))
    maxbeta_rs[count] = as.character(res[[j]]$rs[maxbeta_index])
    maxbeta_pval[count] = res[[j]]$p_wald[maxbeta_index]  
   }
   }
}

      
gene_AA_anno_order$SNPnum = SNPnum
gene_AA_anno_order$topSNP_rs = topSNP_rs
gene_AA_anno_order$topSNP_ps = topSNP_ps
gene_AA_anno_order$topSNP_allele1 = topSNP_allele1
gene_AA_anno_order$topSNP_allele0 = topSNP_allele0
gene_AA_anno_order$topSNP_af = topSNP_af
gene_AA_anno_order$topSNP_pval = topSNP_pval
gene_AA_anno_order$topSNP_beta = topSNP_beta

gene_AA_anno_order$largest_maf_rs = largest_maf_rs
gene_AA_anno_order$largest_maf = largest_maf
gene_AA_anno_order$largest_maf_pval = largest_maf_pval
gene_AA_anno_order$maxbeta_rs = maxbeta_rs
gene_AA_anno_order$maxbeta_pval = maxbeta_pval

	  
save(gene_AA_anno_order,file = "/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order_correct.RData")

#---------------------------------------
# add conservation scores
#---------------------------------------


gene_anno_AA_conservation = gene_anno_AA[match(gene_AA_anno_order$GENE, gene_anno_AA$ensembl_gene_id),]

# check
#sum(gene_anno_AA_conservation[,1]==gene_AA_anno_order[,1])
#> sum(gene_anno_AA_conservation[,1]==gene_AA_anno_order[,1])
#[1] 17616

gene_AA_anno_order$strand = gene_anno_AA_conservation$strand
gene_AA_anno_order$dN = gene_anno_AA_conservation$dN
gene_AA_anno_order$dS = gene_anno_AA_conservation$dS
gene_AA_anno_order$dNdS_ratio = gene_anno_AA_conservation$dNdS_ratio
gene_AA_anno_order$AA_phylop = gene_anno_AA_conservation$AA_phylop
gene_AA_anno_order$AA_phastcon = gene_anno_AA_conservation$AA_phastcon

save(gene_AA_anno_order,file = "/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order.RData")

#---------------------------------------
# add bslmm results
#---------------------------------------

load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order.RData")
load("/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping/res_combinedAA.RData")
load("/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping/res_outchrAA.RData")



PVE_bslmm_combined = unlist(res_combined$AA_bslmm_pve)
PVEsd_bslmm_combined = unlist(res_combined$AA_bslmm_pve_sd)
PGE_bslmm_combined = unlist(res_combined$AA_bslmm_pge)
PGEsd_bslmm_combined = unlist(res_combined$AA_bslmm_pge_sd)

ind = which(unlist(lapply(res_combined$AA_bslmm_pve,length))==0)
gene_AA_anno_order$PVE_bslmm_combined = NA
gene_AA_anno_order[-ind,]$PVE_bslmm_combined = PVE_bslmm_combined

gene_AA_anno_order$PVEsd_bslmm_combined = NA
gene_AA_anno_order[-ind,]$PVEsd_bslmm_combined = PVEsd_bslmm_combined

gene_AA_anno_order$PGE_bslmm_combined = NA
gene_AA_anno_order[-ind,]$PGE_bslmm_combined = PGE_bslmm_combined

gene_AA_anno_order$PGEsd_bslmm_combined = NA
gene_AA_anno_order[-ind,]$PGEsd_bslmm_combined = PGEsd_bslmm_combined




PVE_bslmm_outchr = unlist(res_outchr$AA_bslmm_pve)
PVEsd_bslmm_outchr = unlist(res_outchr$AA_bslmm_pve_sd)
PGE_bslmm_outchr = unlist(res_outchr$AA_bslmm_pge)
PGEsd_bslmm_outchr = unlist(res_outchr$AA_bslmm_pge_sd)

ind = which(unlist(lapply(res_outchr$AA_bslmm_pve,length))==0)
gene_AA_anno_order$PVE_bslmm_outchr = NA
gene_AA_anno_order[-ind,]$PVE_bslmm_outchr = PVE_bslmm_outchr

gene_AA_anno_order$PVEsd_bslmm_outchr = NA
gene_AA_anno_order[-ind,]$PVEsd_bslmm_outchr = PVEsd_bslmm_outchr

gene_AA_anno_order$PGE_bslmm_outchr = NA
gene_AA_anno_order[-ind,]$PGE_bslmm_outchr = PGE_bslmm_outchr

gene_AA_anno_order$PGEsd_bslmm_outchr = NA
gene_AA_anno_order[-ind,]$PGEsd_bslmm_outchr = PGEsd_bslmm_outchr

save(gene_AA_anno_order,file = "/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order_correct.RData")





```

## EA gene annotation
```
path_EA="/net/mulan/home/shanglu/GENOA/data/EA"
path_EA_result="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping/result"
path_EA_output="/net/mulan/home/shanglu/GENOA/data/EA/output/"
load(paste0(path_EA,"/gene_EA_anno_protein.RData"))

#------------------------------------------------------------------------------
# add SNP number, top SNP's position, allele1 and allele0, MAF, pvalue to annotation
#------------------------------------------------------------------------------

load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order_correct.RData")
SNPnum = c()
topSNP_ps = c()
topSNP_allele1 = c()
topSNP_allele0 = c()
topSNP_af = c()
topSNP_rs = c()
topSNP_pval = c()
topSNP_beta = c()
largest_maf = c()
largest_maf_rs = c()
largest_maf_pval = c()
maxbeta_index = c()
maxbeta_rs = c()
maxbeta_pval = c()
count = 0
for(i in 1:22){
  print(i)
  load(paste0(path_EA_result,"/EA_eqtl_res_chr_",i,"pheno2.RData"))
  for(j in 1:length(res)){
    count = count + 1
    if(sum(dim(res[[j]])[1])!=0){
    top_index = which.min(res[[j]]$p_wald)
    SNPnum[count] = dim(res[[j]])[1]
    topSNP_ps[count] = res[[j]]$ps[top_index]
    topSNP_allele1[count] = as.character(res[[j]]$allele1[top_index])
    topSNP_allele0[count] = as.character(res[[j]]$allele0[top_index])
    topSNP_af[count] = res[[j]]$af[top_index]
    topSNP_rs[count] = as.character(res[[j]]$rs[top_index])
    topSNP_pval[count] = res[[j]]$p_wald[top_index]
    topSNP_beta[count] = res[[j]]$beta[top_index]
    
    largest_maf_index = which.max(res[[j]]$af)
    largest_maf_rs[count] = as.character(res[[j]]$rs[largest_maf_index])
    largest_maf[count] = res[[j]]$af[largest_maf_index]
    largest_maf_pval[count] = res[[j]]$p_wald[largest_maf_index]
    
    maxbeta_index = which.max(abs(res[[j]]$beta))
    maxbeta_rs[count] = as.character(res[[j]]$rs[maxbeta_index])
    maxbeta_pval[count] = res[[j]]$p_wald[maxbeta_index]  
   }
   }
}

      
gene_EA_anno_order$SNPnum = SNPnum
gene_EA_anno_order$topSNP_rs = topSNP_rs
gene_EA_anno_order$topSNP_ps = topSNP_ps
gene_EA_anno_order$topSNP_allele1 = topSNP_allele1
gene_EA_anno_order$topSNP_allele0 = topSNP_allele0
gene_EA_anno_order$topSNP_af = topSNP_af
gene_EA_anno_order$topSNP_pval = topSNP_pval
gene_EA_anno_order$topSNP_beta = topSNP_beta

gene_EA_anno_order$largest_maf_rs = largest_maf_rs
gene_EA_anno_order$largest_maf = largest_maf
gene_EA_anno_order$largest_maf_pval = largest_maf_pval
gene_EA_anno_order$maxbeta_rs = maxbeta_rs
gene_EA_anno_order$maxbeta_pval = maxbeta_pval

	  
save(gene_EA_anno_order,file = "/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order_correct.RData")


> SNPnum[1:10]
 [1]  61  58  50  49  49  53  93 146 133 210
> summary(SNPnum)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    1.0   282.5   418.0   491.4   589.0 13379.0      17 

 
#---------------------------------------
# add conservation scores
#---------------------------------------


gene_anno_EA_conservation = gene_anno_EU[match(gene_EA_anno_order$GENE, gene_anno_EU$ensembl_gene_id),]

# check
#sum(gene_anno_EA_conservation[,1]==gene_EA_anno_order[,1])
#> sum(gene_anno_EA_conservation[,1]==gene_EA_anno_order[,1])
#[1] 17360

gene_EA_anno_order$strand = gene_anno_EA_conservation$strand
gene_EA_anno_order$dN = gene_anno_EA_conservation$dN
gene_EA_anno_order$dS = gene_anno_EA_conservation$dS
gene_EA_anno_order$dNdS_ratio = gene_anno_EA_conservation$dNdS_ratio
gene_EA_anno_order$EA_phylop = gene_anno_EA_conservation$EU_phylop
gene_EA_anno_order$EA_phastcon = gene_anno_EA_conservation$EU_phastcon

save(gene_EA_anno_order,file = "/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order.RData")

#---------------------------------------
# add bslmm results
#---------------------------------------
load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order.RData")
load("/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping/res_combinedEA.RData")
load("/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping/res_outchrEA.RData")



PVE_bslmm_combined = unlist(res_combined$EA_bslmm_pve)
PVEsd_bslmm_combined = unlist(res_combined$EA_bslmm_pve_sd)
PGE_bslmm_combined = unlist(res_combined$EA_bslmm_pge)
PGEsd_bslmm_combined = unlist(res_combined$EA_bslmm_pge_sd)

ind = which(unlist(lapply(res_combined$EA_bslmm_pve,length))==0)
gene_EA_anno_order$PVE_bslmm_combined = NA
gene_EA_anno_order[-ind,]$PVE_bslmm_combined = PVE_bslmm_combined

gene_EA_anno_order$PVEsd_bslmm_combined = NA
gene_EA_anno_order[-ind,]$PVEsd_bslmm_combined = PVEsd_bslmm_combined

gene_EA_anno_order$PGE_bslmm_combined = NA
gene_EA_anno_order[-ind,]$PGE_bslmm_combined = PGE_bslmm_combined

gene_EA_anno_order$PGEsd_bslmm_combined = NA
gene_EA_anno_order[-ind,]$PGEsd_bslmm_combined = PGEsd_bslmm_combined




PVE_bslmm_outchr = unlist(res_outchr$EA_bslmm_pve)
PVEsd_bslmm_outchr = unlist(res_outchr$EA_bslmm_pve_sd)
PGE_bslmm_outchr = unlist(res_outchr$EA_bslmm_pge)
PGEsd_bslmm_outchr = unlist(res_outchr$EA_bslmm_pge_sd)

ind = which(unlist(lapply(res_outchr$EA_bslmm_pve,length))==0)
gene_EA_anno_order$PVE_bslmm_outchr = NA
gene_EA_anno_order[-ind,]$PVE_bslmm_outchr = PVE_bslmm_outchr

gene_EA_anno_order$PVEsd_bslmm_outchr = NA
gene_EA_anno_order[-ind,]$PVEsd_bslmm_outchr = PVEsd_bslmm_outchr

gene_EA_anno_order$PGE_bslmm_outchr = NA
gene_EA_anno_order[-ind,]$PGE_bslmm_outchr = PGE_bslmm_outchr

gene_EA_anno_order$PGEsd_bslmm_outchr = NA
gene_EA_anno_order[-ind,]$PGEsd_bslmm_outchr = PGEsd_bslmm_outchr

save(gene_EA_anno_order,file = "/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order_correct.RData")





```

# add eGene annotation

```
load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order_correct.RData")
load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order_correct.RData")

AA_thr = 6.245907e-05
EA_thr = 0.0001385504



gene_AA_anno_order$eGene = NA
gene_AA_anno_order$eGene[which(gene_AA_anno_order$topSNP_pval <= AA_thr)]=1
gene_AA_anno_order$eGene[which(gene_AA_anno_order$topSNP_pval > AA_thr)]=0

save(gene_AA_anno_order, file = "/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order_correct.RData")

> max(which(FDR_EA<0.05))
[1] 4490
> order_p_EA[ max(which(FDR_EA<0.05))]
[1] 0.0002131987


gene_EA_anno_order$eGene = NA
gene_EA_anno_order$eGene[which(gene_EA_anno_order$topSNP_pval <= EA_thr)]=1
gene_EA_anno_order$eGene[which(gene_EA_anno_order$topSNP_pval > EA_thr)]=0

save(gene_EA_anno_order, file = "/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order_correct.RData")


```


# cis- and trans- component of PVE, AA
```
# AA

library(stats)
load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order_correct.RData")

gene_AA_anno_use = gene_AA_anno_order[-which(is.na(gene_AA_anno_order$SNPnum)),]

> dim(gene_AA_anno_use)
[1] 17572    38

#-----------------------------------
# version 1, kick one chr at a time
#-----------------------------------

gene_AA_anno_use$por_cis_outchr = gene_AA_anno_use$PVE_bslmm_outchr * gene_AA_anno_use$PGE_bslmm_outchr
gene_AA_anno_use$por_trans_outchr = gene_AA_anno_use$PVE_bslmm_outchr * (1 - gene_AA_anno_use$PGE_bslmm_outchr)

cis = c(gene_AA_anno_use$por_cis_outchr, gene_AA_anno_use$por_cis_outchr[gene_AA_anno_use$eGene==0],gene_AA_anno_use$por_cis_outchr[gene_AA_anno_use$eGene==1])

trans = c(gene_AA_anno_use$por_trans_outchr, gene_AA_anno_use$por_trans_outchr[gene_AA_anno_use$eGene==0],gene_AA_anno_use$por_trans_outchr[gene_AA_anno_use$eGene==1])

additive = c(gene_AA_anno_use$PVE_bslmm_outchr, gene_AA_anno_use$PVE_bslmm_outchr[gene_AA_anno_use$eGene==0],gene_AA_anno_use$PVE_bslmm_outchr[gene_AA_anno_use$eGene==1])


effect = c(cis, trans, additive)

effect_type = c(rep("cis-PVE",length(cis)), rep("trans-PVE",length(trans)), rep("total PVE",length(additive)))
effect_type <- factor(effect_type,levels = c('cis-PVE','trans-PVE','total PVE'),ordered = TRUE)
tmp = c(rep("All genes", length(gene_AA_anno_use$por_cis_outchr)), rep("non eQTL genes", sum(gene_AA_anno_use$eGene==0)),rep("eQTL genes", sum(gene_AA_anno_use$eGene==1)))

effect_class = factor(rep(tmp, 3),levels = c('All genes','non eQTL genes','eQTL genes'),ordered = TRUE)

dat = data.frame(effect,effect_type,effect_class)

library(ggplot2)
pdf(paste0("Fig3_PVE_outchr_AA.pdf"),width=18, height=12)
ggplot(dat, aes(x = effect_class, y = effect,fill = effect_type)) +
scale_y_continuous(name = "Percent variance explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
scale_x_discrete(name = "") + 
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) + 
geom_boxplot( aes(x = effect_class, y = effect,fill = effect_type),width = 0.2,position=position_dodge(0.9)) +
ggtitle("African American") +
theme_bw(base_size = 30) +
theme(plot.title = element_text(size = 40,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()

#-----------------------------------
# version 2, combined chr
#-----------------------------------

gene_AA_anno_use$por_cis_combined = gene_AA_anno_use$PVE_bslmm_combined * gene_AA_anno_use$PGE_bslmm_combined
gene_AA_anno_use$por_trans_combined = gene_AA_anno_use$PVE_bslmm_combined * (1 - gene_AA_anno_use$PGE_bslmm_combined)

cis = c(gene_AA_anno_use$por_cis_combined, gene_AA_anno_use$por_cis_combined[gene_AA_anno_use$eGene==0],gene_AA_anno_use$por_cis_combined[gene_AA_anno_use$eGene==1])

trans = c(gene_AA_anno_use$por_trans_combined, gene_AA_anno_use$por_trans_combined[gene_AA_anno_use$eGene==0],gene_AA_anno_use$por_trans_combined[gene_AA_anno_use$eGene==1])

additive = c(gene_AA_anno_use$PVE_bslmm_combined, gene_AA_anno_use$PVE_bslmm_combined[gene_AA_anno_use$eGene==0],gene_AA_anno_use$PVE_bslmm_combined[gene_AA_anno_use$eGene==1])


save(gene_AA_anno_use, file = "/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_use.RData")

effect = c(cis, trans, additive)

effect_type = c(rep("cis-PVE",length(cis)), rep("trans-PVE",length(trans)), rep("total PVE",length(additive)))
effect_type <- factor(effect_type,levels = c('cis-PVE','trans-PVE','total PVE'),ordered = TRUE)
tmp = c(rep("All genes", length(gene_AA_anno_use$por_cis_combined)), rep("non eQTL genes", sum(gene_AA_anno_use$eGene==0)),rep("eQTL genes", sum(gene_AA_anno_use$eGene==1)))

effect_class = factor(rep(tmp, 3),levels = c('All genes','non eQTL genes','eQTL genes'),ordered = TRUE)

dat = data.frame(effect,effect_type,effect_class)

library(ggplot2)
pdf(paste0("Fig3_PVE_combined_AA.pdf"),width=18, height=12)
ggplot(dat, aes(x = effect_class, y = effect,fill = effect_type)) +
scale_y_continuous(name = "Percent variance explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
scale_x_discrete(name = "") + 
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) + 
geom_boxplot( aes(x = effect_class, y = effect,fill = effect_type),width = 0.2,position=position_dodge(0.9)) +
ggtitle("African American") +
theme_bw(base_size = 30) +
theme(plot.title = element_text(size = 40,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()

save(gene_AA_anno_use,file = "/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_use.RData")

```
# cis- and trans- component of PVE, EA
```
# EA

library(stats)
load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order_correct.RData")

gene_EA_anno_use = gene_EA_anno_order[-which(is.na(gene_EA_anno_order$SNPnum)),]

> dim(gene_EA_anno_use)
[1] 17343    38
#-----------------------------------
# version 1, kick one chr at a time
#-----------------------------------

gene_EA_anno_use$por_cis_outchr = gene_EA_anno_use$PVE_bslmm_outchr * gene_EA_anno_use$PGE_bslmm_outchr
gene_EA_anno_use$por_trans_outchr = gene_EA_anno_use$PVE_bslmm_outchr * (1 - gene_EA_anno_use$PGE_bslmm_outchr)

cis = c(gene_EA_anno_use$por_cis_outchr, gene_EA_anno_use$por_cis_outchr[gene_EA_anno_use$eGene==0],gene_EA_anno_use$por_cis_outchr[gene_EA_anno_use$eGene==1])

trans = c(gene_EA_anno_use$por_trans_outchr, gene_EA_anno_use$por_trans_outchr[gene_EA_anno_use$eGene==0],gene_EA_anno_use$por_trans_outchr[gene_EA_anno_use$eGene==1])

additive = c(gene_EA_anno_use$PVE_bslmm_outchr, gene_EA_anno_use$PVE_bslmm_outchr[gene_EA_anno_use$eGene==0],gene_EA_anno_use$PVE_bslmm_outchr[gene_EA_anno_use$eGene==1])



effect = c(cis, trans, additive)

effect_type = c(rep("cis-PVE",length(cis)), rep("trans-PVE",length(trans)), rep("total PVE",length(additive)))
effect_type <- factor(effect_type,levels = c('cis-PVE','trans-PVE','total PVE'),ordered = TRUE)
tmp = c(rep("All genes", length(gene_EA_anno_use$por_cis_outchr)), rep("non eQTL genes", sum(gene_EA_anno_use$eGene==0)),rep("eQTL genes", sum(gene_EA_anno_use$eGene==1)))

effect_class = factor(rep(tmp, 3),levels = c('All genes','non eQTL genes','eQTL genes'),ordered = TRUE)

dat = data.frame(effect,effect_type,effect_class)

library(ggplot2)
pdf(paste0("Fig3_PVE_outchr_EA.pdf"),width=18, height=12)
ggplot(dat, aes(x = effect_class, y = effect,fill = effect_type)) +
scale_y_continuous(name = "Percent variance explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
scale_x_discrete(name = "") + 
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) + 
geom_boxplot( aes(x = effect_class, y = effect,fill = effect_type),width = 0.2,position=position_dodge(0.9)) +
ggtitle("African American") +
theme_bw(base_size = 30) +
theme(plot.title = element_text(size = 40,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()

#-----------------------------------
# version 2, combined chr
#-----------------------------------

gene_EA_anno_use$por_cis_combined = gene_EA_anno_use$PVE_bslmm_combined * gene_EA_anno_use$PGE_bslmm_combined
gene_EA_anno_use$por_trans_combined = gene_EA_anno_use$PVE_bslmm_combined * (1 - gene_EA_anno_use$PGE_bslmm_combined)

cis = c(gene_EA_anno_use$por_cis_combined, gene_EA_anno_use$por_cis_combined[gene_EA_anno_use$eGene==0],gene_EA_anno_use$por_cis_combined[gene_EA_anno_use$eGene==1])

trans = c(gene_EA_anno_use$por_trans_combined, gene_EA_anno_use$por_trans_combined[gene_EA_anno_use$eGene==0],gene_EA_anno_use$por_trans_combined[gene_EA_anno_use$eGene==1])

additive = c(gene_EA_anno_use$PVE_bslmm_combined, gene_EA_anno_use$PVE_bslmm_combined[gene_EA_anno_use$eGene==0],gene_EA_anno_use$PVE_bslmm_combined[gene_EA_anno_use$eGene==1])


effect = c(cis, trans, additive)

effect_type = c(rep("cis-PVE",length(cis)), rep("trans-PVE",length(trans)), rep("total PVE",length(additive)))
effect_type <- factor(effect_type,levels = c('cis-PVE','trans-PVE','total PVE'),ordered = TRUE)
tmp = c(rep("All genes", length(gene_EA_anno_use$por_cis_combined)), rep("non eQTL genes", sum(gene_EA_anno_use$eGene==0)),rep("eQTL genes", sum(gene_EA_anno_use$eGene==1)))

effect_class = factor(rep(tmp, 3),levels = c('All genes','non eQTL genes','eQTL genes'),ordered = TRUE)

dat = data.frame(effect,effect_type,effect_class)

library(ggplot2)
pdf(paste0("Fig3_PVE_combined_EA.pdf"),width=18, height=12)
ggplot(dat, aes(x = effect_class, y = effect,fill = effect_type)) +
scale_y_continuous(name = "Percent variance explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
scale_x_discrete(name = "") + 
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) + 
geom_boxplot( aes(x = effect_class, y = effect,fill = effect_type),width = 0.2,position=position_dodge(0.9)) +
ggtitle("African American") +
theme_bw(base_size = 30) +
theme(plot.title = element_text(size = 40,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()


save(gene_EA_anno_use,file = "/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_use.RData")


```

Next, we estimated the genetic architecture underlying gene expression variation through heritability estimation and partitioning. For each gene in turn, we estimated the proportion of variance (PVE) in gene expression levels that are accounted for by all SNPs using the Bayesian sparse linear mixed model (BSLMM). This quantity is commonly referred to as SNP heritability. We used Benjamini-Hochberg false discovery rates (FDR) to correct for multiple testing 59, with an FDR <0.1 used as the threshold for significant PVE. 

```
################################################
# AA
################################################

#-------------------------------
# AA, combined version
#-------------------------------

gene_AA_anno_use_tmp = gene_AA_anno_use[gene_AA_anno_use$PVEsd_bslmm_combined != 0,]
gene_AA_anno_use_tmp$PVE_bslmm_combined_z = gene_AA_anno_use_tmp$PVE_bslmm_combined/gene_AA_anno_use_tmp$PVEsd_bslmm_combined

# FDR <- p.adjust(p, method="BH")
# example text: Nominal P-value thresholds for FDR<0.05 were: xxx

gene_AA_anno_use_tmp$PVE_bslmm_combined_p  = 2*pnorm(abs(gene_AA_anno_use_tmp$PVE_bslmm_combined_z), lower.tail = F)

porder = gene_AA_anno_use_tmp$PVE_bslmm_combined_p[order(gene_AA_anno_use_tmp$PVE_bslmm_combined_p)]

PK = c(1:length(porder))*0.1/length(porder)
tail(porder[which(porder<=PK)])
combined_pve_FDR <- p.adjust(gene_AA_anno_use_tmp$PVE_bslmm_combined_p, method="BH")
gene_AA_anno_use_tmp$combined_pve_FDR = combined_pve_FDR

summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR<0.1])
summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR>=0.1])

summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR<0.05])
summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR>=0.05])

> summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR<0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.07265 0.17679 0.22013 0.24132 0.28014 0.79457 
> summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR>=0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01418 0.04129 0.05995 0.07143 0.09370 0.22952 
> 
> 
> summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR<0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.08485 0.20108 0.24596 0.26676 0.30759 0.79457 
> summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR>=0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01418 0.04209 0.06257 0.07616 0.09977 0.27725 

#-------------------------------
# AA, out chr version
#-------------------------------


gene_AA_anno_use_tmp$PVE_bslmm_outchr_z = gene_AA_anno_use_tmp$PVE_bslmm_outchr/gene_AA_anno_use_tmp$PVEsd_bslmm_outchr

# FDR <- p.adjust(p, method="BH")
# example text: Nominal P-value thresholds for FDR<0.05 were: xxx

gene_AA_anno_use_tmp$PVE_bslmm_outchr_p  = 2*pnorm(abs(gene_AA_anno_use_tmp$PVE_bslmm_outchr_z), lower.tail = F)

porder = gene_AA_anno_use_tmp$PVE_bslmm_outchr_p[order(gene_AA_anno_use_tmp$PVE_bslmm_outchr_p)]

PK = c(1:length(porder))*0.1/length(porder)
tail(porder[which(porder<=PK)])
outchr_pve_FDR <- p.adjust(gene_AA_anno_use_tmp$PVE_bslmm_outchr_p, method="BH")
gene_AA_anno_use_tmp$outchr_pve_FDR = outchr_pve_FDR

summary(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$outchr_pve_FDR<0.1])
summary(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$outchr_pve_FDR>=0.1])

summary(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$outchr_pve_FDR<0.05])
summary(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$outchr_pve_FDR>=0.05])

> summary(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$outchr_pve_FDR<0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.06384 0.17468 0.21990 0.24033 0.27820 0.79699 
> summary(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$outchr_pve_FDR>=0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01414 0.04095 0.05990 0.07127 0.09331 0.22511 
> 
> summary(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$outchr_pve_FDR<0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.08106 0.19915 0.24634 0.26564 0.30774 0.79699 
> summary(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$outchr_pve_FDR>=0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01414 0.04172 0.06244 0.07605 0.09974 0.25060 


sum(gene_AA_anno_use_tmp$combined_pve_FDR < 0.1)
sum(gene_AA_anno_use_tmp$outchr_pve_FDR < 0.1)

sum(gene_AA_anno_use_tmp$combined_pve_FDR < 0.05)
sum(gene_AA_anno_use_tmp$outchr_pve_FDR < 0.05)

> sum(gene_AA_anno_use_tmp$combined_pve_FDR < 0.1)
[1] 2717
> sum(gene_AA_anno_use_tmp$outchr_pve_FDR < 0.1)
[1] 2734

> sum(gene_AA_anno_use_tmp$combined_pve_FDR < 0.05)
[1] 1986
> sum(gene_AA_anno_use_tmp$outchr_pve_FDR < 0.05)
[1] 1995
> 


save(gene_AA_anno_use_tmp,file = "/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_use_tmp.RData")

################################################
# EA
################################################

#-------------------------------
# EA, combined version
#-------------------------------

gene_EA_anno_use_tmp = gene_EA_anno_use[gene_EA_anno_use$PVEsd_bslmm_combined != 0,]
gene_EA_anno_use_tmp$PVE_bslmm_combined_z = gene_EA_anno_use_tmp$PVE_bslmm_combined/gene_EA_anno_use_tmp$PVEsd_bslmm_combined

# FDR <- p.adjust(p, method="BH")
# example text: Nominal P-value thresholds for FDR<0.05 were: xxx

gene_EA_anno_use_tmp$PVE_bslmm_combined_p  = 2*pnorm(abs(gene_EA_anno_use_tmp$PVE_bslmm_combined_z), lower.tail = F)

porder = gene_EA_anno_use_tmp$PVE_bslmm_combined_p[order(gene_EA_anno_use_tmp$PVE_bslmm_combined_p)]

PK = c(1:length(porder))*0.1/length(porder)
tail(porder[which(porder<=PK)])
combined_pve_FDR <- p.adjust(gene_EA_anno_use_tmp$PVE_bslmm_combined_p, method="BH")
gene_EA_anno_use_tmp$combined_pve_FDR = combined_pve_FDR

summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR<0.1])
summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR>=0.1])

summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR<0.05])
summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR>=0.05])

> summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR<0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.08774 0.19601 0.24024 0.25468 0.29191 0.78694 
> summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR>=0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01433 0.04597 0.06656 0.07899 0.10224 0.24806 
> 
> summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR<0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.09813 0.22104 0.26780 0.28246 0.32222 0.78694 
> summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR>=0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01433 0.04669 0.06896 0.08395 0.10919 0.28467 


#-------------------------------
# EA, out chr version
#-------------------------------


gene_EA_anno_use_tmp$PVE_bslmm_outchr_z = gene_EA_anno_use_tmp$PVE_bslmm_outchr/gene_EA_anno_use_tmp$PVEsd_bslmm_outchr

# FDR <- p.adjust(p, method="BH")
# example text: Nominal P-value thresholds for FDR<0.05 were: xxx

gene_EA_anno_use_tmp$PVE_bslmm_outchr_p  = 2*pnorm(abs(gene_EA_anno_use_tmp$PVE_bslmm_outchr_z), lower.tail = F)

porder = gene_EA_anno_use_tmp$PVE_bslmm_outchr_p[order(gene_EA_anno_use_tmp$PVE_bslmm_outchr_p)]

PK = c(1:length(porder))*0.1/length(porder)
tail(porder[which(porder<=PK)])
outchr_pve_FDR <- p.adjust(gene_EA_anno_use_tmp$PVE_bslmm_outchr_p, method="BH")
gene_EA_anno_use_tmp$outchr_pve_FDR = outchr_pve_FDR

summary(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR<0.1])
summary(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR>=0.1])

summary(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR<0.05])
summary(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR>=0.05])

> summary(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR<0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.08571 0.19683 0.23739 0.25391 0.29114 0.78826 
> summary(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR>=0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01391 0.04594 0.06667 0.07889 0.10208 0.25270 
> 
> summary(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR<0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.1133  0.2184  0.2675  0.2813  0.3221  0.7883 
> summary(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR>=0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01391 0.04676 0.06880 0.08382 0.10887 0.29267 



sum(gene_EA_anno_use_tmp$combined_pve_FDR < 0.1)
sum(gene_EA_anno_use_tmp$outchr_pve_FDR < 0.1)

sum(gene_EA_anno_use_tmp$combined_pve_FDR < 0.05)
sum(gene_EA_anno_use_tmp$outchr_pve_FDR < 0.05)

> sum(gene_EA_anno_use_tmp$combined_pve_FDR < 0.1)
[1] 2116
> sum(gene_EA_anno_use_tmp$outchr_pve_FDR < 0.1)
[1] 2121
> 
> sum(gene_EA_anno_use_tmp$combined_pve_FDR < 0.05)
[1] 1440
> sum(gene_EA_anno_use_tmp$outchr_pve_FDR < 0.05)
[1] 1446
> 


save(gene_EA_anno_use_tmp,file = "/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_use_tmp.RData")

```

We used Benjamini-Hochberg false discovery rates (FDR) to correct for multiple testing 59, with an FDR <0.1 used as the threshold for significant PVE. In the analysis, we found that 15.61% of genes in AA and 12.16% of genes in EA have a PVE that significantly deviates from zero at FDR < 0.1. In AA, the median PVE is 22.00% across these significant genes (2,750 genes, mean=24.06%; sd=10.14%), with PVE estimates ranging from 7.07% to 79.59% (Figure S4A). In EA, the median PVE is 23.96% across these significant genes (2,111 genes, mean=25.47%; sd=9.37%), with PVE estimates ranging from 8.63% to 79.25% (Figure S4B). 

```
> dim(gene_AA_anno_order)
[1] 17616    31
> sum(gene_AA_anno_use_tmp$combined_pve_FDR < 0.1)
[1] 2715
> sum(gene_AA_anno_use_tmp$outchr_pve_FDR < 0.1)
[1] 2750

2715/17616
2750/17616
> 2715/17616
[1] 0.1541213
> 2750/17616
[1] 0.1561081

> summary(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$outchr_pve_FDR<0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.07067 0.17522 0.21999 0.24063 0.27806 0.79587 
> summary(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$outchr_pve_FDR>=0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01504 0.04122 0.06030 0.07152 0.09351 0.22169 
sd(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$outchr_pve_FDR<0.1])
[1] 0.1014059


sd(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$combined_pve_FDR<0.1])
sd(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$combined_pve_FDR>=0.1])
> sd(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$combined_pve_FDR<0.1])
[1] 0.101844
> sd(gene_AA_anno_use_tmp$PVE_bslmm_outchr[gene_AA_anno_use_tmp$combined_pve_FDR>=0.1])
[1] 0.0392845

> summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR<0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.07358 0.17643 0.22020 0.24208 0.28136 0.79455 
> summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR>=0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01408 0.04143 0.06049 0.07189 0.09407 0.22491 

sd(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR<0.1])
sd(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR>=0.1])
> sd(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR<0.1])
[1] 0.1021632
> sd(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR>=0.1])
[1] 0.03929911


dim(gene_EA_anno_order)
sum(gene_EA_anno_use_tmp$combined_pve_FDR < 0.1)
sum(gene_EA_anno_use_tmp$outchr_pve_FDR < 0.1)
> dim(gene_EA_anno_order)
[1] 17360    31
> sum(gene_EA_anno_use_tmp$combined_pve_FDR < 0.1)
[1] 2103
> sum(gene_EA_anno_use_tmp$outchr_pve_FDR < 0.1)
[1] 2111

2103/17360
2111/17360
> 2103/17360
[1] 0.1211406
> 2111/17360
[1] 0.1216014

> summary(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR<0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.08626 0.19628 0.23955 0.25473 0.29281 0.79250 
> summary(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR>=0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01347 0.04608 0.06699 0.07923 0.10255 0.25213 
sd(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR<0.1])
> sd(gene_EA_anno_use_tmp$PVE_bslmm_outchr[gene_EA_anno_use_tmp$outchr_pve_FDR<0.1])
[1] 0.09365229

> summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR<0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.08254 0.19688 0.23990 0.25549 0.29336 0.78970 
> summary(gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$combined_pve_FDR>=0.1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01387 0.04630 0.06703 0.07932 0.10269 0.25026 


```

The PVE of tested common genes is generally consistent between AA and EA (Pearson’s correlation = 0.567, p-value < 2.2e-16), and the PVE of common eGenes is also consistent between AA and EA (Pearson’s correlation = 0.506, p-value < 2.2e-16) (Figure S5A-S5B). 
```
load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_use.RData")
load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_use.RData")

gene_common_use = merge(gene_EA_anno_use, gene_AA_anno_use, by="GENE")
> dim(gene_common_use_tmp)
[1] 17230    83

gene_common_use_tmp = gene_common_use
> dim(gene_common_use_tmp)
[1] 17230    83


> cor.test(gene_common_use_tmp$PVE_bslmm_outchr.x, gene_common_use_tmp$PVE_bslmm_outchr.y, method="pearson")$p.value
[1] 0

> a=cor.test(gene_common_use_tmp$PVE_bslmm_outchr.x, gene_common_use_tmp$PVE_bslmm_outchr.y, method="pearson")
> a

	Pearson's product-moment correlation

data:  gene_common_use_tmp$PVE_bslmm_outchr.x and gene_common_use_tmp$PVE_bslmm_outchr.y
t = 90.445, df = 17228, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.5571972 0.5774476
sample estimates:
      cor 
0.5674082 

> a$p.value
[1] 0

2*pt(-abs(a$statistic),df=a$parameter)
2*pt(a$statistic, a$parameter, lower=FALSE)
> 2*pt(-abs(a$statistic),df=a$parameter)
t 
0 
> 2*pt(a$statistic, a$parameter, lower=FALSE)
t 
0 

eGene_common_use_tmp = gene_common_use_tmp[gene_common_use_tmp$eGene.x+gene_common_use_tmp$eGene.y==2,]
> dim(eGene_common_use_tmp)
[1] 2997   83

> cor.test(eGene_common_use_tmp$PVE_bslmm_outchr.x, eGene_common_use_tmp$PVE_bslmm_outchr.y, method="pearson")$p.value
[1] 3.281339e-194

> cor.test(eGene_common_use_tmp$PVE_bslmm_outchr.x, eGene_common_use_tmp$PVE_bslmm_outchr.y, method="pearson")

	Pearson's product-moment correlation

data:  eGene_common_use_tmp$PVE_bslmm_outchr.x and eGene_common_use_tmp$PVE_bslmm_outchr.y
t = 32.066, df = 2995, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.4783938 0.5317190
sample estimates:
     cor 
0.505539 


save(eGene_common_use_tmp, file = "/net/mulan/home/shanglu/GENOA/analysis/compare/eGene_common_use_tmp.RData")
save(gene_common_use_tmp, file = "/net/mulan/home/shanglu/GENOA/analysis/compare/gene_common_use_tmp.RData")


```

As one might expect, eGenes tend to have a higher PVE than non-eGenes (p<2e-16) (Figure 3A-3B): the median PVE is 14.32% across eGenes and 5.33% across non-eGenes in AA, and is 14.82% across eGenes and 6.08% across non-eGenes in EA. 
```

summary(gene_AA_anno_use_tmp[gene_AA_anno_use_tmp$eGene==1,]$PVE_bslmm_outchr)
> summary(gene_AA_anno_use_tmp[gene_AA_anno_use_tmp$eGene==1,]$PVE_bslmm_outchr)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01981 0.09989 0.14316 0.16884 0.21139 0.79587 
summary(gene_AA_anno_use_tmp[gene_AA_anno_use_tmp$eGene==0,]$PVE_bslmm_outchr)
> summary(gene_AA_anno_use_tmp[gene_AA_anno_use_tmp$eGene==0,]$PVE_bslmm_outchr)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01504 0.03856 0.05326 0.06645 0.07904 0.75943 



summary(gene_EA_anno_use_tmp[gene_EA_anno_use_tmp$eGene==1,]$PVE_bslmm_outchr)
> summary(gene_EA_anno_use_tmp[gene_EA_anno_use_tmp$eGene==1,]$PVE_bslmm_outchr)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.03029 0.10411 0.14842 0.17243 0.21386 0.79250 
summary(gene_EA_anno_use_tmp[gene_EA_anno_use_tmp$eGene==0,]$PVE_bslmm_outchr)
> summary(gene_EA_anno_use_tmp[gene_EA_anno_use_tmp$eGene==0,]$PVE_bslmm_outchr)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01347 0.04359 0.06081 0.07666 0.09169 0.48236 

```

With BSLMM, we partitioned the PVE of each gene into two parts: one that is explained by cis-SNPs and the other that is explained by trans-SNPs. Consistent with previous studies 37; 55, we found that the majority of PVE is explained by trans-SNPs, with only a fraction explained by cis-SNPs: the median proportion of PVE explained by cis-SNPs is only 0.83% (mean=2.39%; sd=4.82%) across all genes in AA, and is 0.86% (mean=2.04%; sd=4.22%) across all genes in EA. 


```
summary(gene_AA_anno_use$por_cis_outchr)
sd(gene_AA_anno_use$por_cis_outchr)
> summary(gene_AA_anno_use$por_cis_outchr)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00000 0.00560 0.01012 0.03081 0.03180 0.72155 
> sd(gene_AA_anno_use$por_cis_outchr)
[1] 0.05459647

summary(gene_EA_anno_use$por_cis_outchr)
sd(gene_EA_anno_use$por_cis_outchr)
> summary(gene_EA_anno_use$por_cis_outchr)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.000000 0.005951 0.009896 0.026227 0.024112 0.707872 
> sd(gene_EA_anno_use$por_cis_outchr)
[1] 0.0485165


```
As one might expect, cis-SNPs explain a higher proportion of PVE in eGenes than in non-eGenes. Specifically, the median proportion of PVE explained by cis-SNPs is 5.4% (mean=8.0%; sd= 7.81%) across eGenes in AA, and is 5.15% (mean= 7.69%; sd= 7.68%) across eGenes in EA. The proportion of PVE in eGenes explained by cis-SNPs is correlated between AA and EA (Pearson’s correlation coefficient = 0.62, p-value= 1.482197e-322 and the proportion of PVE in eGenes explained by trans-SNPs is also correlated between AA and EA (Pearson’s correlation coefficient = 0.17, p-value= 4.898539e-22). 

```
summary(gene_AA_anno_use$por_cis_outchr[gene_AA_anno_use$eGene==1])
sd(gene_AA_anno_use$por_cis_outchr[gene_AA_anno_use$eGene==1])
> summary(gene_AA_anno_use$por_cis_outchr[gene_AA_anno_use$eGene==1])
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.004466 0.034515 0.054049 0.080029 0.094053 0.721547 
> sd(gene_AA_anno_use$por_cis_outchr[gene_AA_anno_use$eGene==1])
[1] 0.07809626

summary(gene_EA_anno_use$por_cis_outchr[gene_EA_anno_use$eGene==1])
sd(gene_EA_anno_use$por_cis_outchr[gene_EA_anno_use$eGene==1])
> summary(gene_EA_anno_use$por_cis_outchr[gene_EA_anno_use$eGene==1])
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.002903 0.033938 0.051485 0.076944 0.087626 0.707872 
> sd(gene_EA_anno_use$por_cis_outchr[gene_EA_anno_use$eGene==1])
[1] 0.07678144
 

cor.test(eGene_common_use_tmp$por_cis_outchr.x, eGene_common_use_tmp$por_cis_outchr.y)$p.value
cor.test(eGene_common_use_tmp$por_trans_outchr.x, eGene_common_use_tmp$por_trans_outchr.y)$p.value

> cor.test(eGene_common_use_tmp$por_cis_outchr.x, eGene_common_use_tmp$por_cis_outchr.y)$p.value
[1] 1.482197e-322
> cor.test(eGene_common_use_tmp$por_trans_outchr.x, eGene_common_use_tmp$por_trans_outchr.y)$p.value
[1] 4.898539e-22


> cor.test(eGene_common_use_tmp$por_cis_outchr.x, eGene_common_use_tmp$por_cis_outchr.y)

	Pearson's product-moment correlation

data:  eGene_common_use_tmp$por_cis_outchr.x and eGene_common_use_tmp$por_cis_outchr.y
t = 43.649, df = 2995, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.6011577 0.6449466
sample estimates:
      cor 
0.6235409 

> cor.test(eGene_common_use_tmp$por_trans_outchr.x, eGene_common_use_tmp$por_trans_outchr.y)

	Pearson's product-moment correlation

data:  eGene_common_use_tmp$por_trans_outchr.x and eGene_common_use_tmp$por_trans_outchr.y
t = 9.7267, df = 2995, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.1400631 0.2094819
sample estimates:
    cor 
0.17499 





```



prepare for conditional analysis
```R
gene_anno_AA_eGene=gene_AA_anno_order[which(gene_AA_anno_order$eGene==1),]
save(gene_anno_AA_eGene, file = "/net/mulan/home/shanglu/GENOA/analysis/conditional/gene_anno_AA_eGene.RData")
gene_anno_EA_eGene=gene_EA_anno_order[which(gene_EA_anno_order$eGene==1),]
save(gene_anno_EA_eGene, file = "/net/mulan/home/shanglu/GENOA/analysis/conditional/gene_anno_EA_eGene.RData")

```





