
# Figure 1 A&B
```R
setwd("/net/mulan/home/shanglu/GENOA/analysis/figure/")
################################################################
# TSS TES plot
################################################################

### scale to 100kb in gene body
			
			get_xaxis = function(a,b,c,genelen,strand){
			pos = NULL
			if(strand == 1){
			dist = a - b
			}else if(strand == -1){
			dist = -(a - c) 
			}
			if(dist < 0 ){
			pos = dist/1000
			}else if(dist > 0 &  dist <= genelen){
			pos = dist/genelen*100
			}else if(dist > 0 &  dist > genelen){
			pos = 100 + (dist - genelen)/1000
			}
			return(pos)
			}

# AA		
    gene_AA_anno_eGene = gene_AA_anno_use[gene_AA_anno_use$eGene==1,]
		esnp_gene = data.frame(gene_AA_anno_eGene$topSNP_ps, gene_AA_anno_eGene$low, gene_AA_anno_eGene$up,gene_AA_anno_eGene$distance,gene_AA_anno_eGene$strand)
		e_dist_to_start = unlist(apply(esnp_gene, 1, function(x) get_xaxis(x[1],x[2],x[3],x[4],x[5])))+100
		pdf("TSS_AA_abline.pdf")
	    hist(e_dist_to_start,breaks=300, xaxt="n",ylab="Density", main="African American", prob=TRUE, xlab = "Position relative to gene (kb)",
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylim=c(0,0.03))
        axis(1, at=seq(0, 300, by=100), labels=c("-100kb", "TSS", "TES","+100kb"),cex.axis=1.5)
		abline(v=100, col="red")
		abline(v=200, col="red")
		dev.off()
		
     gene_AA_anno_eGene = gene_AA_anno_eGene[gene_AA_anno_eGene$eGene==1,]
		esnp_gene = data.frame(gene_AA_anno_eGene$topSNP_ps, gene_AA_anno_eGene$low, gene_AA_anno_eGene$up,gene_AA_anno_eGene$distance,gene_AA_anno_eGene$strand)
		e_dist_to_start = unlist(apply(esnp_gene, 1, function(x) get_xaxis(x[1],x[2],x[3],x[4],x[5])))+100
		pdf("TSS_AA_abline_cornflowerblue.pdf")
	    hist(e_dist_to_start,breaks=150, xaxt="n",ylab="Density", main="African American", prob=TRUE, xlab = "Position relative to gene (kb)",
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylim=c(0,0.03),col = "cornflowerblue")
        axis(1, at=seq(0, 300, by=100), labels=c("-100kb", "TSS", "TES","+100kb"),cex.axis=1.5)
		abline(v=100, col="red")
		abline(v=200, col="red")
		dev.off()
    
    
# EA		
    gene_EA_anno_eGene = gene_EA_anno_use[gene_EA_anno_use$eGene==1,]
		esnp_gene = data.frame(gene_EA_anno_eGene$topSNP_ps, gene_EA_anno_eGene$low, gene_EA_anno_eGene$up,gene_EA_anno_eGene$distance,gene_EA_anno_eGene$strand)
		e_dist_to_start = unlist(apply(esnp_gene, 1, function(x) get_xaxis(x[1],x[2],x[3],x[4],x[5])))+100
		pdf("TSS_EA_abline.pdf")
	    hist(e_dist_to_start,breaks=300, xaxt="n",ylab="Density", main="European American", prob=TRUE, xlab = "Position relative to gene (kb)",
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylim=c(0,0.03))
        axis(1, at=seq(0, 300, by=100), labels=c("-100kb", "TSS", "TES","+100kb"),cex.axis=1.5)
		abline(v=100, col="red")
		abline(v=200, col="red")
		dev.off()
		
    
    
    gene_EA_anno_eGene = gene_EA_anno_use[gene_EA_anno_use$eGene==1,]
		esnp_gene = data.frame(gene_EA_anno_eGene$topSNP_ps, gene_EA_anno_eGene$low, gene_EA_anno_eGene$up,gene_EA_anno_eGene$distance,gene_EA_anno_eGene$strand)
		e_dist_to_start = unlist(apply(esnp_gene, 1, function(x) get_xaxis(x[1],x[2],x[3],x[4],x[5])))+100
		pdf("TSS_EA_abline_cornflowerblue.pdf")
	    hist(e_dist_to_start,breaks=150, xaxt="n",ylab="Density", main="European American", prob=TRUE, xlab = "Position relative to gene (kb)",
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylim=c(0,0.03),col = "cornflowerblue")
        axis(1, at=seq(0, 300, by=100), labels=c("-100kb", "TSS", "TES","+100kb"),cex.axis=1.5)
		abline(v=100, col="red")
		abline(v=200, col="red")
		dev.off()
    




```


# Figure 1 C
```R
library(ggplot2)
 
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
common_table_eSNP = merge(AA_table_subset, EA_table_subset, by = c("GENE","rs"))
common_table_eSNP$rs = as.character(common_table_eSNP$rs)
common_table_eSNP$allele1.x = as.character(common_table_eSNP$allele1.x)
common_table_eSNP$allele0.x = as.character(common_table_eSNP$allele0.x)
common_table_eSNP$allele1.y = as.character(common_table_eSNP$allele1.y)
common_table_eSNP$allele0.y = as.character(common_table_eSNP$allele0.y)
ind = which(common_table_eSNP$allele1.x != common_table_eSNP$allele1.y)
common_table_eSNP$beta.x[ind] = -common_table_eSNP$beta.x[ind]


tiff(paste0("F1C_cornflowerblue.tiff"), units="in", width=6, height=6, res=150)
ggplot(common_table_eSNP, aes(beta.x, beta.y)) +
geom_point(size=0.5,color = "cornflowerblue")+
labs(x = "Effect size in African American", y = "Effect size in European American")+
geom_smooth(method='lm',formula=y~x,color="red")+
theme_bw(base_size = 22)+
ylim(-2,2)+
xlim(-2,2)
dev.off()
```

Figure 1D
```

> dim(AA_table_subset)
[1] 354931     12
> dim(common_table_eSNP)
[1] 112316     22
> dim(EA_table_subset)
[1] 371309     12
> gene_AAEA_eGene = merge(gene_AA_anno_eGene, gene_EA_anno_eGene,by="GENE")
> dim(gene_AAEA_eGene)
[1] 3048   83
> dim(gene_AA_anno_eGene)
[1] 5475   42
> dim(gene_EA_anno_eGene)
[1] 4402   42

```


# Figure 2: conservation score
```R
################################################################
# Conservation score violin
################################################################
load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_use.RData")
load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_use.RData")
load("/net/mulan/home/shanglu/GENOA/analysis/conditional/gene_anno_AA_eGene.RData")
load("/net/mulan/home/shanglu/GENOA/analysis/conditional/gene_anno_EA_eGene.RData")


gene_anno_common_eGene = merge(gene_anno_AA_eGene,gene_anno_EA_eGene,by = "GENE")
gene_anno_background = merge(gene_AA_anno_use,gene_EA_anno_use,by = "GENE",all = TRUE)

gene_anno_background$phylop = gene_anno_background$AA_phylop
gene_anno_background$phylop[which(is.na(gene_anno_background$AA_phylop))]= gene_anno_background$EA_phylop[which(is.na(gene_anno_background$AA_phylop))]

gene_anno_background$phastcon = gene_anno_background$AA_phastcon
gene_anno_background$phastcon[which(is.na(gene_anno_background$AA_phastcon))]= gene_anno_background$EA_phastcon[which(is.na(gene_anno_background$AA_phastcon))]

gene_anno_background$dNdS_ratio = gene_anno_background$dNdS_ratio.x
gene_anno_background$dNdS_ratio[which(is.na(gene_anno_background$dNdS_ratio.x))]= gene_anno_background$dNdS_ratio.y[which(is.na(gene_anno_background$dNdS_ratio.x))]


common_in_AA_index = which(gene_AA_anno_use$GENE %in% gene_anno_common_eGene$GENE)
gene_anno_AA_nocommon = gene_AA_anno_use[-common_in_AA_index,]
gene_anno_AA_unique = gene_anno_AA_nocommon[gene_anno_AA_nocommon$eGene==1,]

common_in_EA_index = which(gene_EA_anno_use$GENE %in% gene_anno_common_eGene$GENE)
gene_anno_EA_nocommon = gene_EA_anno_use[-common_in_EA_index,]
gene_anno_EA_unique = gene_anno_EA_nocommon[gene_anno_EA_nocommon$eGene==1,]


# phylop
Background = gene_anno_background$phylop
Common = gene_anno_common_eGene$AA_phylop
AA_unique = gene_anno_AA_unique$AA_phylop
EA_unique = gene_anno_EA_unique$EA_phylop
Phylop_score= c(Background,Common,EA_unique,AA_unique)
Class = factor(c(rep("Background",length(Background)),rep("Common",length(Common)),rep("EA unique",length(EA_unique)),rep("AA unique",length(AA_unique))),levels = c("Background","EA unique","AA unique","Common"),order=T)
dat = data.frame(Phylop_score,Class)

library(ggplot2)



# box
pdf("Fig2_phylop_box_update.pdf")
ggplot(dat, aes(x=Class, y=Phylop_score)) +
    geom_boxplot(alpha=0.4) +
    stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme(legend.position="none") +theme_bw(base_size = 22)+
    scale_fill_brewer(palette="Set3")+labs(title="Phylop score",x="", y = "Phylop score")
dev.off()


# phastcon
Background = gene_anno_background$phastcon
Common = gene_anno_common_eGene$AA_phastcon
AA_unique = gene_anno_AA_unique$AA_phastcon
EA_unique = gene_anno_EA_unique$EA_phastcon
phastcon_score= c(Background,Common,EA_unique,AA_unique)
Class = factor(c(rep("Background",length(Background)),rep("Common",length(Common)),rep("EA unique",length(EA_unique)),rep("AA unique",length(AA_unique))),levels = c("Background","EA unique","AA unique","Common"),order=T)
dat = data.frame(phastcon_score,Class)

pdf("Fig2_phastcon_box_update.pdf")
ggplot(dat, aes(x=Class, y=phastcon_score)) +
    geom_boxplot(alpha=0.4) +
    stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme(legend.position="none") +theme_bw(base_size = 22)+
    scale_fill_brewer(palette="Set3")+labs(title="Phastcon score",x="", y = "Phastcon score")
dev.off()


# dNdS
Background = na.omit(gene_anno_background$dNdS_ratio)
Common = na.omit(gene_anno_common_eGene$dNdS_ratio.x)
AA_unique = na.omit(gene_anno_AA_unique$dNdS_ratio)
EA_unique = na.omit(gene_anno_EA_unique$dNdS_ratio)
dNdS_ratio_score= c(Background,Common,EA_unique,AA_unique)
Class = factor(c(rep("Background",length(Background)),rep("Common",length(Common)),rep("EA unique",length(EA_unique)),rep("AA unique",length(AA_unique))),levels = c("Background","EA unique","AA unique","Common"),order=T)
dat = data.frame(dNdS_ratio_score,Class)


pdf("Fig2_dNdS_box_update.pdf")
ggplot(dat, aes(x=Class, y=dNdS_ratio_score)) +
    geom_boxplot(alpha=0.4) +
    stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme(legend.position="none") +theme_bw(base_size = 22)+
    scale_fill_brewer(palette="Set3")+labs(title="dN/dS score",x="", y = "dN/dS score")
dev.off()



# Fig 2D
library(qvalue)
library(data.table)
library(dplyr)

 pathAA="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
 load(paste0(pathAA,"/AA_table.RData"))
 pathEA="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"
 load(paste0(pathEA,"/EA_table.RData"))
 load("/net/mulan/home/shanglu/GENOA/analysis/compare/Fst_result.RData")

AA_thr = 6.245907e-05
EA_thr = 0.0001385504
AA_table_subset = AA_table[AA_table$p_wald <= AA_thr,]
EA_table_subset = EA_table[EA_table$p_wald <= EA_thr,]
AA_table$signif = (AA_table$p_wald<=AA_thr)*1
EA_table$signif = (EA_table$p_wald<=EA_thr)*1
common_signif = merge(AA_table_subset,EA_table_subset, by=c("GENE","rs") )
common_signif_rs = common_signif$rs
AA_specific = unique(AA_table_subset[-which(AA_table_subset$rs %in% common_signif_rs),]$rs)
EA_specific = unique(EA_table_subset[-which(EA_table_subset$rs %in% common_signif_rs),]$rs)
noeSNP = intersect(AA_table[AA_table$signif==0,]$rs,EA_table[EA_table$signif==0,]$rs)

Fst_result$TYPE = "non eSNPs"
Fst_result$TYPE[which(Fst_result$SNP  %in% AA_specific)] = "AA unique"
Fst_result$TYPE[which(Fst_result$SNP  %in% EA_specific)] = "EA unique"
Fst_result$TYPE[which(Fst_result$SNP  %in% common_signif_rs)] = "Common"

Fst_result_result = Fst_result %>%  group_by(TYPE) %>%
summarise(meanFst = mean(Fst), 
          medianFst = median(Fst),
          minFst = min(Fst),
          maxFst = max(Fst))
> table(Fst_result$TYPE)

AA unique    Common EA unique non eSNPs 
   197188    101285    216653  24391818 
      
> Fst_result_result
# A tibble: 4 x 5
  TYPE      meanFst medianFst minFst maxFst
  <chr>       <dbl>     <dbl>  <dbl>  <dbl>
1 AA unique  0.0895   0.0559       0  0.750
2 Common     0.0888   0.0472       0  0.654
3 EA unique  0.0790   0.0440       0  0.694
4 non eSNPs  0.0338   0.00902      0  0.807

Fst = Fst_result$Fst
Class = factor(Fst_result$TYPE,levels = c("non eSNPs","AA unique","EA unique","Common"),order=T)

dat = data.frame(Fst,Class)

library(ggplot2)

# change order 
Fst = Fst_result$Fst
Class = factor(Fst_result$TYPE,levels = c("non eSNPs","EA unique","AA unique","Common"),order=T)
dat = data.frame(Fst,Class)
library(ggplot2)
tiff(paste0("/net/mulan/home/shanglu/GENOA/analysis/figure/Fig2D_Fst_update_order.tiff"), units="in", width=5, height=5, res=150)
ggplot(dat, aes(x=Class, y=Fst)) +
    geom_boxplot(alpha=0.5,fill="white") +
    ylim(0,1)+
    stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme(legend.position="none") +theme_bw(base_size = 15)+
    scale_fill_brewer(palette="Set3")+labs(title="Fst",x="", y = "Fst")
dev.off()


```

# Figure 3: PVE
```R
library(ggplot2)

gene_AA_anno_use$por_cis_combined = gene_AA_anno_use$PVE_bslmm_combined * gene_AA_anno_use$PGE_bslmm_combined
gene_AA_anno_use$por_trans_combined = gene_AA_anno_use$PVE_bslmm_combined * (1 - gene_AA_anno_use$PGE_bslmm_combined)
cis = c(gene_AA_anno_use$por_cis_combined, gene_AA_anno_use$por_cis_combined[gene_AA_anno_use$eGene==0],gene_AA_anno_use$por_cis_combined[gene_AA_anno_use$eGene==1])
trans = c(gene_AA_anno_use$por_trans_combined, gene_AA_anno_use$por_trans_combined[gene_AA_anno_use$eGene==0],gene_AA_anno_use$por_trans_combined[gene_AA_anno_use$eGene==1])
additive = c(gene_AA_anno_use$PVE_bslmm_combined, gene_AA_anno_use$PVE_bslmm_combined[gene_AA_anno_use$eGene==0],gene_AA_anno_use$PVE_bslmm_combined[gene_AA_anno_use$eGene==1])
effect = c(cis, trans, additive)
effect_type = c(rep("cis-PVE",length(cis)), rep("trans-PVE",length(trans)), rep("total PVE",length(additive)))
effect_type <- factor(effect_type,levels = c('cis-PVE','trans-PVE','total PVE'),ordered = TRUE)
tmp = c(rep("All genes", length(gene_AA_anno_use$por_cis_combined)), rep("non cis-eQTL genes", sum(gene_AA_anno_use$eGene==0)),rep("cis-eQTL genes", sum(gene_AA_anno_use$eGene==1)))
effect_class = factor(rep(tmp, 3),levels = c('All genes','non cis-eQTL genes','cis-eQTL genes'),ordered = TRUE)
dat = data.frame(effect,effect_type,effect_class)


pdf(paste0("Fig3_PVE_combined_AA_update.pdf"),width=18, height=12)
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


gene_EA_anno_use$por_cis_combined = gene_EA_anno_use$PVE_bslmm_combined * gene_EA_anno_use$PGE_bslmm_combined
gene_EA_anno_use$por_trans_combined = gene_EA_anno_use$PVE_bslmm_combined * (1 - gene_EA_anno_use$PGE_bslmm_combined)
cis = c(gene_EA_anno_use$por_cis_combined, gene_EA_anno_use$por_cis_combined[gene_EA_anno_use$eGene==0],gene_EA_anno_use$por_cis_combined[gene_EA_anno_use$eGene==1])
trans = c(gene_EA_anno_use$por_trans_combined, gene_EA_anno_use$por_trans_combined[gene_EA_anno_use$eGene==0],gene_EA_anno_use$por_trans_combined[gene_EA_anno_use$eGene==1])
additive = c(gene_EA_anno_use$PVE_bslmm_combined, gene_EA_anno_use$PVE_bslmm_combined[gene_EA_anno_use$eGene==0],gene_EA_anno_use$PVE_bslmm_combined[gene_EA_anno_use$eGene==1])
effect = c(cis, trans, additive)
effect_type = c(rep("cis-PVE",length(cis)), rep("trans-PVE",length(trans)), rep("total PVE",length(additive)))
effect_type <- factor(effect_type,levels = c('cis-PVE','trans-PVE','total PVE'),ordered = TRUE)
tmp = c(rep("All genes", length(gene_EA_anno_use$por_cis_combined)), rep("non cis-eQTL genes", sum(gene_EA_anno_use$eGene==0)),rep("cis-eQTL genes", sum(gene_EA_anno_use$eGene==1)))
effect_class = factor(rep(tmp, 3),levels = c('All genes','non cis-eQTL genes','cis-eQTL genes'),ordered = TRUE)
dat = data.frame(effect,effect_type,effect_class)

library(ggplot2)
pdf(paste0("Fig3_PVE_combined_EA_update.pdf"),width=18, height=12)
ggplot(dat, aes(x = effect_class, y = effect,fill = effect_type)) +
scale_y_continuous(name = "Percent variance explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
scale_x_discrete(name = "") + 
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) + 
geom_boxplot( aes(x = effect_class, y = effect,fill = effect_type),width = 0.2,position=position_dodge(0.9)) +
ggtitle("European American") +
theme_bw(base_size = 30) +
theme(plot.title = element_text(size = 40,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()
```
# Figure 4
```R
#---------------
# 4A & 4C
#---------------

		library(ggplot2);library(reshape2)

		pdf("Fig4_hist_indep_eqtl_AA.pdf")
		num = c( table((gene_anno_AA_eGene$indep_snp)))
		xaxis = c(1:9)
		dat = data.frame(xaxis, num)
		dat$xaxis = as.factor(dat$xaxis)
		ggplot(dat, aes(x=xaxis,y=num)) + 
		geom_bar(alpha=0.8, fill = "cornflowerblue",stat="identity", position=position_dodge()) + 
		theme_bw(base_size = 22) + 
		labs(title="African American",x="Number of independent eQTLs", y = "Number of eGenes")
		dev.off()
		
		pdf("Fig4_hist_indep_eqtl_EA.pdf")
		num = c( table((gene_anno_EA_eGene$indep_snp)),0,0)
		xaxis = c(1:9)
		dat = data.frame(xaxis, num)
		dat$xaxis = as.factor(dat$xaxis)
		ggplot(dat, aes(x=xaxis,y=num)) + 
		geom_bar(alpha=0.8, fill = "cornflowerblue",stat="identity", position=position_dodge()) + 
		theme_bw(base_size = 22) + 
		labs(title="European American",x="Number of independent eQTLs", y = "Number of eGenes")
		dev.off()

load("/net/mulan/home/shanglu/GENOA/analysis/conditional/gene_anno_AA_eGene.RData")
load("/net/mulan/home/shanglu/GENOA/analysis/conditional/gene_anno_EA_eGene.RData")
		pdf("Fig4_hist_indep_eqtl_AA_update.pdf")
		num = c( table((gene_anno_AA_eGene$indep_snp)))
		xaxis = c(1:9)
		dat = data.frame(xaxis, num)
		dat$xaxis = as.factor(dat$xaxis)
		ggplot(dat, aes(x=xaxis,y=num)) + 
		geom_bar(alpha=0.8, fill = "cornflowerblue",stat="identity", position=position_dodge()) + 
		geom_text(aes(label=num), vjust=-0.8, color="black",
            position = position_dodge(0.9), size=6)+
		theme_bw(base_size = 22) + 
		ylim(0,3800)+
		labs(title="African American",x="Number of independent eQTLs", y = "Number of eGenes")
		dev.off()
		
		pdf("Fig4_hist_indep_eqtl_EA_update.pdf")
		num = c( table((gene_anno_EA_eGene$indep_snp)),0,0)
		xaxis = c(1:9)
		dat = data.frame(xaxis, num)
		dat$xaxis = as.factor(dat$xaxis)
		ggplot(dat, aes(x=xaxis,y=num)) + 
		geom_bar(alpha=0.8, fill = "cornflowerblue",stat="identity", position=position_dodge()) + 
		geom_text(aes(label=num), vjust=-0.8, color="black",
            position = position_dodge(0.9), size=6)+
		theme_bw(base_size = 22) + 
		ylim(0,3800)+
		labs(title="European American",x="Number of independent eQTLs", y = "Number of eGenes")
		dev.off()
		

	    
	    
	    
#---------------
# 4B & 4D
#---------------

# AA
pve = gene_anno_AA_eGene$PVE_bslmm_combined
cispve = gene_anno_AA_eGene$por_cis_combined
indsnp = as.factor(gene_anno_AA_eGene$indep_snp)
dat = data.frame(pve,cispve,indsnp)

pdf("Fig4_conditional_pve_AA.pdf")
ggplot(dat, aes(x=indsnp, y=pve)) +
    geom_boxplot(alpha=0.8, fill="cornflowerblue") +
    #stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme(legend.position="none") +
    theme_bw(base_size = 22)+
    ylim(-0.1,1)+
    scale_fill_brewer(palette="Set3")+
    labs(title="African American",x="Number of independent eQTLs", y = "PVE")
dev.off()



pdf("Fig4_conditional_cispve_AA.pdf")
ggplot(dat, aes(x=indsnp, y=cispve)) +
    geom_boxplot(alpha=0.8, fill="cornflowerblue") +
    #stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme(legend.position="none") +
    theme_bw(base_size = 22)+
    ylim(-0.1,1)+
    scale_fill_brewer(palette="Set3")+
    labs(title="African American",x="Number of independent eQTLs", y = "cis-PVE")
dev.off()

# EA
pve = c(gene_anno_EA_eGene$PVE_bslmm_combined,-2,-2)
cispve = c(gene_anno_EA_eGene$por_cis_combined,-2,-2)
indsnp = as.factor(c(gene_anno_EA_eGene$indep_snp,8,9))
dat = data.frame(pve,cispve,indsnp)

pdf("Fig4_conditional_pve_EA.pdf")
ggplot(dat, aes(x=indsnp, y=pve)) +
    geom_boxplot(alpha=0.8, fill="cornflowerblue") +
    #stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme(legend.position="none") +
    theme_bw(base_size = 22)+
    ylim(-0.1,1)+
    scale_fill_brewer(palette="Set3")+
    labs(title="European American",x="Number of independent eQTLs", y = "PVE")
dev.off()



pdf("Fig4_conditional_cispve_EA.pdf")
ggplot(dat, aes(x=indsnp, y=cispve)) +
    geom_boxplot(alpha=0.8, fill="cornflowerblue") +
    #stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme(legend.position="none") +
    theme_bw(base_size = 22)+
    ylim(-0.1,1)+
    scale_fill_brewer(palette="Set3")+
    labs(title="European American",x="Number of independent eQTLs", y = "cis-PVE")
dev.off()


# 4C & 4F



#---------
# density plot, dist to tss

# table: 
# col 1: eQTL
# col 2: E1, E2, E3, >=E4
# col 3: dist to tss
#---------

# in AA:

path_AAoutcond = "/net/mulan/home/shanglu/GENOA/analysis/conditional/output/AA"

eqtl_table_AA = data.frame()
eqtl_table_AA_rs = NULL
eqtl_table_AA_rs_label = NULL
eqtl_table_AA_rs_dist_tss = NULL

for(line in 1:dim(gene_AA_anno_eGene)[1]){
print(line)
i_index = gene_AA_anno_eGene$chr[line]
j_index = gene_AA_anno_eGene$index_j[line]
indep_table = read.table(paste0(path_AAoutcond,"/chr_",i_index,"_gene_",j_index,"/Indep_eQTL.bim"))
eqtl_table_AA_rs = c(eqtl_table_AA_rs, as.character(indep_table$V2))
eqtl_table_AA_rs_label = c(eqtl_table_AA_rs_label, c(1:dim(indep_table)[1]))
eqtl_table_AA_rs_dist_tss = c(eqtl_table_AA_rs_dist_tss, abs(indep_table$V4-gene_AA_anno_eGene$low[line]))
}
eqtl_table_AA = data.frame(eqtl_table_AA_rs,eqtl_table_AA_rs_label,eqtl_table_AA_rs_dist_tss)

eqtl_table_AA$labelnew = eqtl_table_AA$eqtl_table_AA_rs_label
eqtl_table_AA$labelnew = as.character(eqtl_table_AA$labelnew)
eqtl_table_AA$labelnew[eqtl_table_AA$eqtl_table_AA_rs_label>=4] = ">=4"
eqtl_table_AA$eqtl_table_AA_rs_dist_tss = eqtl_table_AA$eqtl_table_AA_rs_dist_tss/1000
eqtl_table_AA$eQTL_order = eqtl_table_AA$labelnew
library(dplyr)
mu <- eqtl_table_AA %>% group_by(eQTL_order) %>% summarize(grp.median = median(eqtl_table_AA_rs_dist_tss))

pdf("Fig4_density_eqtl_AA.pdf",width=10, height=8)
eqtl_table_AA2 = eqtl_table_AA[eqtl_table_AA$eqtl_table_AA_rs_dist_tss<=200,]
ggplot(eqtl_table_AA2, aes(x = eqtl_table_AA_rs_dist_tss)) + 
geom_density(aes(fill = eQTL_order), alpha = 0.4) + 
geom_vline(data = mu, aes(xintercept = grp.median, color = eQTL_order), linetype = "dashed") + 
labs(title="African American",x ="Distance to TSS(kb)")+
theme(legend.position="bottom")+
theme_bw(base_size = 22) 
dev.off()


# in EA:

path_EAoutcond = "/net/mulan/home/shanglu/GENOA/analysis/conditional/output/EA"

eqtl_table_EA = data.frame()
eqtl_table_EA_rs = NULL
eqtl_table_EA_rs_label = NULL
eqtl_table_EA_rs_dist_tss = NULL

for(line in 1:dim(gene_EA_anno_eGene)[1]){
print(line)
i_index = gene_EA_anno_eGene$chr[line]
j_index = gene_EA_anno_eGene$index_j[line]
indep_table = read.table(paste0(path_EAoutcond,"/chr_",i_index,"_gene_",j_index,"/Indep_eQTL.bim"))
eqtl_table_EA_rs = c(eqtl_table_EA_rs, as.character(indep_table$V2))
eqtl_table_EA_rs_label = c(eqtl_table_EA_rs_label, c(1:dim(indep_table)[1]))
eqtl_table_EA_rs_dist_tss = c(eqtl_table_EA_rs_dist_tss, abs(indep_table$V4-gene_EA_anno_eGene$low[line]))
}
eqtl_table_EA = data.frame(eqtl_table_EA_rs,eqtl_table_EA_rs_label,eqtl_table_EA_rs_dist_tss)

eqtl_table_EA$labelnew = eqtl_table_EA$eqtl_table_EA_rs_label
eqtl_table_EA$labelnew = as.character(eqtl_table_EA$labelnew)
eqtl_table_EA$labelnew[eqtl_table_EA$eqtl_table_EA_rs_label>=4] = ">=4"
eqtl_table_EA$eqtl_table_EA_rs_dist_tss = eqtl_table_EA$eqtl_table_EA_rs_dist_tss/1000
eqtl_table_EA$eQTL_order = eqtl_table_EA$labelnew
library(dplyr)
mu <- eqtl_table_EA %>% group_by(eQTL_order) %>% summarize(grp.median = median(eqtl_table_EA_rs_dist_tss))

pdf("Fig4_density_eqtl_EA.pdf",width=10, height=8)
eqtl_table_EA2 = eqtl_table_EA[eqtl_table_EA$eqtl_table_EA_rs_dist_tss<=200,]
ggplot(eqtl_table_EA2, aes(x = eqtl_table_EA_rs_dist_tss)) + 
geom_density(aes(fill = eQTL_order), alpha = 0.4) + 
geom_vline(data = mu, aes(xintercept = grp.median, color = eQTL_order), linetype = "dashed") + 
labs(title="European American",x ="Distance to TSS(kb)")+
theme(legend.position="bottom")+
theme_bw(base_size = 22) 
dev.off()



```

# Figure S1: PCs vs eGene number and eSNP number
```R
# AA

thr = c(6.312936E-05, 6.245907E-05,	6.262234E-05,	6.191039E-05,	6.228262E-05)
pheno = c("_pheno14","","_pheno26","_pheno38","_pheno50")

for(k in 1:5){
print(k)
numbers = c()
for(i in 1:22){
load(paste0("/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping/AA_full_table_chr_",i,pheno[k],".RData"))
numbers[i] = sum(AA_full_table$p_wald<=thr[k])
}
print(sum(numbers))
}

PC_eGene = c(5527, 5475,5432,5398,5376)
PC_eSNP = c(360856,354931,351626,348148,346208)


library(RColorBrewer)
COL = brewer.pal(9,"Set1");
PC = c(0,  5, 10, 15, 20)
nG = PC_eGene
nPC = length(PC)
pdf("FigS1_eGene_vs_PC_AA.pdf", 4, 3)
plot(0, 0, xlim=c(0, 20), ylim=c(0, 10000), type="n", xlab="Number of Genotype PCs", ylab="Number of eGenes", main="", cex.lab=1.2)
abline(v=seq(0, 20, 5), lty=2, col="lightgrey")
abline(h=seq(0, 10000, 1000), lty=2, col="lightgrey")
points(PC, nG, type="o", col=COL[1], pch=20, lwd=2)
#legend("bottomright", legend=c("eGenes vs PCs at 5% FDR in African American"), fill=c(COL[1], "grey"), bg="white", cex=1.2)
dev.off()

#  y: #eSNPs, x: #PCs
PC = c(0,  5, 10, 15, 20)
nG = PC_eSNP
nPC = length(PC)
pdf("eSNP_vs_PC_AA.pdf", 4,3)
plot(0, 0, xlim=c(0, 20), ylim=c(0, 600000), type="n", xlab="Number of Genotype PCs", ylab="Number of eSNPs", main="",cex.lab=1.2)
abline(v=seq(0, 20, 5), lty=2, col="lightgrey")
abline(h=seq(0, 600000, 50000), lty=2, col="lightgrey")
points(PC, nG, type="o", col=COL[2], pch=20, lwd=2)
#legend("bottomright", legend=c("eSNPs vs PCs at 5% FDR in African American"), fill=c(COL[1], "grey"), bg="white", cex=1.2)
dev.off()


# EA
thr = c(0.0001436664,	0.0001385504	,0.0001374437,	0.0001358554	,0.0001345381)
pheno = c("_pheno14","","_pheno26","_pheno38","_pheno50")
for(k in 1:5){
print(k)
numbers = c()
for(i in 1:22){
load(paste0("/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping/EA_full_table_chr_",i,pheno[k],".RData"))
numbers[i] = sum(EA_full_table$p_wald<=thr[k])
}
print(sum(numbers))
}
PC_eGene = c(4481,	4402,	4365,	4333,	4281)
PC_eSNP = c(375708,371309,367525,363188,359364)
library(RColorBrewer)
COL = brewer.pal(9,"Set1");
PC = c(0,  5, 10, 15, 20)
nG = PC_eGene
nPC = length(PC)
pdf("FigS1_eGene_vs_PC_EA.pdf", 4, 3)
plot(0, 0, xlim=c(0, 20), ylim=c(0, 10000), type="n", xlab="Number of Genotype PCs", ylab="Number of eGenes", main="",cex.lab=1.2)
abline(v=seq(0, 20, 5), lty=2, col="lightgrey")
abline(h=seq(0, 10000, 1000), lty=2, col="lightgrey")
points(PC, nG, type="o", col=COL[1], pch=20, lwd=2)
#legend("bottomright", legend=c("eGenes vs PCs at 5% FDR in European American"), fill=c(COL[1], "grey"), bg="white", cex=1.2)
dev.off()


#  y: #eQTLs, x: #PCs
PC = c(0,  5, 10, 15, 20)
nG = PC_eSNP
nPC = length(PC)
pdf("FigS1eSNP_vs_PC_EA.pdf", 4, 3)
plot(0, 0, xlim=c(0, 20), ylim=c(0, 600000), type="n", xlab="Number of Genotype PCs", ylab="Number of eSNPs", main="",cex.lab=1.2)
abline(v=seq(0, 20, 5), lty=2, col="lightgrey")
abline(h=seq(0, 600000, 50000), lty=2, col="lightgrey")
points(PC, nG, type="o", col=COL[2], pch=20, lwd=2)
#legend("bottomright", legend=c("eQTLs vs PCs at 5% FDR in European American"), fill=c(COL[1], "grey"), bg="white", cex=1.2)
dev.off()


```


# Figure S2: MAF vs beta
```R
# AA
df = data.frame(gene_anno_AA_eGene$topSNP_af, abs(gene_anno_AA_eGene$topSNP_beta))
colnames(df) = c("MAF","Abs_effect_size")
pdf("FigS2_AA.pdf")
p <- ggplot(df, aes(x=MAF, y=Abs_effect_size))+  
labs(title="African American",
         x="Minor Allele Frequency",
         y="Absolute effect size (beta)")+
  geom_point(size=1,color = "cornflowerblue") + ylim(0,2.5)+
  geom_smooth(method="lm",size=1,colour="red")+
  theme_bw(base_size = 22)
p
dev.off() 

# EA
df = data.frame(gene_anno_EA_eGene$topSNP_af, abs(gene_anno_EA_eGene$topSNP_beta))
colnames(df) = c("MAF","Abs_effect_size")
pdf("FigS2_EA.pdf")
p <- ggplot(df, aes(x=MAF, y=Abs_effect_size))+  
labs(title="European American",
         x="Minor Allele Frequency",
         y="Absolute effect size (beta)")+
  geom_point(size=1,color = "cornflowerblue") + ylim(0,2.5)+
  geom_smooth(method="lm",size=1,colour="red")+
  theme_bw(base_size = 22)
p
dev.off() 
 


```

# Figure S3: Fst in AA and EA, eSNPs vs non-eSNPs
```R
pathAA="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
load(paste0(pathAA,"/AA_table.RData"))
pathEA="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"
load(paste0(pathEA,"/EA_table.RData"))

AA_table_subset = AA_table[AA_table$p_wald <= AA_thr,]
EA_table_subset = EA_table[EA_table$p_wald <= EA_thr,]

load("/net/mulan/home/shanglu/GENOA/analysis/compare/Fst_result.RData")
Fst_result$SNP = as.character(Fst_result$SNP)

AA_table$SNP = as.character(AA_table$rs)
AA_table_Fst = merge(AA_table, Fst_result, by="SNP")
AA_table_Fst$signif = "non-eSNP"
AA_table_Fst$signif[which(AA_table_Fst$p_wald <=AA_thr)]="eSNP"

EA_table$SNP = as.character(EA_table$rs)
EA_table_Fst = merge(EA_table, Fst_result, by="SNP")
EA_table_Fst$signif = "non-eSNP"
EA_table_Fst$signif[which(EA_table_Fst$p_wald <=EA_thr)]="eSNP"

fst = AA_table_Fst$Fst
class = AA_table_Fst$signif

dat = data.frame(fst, class)
p1 = ggplot(dat, aes(x=fst, fill=class)) + geom_density(alpha=0.7) + theme_minimal(base_size = 22)+
ggtitle("African American") +labs(x="Fst")


fst = EA_table_Fst$Fst
class = EA_table_Fst$signif

dat = data.frame(fst, class)
p3 = ggplot(dat, aes(x=fst, fill=class)) + geom_density(alpha=0.7) + theme_minimal(base_size = 22)+
ggtitle("European American") +labs(x="Fst")


pdf(paste0("FS3_Fst_AA.pdf"))
p1
dev.off()

pdf(paste0("FS3_Fst_EA.pdf"))
p3
dev.off()


# use common_table_all_Fst:
common_table_all_Fst$signif_AA = "non-eSNP"
common_table_all_Fst$signif_AA[which(common_table_all_Fst$p_wald.x <=AA_thr)]="eSNP"
common_table_all_Fst$signif_EA = "non-eSNP"
common_table_all_Fst$signif_EA[which(common_table_all_Fst$p_wald.y <=EA_thr)]="eSNP"

fst = common_table_all_Fst$Fst
class = common_table_all_Fst$signif_AA
dat = data.frame(fst, class)
p2 = ggplot(dat, aes(x=fst, fill=class)) + geom_density(alpha=0.5) + theme_minimal(base_size = 22)+
ggtitle("African American") +labs(x="Fst")

fst = common_table_all_Fst$Fst
class = common_table_all_Fst$signif_EA
dat = data.frame(fst, class)
p4 = ggplot(dat, aes(x=fst, fill=class)) + geom_density(alpha=0.5) + theme_minimal(base_size = 22)+
ggtitle("European American") +labs(x="Fst")


pdf(paste0("FS3_Fst_AA_commonpairs.pdf"))
p2
dev.off()

pdf(paste0("FS3_Fst_EA_commonpairs.pdf"))
p4
dev.off()

 ```        
	  


# Fig S4: significant PVE genes 
```R
#------------------------------
# only use significant genes
#------------------------------

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


gene_AA_anno_use_tmp = gene_AA_anno_use_tmp[gene_AA_anno_use_tmp$combined_pve_FDR<=0.1,]
cis = c(gene_AA_anno_use_tmp$por_cis_combined, gene_AA_anno_use_tmp$por_cis_combined[gene_AA_anno_use_tmp$eGene==0],gene_AA_anno_use_tmp$por_cis_combined[gene_AA_anno_use_tmp$eGene==1])
trans = c(gene_AA_anno_use_tmp$por_trans_combined, gene_AA_anno_use_tmp$por_trans_combined[gene_AA_anno_use_tmp$eGene==0],gene_AA_anno_use_tmp$por_trans_combined[gene_AA_anno_use_tmp$eGene==1])
additive = c(gene_AA_anno_use_tmp$PVE_bslmm_combined, gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$eGene==0],gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$eGene==1])
effect = c(cis, trans, additive)
effect_type = c(rep("cis-PVE",length(cis)), rep("trans-PVE",length(trans)), rep("total PVE",length(additive)))
effect_type <- factor(effect_type,levels = c('cis-PVE','trans-PVE','total PVE'),ordered = TRUE)
tmp = c(rep("All genes", length(gene_AA_anno_use_tmp$por_cis_combined)), rep("non cis-eQTL genes", sum(gene_AA_anno_use_tmp$eGene==0)),rep("cis-eQTL genes", sum(gene_AA_anno_use_tmp$eGene==1)))
effect_class = factor(rep(tmp, 3),levels = c('All genes','non cis-eQTL genes','cis-eQTL genes'),ordered = TRUE)
dat = data.frame(effect,effect_type,effect_class)
pdf(paste0("Fig3_PVE_combined_AA_signif_FDR0.1.pdf"),width=18, height=12)
ggplot(dat, aes(x = effect_class, y = effect,fill = effect_type)) +
scale_y_continuous(name = "Percent variance explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
scale_x_discrete(name = "") + 
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) + 
geom_boxplot( aes(x = effect_class, y = effect,fill = effect_type),width = 0.2,position=position_dodge(0.9)) +
ggtitle("African American (FDR<=0.1)") +
theme_bw(base_size = 30) +
theme(plot.title = element_text(size = 40,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()

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
gene_AA_anno_use_tmp = gene_AA_anno_use_tmp[gene_AA_anno_use_tmp$combined_pve_FDR<=0.05,]
cis = c(gene_AA_anno_use_tmp$por_cis_combined, gene_AA_anno_use_tmp$por_cis_combined[gene_AA_anno_use_tmp$eGene==0],gene_AA_anno_use_tmp$por_cis_combined[gene_AA_anno_use_tmp$eGene==1])
trans = c(gene_AA_anno_use_tmp$por_trans_combined, gene_AA_anno_use_tmp$por_trans_combined[gene_AA_anno_use_tmp$eGene==0],gene_AA_anno_use_tmp$por_trans_combined[gene_AA_anno_use_tmp$eGene==1])
additive = c(gene_AA_anno_use_tmp$PVE_bslmm_combined, gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$eGene==0],gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$eGene==1])
effect = c(cis, trans, additive)
effect_type = c(rep("cis-PVE",length(cis)), rep("trans-PVE",length(trans)), rep("total PVE",length(additive)))
effect_type <- factor(effect_type,levels = c('cis-PVE','trans-PVE','total PVE'),ordered = TRUE)
tmp = c(rep("All genes", length(gene_AA_anno_use_tmp$por_cis_combined)), rep("non cis-eQTL genes", sum(gene_AA_anno_use_tmp$eGene==0)),rep("cis-eQTL genes", sum(gene_AA_anno_use_tmp$eGene==1)))
effect_class = factor(rep(tmp, 3),levels = c('All genes','non cis-eQTL genes','cis-eQTL genes'),ordered = TRUE)
dat = data.frame(effect,effect_type,effect_class)
pdf(paste0("Fig3_PVE_combined_AA_signif_FDR0.05.pdf"),width=18, height=12)
ggplot(dat, aes(x = effect_class, y = effect,fill = effect_type)) +
scale_y_continuous(name = "Percent variance explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
scale_x_discrete(name = "") + 
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) + 
geom_boxplot( aes(x = effect_class, y = effect,fill = effect_type),width = 0.2,position=position_dodge(0.9)) +
ggtitle("African American (FDR<=0.05)") +
theme_bw(base_size = 30) +
theme(plot.title = element_text(size = 40,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()


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
> summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR<0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.08485 0.20108 0.24596 0.26676 0.30759 0.79457 
> summary(gene_AA_anno_use_tmp$PVE_bslmm_combined[gene_AA_anno_use_tmp$combined_pve_FDR>=0.05])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01418 0.04209 0.06257 0.07616 0.09977 0.27725 




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


gene_EA_anno_use_tmp = gene_EA_anno_use_tmp[gene_EA_anno_use_tmp$combined_pve_FDR<=0.1,]
cis = c(gene_EA_anno_use_tmp$por_cis_combined, gene_EA_anno_use_tmp$por_cis_combined[gene_EA_anno_use_tmp$eGene==0],gene_EA_anno_use_tmp$por_cis_combined[gene_EA_anno_use_tmp$eGene==1])
trans = c(gene_EA_anno_use_tmp$por_trans_combined, gene_EA_anno_use_tmp$por_trans_combined[gene_EA_anno_use_tmp$eGene==0],gene_EA_anno_use_tmp$por_trans_combined[gene_EA_anno_use_tmp$eGene==1])
additive = c(gene_EA_anno_use_tmp$PVE_bslmm_combined, gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$eGene==0],gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$eGene==1])
effect = c(cis, trans, additive)
effect_type = c(rep("cis-PVE",length(cis)), rep("trans-PVE",length(trans)), rep("total PVE",length(additive)))
effect_type <- factor(effect_type,levels = c('cis-PVE','trans-PVE','total PVE'),ordered = TRUE)
tmp = c(rep("All genes", length(gene_EA_anno_use_tmp$por_cis_combined)), rep("non eQTL genes", sum(gene_EA_anno_use_tmp$eGene==0)),rep("eQTL genes", sum(gene_EA_anno_use_tmp$eGene==1)))
effect_class = factor(rep(tmp, 3),levels = c('All genes','non eQTL genes','eQTL genes'),ordered = TRUE)
dat = data.frame(effect,effect_type,effect_class)
pdf(paste0("Fig3_PVE_combined_EA_signif_FDR0.1.pdf"),width=18, height=12)
ggplot(dat, aes(x = effect_class, y = effect,fill = effect_type)) +
scale_y_continuous(name = "Percent variance explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
scale_x_discrete(name = "") + 
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) + 
geom_boxplot( aes(x = effect_class, y = effect,fill = effect_type),width = 0.2,position=position_dodge(0.9)) +
ggtitle("European American (FDR<=0.1)") +
theme_bw(base_size = 30) +
theme(plot.title = element_text(size = 40,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()

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

gene_EA_anno_use_tmp = gene_EA_anno_use_tmp[gene_EA_anno_use_tmp$combined_pve_FDR<=0.05,]
cis = c(gene_EA_anno_use_tmp$por_cis_combined, gene_EA_anno_use_tmp$por_cis_combined[gene_EA_anno_use_tmp$eGene==0],gene_EA_anno_use_tmp$por_cis_combined[gene_EA_anno_use_tmp$eGene==1])
trans = c(gene_EA_anno_use_tmp$por_trans_combined, gene_EA_anno_use_tmp$por_trans_combined[gene_EA_anno_use_tmp$eGene==0],gene_EA_anno_use_tmp$por_trans_combined[gene_EA_anno_use_tmp$eGene==1])
additive = c(gene_EA_anno_use_tmp$PVE_bslmm_combined, gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$eGene==0],gene_EA_anno_use_tmp$PVE_bslmm_combined[gene_EA_anno_use_tmp$eGene==1])
effect = c(cis, trans, additive)
effect_type = c(rep("cis-PVE",length(cis)), rep("trans-PVE",length(trans)), rep("total PVE",length(additive)))
effect_type <- factor(effect_type,levels = c('cis-PVE','trans-PVE','total PVE'),ordered = TRUE)
tmp = c(rep("All genes", length(gene_EA_anno_use_tmp$por_cis_combined)), rep("non eQTL genes", sum(gene_EA_anno_use_tmp$eGene==0)),rep("eQTL genes", sum(gene_EA_anno_use_tmp$eGene==1)))
effect_class = factor(rep(tmp, 3),levels = c('All genes','non eQTL genes','eQTL genes'),ordered = TRUE)
dat = data.frame(effect,effect_type,effect_class)
pdf(paste0("Fig3_PVE_combined_EA_signif_FDR0.05.pdf"),width=18, height=12)
ggplot(dat, aes(x = effect_class, y = effect,fill = effect_type)) +
scale_y_continuous(name = "Percent variance explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
scale_x_discrete(name = "") + 
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) + 
geom_boxplot( aes(x = effect_class, y = effect,fill = effect_type),width = 0.2,position=position_dodge(0.9)) +
ggtitle("European American (FDR<=0.05)") +
theme_bw(base_size = 30) +
theme(plot.title = element_text(size = 40,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()

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

```


# Figure S5: PVE compare AA vs EA
```R

gene_anno_AAEA = merge(gene_AA_anno_use, gene_EA_anno_use, by = "GENE")
#> dim(gene_anno_AAEA)
#[1] 17220    83

#################################################
# all pve
#################################################

pdf("FigS5_PVE_AA_EA_compare.pdf")
ggplot(gene_anno_AAEA, aes(PVE_bslmm_combined.x, PVE_bslmm_combined.y)) +
geom_point(size=0.5,color="cornflowerblue")+
labs(x = "PVE in African American", y = "PVE in European American")+
geom_smooth(method='lm',size=1,color="red")+
theme_bw(base_size = 22)+
ylim(0,1)+
xlim(0,1)
dev.off()

pdf("FigS5_PVE_AA_EA_compare_loess.pdf")
ggplot(gene_anno_AAEA, aes(PVE_bslmm_combined.x, PVE_bslmm_combined.y)) +
geom_point(size=0.5,color="cornflowerblue")+
labs(x = "PVE in African American", y = "PVE in European American")+
geom_smooth(method='loess',size=1,color="red")+
theme_bw(base_size = 22)+
ylim(0,1)+
xlim(0,1)
dev.off()

######################
# only eGenes
######################

ind1 = which(gene_anno_AAEA$topSNP_pval.x <=AA_thr)
ind2 = which(gene_anno_AAEA$topSNP_pval.y <=EA_thr)
ind = intersect(ind1, ind2)
gene_anno_AAEA_eGene = gene_anno_AAEA[ind,]

#> dim(gene_anno_AAEA_eGene)
#[1] 3048   83

pdf("FigS5_PVE_AA_EA_compare_eGene.pdf")
ggplot(gene_anno_AAEA_eGene, aes(PVE_bslmm_combined.x, PVE_bslmm_combined.y)) +
geom_point(size=0.5,color="cornflowerblue")+
labs(x = "PVE in African American", y = "PVE in European American")+
geom_smooth(method='lm',size=1,color="red")+
theme_bw(base_size = 22)+
ylim(0,1)+
xlim(0,1)
dev.off()

pdf("FigS5_PVE_AA_EA_compare_eGene_loess.pdf")
ggplot(gene_anno_AAEA_eGene, aes(PVE_bslmm_combined.x, PVE_bslmm_combined.y)) +
geom_point(size=0.5,color="cornflowerblue")+
labs(x = "PVE in African American", y = "PVE in European American")+
geom_smooth(method='loess',size=1,color="red")+
theme_bw(base_size = 22)+
ylim(0,1)+
xlim(0,1)
dev.off()

#################################################
# cis-PVE
#################################################

pdf("FigS5_cisPVE_AA_EA_compare.pdf")
ggplot(gene_anno_AAEA, aes(por_cis_combined.x, por_cis_combined.y)) +
geom_point(size=0.5,color="cornflowerblue")+
labs(x = "cis-PVE in African American", y = "cis-PVE in European American")+
geom_smooth(method='lm',size=1,color="red")+
theme_bw(base_size = 22)+
ylim(0,1)+
xlim(0,1)
dev.off()

pdf("FigS5_cisPVE_AA_EA_compare_loess.pdf")
ggplot(gene_anno_AAEA, aes(por_cis_combined.x, por_cis_combined.y)) +
geom_point(size=0.5,color="cornflowerblue")+
labs(x = "cis-PVE in African American", y = "cis-PVE in European American")+
geom_smooth(method='loess',size=1,color="red")+
theme_bw(base_size = 22)+
ylim(0,1)+
xlim(0,1)
dev.off()

######################
# cis-PVE, only eGenes
######################

ind1 = which(gene_anno_AAEA$topSNP_pval.x <=AA_thr)
ind2 = which(gene_anno_AAEA$topSNP_pval.y <=EA_thr)
ind = intersect(ind1, ind2)
gene_anno_AAEA_eGene = gene_anno_AAEA[ind,]

#> dim(gene_anno_AAEA_eGene)
#[1] 3048   83

pdf("FigS5_cisPVE_AA_EA_compare_eGene.pdf")
ggplot(gene_anno_AAEA_eGene, aes(por_cis_combined.x, por_cis_combined.y)) +
geom_point(size=0.5,color="cornflowerblue")+
labs(x = "cis-PVE in African American", y = "cis-PVE in European American")+
geom_smooth(method='lm',size=1,color="red")+
theme_bw(base_size = 22)+
ylim(0,1)+
xlim(0,1)
dev.off()

pdf("FigS5_cisPVE_AA_EA_compare_eGene_loess.pdf")
ggplot(gene_anno_AAEA_eGene, aes(por_cis_combined.x, por_cis_combined.y)) +
geom_point(size=0.5,color="cornflowerblue")+
labs(x = "cis-PVE in African American", y = "cis-PVE in European American")+
geom_smooth(method='loess',size=1,color="red")+
theme_bw(base_size = 22)+
ylim(0,1)+
xlim(0,1)
dev.off()

```

# Fig S6: cis-PVE by indep eQTLs
```R

load("/net/mulan/home/shanglu/GENOA/analysis/conditional/gene_anno_AA_eGene.RData")
load("/net/mulan/home/shanglu/GENOA/analysis/conditional/gene_anno_EA_eGene.RData")


pve_ratio_primary_aa = gene_anno_AA_eGene$pve_primary/gene_anno_AA_eGene$por_cis_combined
pve_ratio_indep_aa = gene_anno_AA_eGene$pve_indep/gene_anno_AA_eGene$por_cis_combined
type_aa = c(rep("Primary eQTL", length(pve_ratio_primary_aa)), rep("All independent eQTLs", length(pve_ratio_indep_aa)))
pve_ratio_aa = c(pve_ratio_primary_aa, pve_ratio_indep_aa)

pve_ratio_primary_ea = gene_anno_EA_eGene$pve_primary/gene_anno_EA_eGene$por_cis_combined
pve_ratio_indep_ea = gene_anno_EA_eGene$pve_indep/gene_anno_EA_eGene$por_cis_combined
type_ea = c(rep("Primary eQTL", length(pve_ratio_primary_ea)), rep("All independent eQTLs", length(pve_ratio_indep_ea)))
pve_ratio_ea = c(pve_ratio_primary_ea, pve_ratio_indep_ea)

type = c(type_aa,type_ea)
type = factor(type,levels = c("Primary eQTL","All independent eQTLs"),ordered = TRUE)

pve_ratio = c(pve_ratio_aa, pve_ratio_ea)
class = c(rep("African Amercan", length(type_aa)),rep("European Amercan", length(type_ea)))
dat = data.frame(type,pve_ratio,class)


pdf(paste0("FigS6_pve_indep_AAEA.pdf"),width=10, height=10)
ggplot(dat, aes(x = type, y = pve_ratio,fill=class)) +
scale_y_continuous(name = "Proportion of cis-PVE explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) +
geom_boxplot( aes(x = type, y = pve_ratio,fill=class),width = 0.2,position=position_dodge(0.9)) +
theme_bw(base_size = 30) +
#theme(plot.title = element_text(size = 20,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
#scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()


dat$pve_ratio[dat$pve_ratio>1]=1
pdf(paste0("FigS6_pve_indep_AAEA_greaterthan1settobe1.pdf"),width=10, height=10)
ggplot(dat, aes(x = type, y = pve_ratio,fill=class)) +
scale_y_continuous(name = "Proportion of cis-PVE explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) +
geom_boxplot( aes(x = type, y = pve_ratio,fill=class),width = 0.2,position=position_dodge(0.9)) +
theme_bw(base_size = 30) +
#theme(plot.title = element_text(size = 20,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
#scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()

```

# Figure S7: 
The proportion of cis-SNP heritability for gene expression explained by identified primary eQTLs and indepdent eQTLs depends on the number of eQTLs identifed
```R
# AA

dat = gene_anno_AA_eGene
dat$indep_snp = as.factor(dat$indep_snp)
dat$por_cispve_by_primary = dat$pve_primary/dat$por_cis_combined
dat$por_cispve_by_indep = dat$pve_indep/dat$por_cis_combined

pdf(paste0("FigS7_cispve_AA_byprimary.pdf"),width=10, height=10)
ggplot(dat, aes(x = indep_snp, y = por_cispve_by_primary)) +
scale_y_continuous(name = "Proportion of cis-PVE explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE,fill="cornflowerblue",alpha=0.7) +
geom_boxplot( aes(x = indep_snp, y = por_cispve_by_primary),width = 0.2,position=position_dodge(0.9),fill="plum2",alpha=0.7) +
theme_bw(base_size = 30) +
#theme(plot.title = element_text(size = 20,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette="RdBu")+
theme(legend.position="bottom")+
labs(title="African American",x = "Number of independent eQTLs", y = "Proportion of cis-PVE explained")
dev.off()

pdf(paste0("FigS7_cispve_AA_byindep.pdf"),width=10, height=10)
ggplot(dat, aes(x = indep_snp, y = por_cispve_by_indep)) +
scale_y_continuous(name = "Proportion of cis-PVE explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE,fill="cornflowerblue",alpha=0.7) +
geom_boxplot( aes(x = indep_snp, y = por_cispve_by_indep),width = 0.2,position=position_dodge(0.9),fill="plum2",alpha=0.7) +
theme_bw(base_size = 30) +
#theme(plot.title = element_text(size = 20,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette="RdBu")+
theme(legend.position="bottom")+
labs(title="African American",x = "Number of independent eQTLs", y = "Proportion of cis-PVE explained")
dev.off()

# EA

dat = gene_anno_EA_eGene
dat$indep_snp = as.factor(dat$indep_snp)
dat$por_cispve_by_primary = dat$pve_primary/dat$por_cis_combined
dat$por_cispve_by_indep = dat$pve_indep/dat$por_cis_combined

pdf(paste0("FigS7_cispve_EA_byprimary.pdf"),width=10, height=10)
ggplot(dat, aes(x = indep_snp, y = por_cispve_by_primary)) +
scale_y_continuous(name = "Proportion of cis-PVE explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE,fill="cornflowerblue",alpha=0.7) +
geom_boxplot( aes(x = indep_snp, y = por_cispve_by_primary),width = 0.2,position=position_dodge(0.9),fill="plum2",alpha=0.7) +
theme_bw(base_size = 30) +
#theme(plot.title = element_text(size = 20,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette="RdBu")+
theme(legend.position="bottom")+
labs(title="European American",x = "Number of independent eQTLs", y = "Proportion of cis-PVE explained")
dev.off()

pdf(paste0("FigS7_cispve_EA_byindep.pdf"),width=10, height=10)
ggplot(dat, aes(x = indep_snp, y = por_cispve_by_indep)) +
scale_y_continuous(name = "Proportion of cis-PVE explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE,fill="cornflowerblue",alpha=0.7) +
geom_boxplot( aes(x = indep_snp, y = por_cispve_by_indep),width = 0.2,position=position_dodge(0.9),fill="plum2",alpha=0.7) +
theme_bw(base_size = 30) +
#theme(plot.title = element_text(size = 20,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette="RdBu")+
theme(legend.position="bottom")+
labs(title="European American",x = "Number of independent eQTLs", y = "Proportion of cis-PVE explained")
dev.off()


```

# Figure S8, number of eQTL vs gene length
```R
cor.test(gene_anno_AA_eGene$indep_snp, gene_anno_AA_eGene$distance)$estimate
cor.test(gene_anno_AA_eGene$indep_snp, gene_anno_AA_eGene$distance)$p.value
cor.test(gene_anno_EA_eGene$indep_snp, gene_anno_EA_eGene$distance)$estimate
cor.test(gene_anno_EA_eGene$indep_snp, gene_anno_EA_eGene$distance)$p.value

> cor.test(gene_anno_AA_eGene$indep_snp, gene_anno_AA_eGene$distance)$estimate
        cor 
0.001300337 
> cor.test(gene_anno_AA_eGene$indep_snp, gene_anno_AA_eGene$distance)$p.value
[1] 0.9233663
> cor.test(gene_anno_EA_eGene$indep_snp, gene_anno_EA_eGene$distance)$estimate
       cor 
0.05053222 
> cor.test(gene_anno_EA_eGene$indep_snp, gene_anno_EA_eGene$distance)$p.value
[1] 0.0007968679


ind = which(gene_anno_EA_eGene$indep_snp>1)
cor.test(gene_anno_EA_eGene$indep_snp[ind], gene_anno_EA_eGene$distance[ind])$estimate
cor.test(gene_anno_EA_eGene$indep_snp[ind], gene_anno_EA_eGene$distance[ind])$p.value

ind = which(gene_anno_AA_eGene$indep_snp>1)
cor.test(gene_anno_AA_eGene$indep_snp[ind], gene_anno_AA_eGene$distance[ind])$estimate
cor.test(gene_anno_AA_eGene$indep_snp[ind], gene_anno_AA_eGene$distance[ind])$p.value

> ind = which(gene_anno_EA_eGene$indep_snp>1)
> cor.test(gene_anno_EA_eGene$indep_snp[ind], gene_anno_EA_eGene$distance[ind])$estimate
        cor 
0.003970856 
> cor.test(gene_anno_EA_eGene$indep_snp[ind], gene_anno_EA_eGene$distance[ind])$p.value
[1] 0.9093316
> 
> ind = which(gene_anno_AA_eGene$indep_snp>1)
> cor.test(gene_anno_AA_eGene$indep_snp[ind], gene_anno_AA_eGene$distance[ind])$estimate
        cor 
0.005035502 
> cor.test(gene_anno_AA_eGene$indep_snp[ind], gene_anno_AA_eGene$distance[ind])$p.value
[1] 0.8332767


pdf("FigS8_genelength_num_eQTL_AA.pdf")
dp <- ggplot(dat, aes(x=indep_snp, y=distance)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="Number of independent eQTLs", y = "Gene length")
dp + scale_fill_brewer(palette="RdBu") + theme_bw(base_size = 22)
dev.off()

dat$log10gene = log10(dat$distance+1)
pdf("FigS8_genelength_num_eQTL_AA_log10.pdf")
dp <- ggplot(dat, aes(x=indep_snp, y=log10gene)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="Number of independent eQTLs", y = "log10(Gene length+1)")
dp + scale_fill_brewer(palette="RdBu") + theme_bw(base_size = 22)
dev.off()

dat = gene_anno_EA_eGene
dat$indep_snp = as.factor(dat$indep_snp)

pdf("FigS8_genelength_num_eQTL_EA.pdf")
dp <- ggplot(dat, aes(x=indep_snp, y=distance)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="Number of independent eQTLs", y = "Gene length")
dp + scale_fill_brewer(palette="RdBu") + theme_bw(base_size = 22)
dev.off()

dat$log10gene = log10(dat$distance+1)
pdf("FigS8_genelength_num_eQTL_EA_log10.pdf")
dp <- ggplot(dat, aes(x=indep_snp, y=log10gene)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="Number of independent eQTLs", y = "log10(Gene length+1)")
dp + scale_fill_brewer(palette="RdBu") + theme_bw(base_size = 22)
dev.off()

###################
# 3 groups
###################

dat = gene_anno_AA_eGene
dat$eqtltype = ">=3"
dat$eqtltype[dat$indep_snp==1] = "1"
dat$eqtltype[dat$indep_snp==2] = "2"
dat$eqtltype = factor(dat$eqtltype, levels = c("1","2",">=3"),order=T)
dat$log10gene = log10(dat$distance)

pdf("FigS8_genelength_num_eQTL_AA_log10_3groups.pdf")
ggplot(dat, aes(x=eqtltype, y=log10gene)) + 
  geom_boxplot(alpha=0.8, fill="cornflowerblue")+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  labs(title="African American",x="Number of independent eQTLs", y = "Gene length(log10)")+
  theme(legend.position="none") +theme_bw(base_size = 22)+
  scale_fill_brewer(palette="Set3") + theme_bw(base_size = 22)
dev.off()

dat = gene_anno_EA_eGene
dat$eqtltype = ">=3"
dat$eqtltype[dat$indep_snp==1] = "1"
dat$eqtltype[dat$indep_snp==2] = "2"
dat$eqtltype = factor(dat$eqtltype, levels = c("1","2",">=3"),order=T)
dat$log10gene = log10(dat$distance)

pdf("FigS8_genelength_num_eQTL_EA_log10_3groups.pdf")
ggplot(dat, aes(x=eqtltype, y=log10gene)) + 
  geom_boxplot(alpha=0.8, fill="cornflowerblue")+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  labs(title="European American",x="Number of independent eQTLs", y = "Gene length(log10)")+
  theme(legend.position="none") +theme_bw(base_size = 22)+
  scale_fill_brewer(palette="Set3") + theme_bw(base_size = 22)
dev.off()

    

```

# Figure S9
Scatter plots show the prediction performance measured by R2 (y-axis) versus the SNP heritability of gene expression (x-axis). 

```R
path_pred = "/net/mulan/home/shanglu/GENOA/analysis/prediction"
library(dplyr)
library(ggplot2)
library(patchwork)




# remove predicted genes with negative correlation with the truth
#------------------
# for AA
#------------------
load(paste0(path_pred,"/result/result_","AA","_bslmm_allpairs.RData"))
rho = result$pearson_rho_v2
p = list()
k=0

k = k + 1
rho_pop = rho[,k][which(rho[,k]>=0)]
rho_pop_square = rho_pop^2
gene_pop = names(rho_pop_square)
pve_pop = gene_AA_anno_order$PVE_bslmm_combined[match(gene_pop,gene_AA_anno_order$GENE)]	
df = data.frame(rho_pop_square, pve_pop)
p[[k]] <- ggplot(df, aes(x=pve_pop, y=rho_pop_square))+  
labs(title=paste0("AA->",colnames(rho)[k]),x="PVE",y="R2")+
  geom_point(size=1,color = "cornflowerblue") + 
  ylim(0,1)+
  xlim(0,1)+
  geom_smooth(method="lm",size=1,colour="red")+
  theme_bw(base_size = 22)
	
pdf(paste0("FigS9_AA.pdf"), width=25,height=5)
(p[[1]] | p[[2]] |p[[3]]|p[[4]]|p[[5]]) 
dev.off()

#------------------
# for EA
#------------------
load(paste0(path_pred,"/result/result_","EA","_bslmm_allpairs.RData"))
rho = result$pearson_rho_v2
p = list()
k=0

k = k + 1
rho_pop = rho[,k][which(rho[,k]>=0)]
rho_pop_square = rho_pop^2
gene_pop = names(rho_pop_square)
pve_pop = gene_EA_anno_order$PVE_bslmm_combined[match(gene_pop,gene_EA_anno_order$GENE)]	
df = data.frame(rho_pop_square, pve_pop)
p[[k]] <- ggplot(df, aes(x=pve_pop, y=rho_pop_square))+  
labs(title=paste0("EA->",colnames(rho)[k]),x="PVE",y="R2")+
  geom_point(size=1,color = "cornflowerblue") + 
  ylim(0,1)+
  xlim(0,1)+
  geom_smooth(method="lm",size=1,colour="red")+
  theme_bw(base_size = 22)
	
pdf(paste0("FigS9_EA.pdf"), width=25,height=5)
(p[[1]] | p[[2]] |p[[3]]|p[[4]]|p[[5]]) 
dev.off()



load(paste0(path_pred,"/result/result_","EA","_bslmm_allpairs.RData"))
summary(result$pearson_rho_v2)



```






