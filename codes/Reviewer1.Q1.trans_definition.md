
You chose a cis-window of 100kb up- and down-stream of the TSS and TES, respectively. 
Other studies, including MESA and GEUVADIS used a 1Mb window. 
I think it's fine to use a smaller window, as most eQTLs are near the gene (previous work and Fig 1A-B), 
but this smaller window may inflate the trans component in your heritability partitioning analysis (Figs 3, S4). 
I suggest limiting trans-SNPs to those on different chromosomes, rather than defining trans as those SNPs that are not cis, 
to exclude any potential long-distance cis effects. 
Also, and perhaps more importantly, checking microarray probes for cross-mappability and filtering out bad probes will help 
avoid inflated trans heritability in your partitioning (see Saha and Battle, F1000 2018). 

```
# 1. compare PVE estimate in two versions, scatter plot

pathAAeqtl = "/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
pathEAeqtl = "/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"
pathcompare = "/net/mulan/home/shanglu/GENOA/analysis/compare"

load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_order_correct.RData")
load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order_correct.RData")


gene_AA_anno_order$por_cis_combined = gene_AA_anno_order$PVE_bslmm_combined * gene_AA_anno_order$PGE_bslmm_combined
por_cis_combined = gene_AA_anno_order$PVE_bslmm_combined * gene_AA_anno_order$PGE_bslmm_combined
gene_AA_anno_order$por_trans_combined = gene_AA_anno_order$PVE_bslmm_combined * (1 - gene_AA_anno_order$PGE_bslmm_combined)
por_trans_combined = gene_AA_anno_order$PVE_bslmm_combined * (1 - gene_AA_anno_order$PGE_bslmm_combined)
gene_AA_anno_order$por_cis_outchr = gene_AA_anno_order$PVE_bslmm_outchr * gene_AA_anno_order$PGE_bslmm_outchr
por_cis_outchr = gene_AA_anno_order$PVE_bslmm_outchr * gene_AA_anno_order$PGE_bslmm_outchr
gene_AA_anno_order$por_trans_outchr = gene_AA_anno_order$PVE_bslmm_outchr * (1 - gene_AA_anno_order$PGE_bslmm_outchr)
por_trans_outchr = gene_AA_anno_order$PVE_bslmm_outchr * (1 - gene_AA_anno_order$PGE_bslmm_outchr)

mydata = data.frame(por_cis_combined, por_cis_outchr,por_trans_combined,por_trans_outchr)

pdf(paste0("PVE_compare_AA_cis.pdf"),width=8,height=8)
ggplot(mydata, aes(x=por_cis_combined, y=por_cis_outchr)) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  # geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=25)+
  geom_abline(intercept = 0)+
  ylim(0,1)+
  xlim(0,1)+
  labs(title=paste0("African Americans"),
       x="cis-PVE w/ trans SNPs on the same chr", y = "cis-PVE w/o trans SNPs on the same chr")
dev.off()


pdf(paste0("PVE_compare_AA_trans.pdf"),width=8,height=8)
ggplot(mydata, aes(x=por_trans_combined, y=por_trans_outchr)) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  # geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=25)+
  geom_abline(intercept = 0)+
  ylim(0,1)+
  xlim(0,1)+
  labs(title=paste0("African Americans"),
       x="trans-PVE w/ trans SNPs on the same chr", y = "trans-PVE w/o trans SNPs on the same chr")
dev.off()

/*
cor.test(por_cis_combined, por_cis_outchr)$estimate
cor.test(por_cis_combined, por_cis_outchr)$p.value
cor.test(por_trans_combined, por_trans_outchr)$estimate
cor.test(por_trans_combined, por_trans_outchr)$p.value
> cor.test(por_cis_combined, por_cis_outchr)$estimate
      cor 
0.9988295 
> cor.test(por_cis_combined, por_cis_outchr)$p.value
[1] 0
> cor.test(por_trans_combined, por_trans_outchr)$estimate
      cor 
0.9937301 
> cor.test(por_trans_combined, por_trans_outchr)$p.value
[1] 0

*/

# EA

gene_EA_anno_order$por_cis_combined = gene_EA_anno_order$PVE_bslmm_combined * gene_EA_anno_order$PGE_bslmm_combined
por_cis_combined = gene_EA_anno_order$PVE_bslmm_combined * gene_EA_anno_order$PGE_bslmm_combined
gene_EA_anno_order$por_trans_combined = gene_EA_anno_order$PVE_bslmm_combined * (1 - gene_EA_anno_order$PGE_bslmm_combined)
por_trans_combined = gene_EA_anno_order$PVE_bslmm_combined * (1 - gene_EA_anno_order$PGE_bslmm_combined)
gene_EA_anno_order$por_cis_outchr = gene_EA_anno_order$PVE_bslmm_outchr * gene_EA_anno_order$PGE_bslmm_outchr
por_cis_outchr = gene_EA_anno_order$PVE_bslmm_outchr * gene_EA_anno_order$PGE_bslmm_outchr
gene_EA_anno_order$por_trans_outchr = gene_EA_anno_order$PVE_bslmm_outchr * (1 - gene_EA_anno_order$PGE_bslmm_outchr)
por_trans_outchr = gene_EA_anno_order$PVE_bslmm_outchr * (1 - gene_EA_anno_order$PGE_bslmm_outchr)

mydata = data.frame(por_cis_combined, por_cis_outchr,por_trans_combined,por_trans_outchr)

pdf(paste0("PVE_compare_EA_cis.pdf"),width=8,height=8)
ggplot(mydata, aes(x=por_cis_combined, y=por_cis_outchr)) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  # geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=25)+
  geom_abline(intercept = 0)+
  ylim(0,1)+
  xlim(0,1)+
  labs(title=paste0("European Americans"),
       x="cis-PVE w/ trans SNPs on the same chr", y = "cis-PVE w/o trans SNPs on the same chr")
dev.off()


pdf(paste0("PVE_compare_EA_trans.pdf"),width=8,height=8)
ggplot(mydata, aes(x=por_trans_combined, y=por_trans_outchr)) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  # geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=25)+
  geom_abline(intercept = 0)+
  ylim(0,1)+
  xlim(0,1)+
  labs(title=paste0("European Americans"),
       x="trans-PVE w/ trans SNPs on the same chr", y = "trans-PVE w/o trans SNPs on the same chr")
dev.off()

/*
cor.test(por_cis_combined, por_cis_outchr)$estimate
cor.test(por_cis_combined, por_cis_outchr)$p.value
cor.test(por_trans_combined, por_trans_outchr)$estimate
cor.test(por_trans_combined, por_trans_outchr)$p.value
> cor.test(por_cis_combined, por_cis_outchr)$estimate
      cor 
0.9988572 
> cor.test(por_cis_combined, por_cis_outchr)$p.value
[1] 0
> cor.test(por_trans_combined, por_trans_outchr)$estimate
      cor 
0.9933351 
> cor.test(por_trans_combined, por_trans_outchr)$p.value
[1] 0
> 


*/



```
