
In your gene set enrichment analysis, you used all annotated genes as background, 
which means the pathways that appear to be enriched may simply reflect genes that are expressed in LCLs 
compared to those that aren't. 
I suggest using expressed genes as your background and then seeing which pathways may be enriched in eGenes. 

```

gene_AA_anno_order$eGene = NA
gene_AA_anno_order$eGene[gene_AA_anno_order$topSNP_pval<=5.471313e-05] = 1
gene_AA_anno_order$eGene[gene_AA_anno_order$topSNP_pval>5.471313e-05] = 0
gene_EA_anno_order$eGene = NA
gene_EA_anno_order$eGene[gene_EA_anno_order$topSNP_pval<=0.0001199845] = 1
gene_EA_anno_order$eGene[gene_EA_anno_order$topSNP_pval>0.0001199845] = 0

AA_all_genes = gene_AA_anno_order[,1]
EA_all_genes = gene_EA_anno_order[,1]
AA_EA_all_genes  = unique(c(AA_all_genes, EA_all_genes))
AA_eGenes = gene_AA_anno_order[gene_AA_anno_order$eGene==1,1]
AA_non_eGenes = gene_AA_anno_order[gene_AA_anno_order$eGene==0,1]
EA_eGenes = gene_EA_anno_order[gene_EA_anno_order$eGene==1,1]
EA_non_eGenes = gene_EA_anno_order[gene_EA_anno_order$eGene==0,1]
AA_EA_common_egenes = intersect(AA_eGenes, EA_eGenes)

write.table(AA_all_genes, "AA_all_genes.txt", quote=F, col.names=F, row.names=F)
write.table(EA_all_genes, "EA_all_genes.txt", quote=F, col.names=F, row.names=F)
write.table(AA_EA_all_genes, "AA_EA_all_genes.txt", quote=F, col.names=F, row.names=F)
write.table(AA_eGenes, "AA_eGenes.txt", quote=F, col.names=F, row.names=F)
write.table(AA_non_eGenes, "AA_non_eGenes.txt", quote=F, col.names=F, row.names=F)
write.table(EA_eGenes, "EA_eGenes.txt", quote=F, col.names=F, row.names=F)
write.table(EA_non_eGenes, "EA_non_eGenes.txt", quote=F, col.names=F, row.names=F)
write.table(AA_EA_common_egenes, "AA_EA_common_egenes.txt", quote=F, col.names=F, row.names=F)





```


