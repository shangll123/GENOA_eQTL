
```
cd /home/skardia_lab/clubhouse/research/projects/GENOA_AA_GE/preprocessing_zxu_March2017/Lulu_Mar2018/AJHG_revision/GOfuncR
/garage/skardia_lab/R-3.6.1/bin/R

 load("gene_anno_AA_eGene.RData")
 load("gene_anno_EA_eGene.RData")
 load("gene_AA_anno_order_correct.RData")
 load("gene_EA_anno_order_correct.RData")

# BiocManager::install('Homo.sapiens')
library(GOfuncR)


AA_egene = as.character(gene_anno_AA_eGene$gene)
AA_bg = as.character(gene_AA_anno_order$gene[-which(gene_AA_anno_order$gene %in% AA_egene)])

EA_egene = as.character(gene_anno_EA_eGene$gene)
EA_bg = as.character(gene_EA_anno_order$gene[-which(gene_EA_anno_order$gene %in% EA_egene)])

AAEA_common_egene = intersect(AA_egene, EA_egene)
AAEA_background = unique(c(AA_bg, EA_bg))

/*
> dim(gene_anno_AA_eGene)
[1] 5475   47
> dim(gene_anno_EA_eGene)
[1] 4402   47
> length(AA_bg)
[1] 12141
> length(EA_bg)
[1] 12957
> length(AAEA_background)
[1] 14619
> length(AAEA_common_egene)
[1] 3048
> 


*/

# AA eGene vs AA background
candi_gene_ids = AA_egene
bg_gene_ids = AA_bg
is_candidate = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
input_hyper_bg = data.frame(gene_ids = c(candi_gene_ids, bg_gene_ids),is_candidate)
head(input_hyper_bg)
res_hyper_eGene = go_enrich(input_hyper_bg, n_randsets=100,gene_len=TRUE)
length(res_hyper_eGene)
names(res_hyper_eGene)


result_eGene_AA = res_hyper_eGene[[1]][res_hyper_eGene[[1]]$FWER_overrep<0.1,]
write.csv(result_eGene_AA, "result_eGene_AA.csv",quote=F)




# AA noneGene vs AA background
candi_gene_ids = AA_bg
bg_gene_ids = AA_egene 
is_candidate = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
input_hyper_bg = data.frame(gene_ids = c(candi_gene_ids, bg_gene_ids),is_candidate)
head(input_hyper_bg)
res_hyper_noneGene = go_enrich(input_hyper_bg, n_randsets=100,gene_len=TRUE)
length(res_hyper_noneGene)
names(res_hyper_noneGene)
head(res_hyper_noneGene[[1]])
refined_noneGene = refine(res_hyper_noneGene, fwer=0.1)


result_noneGene_AA = res_hyper_noneGene[[1]][res_hyper_noneGene[[1]]$FWER_overrep<0.1,]
write.csv(result_noneGene_AA, "result_noneGene_AA.csv",quote=F)



    
 
# EA eGene vs EA background
candi_gene_ids = EA_egene
bg_gene_ids = EA_bg
is_candidate = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
input_hyper_bg = data.frame(gene_ids = c(candi_gene_ids, bg_gene_ids),is_candidate)
head(input_hyper_bg)
res_hyper_eGene = go_enrich(input_hyper_bg, n_randsets=100,gene_len=TRUE)
length(res_hyper_eGene)
names(res_hyper_eGene)
head(res_hyper_eGene[[1]])
refined_eGene = refine(res_hyper_eGene, fwer=0.1)

result_eGene_EA = res_hyper_eGene[[1]][res_hyper_eGene[[1]]$FWER_overrep<0.05,]
write.csv(result_eGene_EA, "result_eGene_EA.csv",quote=F)




# EA noneGene vs EA background
candi_gene_ids = EA_bg
bg_gene_ids = EA_egene 
is_candidate = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
input_hyper_bg = data.frame(gene_ids = c(candi_gene_ids, bg_gene_ids),is_candidate)
head(input_hyper_bg)
res_hyper_noneGene = go_enrich(input_hyper_bg, n_randsets=100,gene_len=TRUE)
length(res_hyper_noneGene)
names(res_hyper_noneGene)
head(res_hyper_noneGene[[1]])
refined_noneGene = refine(res_hyper_noneGene, fwer=0.1)


result_noneGene_EA = res_hyper_noneGene[[1]][res_hyper_noneGene[[1]]$FWER_overrep<0.1,]
write.csv(result_noneGene_EA, "result_noneGene_EA.csv",quote=F)

   
    
    
```

```
> result_eGene_AA
            ontology    node_id                                node_name
1 cellular_component GO:0005622                            intracellular
2 cellular_component GO:0043231 intracellular membrane-bounded organelle
3 cellular_component GO:0043229                  intracellular organelle
4 cellular_component GO:0043226                                organelle
5 cellular_component GO:0043227               membrane-bounded organelle
6 cellular_component GO:0005737                                cytoplasm
7 molecular_function GO:0003824                       catalytic activity
8 cellular_component GO:0005739                            mitochondrion
  raw_p_underrep raw_p_overrep FWER_underrep FWER_overrep
1              1  2.076929e-61             1         0.00
2              1  1.221447e-46             1         0.00
3              1  3.988981e-45             1         0.00
4              1  1.955034e-44             1         0.00
5              1  1.870330e-42             1         0.00
6              1  8.572787e-36             1         0.00
7              1  3.051325e-30             1         0.00
8              1  5.646730e-34             1         0.01

> result_eGene_EA
            ontology    node_id          node_name raw_p_underrep raw_p_overrep
1 cellular_component GO:0005622      intracellular              1  1.807875e-32
2 molecular_function GO:0003824 catalytic activity              1  4.970776e-24
  FWER_underrep FWER_overrep
1             1         0.01
2             1         0.01

> result_noneGene_AA
            ontology    node_id                                    node_name
1 biological_process GO:0032501             multicellular organismal process
2 biological_process GO:0007186 G protein-coupled receptor signaling pathway
  raw_p_underrep raw_p_overrep FWER_underrep FWER_overrep
1              1  1.304832e-42             1         0.00
2              1  8.544452e-31             1         0.01


> result_noneGene_EA
[1] ontology       node_id        node_name      raw_p_underrep raw_p_overrep 
[6] FWER_underrep  FWER_overrep  
<0 rows> (or 0-length row.names)

> summary(res_hyper_noneGene[[1]]$FWER_overrep)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.45    1.00    1.00    1.00    1.00    1.00 
> length(res_hyper_noneGene[[1]]$FWER_overrep)
[1] 22030




```


# g:Profiler
```
AAnoneGene = gene_AA_anno_order$GENE[which(gene_AA_anno_order$eGene==0)]
AAeGene = gene_AA_anno_order$GENE[which(gene_AA_anno_order$eGene==1)]
EAnoneGene = gene_EA_anno_order$GENE[which(gene_EA_anno_order$eGene==0)]
EAeGene = gene_EA_anno_order$GENE[which(gene_EA_anno_order$eGene==1)]


write.table(EAeGene, file = "EAeGene.txt", col.names=F, row.names=F, quote=F)
write.table(AAeGene, file = "AAeGene.txt", col.names=F, row.names=F, quote=F)
write.table(AAnoneGene, file = "AAnoneGene.txt", col.names=F, row.names=F, quote=F)
write.table(EAnoneGene, file = "EAnoneGene.txt", col.names=F, row.names=F, quote=F)






```


