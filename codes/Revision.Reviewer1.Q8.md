
```
cd /net/mulan/home/shanglu/GENOA/analysis/compare


> AA_table[1,]
  chr        rs     ps n_miss allele1 allele0    af       beta         se
1   1 rs2286139 761732     14       T       C 0.233 0.02591264 0.05346187
     l_remle    p_wald            GENE
1 0.03708061 0.6279965 ENSG00000187634
> colnames(AA_table)
 [1] "chr"     "rs"      "ps"      "n_miss"  "allele1" "allele0" "af"     
 [8] "beta"    "se"      "l_remle" "p_wald"  "GENE"  
 
 
AA_res = AA_table[,c(1,12,2,3,5,6,7,8,9,11)]
write.table(AA_res, "AA_summary_statistics.txt",row.names=F,quote=F)

EA_res = EA_table[,c(1,12,2,3,5,6,7,8,9,11)]
write.table(EA_res, "EA_summary_statistics.txt",row.names=F,quote=F)


```
