
# use Franc to calculayte local ancestry
tried different softwares


```
# download genetic_map_chr files
for ((i=1;i<=22;i++)); do
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140404_data_for_phase3_paper/shapeit2_scaffolds/genetic_map_chr${i}_combined_b37.20140701.txt
cp genetic_map_chr${i}_combined_b37.20140701.txt genetic_map_chr${i}.txt
done
```

```
# extract rsID for analyzed SNPs

> pathAA="/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping"
> pathEA="/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping"
> 
> load(paste0(pathEA,"/EA_table.RData"))
> load(paste0(pathAA,"/AA_table.RData"))

for(i in 1:22){
	print(i)
	rslist = unique(as.character(AA_table$rs)[which(AA_table$chr==i)])
	write.table(rslist, paste0("/net/mulan/home/shanglu/GENOA/data/AA/franc/SNPlist/AA_rslist_chr",i,".txt"),col.names=F,row.names=F,quote=F)
}



for(i in 1:22){
	print(i)
	rslist = unique(as.character(EA_table$rs)[which(EA_table$chr==i)])
	write.table(rslist, paste0("/net/mulan/home/shanglu/GENOA/data/EA/franc/SNPlist/EA_rslist_chr",i,".txt"),col.names=F,row.names=F,quote=F)
}


patheqtlAA=/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping
pathdataAA=/net/mulan/home/shanglu/GENOA/data/AA
pathSNPlist=/net/mulan/home/shanglu/GENOA/data/AA/franc/SNPlist

for ((i=1;i<=22;i++)); do
echo ${i}
Gfilename=/net/mulan/home/shanglu/GENOA/data/AA/Gene_file/chr_${i}
Gfile=/net/mulan/home/shanglu/GENOA/data/AA/VCFtoBED/chr_${i}
plink --noweb --bfile ${Gfile}  --extract ${pathSNPlist}/AA_rslist_chr${i}.txt --keep ${pathdataAA}/PIN_keep1032.txt --make-bed --out GENOA_AA.${i}
done

patheqtlEA=/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping
pathdataEA=/net/mulan/home/shanglu/GENOA/data/EA
pathSNPlist=/net/mulan/home/shanglu/GENOA/data/EA/franc/SNPlist

for ((i=1;i<=22;i++)); do
echo ${i}
Gfilename=/net/mulan/home/shanglu/GENOA/data/EA/Gene_file/chr_${i}
Gfile=/net/mulan/home/shanglu/GENOA/data/EA/VCFtoBED/chr_${i}
plink --noweb --bfile ${Gfile}  --extract ${pathSNPlist}/EA_rslist_chr${i}.txt --keep ${pathdataEA}/PIN_keep801_EU.txt --make-bed --out GENOA_EA.${i}
done


# separate by chromosome
cd /net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data/HapMap
tar -xvf check_ethnicity_HapMap_reference_populations_b37.tar.bz2
cd /net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data/1000geno
tar -xvf check_ethnicity_1000G_reference_populations_b37.tar.bz2

for chr in {1..22}; do 
plink --noweb --bfile 1000G_YRI_founders --chr $chr --make-bed --out YRI.${chr}
done

for chr in {1..22}; do 
plink --noweb --bfile 1000G_CEU_founders --chr $chr --make-bed --out CEU.${chr}
done



```


```

# Error in rbind(GENOA_AA_testT, CEUT, YRIT) : 
#  number of columns of matrices must match (see arg 2)
# yri.BIM adn CEU.bim must have same set of SNPs
# first need to find common set of SNPs in three panels
=====================
# AA
=====================

pathSNPlist="/net/mulan/home/shanglu/GENOA/data/AA/franc/SNPlist"
pathRef="/net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data/1000geno/check_ethnicity_1000G_ref_pops"
for(i in 1:22){
	print(i)
rs_AA = as.character(read.table(paste0(pathSNPlist,"/AA_rslist_chr",i,".txt"))$V1)
rs_CEU = as.character(read.table(paste0(pathRef,"/CEU.",i,".bim"))$V2)
rs_YRI = as.character(read.table(paste0(pathRef,"/YRI.",i,".bim"))$V2)
rs_common_AA_CEU = intersect(rs_AA,rs_CEU)
rs_common_AA_CEU_YRI = intersect(rs_common_AA_CEU,rs_YRI)
print(length(rs_common_AA_CEU_YRI))
write.table(rs_common_AA_CEU_YRI, paste0("/net/mulan/home/shanglu/GENOA/data/AA/franc/SNPlist/rs_common_AA_CEU_YRI_chr",i,".txt"),col.names=F,row.names=F,quote=F)
}


for i in {1..22}; do 
echo ${i}
pathSNPlist=/net/mulan/home/shanglu/GENOA/data/AA/franc/SNPlist
Gfile=/net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data/1000geno/check_ethnicity_1000G_ref_pops/CEU.${i}
plink --noweb --bfile ${Gfile}  --extract ${pathSNPlist}/rs_common_AA_CEU_YRI_chr${i}.txt  --make-bed --out CEU.${i}
done


for i in {1..22}; do 
echo ${i}
pathSNPlist=/net/mulan/home/shanglu/GENOA/data/AA/franc/SNPlist
Gfile=/net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data/1000geno/check_ethnicity_1000G_ref_pops/YRI.${i}
plink --noweb --bfile ${Gfile}  --extract ${pathSNPlist}/rs_common_AA_CEU_YRI_chr${i}.txt  --make-bed --out YRI.${i}
done

for i in {1..22}; do 
echo ${i}
pathSNPlist=/net/mulan/home/shanglu/GENOA/data/AA/franc/SNPlist
Gfile=/net/mulan/home/shanglu/GENOA/data/AA/franc/GENOA_AA.${i}
plink --noweb --bfile ${Gfile}  --extract ${pathSNPlist}/rs_common_AA_CEU_YRI_chr${i}.txt  --make-bed --out GENOA_AA_use.${i}
done


cd /net/mulan/home/shanglu/GENOA/src/FRANC/franc/franc_interface
python franc.py -p AApar.txt -d /net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data -a GENOA_AA_use -t rfmix -f rfmix -o GENOA_AA_use -m local
# worked, but looks rfmix will take a long time, try other softwares instead of rfmix

chr=1
cd /net/mulan/home/shanglu/GENOA/src/FRANC/franc/franc_interface
python myfranc.py -p AApar${chr}.txt -d /net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data -a GENOA_AA_use -t elai -f rfmix -o GENOA_AA_use -m local

# elai is faster


=====================
# EA
=====================

pathSNPlist="/net/mulan/home/shanglu/GENOA/data/EA/franc/SNPlist"
pathRef="/net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data/1000geno/check_ethnicity_1000G_ref_pops"
for(i in 1:22){
	print(i)
rs_EA = as.character(read.table(paste0(pathSNPlist,"/EA_rslist_chr",i,".txt"))$V1)
rs_CEU = as.character(read.table(paste0(pathRef,"/CEU.",i,".bim"))$V2)
rs_YRI = as.character(read.table(paste0(pathRef,"/YRI.",i,".bim"))$V2)
rs_common_EA_CEU = intersect(rs_EA,rs_CEU)
rs_common_EA_CEU_YRI = intersect(rs_common_EA_CEU,rs_YRI)
print(length(rs_common_EA_CEU_YRI))
write.table(rs_common_EA_CEU_YRI, paste0("/net/mulan/home/shanglu/GENOA/data/EA/franc/SNPlist/rs_common_EA_CEU_YRI_chr",i,".txt"),col.names=F,row.names=F,quote=F)
}


for i in {1..22}; do 
echo ${i}
pathSNPlist=/net/mulan/home/shanglu/GENOA/data/EA/franc/SNPlist
Gfile=/net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data/1000geno/check_ethnicity_1000G_ref_pops/CEU.${i}
plink --noweb --bfile ${Gfile}  --extract ${pathSNPlist}/rs_common_EA_CEU_YRI_chr${i}.txt  --make-bed --out CEU.${i}
done


for i in {1..22}; do 
echo ${i}
pathSNPlist=/net/mulan/home/shanglu/GENOA/data/EA/franc/SNPlist
Gfile=/net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data/1000geno/check_ethnicity_1000G_ref_pops/YRI.${i}
plink --noweb --bfile ${Gfile}  --extract ${pathSNPlist}/rs_common_EA_CEU_YRI_chr${i}.txt  --make-bed --out YRI.${i}
done

for i in {1..22}; do 
echo ${i}
pathSNPlist=/net/mulan/home/shanglu/GENOA/data/EA/franc/SNPlist
Gfile=/net/mulan/home/shanglu/GENOA/data/EA/franc/GENOA_EA.${i}
plink --noweb --bfile ${Gfile}  --extract ${pathSNPlist}/rs_common_EA_CEU_YRI_chr${i}.txt  --make-bed --out GENOA_EA_use.${i}
done


cd /net/mulan/home/shanglu/GENOA/src/FRANC/franc/franc_interface
python franc.py -p EApar.txt -d /net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data -a GENOA_EA_use -t rfmix -f rfmix -o GENOA_EA_use -m local
# nice


# make EA backups
cd /net/mulan/home/shanglu/GENOA/src/FRANC/franc/franc_interface
chr=5
python myfranc.py -p EApar${chr}.txt -d /net/mulan/home/shanglu/GENOA/src/FRANC/franc/EA_data -a GENOA_EA_use -t elai -f rfmix -o GENOA_EA_use -m local



```




## rfmix output read example
```R
===============================
@ read rfmix output

# devtools::install_github("keurcien/simulate") 
rfmix.local.ancestry(rfmix.output)

rfmix.local.ancestry <- function(rfmix.output){
  output <- read.table(rfmix.output, stringsAsFactors = F)
  numAdm <- ncol(output) / 2
  numSnps <- nrow(output)
  mat <- matrix(0, numSnps, numAdm)
  for(adm in 1:numAdm){
    mat[, adm] = paste0(output[, 2 * adm - 1], output[, 2 * adm])
  }
  return(mat)
}


display.ancestry <- function(ancestry.matrix){
    im <- matrix(0, nrow(ancestry.matrix), ncol(ancestry.matrix))
    im[ancestry.matrix == "22"] <- 2
    im[ancestry.matrix == "11"] <- 0
    im[ancestry.matrix == "12"] <- 1
    im[ancestry.matrix == "21"] <- 1
    return(im)
    #image(im, axes = FALSE, xlab = "", ylab = "")
}

rfmixout="GENOA_AA_test.19.0.Viterbi.txt"
res = rfmix.local.ancestry(rfmixout)
disp = display.ancestry(res)

```


## LAMatrix example
```R
===============================
LAMatrix
===============================

library(devtools)
install_github("yizhenzhong/Local_ancestry/LAMatrix")

https://github.com/yizhenzhong/Local_ancestry/blob/0bcf8beecb3ec753830ac01a52f0b23a4df701ad/LAMatrix.Rcheck/LAMatrix-Ex.R

n = 200;# Number of columns (samples)
nc = 10;# Number ofs covariates

cvrt.mat = 2 + matrix(rnorm(n*nc), ncol = nc);# Generate the covariates
snps.mat = floor(runif(n, min = 0, max = 3)); # Generate the genotype
local.mat = floor(runif(n, min = 0, max = 3)); # Generate the local ancestry

# Generate the expression vector
gene.mat = cvrt.mat %*% rnorm(nc) + rnorm(n) + 0.5 * snps.mat + 1 + 0.3 *local.mat;


#Create SlicedData objects
snps = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
genes = SlicedData$new( matrix( gene.mat, nrow = 1 ) );
cvrts = SlicedData$new( t(cvrt.mat) );
locals = SlicedData$new( matrix(local.mat,nrow=1))

# Produce no output files
filename = NULL; # tempfile()

modelLOCAL = 930507L;
me = LAMatrix_main(
 snps = snps,
 gene = genes,
 cvrt = cvrts,
 local = locals,
 output_file_name = filename,
 pvOutputThreshold = 1,
 useModel = modelLOCAL,
 verbose = TRUE,
 pvalue.hist = FALSE);


# Pull results - t-statistic and p-value
beta = me$all$eqtls$beta;
tstat = me$all$eqtls$statistic;
pvalue = me$all$eqtls$pvalue;
rez = c(beta = beta, tstat = tstat, pvalue = pvalue);
print(rez)

# Results from linear
lmdl = lm( gene.mat ~ snps.mat + cvrt.mat + local.mat);
lmout = summary(lmdl)$coefficients[2,c("Estimate","t value","Pr(>|t|)")];
print( lmout );

```

## AA prepare to run LAMatrix
```R


pathAA="/nfs/sph-mulan/home/shanglu/GENOA/data/AA"
load(paste0(pathAA,"/meta_AA.RData"))

load(paste0(pathAA,"/expression/expr_AA_progene_combat_qq.RData"))

load(paste0(pathAA,"/expression/expr_AA_proteingene.RData"))

colnames(expr_AA_progene_combat_qq) = rownames(expr_AA_proteingene)
rownames(expr_AA_progene_combat_qq) = colnames(expr_AA_proteingene)


X=expr_AA_progene_combat_qq

PCx_res=matrix(0,dim(X)[1],dim(X)[2])
PCx_res_qq=matrix(0,dim(X)[1],dim(X)[2])
for(j in 1:dim(X)[2]){
print(j)
PCx_res[,j]=lm(X[,j] ~ factor(gender) + age , data = meta)$residual
PCx_res_qq[,j]=qqnorm(PCx_res[,j],plot.it=F)$x
}

save(PCx_res_qq, file = "/nfs/sph-mulan/home/shanglu/GENOA/data/AA/LA/PCx_res_qq.RData")

load(paste0("/nfs/sph-mulan/home/shanglu/GENOA/data/AA/permuteID_AA.RData"))

for(permu in 1:10){
	print(permu)
	PCx_res_qq_permute = PCx_res_qq[new_id[[permu]],]
	colnames(PCx_res_qq_permute) = rownames(expr_AA_proteingene)
	save(PCx_res_qq_permute, file = paste0("/nfs/sph-mulan/home/shanglu/GENOA/data/AA/LA/PCx_res_qq_pheno_",permu+1,".RData"))
	}
	
# rename original expression
PCx_res_qq_permute = PCx_res_qq
colnames(PCx_res_qq_permute) = rownames(expr_AA_proteingene)
save(PCx_res_qq_permute, file = "/nfs/sph-mulan/home/shanglu/GENOA/data/AA/LA/PCx_res_qq_pheno_1.RData")




for(chr in 1:22){
print(chr)
path_rfmix = paste0("/net/mulan/home/shanglu/GENOA/src/FRANC/franc/test_data/GENOA_AA_use/elai/",chr)
rfmixout=paste0(path_rfmix,"/GENOA_AA_use.",chr,".elairfmix.txt")
res = rfmix.local.ancestry(rfmixout)
print("res done")
disp = display.ancestry(res)
save(disp, file = paste0(path_rfmix,"/disp_chr_",chr,".elairfmix.RData"))
print("disp done")
}

```

## EA prepare to run LAMatrix
```R


pathEA="/net/mulan/home/shanglu/GENOA/data/EA"
load(paste0(pathEA,"/EA_meta_new_order.RData"))
meta=meta_new_order
load(paste0(pathEA,"/expression/expr_EA_progene_combat_qq.RData"))

load(paste0(pathEA,"/expression/expr_EA_proteingene.RData"))

colnames(expr_EA_progene_combat_qq) = rownames(expr_EA_proteingene)
rownames(expr_EA_progene_combat_qq) = colnames(expr_EA_proteingene)


X=expr_EA_progene_combat_qq

PCx_res=matrix(0,dim(X)[1],dim(X)[2])
PCx_res_qq=matrix(0,dim(X)[1],dim(X)[2])
for(j in 1:dim(X)[2]){
print(j)
PCx_res[,j]=lm(X[,j] ~ factor(gender) + age , data = meta)$residual
PCx_res_qq[,j]=qqnorm(PCx_res[,j],plot.it=F)$x
}

save(PCx_res_qq, file = "/net/mulan/home/shanglu/GENOA/data/EA/LA/PCx_res_qq.RData")

load(paste0("/net/mulan/home/shanglu/GENOA/data/EA/permuteID_EA.RData"))

for(permu in 1:10){
	print(permu)
	PCx_res_qq_permute = PCx_res_qq[new_id[[permu]],]
	colnames(PCx_res_qq_permute) = rownames(expr_EA_proteingene)
	save(PCx_res_qq_permute, file = paste0("/net/mulan/home/shanglu/GENOA/data/EA/LA/PCx_res_qq_pheno_",permu+1,".RData"))
	}
	
# rename original expression
PCx_res_qq_permute = PCx_res_qq
colnames(PCx_res_qq_permute) = rownames(expr_EA_proteingene)
save(PCx_res_qq_permute, file = "/net/mulan/home/shanglu/GENOA/data/EA/LA/PCx_res_qq_pheno_1.RData")




for(chr in 1:22){
print(chr)
path_rfmix = paste0("/net/mulan/home/shanglu/GENOA/src/FRANC/franc/EA_data/GENOA_EA_use/elai/",chr)
rfmixout=paste0(path_rfmix,"/GENOA_EA_use.",chr,".elairfmix.txt")
res = rfmix.local.ancestry(rfmixout)
print("res done")
disp = display.ancestry(res)
save(disp, file = paste0(path_rfmix,"/disp_chr_",chr,".elairfmix.RData"))
print("disp done")
}

```

# run LAMatrix by gene in AA on greatlakes

```R
args <- as.numeric(commandArgs(TRUE))
gene = args[1]
pheno = args[2]
print(gene)
print(pheno)

library(LAMatrix)

load("/nfs/sph-mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_use.RData")
chr     = gene_AA_anno_use$chr[gene]
j_index = gene_AA_anno_use$index_j[gene]

path_rfmix = paste0("/nfs/sph-mulan/home/shanglu/GENOA/src/FRANC/franc/test_data/GENOA_AA_use/elai/",chr)
path_AAbed = paste0("/nfs/sph-mulan/home/shanglu/GENOA/data/AA/gene_bed")
path_LA = paste0("/nfs/sph-mulan/home/shanglu/GENOA/data/AA/LA/chr_",chr)
load(paste0("/nfs/sph-mulan/home/shanglu/GENOA/data/AA/LA/PCx_res_qq_pheno_",pheno,".RData"))
load(paste0(path_rfmix,"/disp_chr_",chr,".elairfmix.RData"))
snpinfo = paste0(path_rfmix,"/GENOA_AA_use.",chr,".snpinfo.txt")
mysnps = read.table(snpinfo,header=T)
mysnps$rs = as.character(mysnps$rs)
load(paste0(path_AAbed,"/chr_",chr,"/chr_",chr,"_gene_",j_index,".RData"))

print("data loaded")

	commonsnps = intersect(mysnps$rs , colnames(a))
	snps.genotype = a[,match(commonsnps, colnames(a))]
	snps.localances = disp[match(commonsnps,mysnps$rs),]
	gene.pheno = PCx_res_qq_permute[,match(gene_AA_anno_use$GENE[gene],colnames(PCx_res_qq_permute))]
	snpnum = length(commonsnps)
		
		n = 1032;# Number of columns (samples)
		nc = 0;# Number ofs covariates
		
		res_p_mat = c()
		res_beta_mat = c()
		gene.mat = gene.pheno
		genes = SlicedData$new( matrix( gene.mat, nrow = 1 ))
		
	for(j in 1:snpnum){

		snps.mat = snps.genotype[,j]
		local.mat = snps.localances[j,]
		snps = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
		locals = SlicedData$new( matrix(local.mat,nrow=1))

		modelLOCAL = 930507L;
		me = LAMatrix_main(
 			snps = snps,
			gene = genes,
 			local = locals,
 			pvOutputThreshold = 1,
 			useModel = modelLOCAL,
 			verbose = TRUE,
 			output_file_name = paste0(path_LA,"/chr_",chr,"_gene_",j_index,"_pheno_",pheno),
 			pvalue.hist = FALSE)
 			
 		res_p_mat= c(res_p_mat,unlist(me$all$eqtls[4]))
 		res_beta_mat = c(res_beta_mat,unlist(me$all$eqtls[6]))
 		}
 		res_dat = data.frame(res_p_mat, res_beta_mat)
 		colnames(res_dat) = c("pval","beta")
 		res_dat$SNP = colnames(a)[match(commonsnps, colnames(a))]
 		res_dat$GENE = gene_AA_anno_use$GENE[gene]
 		save(res_dat, file = paste0(path_LA,"/res_dat_chr_",chr,"_gene_",j_index,"_pheno_",pheno,".RData"))

print("Done!")

===============================================
#!/bin/bash
#SBATCH --job-name=AArfmixgreat
#SBATCH --array=1-3000%200
#SBATCH --cpus-per-task=1
#SBATCH --output=/nfs/sph-mulan/home/shanglu/GENOA/data/AA/LA/out/AArfmixgreat%a.out
#SBATCH --error=/nfs/sph-mulan/home/shanglu/GENOA/data/AA/LA/out/AArfmixgreat%a.err
#SBATCH --workdir=/nfs/sph-mulan/home/shanglu/GENOA/src/FRANC/franc/franc_interface
#SBATCH --partition=standard
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=8000MB
#SBATCH --account=xzhousph

bash
module load R
let k=0
for ((gene=1;gene<=3000;gene++)); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
for ((pheno=1;pheno<=11;pheno++)); do
Rscript --verbose ./AA_rfmix_greatlakes.r ${gene} ${pheno}
done
fi
done




```

# run LAMatrix by gene in EA on greatlakes
```R
args <- as.numeric(commandArgs(TRUE))
gene = args[1]
pheno = args[2]
print(gene)
print(pheno)

library(LAMatrix)

load("/nfs/sph-mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_use.RData")
chr     = gene_EA_anno_use$chr[gene]
j_index = gene_EA_anno_use$index_j[gene]

path_rfmix = paste0("/nfs/sph-mulan/home/shanglu/GENOA/src/FRANC/franc/EA_data/GENOA_EA_use/elai/",chr)
path_EAbed = paste0("/nfs/sph-mulan/home/shanglu/GENOA/data/EA/gene_bed")
path_LA = paste0("/nfs/sph-mulan/home/shanglu/GENOA/data/EA/LA/chr_",chr)
load(paste0("/nfs/sph-mulan/home/shanglu/GENOA/data/EA/LA/PCx_res_qq_pheno_",pheno,".RData"))
load(paste0(path_rfmix,"/disp_chr_",chr,".elairfmix.RData"))
snpinfo = paste0(path_rfmix,"/GENOA_EA_use.",chr,".snpinfo.txt")
mysnps = read.table(snpinfo,header=T)
mysnps$rs = as.character(mysnps$rs)
load(paste0(path_EAbed,"/chr_",chr,"/chr_",chr,"_gene_",j_index,".RData"))

print("data loaded")

	commonsnps = intersect(mysnps$rs , colnames(a))
	snps.genotype = a[,match(commonsnps, colnames(a))]
	snps.localances = disp[match(commonsnps,mysnps$rs),]
	gene.pheno = PCx_res_qq_permute[,match(gene_EA_anno_use$GENE[gene],colnames(PCx_res_qq_permute))]
	snpnum = length(commonsnps)
		
		n = 801;# Number of columns (samples)
		nc = 0;# Number ofs covariates
		
		res_p_mat = c()
		res_beta_mat = c()
		gene.mat = gene.pheno
		genes = SlicedData$new( matrix( gene.mat, nrow = 1 ))
		
	for(j in 1:snpnum){

		snps.mat = snps.genotype[,j]
		local.mat = snps.localances[j,]
		snps = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
		locals = SlicedData$new( matrix(local.mat,nrow=1))

		modelLOCAL = 930507L;
		me = LAMatrix_main(
 			snps = snps,
			gene = genes,
 			local = locals,
 			pvOutputThreshold = 1,
 			useModel = modelLOCAL,
 			verbose = TRUE,
 			output_file_name = paste0(path_LA,"/chr_",chr,"_gene_",j_index,"_pheno_",pheno),
 			pvalue.hist = FALSE)
 			
 		res_p_mat= c(res_p_mat,unlist(me$all$eqtls[4]))
 		res_beta_mat = c(res_beta_mat,unlist(me$all$eqtls[6]))
 		}
 		res_dat = data.frame(res_p_mat, res_beta_mat)
 		colnames(res_dat) = c("pval","beta")
 		res_dat$SNP = colnames(a)[match(commonsnps, colnames(a))]
 		res_dat$GENE = gene_EA_anno_use$GENE[gene]
 		save(res_dat, file = paste0(path_LA,"/res_dat_chr_",chr,"_gene_",j_index,"_pheno_",pheno,".RData"))

print("Done!")

===========================================

#!/bin/bash
#SBATCH --job-name=EArfmixgreat
#SBATCH --array=1-2000%500
#SBATCH --cpus-per-task=1
#SBATCH --output=/nfs/sph-mulan/home/shanglu/GENOA/data/EA/LA/out/EArfmixgreat%a.out
#SBATCH --error=/nfs/sph-mulan/home/shanglu/GENOA/data/EA/LA/out/EArfmixgreat%a.err
#SBATCH --workdir=/nfs/sph-mulan/home/shanglu/GENOA/src/FRANC/franc/franc_interface
#SBATCH --partition=standard
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=8000MB
#SBATCH --account=xzhousph

bash
module load R
let k=0
for ((gene=4001;gene<=6000;gene++)); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
for ((pheno=1;pheno<=11;pheno++)); do
Rscript --verbose ./EA_rfmix_greatlakes.r ${gene} ${pheno}
done
fi
done



```


# AA colloect results
```R
args <- as.numeric(commandArgs(TRUE))
chr = args[1]
pheno = args[2]
print(chr)
print(pheno)

path_LA = paste0("/net/mulan/home/shanglu/GENOA/data/AA/LA/chr_",chr)
load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_use.RData")
gene_AA_anno_use_chr = gene_AA_anno_use[which(gene_AA_anno_use$chr==chr),]
num = sum(gene_AA_anno_use$chr==chr)

pval = NULL
beta = NULL
SNP = NULL
GENE = NULL
result_dat = data.frame(pval,beta,SNP,GENE)
for(line in 1:num){
print(line)
j_index = gene_AA_anno_use_chr$index_j[line]
if(file.exists(paste0(path_LA,"/res_dat_chr_",chr,"_gene_",j_index,"_pheno_",pheno,".RData"))){
load(paste0(path_LA,"/res_dat_chr_",chr,"_gene_",j_index,"_pheno_",pheno,".RData"))
result_dat = rbind(result_dat, res_dat)
}
}
save(result_dat, file = paste0(path_LA,"/result_dat_chr_",chr,"_pheno_",pheno,".RData"))


============================

#!/bin/bash
#SBATCH --job-name=AAres
#SBATCH --array=1-11
#SBATCH --cpus-per-task=1
#SBATCH --output=/net/mulan/home/shanglu/GENOA/data/AA/LA/out/AAres%a.out
#SBATCH --error=/net/mulan/home/shanglu/GENOA/data/AA/LA/out/AAres%a.err
#SBATCH --partition=mulan
#SBATCH --mem-per-cpu=5000MB


bash

let k=0
chr=22
for ((pheno=1;pheno<=11;pheno++)); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
Rscript --verbose ./result_dat_bychr_bypheno.r ${chr}  ${pheno}
fi
done


====================================

for(pheno in 1:11){
print("pheno")
print(pheno)
pval = NULL
beta = NULL
SNP = NULL
GENE = NULL
result_table = data.frame(pval,beta,SNP,GENE)

for(chr in 1:22){
print(chr)
path_LA = paste0("/net/mulan/home/shanglu/GENOA/data/AA/LA/chr_",chr)
load(paste0(path_LA,"/result_dat_chr_",chr,"_pheno_",pheno,".RData"))
result_table = rbind(result_table, result_dat)
}
save(result_table, file = paste0("/net/mulan/home/shanglu/GENOA/data/AA/LA/AA_LA_table_pheno_",pheno,".RData"))
}


```


# EA collect result 
```R
args <- as.numeric(commandArgs(TRUE))
chr = args[1]
pheno = args[2]
print(chr)
print(pheno)

path_LA = paste0("/net/mulan/home/shanglu/GENOA/data/EA/LA/chr_",chr)
load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_use.RData")
gene_EA_anno_use_chr = gene_EA_anno_use[which(gene_EA_anno_use$chr==chr),]
num = sum(gene_EA_anno_use$chr==chr)

pval = NULL
beta = NULL
SNP = NULL
GENE = NULL
result_dat = data.frame(pval,beta,SNP,GENE)
for(line in 1:num){
print(line)
j_index = gene_EA_anno_use_chr$index_j[line]
if(file.exists(paste0(path_LA,"/res_dat_chr_",chr,"_gene_",j_index,"_pheno_",pheno,".RData"))){
load(paste0(path_LA,"/res_dat_chr_",chr,"_gene_",j_index,"_pheno_",pheno,".RData"))
result_dat = rbind(result_dat, res_dat)
}
}
save(result_dat, file = paste0(path_LA,"/result_dat_chr_",chr,"_pheno_",pheno,".RData"))

=================================

#!/bin/bash
#SBATCH --job-name=EAres
#SBATCH --array=1-300
#SBATCH --cpus-per-task=1
#SBATCH --output=/net/mulan/home/shanglu/GENOA/data/EA/LA/out/EAres%a.out
#SBATCH --error=/net/mulan/home/shanglu/GENOA/data/EA/LA/out/EAres%a.err
#SBATCH --partition=mulan
#SBATCH --mem-per-cpu=5000MB


bash

let k=0
for ((chr=1;chr<=22;chr++)); do
for ((pheno=1;pheno<=11;pheno++)); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
Rscript --verbose ./result_dat_bychr_byphenoEA.r ${chr}  ${pheno}
fi
done
done

=================================

for(pheno in 1:11){
print("pheno")
print(pheno)
pval = NULL
beta = NULL
SNP = NULL
GENE = NULL
result_table = data.frame(pval,beta,SNP,GENE)

for(chr in 1:22){
print(chr)
path_LA = paste0("/net/mulan/home/shanglu/GENOA/data/EA/LA/chr_",chr)
load(paste0(path_LA,"/result_dat_chr_",chr,"_pheno_",pheno,".RData"))
result_table = rbind(result_table, result_dat)
}
save(result_table, file = paste0("/net/mulan/home/shanglu/GENOA/data/EA/LA/EA_LA_table_pheno_",pheno,".RData"))
}






```

# analysis results
```R

# AA
library(qvalue)
library(data.table)
library(dplyr)

LAlist = list()
for(pheno in 1:11){
print(pheno)
load(paste0("/net/mulan/home/shanglu/GENOA/data/AA/LA/AA_LA_table_pheno_",pheno,".RData"))

LA_table = result_table %>%  group_by(GENE) %>%
summarise(minp = min(pval))
LAlist[[pheno]] =    LA_table      
 }
          
for(pheno in 1:11){
print(dim(LAlist[[pheno]]))
}          


    
min_p_AA = LAlist[[1]]$minp
order_p_AA = min_p_AA[order(min_p_AA)]
min_p_permu_AA = list()
for(i in 1:10){ 
min_p_permu_AA[[i]] = LAlist[[i+1]]$minp
}
min_p_AA_permu_mat = matrix(unlist(min_p_permu_AA), ncol = length(order_p_AA), byrow = TRUE) 



count = 0 
num_real = NULL 
FDR_AA = NULL 
for(threshold in order_p_AA){ 
count = count + 1 
num_real[count] = sum(min_p_AA<=threshold) 
FDR_AA[count] = mean(apply(min_p_AA_permu_mat, 1, function(x) sum(x<=threshold)/num_real[count])) 

if(FDR_AA[count]>0.2){
break
}

}


max(which(FDR_AA<0.05))
order_p_AA[ max(which(FDR_AA<0.05))]

> max(which(FDR_AA<0.05))
[1] 5553
> order_p_AA[ max(which(FDR_AA<0.05))]
[1] 5.962974e-05    
          
#--------------------------------------------     
# analysis for comparison between LA and GA
#--------------------------------------------  

pheno=1
load(paste0("/net/mulan/home/shanglu/GENOA/data/AA/LA/AA_LA_table_pheno_",pheno,".RData"))


LA_thr_AA = 5.962974e-05  

LA_table = result_table %>%  group_by(GENE) %>%
summarise(minp = min(pval))

sum(LA_table$minp<=LA_thr_AA)

LA_eGene = LA_table$GENE[which(LA_table$minp<=LA_thr_AA)]

#--------------------------------------------  
# eGene overlap compare
#--------------------------------------------  

eGene_GA = gene_AA_anno_use$GENE[which(gene_AA_anno_use$eGene==1)]
eGene_LA = LA_eGene

# 5312 genes overlap between GA and LA
length(intersect(eGene_GA, eGene_LA))
> length(intersect(eGene_GA, eGene_LA))
[1] 5312

#--------------------------------------------  
# compare beta
#--------------------------------------------  

load("/net/mulan/home/shanglu/GENOA/analysis/AA/eqtlmapping/AA_table.RData")

> dim(AA_table)
[1] 14511338       12
> dim(result_table)
[1] 14352738        4

result_table$rs = result_table$SNP
result_table$LAbeta = result_table$beta

common_table = merge(AA_table, result_table, by = c("GENE","rs"))
save(common_table, file = "/net/mulan/home/shanglu/GENOA/data/AA/LA/common_table_LA_AA.RData")

mydata = common_table
p = ggplot(mydata, aes(x=beta.x, y=beta.y)) +    
  geom_point(shape=19, fill="violet", color="violet", size=0.5)+
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  labs(title=paste0("Effect sizes in LA and GA"),
       x="Effect size in GA", y = "Effect size in LA")+
  xlim(-2.1,2.1)+
  ylim(-2.1,2.1) 

path_LA = "/net/mulan/home/shanglu/GENOA/data/AA/LA"
tiff(paste0(path_LA,"/AA_effectsize_LAvsGA.tiff"), units="in", width=10, height=10, res=150)    
p
dev.off()


> summary(common_table$beta.x)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
-1.986350 -0.077252 -0.001134 -0.001270  0.074448  2.029662 
>  summary(common_table$beta.y)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
-2.003489 -0.077968 -0.001122 -0.001306  0.075160  2.073039 


cor.test(common_table$beta.x, common_table$beta.y)$estimate
cor.test(common_table$beta.x, common_table$beta.y)$p.value
> cor.test(common_table$beta.x, common_table$beta.y)$estimate
      cor 
0.9914589 
> cor.test(common_table$beta.x, common_table$beta.y)$p.value
[1] 0

cor.test(common_table$p_wald, common_table$pval, method = "spearman")$estimate
cor.test(common_table$p_wald, common_table$pval, method = "spearman")$p.value
> cor.test(common_table$p_wald, common_table$pval, method = "spearman")$estimate
      rho 
0.9628653 
> cor.test(common_table$p_wald, common_table$pval, method = "spearman")$p.value
[1] 0

#--------------------------------------------  
# compare eSNPs
#--------------------------------------------

# in LA, AA
AA_thr = 6.245907e-05
AA_table_eSNP = AA_table[which(AA_table$p_wald <= AA_thr),]
LA_table_eSNP = result_table[which(result_table$pval <= LA_thr_AA),]
> dim(LA_table_eSNP)
[1] 357072      4
> dim(AA_table_eSNP)
[1] 354931     12

LA_table_eSNP$rs = LA_table_eSNP$SNP
LA_table_eSNP$LAbeta = LA_table_eSNP$beta

common_eSNP_AA = merge(AA_table_eSNP, LA_table_eSNP, by = c("GENE","rs"))
save(common_eSNP_AA, file = "/net/mulan/home/shanglu/GENOA/data/AA/LA/common_eSNP_table_LA_AA.RData")

> dim(common_eSNP_AA)
[1] 325571     16
> 357072-325571
[1] 31501
> 354931 - 325571
[1] 29360


cor.test(common_eSNP_AA$beta.x, common_eSNP_AA$beta.y)$estimate
cor.test(common_eSNP_AA$beta.x, common_eSNP_AA$beta.y)$p.value

cor.test(common_eSNP_AA$p_wald, common_eSNP_AA$pval, method = "spearman")$estimate
cor.test(common_eSNP_AA$p_wald, common_eSNP_AA$pval, method = "spearman")$p.value

> cor.test(common_eSNP_AA$beta.x, common_eSNP_AA$beta.y)$estimate
      cor
0.9987615
> cor.test(common_eSNP_AA$beta.x, common_eSNP_AA$beta.y)$p.value
[1] 0
>
> cor.test(common_eSNP_AA$p_wald, common_eSNP_AA$pval, method = "spearman")$estimate
      rho
0.9636551
> cor.test(common_eSNP_AA$p_wald, common_eSNP_AA$pval, method = "spearman")$p.value
[1] 0

library(ggplot2)
mydata = common_eSNP_AA
p = ggplot(mydata, aes(x=beta.x, y=beta.y)) +    
  geom_point(shape=19, fill="violet", color="violet", size=0.5)+
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  labs(title=paste0("Effect sizes in LA and GA"),
       x="Effect size in GA", y = "Effect size in LA")+
  xlim(-2.1,2.1)+
  ylim(-2.1,2.1) 

path_LA = "/net/mulan/home/shanglu/GENOA/data/AA/LA"
tiff(paste0(path_LA,"/AA_eSNP_effectsize_LAvsGA.tiff"), units="in", width=10, height=10, res=150)    
p
dev.off()


============================================

# EA


for(pheno in 1:11){
print("pheno")
print(pheno)
pval = NULL
beta = NULL
SNP = NULL
GENE = NULL
result_table = data.frame(pval,beta,SNP,GENE)

for(chr in 1:22){
print(chr)
path_LA = paste0("/net/mulan/home/shanglu/GENOA/data/EA/LA/chr_",chr)
load(paste0(path_LA,"/result_dat_chr_",chr,"_pheno_",pheno,".RData"))
result_table = rbind(result_table, result_dat)
}
save(result_table, file = paste0("/net/mulan/home/shanglu/GENOA/data/EA/LA/EA_LA_table_pheno_",pheno,".RData"))
}



       
library(qvalue)
library(data.table)
library(dplyr)

LAlist = list()
for(pheno in 1:11){
print(pheno)
load(paste0("/net/mulan/home/shanglu/GENOA/data/EA/LA/EA_LA_table_pheno_",pheno,".RData"))

LA_table = result_table %>%  group_by(GENE) %>%
summarise(minp = min(pval))
LAlist[[pheno]] =    LA_table      
 }
          
for(pheno in 1:11){
print(dim(LAlist[[pheno]]))
}          

min_p_EA = LAlist[[1]]$minp
order_p_EA = min_p_EA[order(min_p_EA)]
min_p_permu_EA = list()
for(i in 1:10){ 
min_p_permu_EA[[i]] = LAlist[[i+1]]$minp
}
min_p_EA_permu_mat = matrix(unlist(min_p_permu_EA), ncol = length(order_p_EA), byrow = TRUE) 



count = 0 
num_real = NULL 
FDR_EA = NULL 
for(threshold in order_p_EA){ 
count = count + 1 
num_real[count] = sum(min_p_EA<=threshold) 
FDR_EA[count] = mean(apply(min_p_EA_permu_mat, 1, function(x) sum(x<=threshold)/num_real[count])) 

if(FDR_EA[count]>0.2){
break
}

}


max(which(FDR_EA<0.05))
order_p_EA[ max(which(FDR_EA<0.05))]

> max(which(FDR_EA<0.05))
[1] 4586
> order_p_EA[ max(which(FDR_EA<0.05))]
[1] 0.0001376416
          
#--------------------------------------------     
# analysis for comparison between LA and GA
#--------------------------------------------  
load("/net/mulan/home/shanglu/GENOA/data/EA/gene_EA_anno_use.RData")
pheno=1
load(paste0("/net/mulan/home/shanglu/GENOA/data/EA/LA/EA_LA_table_pheno_",pheno,".RData"))


LA_thr_EA = 0.0001376416

LA_table = result_table %>%  group_by(GENE) %>%
summarise(minp = min(pval))

sum(LA_table$minp<=LA_thr_EA)

LA_eGene = LA_table$GENE[which(LA_table$minp<=LA_thr_EA)]

#--------------------------------------------  
# eGene overlap compare
#--------------------------------------------  


eGene_GA = gene_EA_anno_use$GENE[which(gene_EA_anno_use$eGene==1)]
eGene_LA = LA_eGene


> length(LA_eGene)
[1] 4586
> length(eGene_GA)
[1] 4402
> length(intersect(eGene_GA, eGene_LA))
[1] 4333

#--------------------------------------------  
# compare beta
#--------------------------------------------  

load("/net/mulan/home/shanglu/GENOA/analysis/EA/eqtlmapping/EA_table.RData")

dim(EA_table)
dim(result_table)


result_table$rs = result_table$SNP
result_table$LAbeta = result_table$beta

common_table = merge(EA_table, result_table, by = c("GENE","rs"))
save(common_table, file = "/net/mulan/home/shanglu/GENOA/data/EA/LA/common_table_LA_EA.RData")

mydata = common_table
p = ggplot(mydata, aes(x=beta.x, y=beta.y)) +    
  geom_point(shape=19, fill="cornflowerblue", color="cornflowerblue", size=0.5)+
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  labs(title=paste0("Effect sizes in LA and GA"),
       x="Effect size in GA", y = "Effect size in LA")+
  xlim(-2.1,2.1)+
  ylim(-2.1,2.1) 

path_LA = "/net/mulan/home/shanglu/GENOA/data/EA/LA"
tiff(paste0(path_LA,"/EA_effectsize_LAvsGA.tiff"), units="in", width=10, height=10, res=150)    
p
dev.off()




cor.test(common_table$beta.x, common_table$beta.y)$estimate
cor.test(common_table$beta.x, common_table$beta.y)$p.value

cor.test(common_table$p_wald, common_table$pval, method = "spearman")$estimate
cor.test(common_table$p_wald, common_table$pval, method = "spearman")$p.value

> cor.test(common_table$beta.x, common_table$beta.y)$estimate

      cor 
0.9941167 
> cor.test(common_table$beta.x, common_table$beta.y)$p.value
[1] 0
> 
> cor.test(common_table$p_wald, common_table$pval, method = "spearman")$estimate
      rho 
0.9721395 
> cor.test(common_table$p_wald, common_table$pval, method = "spearman")$p.value
[1] 0


#--------------------------------------------  
# compare eSNPs
#--------------------------------------------

# in LA, EA
EA_thr = 0.0001385504
EA_table_eSNP = EA_table[which(EA_table$p_wald <= EA_thr),]
LA_table_eSNP = result_table[which(result_table$pval <= LA_thr_EA),]
dim(LA_table_eSNP)
dim(EA_table_eSNP)

> dim(LA_table_eSNP)
[1] 372240      4
> dim(EA_table_eSNP)
[1] 371309     12
> 


LA_table_eSNP$rs = LA_table_eSNP$SNP
LA_table_eSNP$LAbeta = LA_table_eSNP$beta

common_eSNP_EA = merge(EA_table_eSNP, LA_table_eSNP, by = c("GENE","rs"))
save(common_eSNP_EA, file = "/net/mulan/home/shanglu/GENOA/data/EA/LA/common_eSNP_table_LA_EA.RData")

> dim(common_eSNP_EA)
[1] 347949     16

372240-347949
371309 - 347949
> 372240-347949
[1] 24291
> 371309 - 347949
[1] 23360



cor.test(common_eSNP_EA$beta.x, common_eSNP_EA$beta.y)$estimate
cor.test(common_eSNP_EA$beta.x, common_eSNP_EA$beta.y)$p.value

cor.test(common_eSNP_EA$p_wald, common_eSNP_EA$pval, method = "spearman")$estimate
cor.test(common_eSNP_EA$p_wald, common_eSNP_EA$pval, method = "spearman")$p.value

> cor.test(common_eSNP_EA$beta.x, common_eSNP_EA$beta.y)$estimate
      cor 
0.9994796 
> cor.test(common_eSNP_EA$beta.x, common_eSNP_EA$beta.y)$p.value
[1] 0
> cor.test(common_eSNP_EA$p_wald, common_eSNP_EA$pval, method = "spearman")$estimate
      rho 
0.9831693 
> cor.test(common_eSNP_EA$p_wald, common_eSNP_EA$pval, method = "spearman")$p.value
[1] 0







mydata = common_table
p = ggplot(mydata, aes(x=beta.x, y=beta.y)) +    
  geom_point(shape=19, fill="cornflowerblue", color="cornflowerblue", size=0.5)+
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  labs(title=paste0("Effect sizes in LA and GA"),
       x="Effect size in GA", y = "Effect size in LA")+
  xlim(-2.1,2.1)+
  ylim(-2.1,2.1) 

path_LA = "/net/mulan/home/shanglu/GENOA/data/EA/LA"
tiff(paste0(path_LA,"/EA_effectsize_LAvsGA.tiff"), units="in", width=10, height=10, res=150)    
p
dev.off()




#--------------------------------------------  
# plot -log10(p)
#--------------------------------------------

# EA
library(ggplot2)
mydata = common_table
common_table$log10pGA = -log10(common_table$p_wald)
common_table$log10pLA = -log10(common_table$pval)

> summary(common_table$log10pGA)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  0.0000   0.1614   0.4099   1.0414   0.9120 198.8397 
> summary(common_table$log10pLA )
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  0.0000   0.1638   0.4171   1.0807   0.9332 207.2681 
  
mydata = common_table
p = ggplot(mydata, aes(x=log10pGA, y=log10pLA)) +    
  geom_point(shape=19, fill="cornflowerblue", color="cornflowerblue", size=0.5)+
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  labs(title=paste0("-log10(p) in LA and GA"),
       x="-log10(p) in GA", y = "-log10(p) in LA")+
  xlim(0,210)+
  ylim(0,210) 

path_LA = "/net/mulan/home/shanglu/GENOA/data/EA/LA"
tiff(paste0(path_LA,"/EA_log10p_LAvsGA.tiff"), units="in", width=10, height=10, res=150)    
p
dev.off()


# AA
library(ggplot2)

common_table$log10pGA = -log10(common_table$p_wald)
common_table$log10pLA = -log10(common_table$pval)

> summary(common_table$log10pGA)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  0.0000   0.1532   0.3834   0.8331   0.8245 264.2144 
> summary(common_table$log10pLA )
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  0.0000   0.1549   0.3885   0.8564   0.8400 279.1684 
  
mydata = common_table
p = ggplot(mydata, aes(x=log10pGA, y=log10pLA)) +    
  geom_point(shape=19, fill="violet", color="violet", size=0.5)+
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=22)+
  labs(title=paste0("-log10(p) in LA and GA"),
       x="-log10(p) in GA", y = "-log10(p) in LA")+
  xlim(0,280)+
  ylim(0,280) 

path_LA = "/net/mulan/home/shanglu/GENOA/data/AA/LA"
tiff(paste0(path_LA,"/AA_log10p_LAvsGAviolet.tiff"), units="in", width=10, height=10, res=150)    
p
dev.off()





```








