# SAIGE PROJECT OVERVIEW

The project revolves around neurodevelopmental diseases which are Autism, ADHD, Schizophrenia and Bipolar Disorder.
The disorders described are characterized by the involvement of both rare and common genetic variants that influence brain function and development. They have a polygenic nature meaning that each disorder is influenced by multiple genes, with no single gene being solely responsible. Instead, numerous genetic variants, often spread across many genes, contribute to the overall risk.

The SAIGE-GENE+ is a statistical tool designed to perform gene-based association tests, particularly in large-scale biobank data sets. 

Data for this porject is coming from dried blood samples and it is stored in a file called iPSYCH2016 which has been updated from iPSYCH2012 and iPSYCH2015.

---

PREPARATION FOR SAIGE
(The work will be done in a closed zone)
1.	Install conda in the closed zone
2.	Create an enviroment to conduct the SAIGE methodology
3.	Conda install SAIGE

https://anaconda.org/bioconda/r-saige

For the data it is also a must to have an enviroment with plink installed since we are going to need it to obtain some input for SAIGE.


## OVERVIEW of STEPS: (https://saigegit.github.io/SAIGE-doc/)

<img width="750" alt="Captura de pantalla 2024-10-29 a las 12 36 17" src="https://github.com/user-attachments/assets/5a31d46f-60d7-4f30-93b3-3d2cdd50d7b9">




# GRM (STEP 0)

Use iPSYCH common variants to construct GRM
1. Keep only MAF>=5% in ipsych.vcf --> filtered_data
```sh
plink --vcf ipsych.vcf --maf 0.05 --make-bed --out filtered_data 
```

2. Do LD pruning for common variants --> pruned_data
```sh
plink --bfile filtered_data --indep-pairwise 50 5 0.2 --out pruned_data
plink --bfile filtered_data --extract pruned_data.prune.in --make-bed --out pruned_filtered_data
```

This will output new .bed, .bim, and .fam files containing only the pruned variants.

3. Generate sparse GRM
```sh
Rscript createSparseGRM.R       \
    --plinkFile=/home/ialbacoto/Alba_PiB_project_2024fall/data/pruned_filtered_data \
    --nThreads=72  \
    --outputPrefix=/home/ialbacoto/Alba_PiB_project_2024fall/people/albacoto/pruned_filtered_sparseGRM      \
    --numRandomMarkerforSparseKin=5000      \
    --relatednessCutoff=0.05
```

   - nThreads: number of threads to use for parallel processing
   - numRandomMarkerforSparseKin=5000: number of random markers for calculating kinship/relatedness
   - Cutoff: cutoff for excluding individuals with high relatedness

High-performance computing cluster or a powerful machine with many cores uses 72 threads for faster execution.
More markers typically lead to more accurate estimates of relatedness.
Relatedness cutoff of 0.05 (stricter), meaning individuals with closer relationships (e.g., 1st or 2nd degree relatives) are excluded.


OBTAIN:
-	sparseGRM.mtx --> file containing the sparse GRM
-	sparseGRM.mtx.sampleIDs.txt --> file containing sample IDs in the sparse GRM.


# STEP 1
Previous step before performing STEP 1 is to separate the raw vcf file depending on the disorder. That way we are able to see each disorder separately.


1. First we will separate the raw vcf file from the asd_adhd_sz_bp_ctl.ped in _Rstudio_. We need to keep individuals with either ASD==1 (positive for ASD) or Control==1 (which are the controls) --> We will obtain a file called ASD_samples.txt
  
```sh
install.packages("data.table")
library(data.table)
```

```sh
ped_file <- "yourfile.ped"  # Replace with your actual .ped file path
ped_data <- fread(ped_file, header = FALSE)
```

```sh
asd_samples <- ped_data[V13 == 1, V1]  # ASD cases (where V12 == 1)
control_samples <- ped_data[V25 == 1, V1]  # Control cases (where V25 == 1)

# Combine and save to a file
all_samples <- c(asd_samples, control_samples)
writeLines(all_samples, "ASD_samples.txt")
```


2. Following, we run in the _terminal, closed zone_ a bcftools comand to create the new vcf file with only ASD & control data.
```sh
bcftools view -S ASD_samples.txt -a -c 1 -Oz -o asd_ctl.vcf.gz asd_adhd_sz_bp_ctl_hg38_inDGCCregions_vepPICK_casecontrol.vcf.gz
```

```-S ASD_samples.txt```: Specifies the file containing the sample IDs we want to keep (ASD samples).

```-a```: This option is set to include all alleles, which is useful to keep all information related to variants.

```-c 1```: Keeps only sites with at least one allele for the specified samples, ensuring that we only get variants that are present in the filtered list.

```-Oz```: This option compresses the output file in the BGZF format, which is suitable for VCF files.

```-o asd_ctl.vcf.gz```: Specifies the name of the output file.


### STEP1 
Obtain plink files from the vcf file
```sh
plink --vcf asd_ctl.vcf.gz --make-bed --out asd_ctl
```

(Optional) get ids for 1000 random markers for each MAC category. Calculate allele counts for each marker in the large plink file.
```sh
plink2 --bfile asd_ctl --freq counts --out asd_ctl
```

Randomly extract IDs for markers falling in the two MAC categories
```sh
cat <(tail -n +2 asd_ctl.acount | awk '((2*$6-$5) < 20 && (2*$6-$5) >= 10) || ($5 < 20 && $5 >= 10) {print $2}' | shuf -n 1000) \
<(tail -n +2 asd_ctl.acount | awk ' $5 >= 20 && (2*$6-$5)>= 20 {print $2}' | shuf -n 1000) > asd_ctl.forCate_vr.markerid.list
```

Extract markers from the large plink file
```sh
plink2 --bfile asd_ctl --extract asd_ctl.forCate_vr.markerid.list --make-bed --out asd_ctl.forCate_vr
```

* Make sure that the chromosomes X & Y in the .bim file are substituted by number 23 & 24 (because SAIGE does not support letters).

After this we are finally able to start with STEP1:

INPUT FILES:
- Sparse GRM files --> obtained in STEP0
- Plink file --> obtained in the steps above
- Pheno file --> We obtain it from the csv file. From *csv --> txt* file.

  It is required that the file contains one column for sample IDs and one column for the phenotype. It may contain columns for covariates. (can use *rstudio*).
  Make sure to have y_binary column: 0 for control + 1 for ASD & make necessary changes to covariates (we can remove columns that are not necessary).

  
  


```sh
 Rscript step1_fitNULLGLMM.R     \
     --plinkFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/asd_ctl.forCate_vr"  \
     --sparseGRMFile="/home/ialbacoto/Alba_PiB_project_2024fall/people/albacoto/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"   \
     --sparseGRMSampleIDFile="/home/ialbacoto/Alba_PiB_project_2024fall/people/albacoto/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"     \
     --phenoFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/phenotype_data.txt" \
     --phenoCol=y_binary \
     --covarColList=x1,x2 \
     --sampleIDColinphenoFile=vcfID \
     --traitType=binary        \
     --outputPrefix="/home/ialbacoto/Alba_PiB_project_2024fall/data" \
     --nThreads=64    \
     --useSparseGRMtoFitNULL=FALSE    \
     --isCateVarianceRatio=TRUE      \
     --useSparseGRMforVarRatio=TRUE  \
     --IsOverwriteVarianceRatioFile=TRUE
```

```qCovarColList```: This list is used for quantitative covariates.


To have all the covariates I generated a list with all covariates so I don't need to enter each name manually (done with *rstudio*) --> covarColList

```sh
Rscript step1_fitNULLGLMM.R     \
    --plinkFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/asd_ctl.forCate_vr"  \
    --sparseGRMFile="/home/ialbacoto/Alba_PiB_project_2024fall/people/albacoto/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"   \
    --sparseGRMSampleIDFile="/home/ialbacoto/Alba_PiB_project_2024fall/people/albacoto/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"     \
    --phenoFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/phenotype_data.txt" \
    --phenoCol=y_binary \
    --covarColList=$covarColList \
    --sampleIDColinphenoFile=vcfID \
    --traitType=binary        \
    --outputPrefix="/home/ialbacoto/Alba_PiB_project_2024fall/data" \
    --nThreads=64    \
    --useSparseGRMtoFitNULL=FALSE    \
    --isCateVarianceRatio=TRUE      \
    --useSparseGRMforVarRatio=TRUE  \
    --IsOverwriteVarianceRatioFile=TRUE
```
    
OUTPUT FILES:
- model file --> *data.Rda*.

  Inside this file: theta --> vector with 2 elements & coefficients --> coefficient estimates for covariates in the model and intercept estimate
  
- variance ratio file --> *data.varianceRatio.txt*

both are going to be inputs for STEP2




### STEP2

INPUT FILES:
- Group file --> with genes, marker IDs and annotations.

To obtain this file we should perform the following steps in Rstudio:
1. Load full ASD .bim file *asd_ctl.bim* and extract marker IDs
2. Load .gene.marker.ann.txt file
3. Filter the .gene.marker.ann.txt for ASD markers
4. Merge the data and create var & anno lines (to obtain the correct formatting)

- Sample file --> sample IDs extracted from .fam file (with Rstudio)

- GMMAT & and varianceRatio are the files generated in STEP1


```sh
Rscript step2_SPAtests.R        \
     --bedFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/asd_ctl.bed"       \
     --bimFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/asd_ctl.bim"       \
     --famFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/asd_ctl.fam"       \
     --SAIGEOutputFile="/home/ialbacoto/Alba_PiB_project_2024fall/genotype_groupTest_out.txt" \
     --chrom=1 \
     --LOCO=TRUE    \
     --AlleleOrder=alt-first \
     --minMAF=0 \
     --minMAC=0.5 \
     --sampleFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/sample_ids.txt" \
     --GMMATmodelFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/output/data.rda" \
     --varianceRatioFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/output/data.varianceRatio.txt"      \
     --sparseGRMFile="/home/ialbacoto/Alba_PiB_project_2024fall/people/albacoto/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"   \
     --sparseGRMSampleIDFile="/home/ialbacoto/Alba_PiB_project_2024fall/people/albacoto/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"    \
     --groupFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/ASD_group_file_formatted.txt"    \
     --annotation_in_groupTest="pLoF,severeMis:pLoF,severeMis:pLoF:moderateMis,severeMis:pLoF:moderateMis:synonymous"       \
     --maxMAF_in_groupTest=0.0001,0.001,0.01
```

```AlleleOrder```: can be alt-first or ref-first. 
- ref-first: This setting treats the first allele listed in the genotype file as the reference allele (typically the more common allele in the population). The second allele is considered the alternative (less common or variant) allele.
- alt-first: This setting treats the first allele as the alternative allele and the second as the reference.



OUPUT FILES:
- genotype_groupTest_out.txt --> file with region- or gene-based association test results
- genotype_groupTest_out.txt.singleAssoc.txt--> file with association test results for single markers in the set-based tests




 

