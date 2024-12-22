# SAIGE PROJECT OVERVIEW

The project revolves around neurodevelopmental diseases which are Autism, ADHD, Schizophrenia and Bipolar Disorder.
The disorders described are characterized by the involvement of both rare and common genetic variants that influence brain function and development. They have a polygenic nature meaning that each disorder is influenced by multiple genes, with no single gene being solely responsible. Instead, numerous genetic variants, often spread across many genes, contribute to the overall risk.

The SAIGE-GENE+ is a statistical tool designed to perform gene-based association tests, particularly in large-scale biobank data sets. 
What SAIGE does is to test vairants with MAF (Minor Allele Frequency) <= 0.1, 0.01, 0.001 % (we use multiple mAF cutoffs to improve POWER).
* To incorporate multiple --> MAF cutoffs --> functional annotations, multiple tests are needed for each gene and results need to be combined using:
    1. minimum p-value
    2. Cauchy combination method

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
ped_file <- "asd_adhd_sz_bp_ctl.ped"  # Replace with .ped file path
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
plink --vcf adhd_ctl.vcf.gz --keep-allele-order --make-bed --out adhd2_ctl
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

  ¡¡¡Change X & Y chromosmes to 23 & 24 in plink bim file!!! 


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

version2
```sh
Rscript step2_SPAtests.R        \
     --bedFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/asd_ctl.bed"       \
     --bimFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/asd_ctl.bim"       \
     --famFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/asd_ctl.fam"       \
     --SAIGEOutputFile="/home/ialbacoto/Alba_PiB_project_2024fall/genotype_groupTest_out.txt" \
     --LOCO=FALSE    \
     --AlleleOrder=ref-first \
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

version3
```sh
Rscript step2_SPAtests.R        \
     --bedFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/asd_ctl.bed"       \
     --bimFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/asd_ctl.bim"       \
     --famFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/asd_ctl.fam"       \
     --SAIGEOutputFile="/home/ialbacoto/Alba_PiB_project_2024fall/genotype_groupTest_out.txt" \
     --LOCO=FALSE    \
     --AlleleOrder=ref-first \
     --minMAF=0 \
     --minMAC=0.5 \
     --sampleFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/sample_ids.txt" \
     --GMMATmodelFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/output/data.rda" \
     --varianceRatioFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/output/data.varianceRatio.txt"      \
     --sparseGRMFile="/home/ialbacoto/Alba_PiB_project_2024fall/people/albacoto/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"   \
     --sparseGRMSampleIDFile="/home/ialbacoto/Alba_PiB_project_2024fall/people/albacoto/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"    \
     --groupFile="/home/ialbacoto/Alba_PiB_project_2024fall/data/ASD/ASD_group_file_formatted.txt"    \
     --annotation_in_groupTest="pLoF,severeMis:pLoF,severeMis:pLoF:moderateMis,severeMis:pLoF:moderateMis:synonymous"       \
     --maxMAF_in_groupTest=0.0001,0.0003,0.00055
```

Next, I will try it changing the last row of the command to: ```--maxMAC_in_groupTest=1,2,5,10```



```AlleleOrder```: can be alt-first or ref-first. 
- ref-first: This setting treats the first allele listed in the genotype file as the reference allele (typically the more common allele in the population). The second allele is considered the alternative (less common or variant) allele.
- alt-first: This setting treats the first allele as the alternative allele and the second as the reference.



OUPUT FILES:
- genotype_groupTest_out.txt --> file with region- or gene-based association test results. The following columns will appear:
    1. Region: set name
    2. Group: annotation mask
    3. max_MAF: maximum MAF cutoff
    4. Pvalue: p value for SKAT-O test
    5. Pvalue_Burden: p value for BURDEN test
    6. Pvalue_SKAT: p value for SKAT test
    7. BETA_Burden: effect size of BURDEN test
    8. SE_Burden: standard error of BETA_Burden
    9. MAC: minor allele count in the set
    10. MAC_case: minor allele count in cases
    11. MAC_control: minor allele count in controls
    12. Number_rare: number of markers that are not ultra-rare with MAC > MACCutoff_to_CollapseUltraRare (=10 by default)
    13. Number_ultra_rare: number of markers that are ultra-rare with MAC <= MACCutoff_to_CollapseUltraRare (=10 by default)

- genotype_groupTest_out.txt.singleAssoc.txt--> file with association test results for single markers in the set-based tests




# ANALYSIS OF THE RESULTS
Perform the analysis with Rstudio

## Comparison of my ASD output vs Paper (NATURE, Table11)
We are going to compare the output of STEP2 with one table of paper: "Rare coding variation provides insight into the genetic architecture and phenotypic context of autism". That table contains some information that is going to help us get a comparison of out otuput with other outputs so we can validate and give some sense into our results. The information of utility that Table 11 gives us is: genes and its respective p-value.

1. Load Table 11 & my data (my results)
2. Filter Table 11 to only include genes present in my data
3. Merge the 2 tables based on common gene column
4. Select only the relevenat columns and calculate the difference between 2 pvalues

What we can obtain from the analysis/comparison:
- Genes with lowest p-values difference (ex: 20 genes)
- 20 genes with smallest p_TADA_ASD (which is the pvalue from the paper)
- 20 genes with smallest pvalue (which is the pvalue that I obtained)
We can overlap these last 2 combination of genes to find tihe overlapping genes

FOR VISUALIZATION: Scatter plot & bar plot

 
## Comparison of my ADHD output vs Paper (Table3)
We will make the same comparison as before but with the ADHD results from SAIGE+ (STEP2) and comparing it with another table from another paper: "Rare de novo damaging DNA variants are enriched in attention-deficit/hyperactivity disorder and implicate risk genes"

We will follow the same steps

What we can obtain from the analysis/comparison:
- Genes with lowest p-values difference (ex: 20 genes)
- 100 genes with smallest p_TADA_ASD (which is the pvalue from the paper)
- 100 genes with smallest pvalue (which is the pvalue that I obtained)
We can overlap these last 2 combination of genes to find tihe overlapping genes (we did it with a larger amount of genes this time because with lower genes the overlpaping genes dataset was empty).

FOR VISUALIZATION: Scatter plot & bar plot


## Comparison of the cross disorder between ASD & ADHD

1. DEFINE RARE VARIANTS
    - Filter the VCF for rare variants based on MAF (<1%) --> use bcftools

        ```bcftools view -i 'INFO/MAF<0.01' yourfile.vcf.gz -Oz -o rare_variants.vcf.gz```

        This creates a new VCF file (rare_variants.vcf.gz) containing only rare variants.

2. EXTRACT VARIANTS BASED ON GENE-MARKER MATCHES
    - Match variants in the VCF to markerID in gene.maker.annotation:

* Structure of gene.maker.annotation:
    gene: Gene names.
    markerID: Variant IDs (e.g., chromosomal positions like chr:pos).
    ann: Annotation information (e.g., pLoF, missense, etc.).


  Extract all markerID values from gene.maker.annotation:
   
    ```awk '{print $2}' gene.maker.annotation > marker_ids.txt```

    Use bcftools to filter the VCF for these markers:

    ```bcftools view -T marker_ids.txt rare_variants.vcf.gz -Oz -o rare_variants_in_genes.vcf.gz```

3. CLASSIFY VARIANTS INTO ANNOTATION GROUPS
    - Separate variants into Class 1 and Class 2 based on the ann column:
        Class 1: pLoF or severeMis.

       ```grep -E 'pLoF|missense' gene.maker.annotation | awk '{print $2}' > class1_marker_ids.txt```
      
       ``` bcftools view -T class1_marker_ids.txt rare_variants_in_genes.vcf.gz -Oz -o class1_variants.vcf.gz ```
      
        Class 2: moderateMis.

        ``` grep 'moderateMis' gene.maker.annotation | awk '{print $2}' > class2_marker_ids.txt```
   
        ```bcftools view -T class2_marker_ids.txt rare_variants_in_genes.vcf.gz -Oz -o class2_variants.vcf.gz```
      
5. IDENTIFY CARRIERS
    - Extract genotypes for all individuals in the filtered VCFs:
        For Class 1:

         ```bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' class1_variants.vcf.gz > class1_genotypes.txt```

        For Class 2:

         ```bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' class2_variants.vcf.gz > class2_genotypes.txt```
      
    - Identify carriers (carriers have genotypes 0/1, 1/0, or 1/1)

         ```awk '{for (i=5; i<=NF; i++) if ($i ~ /1/) print i}' class1_genotypes.txt > class1_carriers.txt```
      
         ```awk '{for (i=5; i<=NF; i++) if ($i ~ /1/) print i}' class2_genotypes.txt > class2_carriers.txt```

6. LINK CARRIERS TO PHENOTYPES
    - Match carrier IDs to the .ped file:
      Extract rows from the .ped file that correspond to the carriers:

        ```grep -Ff class1_carriers.txt input.ped > class1_carrier_phenotypes.txt```
      
        ```grep -Ff class2_carriers.txt input.ped > class2_carrier_phenotypes.txt```
   
7. ANALYZE RESULTS
    - Compare the burden of rare variants across classes and phenotypes:
      Count the number of carriers for Class 1 and Class 2:
    
        ```wc -l class1_carrier_phenotypes.txt```
        ``` wc -l class2_carrier_phenotypes.txt```

      Summarize the phenotypes of carriers:
      
        ```cut -f6 class1_carrier_phenotypes.txt | sort | uniq -c```
      
        ```cut -f6 class2_carrier_phenotypes.txt | sort | uniq -c```


    - Visualize results: Create bar plots or summary tables showing the burden of rare variants for each phenotype.






