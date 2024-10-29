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




# GRM

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



