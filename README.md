### PIPELINE FOR SUPPORTING THE TRACTOR LOCAL ANCESTRY GWAS ANALYSES 
(Based on Atkinson et al. Nat Genet 53, 195-204, 2021;
https://doi.org/10.1038/s41588-020-00766-y)

Developed by Danilo J. Carmona 

### Software Requirements
1. python 3.8
2. R 4.3
3. bcftools 1.10
4. RFMix 2.03
5. SHAPEIT 4.1.2

Pipeline was tested with: Python 3.8.18, R 4.3.2, bcftools 1.10.2, RFMix 2.03, SHAPEIT 4.1.2.

### Directory Structure before running

```bash
tractor/
│
├── *.vcf.gz                # Sample genotype file (input)
├── contigs.txt                # Chromosome harmonization file
├── *.txt              # Phenotype and covariate file
├── b38_sample_map.txt        # Ancestry mapping file (used in LAI step)
│
├── ref_1kGP/                  # Reference panel directory (must exist)
│   ├── ref_1.vcf.gz
│   ├── ref_2.vcf.gz
│   └── ...
│
├── gmap_1kGP/                 # Genetic maps directory (must exist)
│   ├── chr1.b38.gmap
│   ├── chr2.b38.gmap
│   └── ...
│
├── pre_phasing.sh
├── phasing.sh
├── la_inf.sh
├── extract.sh
├── run_tractor.sh
│
├── sel_anc.py
├── extract_tracts.py
├── run_tractor.R
```

### Execution order
1. pre_phasing.sh
2. phasing.sh
3. la_inf.sh
4. extract.sh
5. run_tractor.sh

### Information for the scripts in logic order of use:

i. File containing the script  
ii. Goal of the script  
iii. Usage example of the script  
iv. Input files  
v.  Output files  
vi. What is still missing  
```text
1. Pre-phasing:
    i. File: pre_phasing.sh
    ii. Goal: To generate chromosome-specific sample and reference BCF files containing only shared variants between both datasets, preparing them for downstream phasing.
    iii. Usage example: bash pre_phasing.sh -sample geno.vcf.gz -chr "1-22" -threads 8 -jobs 4 (threads and jobs optional)
    iv. Input:
        - *.vcf.gz sample file for which #CHROM column has the format 1 2 ...
        - reference ref_1kGP/ref_${chr}.vcf.gz
        - contigs.txt file used to harmonize chromosome naming between sample and reference
    v. Output:
        - geno_chr/geno_${chr}.bcf
        - geno_chr/geno_${chr}_com.bcf
        - ref_1kGP/ref_${chr}_com.bcf
    vi. Missing: 
        - flexibility in the extension of the input files
        - changing the format of the column CHR should be conditional

2. Phasing:
    i. File: phasing.sh
    ii. Goal: To use SHAPEIT4 for phasing the sample *.vcf file for specific chromosomes 
    iii. Usage example: bash phasing.sh -chr 1-5,20 -threads 8 -jobs 4 (threads and jobs optional)
    iv. Input:
        - geno_chr/geno_${chr}_com.bcf Sample genotype BCF file for each chromosome (output from pre-phasing step)
        - ref_1kGP/ref_${chr}_com.bcf Reference panel BCF file for each chromosome
        - gmap_1kGP/chr${chr}.b38.gmap Genetic map file with distances in Morgan units
    v. Output:
        - geno_ph/geno_${chr}_ph.bcf Phased genotype BCF file per chromosome
    vi. Missing: 
        - specify a condition for using chunks
        - flexibility in the extension of the input files
        - add an initial step checking the file is phased
        
3. Local Ancestry Inference (LAI):
    i. File: la_inf.sh
    ii. Goal: To infer local ancestries from the reference by using rfmix
    iii. Usage example: bash la_inf.sh -ancs amr,eur -chr 1-5,21 -threads 8 -jobs 4 (threads and jobs optional)
    iv. Input: 
        - sample geno_ph/geno_${chr}_ph.bcf file(s) with phased variants matching with reference data at specified chromosomes (output from phasing step)
        - reference ref_1kGP/ref_${chr}_com.bcf file with variants matching with sample data at specified chromosomes (output from phasing step)
        - ancestry b38_sample_map.txt file with a column of sample_id's and another with ancestries for the reference panel
        - gmap_1kGP/chr${chr}.b38.gmap file mapping the distances between variants in Morgan units
    v. Output:
        - dec_chr/dec_${chr}.msp.tsv file with the tracts at specified chromosomes
    vi. Missing:
        - specify a condition for using chunks
        - flexibility in the extension of the input files
        - add an input with the general ancestries file for building the b38_${ancs}_map.txt file
        
4. Extracting tracts:
    i. File: extract.sh
    ii. Goal: To generate variables for each individual at each variant (hapcounts and dosages) on separated files 
    iii. Usage example: bash extract.sh -ancs amr,eur -chr 1-5,21 -jobs 22 (jobs optional)
    iv. Input:
        - sample geno_ph/geno_${chr}_ph.* (bcf, vcf or vcf.gz) file with phased variants matching with reference data at specified chromosomes (same input as for LAI step)
        - dec_chr/dec_${chr}.msp.tsv file with the tracts at specified chromosomes (output from LAI step)
    v. Output:
        - sample geno_${chr}_ph.anc${anc}hapcount.txt file with haplotype counts for each individual at each snp
        - sample geno_${chr}_ph.anc${anc}dosage.txt file with ancestry dosages for each individual at each snp
    vi. Missing:
        - flag to specify name of input files 
        - flag to specify number of threads to be used per chromosome
        - specify a condition for using chunks
        - flexibility in the extension of the input files
    
5. Association study:
    i. File: run_tractor.sh
    ii. Goal: To perform the GWAS by including hapcount and dosage variables (tractor-GWAS) 
    iii. Usage example: bash run_tractor.sh \
                          -chr 1-10,20 \
                          -jobs 4 \
                          -threads 10 \
                          -phenofile covarlist.txt \
                          -covarcollist "sex,age,rnd_source,pc_1,pc_2,pc_3,pc_4,pc_5" \
                          -method logistic \
                          -sampleidcol id \
                          -phenocol ad
    iv. Input:
        - geno_${chr}_ph.* hapdose prefix files generated from previous steps (per chromosome)
        - covarlist.txt phenotype + covariate file containing: sample IDs, phenotype column, specified covariate columns
    v. Output:
        - geno_out/geno_${chr}_out file containing the statistical summary of the performed GWAS
    vi. Missing:
        - specify a condition for using chunks
        - flexibility in the extension of the input files

=^.^=


