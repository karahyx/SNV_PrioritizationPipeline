# :dna: PriorizationPipeline

An adaptation of the TCAG Small Variant Prioritization Pipeline originally developed to be used for Broad Institute's GATK (Genome Analysis Toolkit) to the Illumina DRAGEN (Dynamic Read Analysis for GENomics) Bio-IT Platform (https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/DRAGEN/GPipelineIntro_fDG.htm). 

## :green_book: Description

The PrioritizationPipeline script is used to annotate small variants called by DRAGEN with features including frequency and quality filters, damage type and ranks, as well as main and secondary findings for the purpose of variant prioritization. The script also generates pre-defined summary statistics such as chromosome-wise zygosity type counts and the number of variants in each category of interest. 

PrioritizationPipeline relies on the TCAG Small Variant Annotate Pipeline rev27.7 developed and maintained by Thomas. The documentation can be found here [TCAG_SMALL VARIANT_ANNOTATION_PIPELINE_rev27.7_hg38_JUN2022.pdf](./TCAG_SMALL_VARIANT_ANNOTATION_PIPELINE_rev27.7_hg38_JUN2022.pdf). In comparison, [the old prioritization script](./R/pipeline_old.R) relies on [TCAG_SMALL VARIANT_ANNOTATION_PIPELINE_rev27.4_hg18_AUG2022.pdf](./TCAG_SMALL_VARIANT_ANNOTATION_PIPELINE_rev27.4_hg18_AUG2022.pdf).

Two versions of PrioritizationPipeline are provided. [Version 1](./R/pipeline_new.R) is a general and more detailed version for users who wish to run the script on one file at a time and on their local machine. The individual sections that add different types of annotations and make up the Main Section are provided. As the file name suggests, [Version 2](./R/pipeline_new_hpf.R) is designed to be called by a shell script that receives user input and usually used to process multiple files in one go on HPF platforms. The shell scripts in the [bash](./bash) folder are used to perform such tasks. Detailed instructions for how to run the file on HPF can be found in the Instructions section below.

The main change from the old script is the utilization of functions that reduce the program to smaller, more manageable chunks and allow for reusability and extension. Detailed changes can be found in the Changes From the Old Script section below.

## :desktop_computer: Instructions

### Working with PrioritizationPipeline Version 1 (General Version)

#### Script Structure Overview

Section | Name | Content | Purpose | Column Definitions | 
--- | --- | --- | --- | --- |
0 | Variables & Cutoffs | 0.1. Input Variables <br /> 0.2. Output Variables <br /> 0.3. Internal Variables <br /> 0.4. Cutoffs <br /> 0.4.1. High-quality Filter <br /> 0.4.2. Define Damage <br /> 0.4.3. Main Findings | Modify file locations and cutoffs here | N/A |
1 | Functions | 1.1. Frequency Filter <br /> 1.2. Quality Filter <br /> 1.3. Coding Tag <br /> 1.4. Define Damage <br /> 1.5. Phenotype Filter <br /> 1.6. Main Findings <br /> 1.7. Secondary Findings <br /> 1.8. Get Final Results <br /> 1.8.1. Rare Variants <br /> 1.8.2. HQ Variants (FILTER == "PASS") <br /> 1.8.3. HQ Rare Variants <br /> 1.9. Stats | All functions used in the script can be found here | N/A |
1.5 | File Import & Pre-processing | N/A | Imports the original variant data and remove those variants with homozygous reference or unknown zygosity | N/A |
2 | Frequency Filter | freq_cutoff = 0.05 <br /> freq_cutoff = 0.01 <br /> freq_cutoff = 0.005 <br /> freq_cutoff = 0.0015 <br /> freq_cutoff = 0 | Add a frequency filter to filter for variants that pass a specific allele frequency cutoff | `F_Rare` = the smallest allele frequency that the variant passes |
3 | Quality Filter | 3.1. Pass tag <br /> 3.2. Quality tag | 3.1. Add a pass tag that indicates whether the variant has a "PASS" FILTER <br /> <br /> 3.2. Add a quality tag that indicates whether variants with "PASS" FILTER pass the DP cutoff | 3.1. `F_Pass` = whether the variant has a "PASS" FILTER <br /> <br /> 3.2. `F_Qual_tag` = <br /> <ul> <li> "OK" if the variant has DP greater than or equal to 2 </li> <li> "LowQuality" otherwise </li> |
4 | Coding Tag | Coding <br /> ncRNA <br /> Other | Add a coding tag that indicates whether the variant's type of sequence overlapped is coding, ncRNA or other types | `F_Coding` = <br /> <ul> <li> "Coding" if the variant's `typeseq_priority` is one of `exonic`, `exonic;splicing`, or `splicing` </li> <li> "ncRNA" if the variant's `typeseq_priority` is one of `ncRNA_exonic`, `ncRNA_splicing`, or `ncRNA_exonic;ncRNA_splicing` </li> <li> "Other" otherwise </li> |
5 | Define Damage | 5.0. Variable Initialization <br /> 5.1. Coding LOF <br /> 5.2. Missense <br /> 5.3. Other Coding <br /> 5.4. Splicing Predictions <br /> 5.5. UTR <br /> 5.6. Non-coding | 5.0. Initialize columns `F_DamageType = "NotDmg"` , `F_DamageRank = 0`, and `F_S_DamageType = "NotDmg"` <br /> <br /> 5.1.-5.5. Add specific damaging type tags | 5.1. The variant is Coding LOF if <br />  |
6 | Phenotype Filter |  | 6.1. HPO dominant <br /> 6.2. CGD dominant <br /> 6.3. Phenotype ranks |
7 | Main Findings | | 7.1. Recessive Homozygous <br /> 7.2. X-linked Haploid <br /> 7.3. Potential Compound Heterozygous <br /> 7.4. Dominant <br /> 7.5. Heterozygous Hotzone | 
8 | Secondary Findings | | 8.0. Pathogenicity flag <br /> 8.1. Rank 1 <br /> 8.1.0. Dominant, Pathogenic <br /> 8.1.1. Recessive, Homozygous, Pathogenic <br /> 8.1.2. Recessive, Potential Compound Heterozygous, Pathogenic <br /> 8.1.3. X-linked, Homozygous/Haploid, Pathogenic <br /> 8.1.4. Complex, Homozygous, Pathogenic <br /> 8.1.5. Complex, Potential Compound Heterozygous, Pathogenic <br /> 8.1.6. Complex, Single Heterozygous, Uncertain <br /> 8.1.7. Recessive, Single Heterozygous <br /> 8.1.8. X-linked, Heterozygous, Carrier <br /> 8.2. Rank 2 <br /> 8.3. Rank 3 <br /> 8.4. ACMG Disease |
9 | Main | | Step 1. File import <br /> Step 2. Re-format column names <br /> Step 3. Process the full data <br /> Step 4. Free up memory <br /> Step 5. Annotate data <br /> Step 6. Get chromosome counts and chromosome-wise zygosity counts <br /> Step 7. Get stats list for each data set <br /> Step 8. Convert stats lists to readable data frames <br /> Step 9. Get all stats in one data frame <br /> Step 10. Output results as .txt files |

#### Running the Script

### Working with PrioritizationPipeline Version 2 (HPF Version)

#### Script Structure Overview

#### Running the Script

### Stats Summary

## Changes From the Old Script

## Credits

## License
