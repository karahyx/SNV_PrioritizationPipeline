# :dna: PriorizationPipeline

An adaptation of the TCAG Small Variant Prioritization Pipeline originally developed to be used for Broad Institute's GATK (Genome Analysis Toolkit) to the Illumina DRAGEN (Dynamic Read Analysis for GENomics) Bio-IT Platform (https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/DRAGEN/GPipelineIntro_fDG.htm). 

## :green_book: Description

The PrioritizationPipeline script is used to annotate small variants called by DRAGEN with features including frequency and quality filters, damage type and ranks, as well as main and secondary findings for the purpose of variant prioritization. The script also generates pre-defined summary statistics such as chromosome-wise zygosity type counts and the number of variants in each category of interest. 

PrioritizationPipeline relies on the TCAG Small Variant Annotate Pipeline rev27.7 developed and maintained by Thomas. The documentation can be found here [TCAG_SMALL VARIANT_ANNOTATION_PIPELINE_rev27.7_hg38_JUN2022.pdf](./TCAG_SMALL_VARIANT_ANNOTATION_PIPELINE_rev27.7_hg38_JUN2022.pdf). In comparison, [the old prioritization script](./R/pipeline_old.R) relies on [TCAG_SMALL VARIANT_ANNOTATION_PIPELINE_rev27.4_hg18_AUG2022.pdf](./TCAG_SMALL_VARIANT_ANNOTATION_PIPELINE_rev27.4_hg18_AUG2022.pdf).

Two versions of PrioritizationPipeline are provided. [Version 1](./R/pipeline_new.R) is a general and more detailed version for users who wish to run the script on one file at a time and on their local machine. The individual sections that add different types of annotations and make up the Main Section are provided. As the file name suggests, [Version 2](./R/pipeline_new_hpf.R) is designed to be called by a shell script that receives user input and usually used to process multiple files in one go on HPF platforms. The shell scripts in the [bash](./bash) folder are used to perform such tasks. Detailed instructions for how to run the file on HPF can be found in the Instructions section below.

The main change from the old script is the utilization of functions that reduce the program to smaller, more manageable chunks and allow for reusability and extension. Detailed changes can be found in the Changes From the Old Script section below.

## :desktop_computer: Instructions

### Working with PrioritizationPipeline Version 1 (General Version)

Section | Name | Purpose | Subsections |
--- | --- | --- | --- |
0 | Variables & Cutoffs |  | N/A |
1 | Functions |  | N/A |
1.5 | File Import & Pre-processing |  | N/A |
2 | Frequency Filter |  | N/A |
3 | Quality Filter |  | N/A |
4 | Coding Tag |  | N/A |
5 | Define Damage |  | 5.0. Variable Initialization <br /> 5.1. Coding LOF <br /> 5.2. Missense <br /> 5.3. Other coding <br /> 5.4. Splicing predictions <br /> 5.5. UTR <br /> 5.6. Non-coding |
6 | Phenotype Filter |  | 6.1. HPO dominant <br /> 6.2. CGD dominant <br /> 6.3. Phenotype ranks |
7 | Main Findings | | 7.1. Recessive Homozygous <br /> 7.2. X-linked Haploid <br /> 7.3. Potential Compound Heterozygous <br /> 7.4. Dominant <br /> 7.5. Heterozygous Hotzone | 
8 | Secondary Findings | | 8.0. Pathogenicity flag <br /> 8.1. Rank 1 <br /> 8.1.0. Dominant, Pathogenic <br /> 8.1.1. Recessive, Homozygous, Pathogenic <br /> 8.1.2. Recessive, Potential Compound Heterozygous, Pathogenic <br /> 8.1.3. X-linked, Homozygous/Haploid, Pathogenic <br /> 8.1.4. Complex, Homozygous, Pathogenic <br /> 8.1.5. Complex, Potential Compound Heterozygous, Pathogenic <br /> 8.1.6. Complex, Single Heterozygous, Uncertain <br /> 8.1.7. Recessive, Single Heterozygous <br /> 8.1.8. X-linked, Heterozygous, Carrier <br /> 8.2. Rank 2 <br /> 8.3. Rank3 <br /> 8.4. ACMG Disease |
9 | Main | | |


### Working with PrioritizationPipeline Version 2 (HPF Version)


## Results
### Column Definitions

## Changes From the Old Script

## Contributions

## License
