# :dna: PrioritizationPipeline

An adaptation of the TCAG Small Variant Prioritization Pipeline originally developed to be used for Broad Institute's [GATK (Genome Analysis Toolkit)](https://gatk.broadinstitute.org/hc/en-us) to the [Illumina DRAGEN (Dynamic Read Analysis for GENomics) Bio-IT Platform](https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/DRAGEN/GPipelineIntro_fDG.htm). 

## :green_book: Description

The PrioritizationPipeline script is used to annotate small variants called by DRAGEN with features including frequency and quality filters, damage type and tiers, as well as main and secondary findings for variant prioritization. The script also generates pre-defined summary statistics such as chromosome-wise zygosity type counts and the number of variants in each category of interest. 

Two versions of Pipeline_v17 are available. [Version 1](./R/Pipeline_v17_Version1_IND) is used to process single samples, while [Version 2](./R/Pipeline_v17_Version2_FAM) provides the workflow for family-based studies. Version 2 adds the same filters as Version 1, but reads in pedigrees to identify compound heterozygotes (Version 1 only identifies potential compound heterozygotes). The shell scripts in the [bash](./bash) folder are used to receive user input and process one or multiple files at once on high-performance computing (HPC) systems using the job scheduler SLURM.

The main changes from the previous version include updates of cutoffs and definitions, minor error fixes, and the utilization of functions that reduce the program to smaller, more manageable chunks that allow for reusability and extension. Detailed changes can be found in the **Changes From the Old Script** section below.

## :desktop_computer: Instructions

### :brain: Working with PrioritizationPipeline v17 Version 1 (Single Samples)

#### :mag_right: Script Structure Overview

##### :seedling: Main Script

<table>
<thead>
  <tr>
    <th>Section</th>
    <th>Content</th>
    <th>Purpose</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td rowspan="3"> (0) Settings </td>
    <td> 0.1. Libraries </td>
    <td> Import all libraries used in the script </td>
  </tr>
  <tr>
    <td> 0.2. Input variables </td>
    <td> <code>args[1]</code>: root path <br>
         <code>args[2]</code>: sub-path and filename of second R code file defining functions <br>
         <code>args[3]</code>: variant input file path <br>
         <code>args[4]</code>: genome name used at columns <br>
         <code>args[5]</code>: sub-path and filename of stats visualization Rmd file <br>
         <code>args[6]</code>: output sub-path </td>
  </tr>
  <tr>
    <td> 0.3. Functions </td>
    <td> Import functions used in the script </td>
  </tr>
  <tr>
    <td rowspan="5">(1) Variables & Cutoffs</td>
    <td> 1.1. Internal variables </td>
    <td rowspan="5"> Modify internal variables and cutoffs here </td>
  </tr>
  <tr>
    <td> 1.2. Cutoffs </td>
  </tr>
  <tr>
    <td> 1.2.1. High-quality filter </td>
  </tr>
  <tr>
    <td> 1.2.2. Define damage </td>
  </tr>
  <tr>
    <td> 1.2.3. Main findings </td>
  </tr>
  <tr>
    <td rowspan="12"> (2) Main </td>
    <td> 2.1. File import </td>
    <td> Imports the original variant data </td>
  </tr>
  <tr>
    <td> 2.1.1. Check if multi-sample </td>
    <td> Check if the imported variant data contains only one sample or multiple samples </td>
  </tr>
  <tr>
    <td> 2.2. Re-format column names</td>
    <td> Remove <code>{genome_name}:</code> from several columns for easier access and processing </td>
  </tr>
  <tr>
    <td> 2.3. Process the original imported variant data </td>
    <td> 
      <ul>
        <li> Remove variants with homozygous reference (hom-ref) or unknown zygosity from the data </li>
        <li> Calculate the alternate allele frequency </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 2.3.1. Add some fields for allele frequency </td>
    <td> 
      <ul>
          <li> <code>FreqMaxSimple_AfrAmrEasNfeSasOth</code> = the maximum of a variant's allele frequencies (exomes and genomes) from the African-American/African, Latino/Admixed American, East Asian, Non-Finnish European, South Asian, and other subsets </li>
          <li> <code>FreqHomCount_AfrAmrEasNfeSasOth</code> = the maximum between the counts of homozygous individuals in samples from the gnomAD exome dataset and from the gnomAD genome dataset </li>
          <li> <code>dbsnp_region_count</code> = a variant's total number of overlap-based match for dbSNP </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 2.3.2. Add some fields for gene constraint </td>
    <td> <code>F_GeneConstr</code> = 
        <ul>
            <li> 1, if a variant's LOF observed/expected (oe) metric from the gnomAD constraint matrix is less than or equal to 0.15 or its missense oe metric from the gnomAD constraint matrix is less than or equal to 0.70 (i.e. <code>gnomad_oe_lof &le; 0.15 | gnomad_oe_mis &le; 0.70</code>) </li>
            <li> 0, otherwise </li>
        </ul>
    </td>
  </tr>
  <tr>
    <td> 2.4. Free up memory</td>
    <td> Remove <code>v_full.temp.df</code> from the current workspace </td>
  </tr>
  <tr>
    <td> 2.5. Add filters </td>
    <td> 
      <ul>
        <li> Remove alternate contigs and unlocalized/unplaced sequence from the data </li>
        <li> Obtain variants with a maximum frequency of 0.05 and annotate them with filtering tags </li>
        <li> Obtain high-quality variants that have a "PASS" FILTER and annotate them with filtering tags </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 2.6. Get chromosome counts and chromosome-wise zygosity counts </td>
    <td> 
      <ul>
        <li> Obtain the number of variants in each chromosome </li> 
        <li> Obtain the number of alt-alt, hom-alt, ref-alt, and the percentage of hom-alt in each chromosome </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 2.7. Get all variant stats </td>
    <td> Obtain the summary statistics for all data sets in one table </td>
  </tr>
    <tr>
    <td> 2.8. Process the results before outputting </td>
    <td> 
        <ul>
            <li> Change tier levels from 0, 1 and 2 to "-", "Low" and "High" in columns <code>F_DamageTier</code> and <code>F_PhenoTier</code> </li>
            <li> Add <code>{genome_name}:</code> back to the sample columns </li>
        </ul>
    </td>
  </tr>
  <tr>
    <td> 2.9. Output desired results as .txt files</td>
    <td> Output the following to the user-defined output directory:
        <ul>
            <li> Rare 5% variants </li>
            <li> Rare 5% damaging variants </li>
            <li> Rare 5% secondary findings variants </li> 
            <li> Chromosome zygosity counts tables </li>
            <li> A summary statistics table </li> 
            <li> A PDF file that visualizes the summary statistics </li>
        </ul>
    </td>
  </tr>
</tbody>
</table>

##### :seedling: Filter Functions

<table>
<thead>
  <tr>
    <th>Section</th>
    <th>Content</th>
    <th>Purpose</th>
    <th>Column Definitions</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td rowspan="11"> (3) Overview <br> </td>
    <td> 3.0. Alternate contigs and unlocalized/unplaced sequence filter </td>
    <td rowspan="11"> All functions used in the script can be found here </td>
    <td rowspan="11"> N/A </td>
  </tr>
  <tr>
    <td> 3.1. Frequency filter </td>
  </tr>
  <tr>
    <td> 3.2. Quality filter </td>
  </tr>
  <tr>
    <td> 3.3. Coding tag </td>
  </tr>
  <tr>
    <td> 3.4. Define damage </td>
  </tr>
  <tr>
    <td> 3.5. Phenotype filter </td>
  </tr>
  <tr>
    <td> 3.6. Main findings </td>
  </tr>
  <tr>
    <td> 3.7. Secondary findings </td>
  </tr>
  <tr>
    <td> 3.8. Get final results </td>
  </tr>
  <tr>
    <td> 3.9. Stats </td>
  </tr>
  <tr>
    <td> 3.10. Change tier names to "Low" and "High" </td>
  </tr>
  <tr>
    <td> 3.0. Alternate Contigs and Unlocalized/Unplaced Sequence Filter</td>
    <td> N/A </td>
    <td> Remove <a href="https://gatk.broadinstitute.org/hc/en-us/articles/360041155232-Reference-Genome-Components"> alternate contigs and unlocalized/unplaced sequence</a> from the variant data </td>
    <td> N/A </td>
  </tr>
  <tr>
    <td rowspan="5"> 3.1. Frequency Filter </td>
    <td> allele frequency cutoff = 0.05 </td>
    <td rowspan="5"> Add a frequency filter that filters for variants that pass a specific allele frequency cutoff </td>
    <td rowspan="5"> <code>F_Rare</code> = <br> the smallest allele frequency cutoff that the variant passes </td>
  </tr>
  <tr>
    <td> allele frequency cutoff = 0.01 </td>
  </tr>
  <tr>
    <td> allele frequency cutoff = 0.001 </td>
  </tr>
  <tr>
    <td> allele frequency cutoff = 0.0001 </td>
  </tr>
  <tr>
    <td> allele frequency cutoff = 0 </td>
  </tr>
  <tr>
    <td rowspan="3"> 3.2. Quality Filter </td>
    <td> 3.2.1. Pass tag </td>
    <td> Add a pass tag that indicates whether the variant has FILTER = "PASS" </td>
    <td> <code>F_Pass</code> = 
      <ul>
        <li> 1, if the variant has a "PASS" FILTER </li>
        <li> 0, otherwise </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td rowspan="2"> 3.2.2. Quality tag </td>
    <td> If single sample, add a quality tag that indicates whether variants with a "PASS" FILTER pass the DP cutoff </td>
    <td> <code>F_Qual</code> = 
      <ul>
        <li> 1, if the variant has a "PASS" FILTER and DP &ge; 2 </li>
        <li> 0, otherwise </li>
      </ul>
    </td>
  </tr>
  <tr>
      <td> If multi-sample, add a quality tag that indicates whether variants with a "PASS" FT (i.e. sample genotype filter indicating if this genotype was "called") pass the DP cutoff  </td>
      <td> <code>F_Qual</code> = 
        <ul>
          <li> 1, if the variant has a "PASS" FT and DP &ge; 2 </li>
          <li> 0, otherwise </li>
        </ul>
    </td>
  </tr>
  <tr>
    <td> 3.3. Coding Tag </td>
    <td> Add coding tags to the variant data </td>
    <td> Add a coding tag that indicates whether a variant's type of sequence overlapped with respect to known genes/transcripts is Coding, ncRNA, or Other </td>
    <td> <code>F_Coding</code> = 
      <ul>
        <li> "Coding" if the variant's <code>typeseq_priority</code> is one of <code>exonic</code>, <code>exonic;splicing</code>, or <code>splicing</code> </li>
        <li> "ncRNA" if the variant's <code>typeseq_priority</code> is one of <code>ncRNA_exonic</code>, <code>ncRNA_splicing</code>, or <code>ncRNA_exonic;ncRNA_splicing</code> </li>
        <li> "Other" otherwise </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td rowspan="7"> 3.4. Define Damage </td>
    <td> 3.4.0. Variable initialization <br><br> :point_up:can be found in the <code>get_rare05_variants()</code> function </td>
    <td> Initialize the following columns:<br><code>F_DamageType = "NotDmg"</code><br><code>F_DamageTier = 0</code> </td>
    <td> <code>F_DamageType</code> = a variant's damage type
      <ul>
        <li> one of <code>LOF</code>, <code>Missense</code>, <code>OtherC</code>, <code>Splc</code>, <code>UTR</code>, <code>DmgNcRNA</code>, or <code>NotDmg</code> </li>
      </ul>
      <br> <code>F_DamageTier</code> = a variant's damage tier &isin; {0, 1, 2}
        <ul>
          <li> Note that the higher the tier, the more damaging a variant is </li>
        </ul>
    </td>
  </tr>
  <tr>
    <td> 3.4.1. Coding LOF </td>
    <td> <code>add_coding_lof_tag()</code> identify variants of type <code>Coding LOF</code> and change their <code>F_DamageTier</code> to 1 based on specific conditions </td>
    <td>The variant is <code>Coding LOF</code> and has <code>F_DamageTier = 1</code> if it
      <ul>
        <li> is coding </li>
        <li> causes frameshift or point mutations in the coding sequence; or its type of sequence overlapped is "splicing" </li>
        <li> has <code>distance_spliceJunction &le; 2</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.4.2. Missense </td>
    <td> Identify variants of type <code>Missense</code> and change their <code>F_DamageTier</code> to either 1 or 2 based on specific conditions
    </td>
    <td> The variant is <code>Missense</code> and has <code>F_DamageTier = 1</code> if it
      <ul>
        <li> is coding </li> 
        <li> is a nonsynonymous SNV </li>
        <li> satiesfies one or more one of the following
          <ul>
            <li> <code>REVEL_score &ge; 0.25 </code> </li>
            <li> <code>phylopMam_avg &ge; 1.3</code> & <code>phylopVert100_avg &ge; 3.9</code> </li>
            <li> <code>CADD_phred &ge; 30</code> or <code>MPC_score &ge; 2</code> </li>
          </ul>
        </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.4.3. Other coding </td>
    <td> Identify variants of type <code>Other coding</code> and change their <code>F_DamageTier</code> to either 1 or 2 based on specific cutoffs </td>
    <td> The variant is <code>Other coding</code> if it
      <ul>
        <li> is coding </li> 
        <li> causes nonframeshift mutations in the coding sequence </li>
        <li> is not an exact match to common dbSNP track UCSC, and <code>phylopMam_avg &ge; phylopMam_cutoff</code> or <code>phylopVert100_avg &ge; phylopVert_cutoff</code> or <code>CADD_phred &ge; CADD_phred_cutoff</code> </li>
      </ul><br>
      <code>F_DamageTier = 1</code> if 
      <ul>
        <li> <code>phylopMam_cutoff = 1.1</code> </li>
        <li> <code>phylopVert_cutoff = 1.6</code> </li>
        <li> <code>CADD_phred_cutoff = 13.7</code> </li>
      </ul><br>
      <code>F_DamageTier</code> = 2 if 
      <ul>
        <li> <code>phylopMam_cutoff = 1.3</code> </li>
        <li> <code>phylopVert_cutoff = 3.9</code> </li>
        <li> <code>CADD_phred_cutoff = 21.1</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.4.4. Splicing predictions </td>
    <td> Predict whether the variant is of type <code>Splicing</code> and change their <code>F_DamageTier</code> to either 1 or 2 based on specific conditions </td>
    <td> The variant is predicted to be <code>Splicing</code> and has <code>F_DamageTier = 1</code> if
      <ul>
        <li> its <code>F_DamageType</code> is not "LOF", "Splc" or "Missense" </li>
        <li> its <code>F_DamageTier &ne; 2</code></li>
        <li> satisfies one or more of the following: </li>
          <ul>
            <li> <code>spliceAI_DS_AG > 0.2</code> & <code>|spliceAI_DP_AG| &le; 50</code> </li>
            <li> <code>spliceAI_DS_AL > 0.2</code> & <code>|spliceAI_DP_AL| &le; 50</code> </li>
            <li> <code>spliceAI_DS_DG > 0.2</code> & <code>|spliceAI_DP_DG| &le; 50</code> </li>
            <li> <code>spliceAI_DS_DL > 0.2</code> & <code>|spliceAI_DP_DL| &le; 50</code> </li>
            <li> has variant type "del" (deletion) or "ins" (insertion) and its type of sequence overlapped is "splicing" or "exonic;splicing" </li>
          </ul>
      </ul> <br>
      The variant is predicted to be <code>Splicing</code> and has <code>F_DamageTier = 2</code> if
      <ul>
        <li> its <code>F_DamageType</code> is not "LOF" </li>
        <li> satisfies one or more of the following: </li>
          <ul>
            <li> <code>spliceAI_DS_AG > 0.5</code> & <code>|spliceAI_DP_AG| &le; 50</code> </li>
            <li> <code>spliceAI_DS_AL > 0.5</code> & <code>|spliceAI_DP_AL| &le; 50</code> </li>
            <li> <code>spliceAI_DS_DG > 0.5</code> & <code>|spliceAI_DP_DG| &le; 50</code> </li>
            <li> <code>spliceAI_DS_DL > 0.5</code> & <code>|spliceAI_DP_DL| &le; 50</code> </li>
            <li> <code>dbscSNV_ADA_SCORE > 0.6</code> & <code>dbscSNV_RF_SCORE > 0.6</code> </li>
          </ul>
      </ul>
    </td>
  </tr>
  <tr>
    <td>3.4.5. UTR</td>
    <td> Identify variants of type <code>UTR</code> and change their <code>F_DamageTier</code> to 1 or 2 based on specific cutoffs </td>
    <td> The variant is <code>UTR</code> if 
      <ul>
        <li> the type of sequence overlapped with respect to known genes/transcripts is "UTR3", "UTR5", or both </li>
        <li> has a PhastCons score for the Placental Mammal genome group </li>
        <li> satisfies <code>phylopMam_avg &ge; phylopMam_cutoff</code> or <code>CADD_phred &ge; CADD_phred_cutoff</code> </li>
      </ul> <br>
      <code>F_DamageTier = 1</code> if 
        <ul>
          <li> <code>phylopMam_cutoff = 1.1</code> </li>
          <li> <code>CADD_phred_cutoff = 13.7</code> </li>
        </ul> <br>
      <code>F_DamageTier = 2</code> if 
        <ul>
          <li> <code>phylopMam_cutoff = 1.3</code> </li>
          <li> <code>CADD_phred_cutoff = 21.1</code> </li>
        </ul>
    </td>
  </tr>
  <tr>
    <td>3.4.6. Non-coding</td>
    <td> Identify variants of type <code>Non-coding</code> and change their <code>F_DamageTier</code> to either 1 or 2 based on specific cutoffs </td>
      <td> The variant is <code>Non-coding</code> if 
        <ul>
          <li> it is identified as "ncRNA" in <code>F_Coding</code> </li>
          <li> its full gene name is not "pseudogene" </li>
          <li> its <code>F_DamageTier = 0</code> (i.e. it's not damaging) </li>
          <li> satisfies one or more of the following: </li>
            <ul>
              <li> has a PhastCons score for the Placental Mammal genome group AND, <code>phylopMam_avg &ge; phylopMam_cutoff</code> OR <code>phylopVert100_avg &ge; phylopVert_cutoff</code> </li>
              <li> <code>CADD_phred &ge; CADD_phred_cutoff</code> </li>
            </ul>
        </ul> <br>
        <code>F_DamageTier = 1</code> if 
          <ul>
            <li> <code>phylopMam_cutoff = 1.1</code> </li>
            <li> <code>phylopVert_cutoff = 1.6</code> </li>
            <li> <code>CADD_phred_cutoff = 13.7</code> </li>
          </ul> <br>
        <code>F_DamageTier = 2</code> if 
        <ul>
            <li> <code>phylopMam_cutoff = 1.3</code> </li>
            <li> <code>phylopVert_cutoff = 3.9</code> </li>
            <li> <code>CADD_phred_cutoff = 21.1</code> </li>
          </ul>
      </td>
  </tr>
  <tr>
    <td rowspan="2"> 3.5. Phenotype Filter </td>
    <td> 3.5.1. OMIM dominant <br> 3.5.2. CGD dominant </td>
    <td> Add a tag <code>F_AXD</code> that indicates whether the variant has autosomal dominant (AD) as their mode of inheritance based on the OMIM and CGD annotations </td>
    <td> <code>F_AXD</code> = 
      <ul>
        <li> 1, if the pattern "i:AD" is found in column <code>omim_phenotype</code> or the pattern "AD" is found in column <code>CGD_inheritance</code> </li>
        <li> 0, otherwise </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.5.3. Phenotype tiers </td>
    <td> Add a tag <code>F_PhenoTier</code> that indicates the phenotype tier of a variant. Here, tier 1 and 2 are used to differentiate mouse and human phenotype annotations, respectively. Human phenotype annotations take priority. </td>
    <td> <code>F_PhenoTier</code> =
      <ul>
        <li> 1, if the variant has an MPO annotation imported from MGI (Mouse Genome Informatics) and mapped from an orthologous mouse gene </li>
        <li> 2, if the variant has one or more annotations from OMIM, HPO, or CGD </li>
        <li> 0, otherwise </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td rowspan="5"> 3.6. Main Findings </td>
    <td> 3.6.1. Recessive homozygous </td>
    <td> Add a tag <code>FM_HOM</code> that indicates whether the variant is recessive homozygous </td>
    <td> <code>FM_HOM = 1</code> if the variant
      <ul>
        <li> is a rare variant with a maximum allele frequency of 0.05 </li>
        <li> is damaging, i.e. <code>F_DamageType</code> != "NotDmg" </li>
        <li> has zygosity <code>hom-alt</code> (homozygous alternative) </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.6.2. X-linked haploid </td>
    <td> Add a tag <code>FM_XHAP</code> that indicates whether the variant is an X-linked haploid </td>
    <td> <code>FM_XHAP = 1</code> if the variant
      <ul>
        <li> is a rare variant with a maximum allele frequency of 1e-4 or 0.0001 </li>
        <li> is damaging, i.e. <code>F_DamageType</code> != "NotDmg" </li>
        <li> has zygosity <code>hap</code> (homozygous alternative) </li>
        <li> is found in chromosome X </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.6.3. Potential compound heterozygous </td>
    <td> <code>add_potential_compound_heterozygous_tag()</code> adds a tag <code>FM_PCHET</code> that indicates whether the variant is a potential compound heterozygote. <br>
    <br> <code>add_potential_dmg_compound_heterozygous_tag()</code> adds a tag <code>FM_PCHET_DMG</code> that indicates whether the variant is a damaging potential compound heterozygote. <br>
    <br> Note that the method used involves looking for multiple mutations on the same gene.  </td>
    <td> <code>FM_PCHET</code> = 
      <ul>
        <li> 1, if the variant has a maximum allele frequency of 0.05 and is damaging </li>
        <li> 2, if the variant has a maximum allele frequency of 0.05, is damaging, and is high quality (i.e. has a "PASS" FILTER and DP &ge; 2) </li>
        <li> 0, otherwise </li>
      </ul>
      <br> A variant has <code>FM_PCHET_DMG = 1</code> if it
      <ul>
        <li> is labelled as potential compound heterozygous, i.e. <code>FM_PCHET = 2</code> </li>
        <li> has <code>F_DamageTier = 2</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.6.4. Dominant </td>
    <td> Add a tag <code>FM_AXDOM</code> that indicates whether the variant is autosomal dominant (AD) </td>
    <td> <code>FM_AXDOM = 1</code> if the variant
      <ul>
        <li> is a rare variant with a maximum allele frequency of 1e-4 or 0.0001 </li>
        <li> is damaging, i.e. <code>F_DamageType</code> != "NotDmg" </li>
        <li> has <code>F_AXD == 1</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.6.5. Predicted dominant </td>
    <td> Add a tag <code>FM_PDDOM</code> that indicates whether the variant is predicted to be haploinsufficient, dominant gain-of-function, or dominant negative </td>
    <td> <code>FM_PDDOM = 1</code> if the variant
      <ul>
        <li> is a rare variant with a maximum allele frequency of 1e-4 or 0.0001 </li>
        <li> has zygosity <code>ref-alt</code> (heterozygous reference) </li>
        <li> has <code>F_GeneConstr == 1</code> </li>
        <li> satisfies one of the following: </li>
          <ul>
            <li> <code>F_DamageType</code> is one of "LOF", "Splc", or "OtherC" </li>
            <li> <code>F_DamageType</code> is "Missense" and <code>F_DamageTier = 2</code> </li>
          </ul>
      </ul>
    </td>
  </tr>
  <tr>
    <td rowspan="14"> 3.7. Secondary Findings </td>
    <td> 3.7.0. Pathogenicity flag </td>
    <td> Add two pathogenicity related flags:
      <ul> 
        <li> flag <code>F_Clinvar_Pathg</code> indicates whether the variant has at least one record submitted with pathogenic or likely pathogenic based on ClinVar </li>
        <li> flag <code>F_Clinvar_notPathg</code> indicates whether the variant has no current value of pathogenicity based on ClinVar </li>
      </ul>
    </td>
    <td> <code>F_Clinvar_Pathg = 1</code> if <code>Clinvar_SIG_Submission</code> is one of <code>Pathogenic</code>, <code>Pathogenic\x2c_low_penetrance</code>, <code>Likely_pathogenic</code>, <code>Likely_pathogenic\x2c_low_penetrance</code>, <code>Established_risk_allele</code>, <code>Likely_risk_allele</code>, or <code>risk_factor</code>
    <br><br> <code>F_Clinvar_notPathg = 1</code> if <code>Clinvar_SIG_Submission</code> is not NA or one of <code>Pathogenic</code>, <code>Pathogenic\x2c_low_penetrance</code>, <code>Likely_pathogenic</code>, <code>Likely_pathogenic\x2c_low_penetrance</code>, <code>Established_risk_allele</code>, <code>Likely_risk_allele</code>, or <code>risk_factor</code>
    </td>
  </tr>
  <tr>
    <td> 3.7.1. Tier 1 </td>
    <td> Add a tag <code>FS1_Select</code> that indicates whether a variant belongs to tier 1 for secondary findings </td>
    <td> <code>FS1_Select = 1</code> if the variant
      <ul>
        <li> has CGD disease annotations </li>
        <li> satisfies one or more of the following: </li>
          <ul>
            <li> is "Coding LOF" with <code>distance_spliceJunction < 3</code>, i.e. <code>F_S_DamageType = "LOF"</code> </li>
            <li> indicated as pathogenic or likely pathogenic by ClinVar </li>
          </ul>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.2. Tier 2 </td>
    <td> Add a tag <code>FS2_Select</code> that indicates whether a variant belongs to tier 2 for secondary findings </td>
    <td> <code>FS2_Select = 1</code> if the variant
      <ul>
        <li> has CGD disease annotations </li>
        <li> is high quality, i.e. has a "PASS" FILTER and DP &ge; 2 </li>
        <li> has a maximum allele frequency of 0.01 </li>
        <li> satisfies one or more of the following: </li>
          <ul>
            <li> is of type Coding LOF, i.e. <code>F_DamageType = "LOF"</code></li>
            <li> indicated as pathogenic or likely pathogenic by ClinVar, i.e. <code>F_Clinvar_Pathg = 1</code> </li>
            <li> not indicated as no current value of pathogenic by ClinVar, i.e. <code>F_Clinvar_notPathg = 0</code> 
          </ul>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.3. Tier 3 </td>
    <td> Add a tag <code>FS3_Select</code> that indicates whether a variant belongs to tier 3 for secondary findings </td>
    <td> <code>FS3_Select = 1</code> if the variant
      <ul>
        <li> has CGD disease annotations </li>
        <li> has a "PASS" FILTER and DP &ge; 2 </li>
        <li> has a maximum allele frequency of 0.01 </li>
        <li> is damaging </li>
        <li> satisfies one or more of the following: </li>
          <ul>
            <li> is of type Coding LOF, i.e. <code>F_DamageType = "LOF"</code></li>
            <li> indicated as pathogenic or likely pathogenic by ClinVar, i.e. <code>F_Clinvar_Pathg = 1</code> </li>
            <li> not indicated as no current value of pathogenic by ClinVar, i.e. <code>F_Clinvar_notPathg = 0</code> 
          </ul>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.4. Dominant, pathogenic </td>
    <td> Add a tag <code>FS1_AD_Pathg_Any</code> that indicates whether a variant is dominant and pathogenic </td>
    <td> A variant has <code>FS1_AD_Pathg_Any = 1</code> if
      <ul>
        <li> it belongs to tier 1 for secondary findings, i.e. <code>FS1_Select = 1</code> </li>
        <li> it has autosomal dominant (AD) as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.5. Recessive, homozygous, pathogenic </td>
    <td> Add a tag <code>FS1_AR_Pathg_Hom</code> that indicates whether a variant is recessive, homozygous and pathogenic </td>
    <td> A variant has <code>FS1_AR_Pathg_Hom</code> if it
      <ul>
        <li> belongs to tier 1 in secondary findings, i.e. <code>FS1_Select = 1</code> </li>
        <li> has autosomal recessive (AR) as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> has zygosity "homozygous alternative" (hom-alt) </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.6. Recessive, potential compound heterozygous, pathogenic </td>
    <td> Two tags are being added here:
      <ul>
        <li> Add a tag <code>F_CmpHet_S1</code> that indicates whether a variant is a potential compound heterozygote </li>
        <li> Add a tag <code>FS1_AR_Pathg_PotCompHet</code> that indicates whether a variant is recessive, potential compound heterozygous, and pathogenic </li>
      </ul>
    </td>
    <td>
      <code>F_CmpHet_S1</code> =
      <ul>
        <li> 1, if it has a "PASS" FILTER and belongs to tier 1 for secondary findings </li>
        <li> 2, if it has a "PASS" FILTER and DP &ge; 2 and belongs to tier 1 for secondary findings </li>
        <li> 0, otherwise </li>
      </ul><br>
      A variant has <code>FS1_AR_Pathg_PotCompHet = 1</code> if it
      <ul>
        <li> belongs to tier 1 for secondary findings </li>
        <li> has autosomal recessive (AR) as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> has <code>F_CmpHet_S1 &ge; 1</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.7. X-linked, homozygous/haploid, pathogenic </td>
    <td> Add a tag <code>FS1_XL_Pathg_Hom</code> that indicates whether a variant is X-linked homozygous and pathogenic <br><br>
    Add a tag <code>FS1_XL_Pathg_Hap</code> that indicates whether a variant is X-linked haploid and pathogenic </td>
    <td> <code>FS1_XL_Pathg_Hom = 1</code> or <code>FS1_XL_Pathg_Hap = 1</code> if a variant
      <ul>
        <li> belongs to tier 1 for secondary findings </li>
        <li> has X-linked (XL) as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.8. Complex, homozygous/haploid, pathogenic </td>
    <td> Add a tag <code>FS1_CX_Pathg_HomHap</code> that indicates whether a variant is complex, homozygous/haploid, and pathogenic </td>
    <td> A variant has <code>FS1_CX_Pathg_HomHap = 1</code> if it
      <ul>
        <li> belongs to tier 1 for secondary findings </li>
        <li> does not have AD, AR, or XL as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> has zygosity <code>hom-alt</code> (homozygous alternative) </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.9. Complex, potential compound heterozygous, pathogenic </td>
    <td> Add a tag <code>FS1_CX_Pathg_PotCompHet</code> that indicates whether a variant is complex, potential compound heterozygous, and pathogenic </td>
    <td> A variant has <code>FS1_CX_Pathg_HomHap = 1</code> if it
      <ul>
        <li> belongs to tier 1 for secondary findings </li>
        <li> does not have AD, AR, or XL as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> has <code>F_CmpHet_S1 &ge; 1</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.10. Complex, single heterozygous, uncertain </td>
    <td> Add a tag <code>FS1_CX_Uncertain</code> that indicates whether a variant is complex, single heterozygous, and pathogenicity uncertain </td>
    <td> A variant has <code>FS1_CX_Uncertain = 1</code> if it
      <ul>
        <li> belongs to tier 1 for secondary findings </li>
        <li> does not have AD, AR, or XL as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> is not potential compound heterozygous, i.e. has <code>F_CmpHet_S1 = 0</code> </li>
        <li> has zygosity <code>ref-alt</code> (heterozygous reference) or <code>alt-alt</code> (heterozygous alternate) </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.11. Recessive, single heterozygous, carrier </td>
    <td> Add a tag <code>FS1_AR_Carrier</code> that indicates whether a variant is recessive, single heterozygous, and a carrier </td>
    <td> A variant has <code>FS1_AR_Carrier = 1</code> if it
      <ul>
        <li> belongs to tier 1 for secondary findings </li>
        <li> has AR as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> is not potential compound heterozygous, i.e. has <code>F_CmpHet_S1 = 0</code> </li>
        <li> has zygosity <code>ref-alt</code> (heterozygous reference) or <code>alt-alt</code> (heterozygous alternate) </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.12. X-linked, heterozygous, carrier </td>
    <td> Add a tag <code>FS1_XL_Carrier</code> that indicates whether a variant is X-linked, heterozygous, and a carrier </td>
    <td> A variant has <code>FS1_XL_Carrier = 1</code> if it
      <ul>
        <li> belongs to tier 1 for secondary findings </li>
        <li> has XL as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> is found in chromosome X </li>
        <li> has zygosity <code>ref-alt</code> (heterozygous reference) or <code>alt-alt</code> (heterozygous alternate) </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 3.7.13. ACMG disease </td>
    <td> Add two tags:
      <ul>
        <li> <code>F_ACMG</code> indicates whether the variant has an ACMG disease annotation </li>
        <li> <code>F_ACMG_Coding</code> indicates whether the variant is coding and has an ACMG disease annotation </li>
      </ul> 
    </td>
  <td> <code>F_ACMG = 1</code> if the variant does not have NA in the <code>ACMG_disease</code> column <br><br>
    <code>F_ACMG_Coding = 1</code> if the variant
      <ul>
        <li> does not have NA in the <code>ACMG_disease</code> column </li>
        <li> whose type of sequence overlapped with respect to known genes/transcripts is one of "exonic", "splicing", or "exonic;splicing" </li>
      </ul>
    </td>
  </tr>
</tbody>
</table>


### :brain: Working with PrioritizationPipeline v17 Version 2 (Family Samples)

#### :mag_right: Script Structure Overview

The family-based pipeline adds filtering tags based on the child columns and identifies true compound heterozygotes in the child using information from the parents. The filter and statistics functions are the same as the ones used in Version 1. The steps in the main script are modified and family-specific functions were added.

##### :seedling: Main Script

<table>
<thead>
  <tr>
    <th>Section</th>
    <th>Content</th>
    <th>Purpose</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td rowspan="3"> (0) Settings </td>
    <td> 0.1. Libraries </td>
    <td> Import all libraries used in the script </td>
  </tr>
  <tr>
    <td> 0.2. Input variables </td>
    <td> <code>args[1]</code>: root path <br>
         <code>args[2]</code>: sub-path and filename of second R code file defining functions <br>
         <code>args[3]</code>: family variant input file path <br>
         <code>args[4]</code>: child genome name used at columns <br>
         <code>args[5]</code>: sub-path and filename of pedigree file <br>
         <code>args[6]</code>: sub-path and filename of stats visualization Rmd file <br>
         <code>args[7]</code>: output sub-path </td>
  </tr>
  <tr>
    <td> 0.3. Functions </td>
    <td> Import functions used in the script </td>
  </tr>
  <tr>
    <td rowspan="5">(1) Variables & Cutoffs</td>
    <td> 1.1. Internal variables </td>
    <td rowspan="5"> Modify internal variables and cutoffs here </td>
  </tr>
  <tr>
    <td> 1.2. Cutoffs </td>
  </tr>
  <tr>
    <td> 1.2.1. High-quality filter </td>
  </tr>
  <tr>
    <td> 1.2.2. Define damage </td>
  </tr>
  <tr>
    <td> 1.2.3. Main findings </td>
  </tr>
  <tr>
    <td rowspan="13">(2) Main</td>
    <td> 2.1. Import pedigree file </td>
    <td> See title </td>
  </tr>
  <tr>
    <td> 2.2. Find the parents </td>
    <td> Identify the child and the parents from the pedigree file and extract the parents' sample IDs </td>
  </tr>
  <tr>
    <td> 2.3. Import variant data </td>
    <td> 
        <ul>
            <li> Imports the original variant data containing family samples </li>
            <li> Check if all samples exist in the data. An error is returned if at least one sample is not found </li>
        </ul> 
    </td>
  </tr>
  <tr>
    <td> 2.4. Re-format column names</td>
    <td> Remove <code>{genome_name}:</code> from several columns using the child's genome name to identify which columns to use to add the filtering tags </td>
  </tr>
  <tr>
    <td> 2.5. Process the original imported variant data </td>
    <td> 
      <ul>
        <li> Remove variants with homozygous reference (hom-ref) or unknown zygosity from the data </li>
        <li> Calculate the alternate allele frequency </li>
      </ul>
    </td>
  </tr>
    <tr>
    <td> 2.5.1. Add some fields for allele frequency </td>
    <td> 
      <ul>
          <li> <code>FreqMaxSimple_AfrAmrEasNfeSasOth</code> = the maximum of a variant's allele frequencies (exomes and genomes) from the African-American/African, Latino/Admixed American, East Asian, Non-Finnish European, South Asian, and other subsets </li>
          <li> <code>FreqHomCount_AfrAmrEasNfeSasOth</code> = the maximum between the counts of homozygous individuals in samples from the gnomAD exome dataset and from the gnomAD genome dataset </li>
          <li> <code>dbsnp_region_count</code> = a variant's total number of overlap-based match for dbSNP </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 2.5.2. Add some fields for gene constraint </td>
    <td> <code>F_GeneConstr</code> = 
        <ul>
            <li> 1, if a variant's LOF observed/expected (oe) metric from the gnomAD constraint matrix is less than or equal to 0.15 or its missense oe metric from the gnomAD constraint matrix is less than or equal to 0.70 (i.e. <code>gnomad_oe_lof &le; 0.15 | gnomad_oe_mis &le; 0.70</code>) </li>
            <li> 0, otherwise </li>
        </ul>
    </td>
  </tr>
  <tr>
    <td> 2.6. Free up memory </td>
    <td> Remove <code>v_full.temp.df</code> from the current workspace </td>
  </tr>
  <tr>
    <td> 2.7. Add filters </td>
    <td> 
      <ul>
        <li> Remove alternate contigs and unlocalized/unplaced sequence from the data </li>
        <li> Obtain variants with a maximum frequency of 0.05 and annotate them with filtering tags </li>
          <ul> 
              <li> Note that a new functionality that identifies true compound heterozygotes can be found here </li>
          </ul>
        <li> Obtain high-quality variants that have a "PASS" FILTER and annotate them with filtering tags </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 2.8. Get chromosome counts and chromosome-wise zygosity counts </td>
    <td> 
      <ul>
        <li> Obtain the number of variants in each chromosome </li> 
        <li> Obtain the number of alt-alt, hom-alt, ref-alt and the percentage of hom-alt in each chromosome </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td> 2.9. Get all variant stats </td>
    <td> Obtain the summary statistics for all data sets in one table </td>
  </tr>
  <tr>
    <td> 2.10. Process the results before outputting </td>
    <td> 
        <ul>
            <li> Change tier levels from 0, 1 and 2 to "-", "Low" and "High" in columns <code>F_DamageTier</code> and <code>F_PhenoTier</code> </li>
            <li> Add <code>{genome_name}:</code> back to the columns </li>
        </ul>
    </td>
  </tr>
  <tr>
    <td> 2.11. Output desired results as .txt files </td>
    <td> Output the following to the user-defined output directory:
        <ul>
            <li> Rare 5% variants </li>
            <li> Rare 5% damaging variants </li>
            <li> Rare 5% secondary findings variants </li> 
            <li> Chromosome zygosity counts tables </li>
            <li> A summary statistics table </li> 
            <li> A PDF file that visualizes the summary statistics </li>
            <li> True compound heterozygotes based on main findings </li>
            <li> True compound heterozygotes based on secondary findings </li>
        </ul>
    </td>
  </tr>
</tbody>
</table>

##### :seedling: Family-based Columns
The family-based script contains two new columns: <code>FM_Fam_CmpHet</code> and <code>FS_Fam_CmpHet</code>. Their definitions are shown below.

<table>
<thead>
  <tr>
    <th>Section</th>
    <th>Purpose</th>
    <th>Column Definitions</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td> (4) Family-based functions </td>
    <td> Identifies true compound heterozygotes in the child based on the parents' genotypes </td>
    <td> <code>FM_Fam_CmpHet = 1</code> if
         <ul>
             <li> the variant was tagged as a potential compound heterozygote based on the criteria defind in Main Findings (i.e. <code>FM_PCHET &ge; 1</code>) </li>
             <li> there are at least two variants that are located at different loci within the same gene </li>
         </ul>
         <br>
         <code>FS_Fam_CmpHet = 1</code> if
         <ul>
             <li> the variant was tagged as a potential compound heterozygote based on the criteria defind in Secondary Findings (i.e. <code>F_CmpHet_S1 &ge; 1</code>) </li>
             <li> there are at least two variants that are located at different loci within the same gene </li>
         </ul>
    </td>
  </tr>
</tbody>
</table>

## :bulb: Changes From the Old Script
* Updated arguments to avoid path repetition
* Added <code>frameshift block substitution</code> and <code>nonframeshift block subsitution</code> to <code>eff_lof.chv</code> and <code>eff_other_sub.chv</code> upon [updates from ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/gene/)
* Used <code>data.table::fread</code> to achieve a faster speed when importing the original variant data
  * When reading in data table, indicated "." as an NA value to maintain numeric as numeric
* Added a line to remove column DP immediately after reading in the original variant data - this is because at one point the script removes the <code>{genome_name}.</code> part in columns that start with it, which includes <code>{genome_name}.DP</code>. After the removal of <code>{genome_name}.</code>, there would be two DP columns and the first DP column would be used by default, which is not desired
* Added a simple definition of max frequency (i.e. column <code>FreqMaxSimple_AfrAmrEasNfeSasOth</code>) based on simple max over major populations (AFR, AMR, EAS, NFE, SAS, OTH) for manual filtering
* Added the column <code>dbsnp_region_count</code> for manual filtering
* Created a definition for constrained genes (i.e. column <code>F_GeneConstr</code>), used for predicted dominant and manual filtering
* Added a filtering step that removes alternate contigs and unlocalized/unplaced sequence from the variants before adding other filters
* Replaced <code>F_Qual_tag</code> with <code>F_Pass</code> and assigned a new definition to <code>F_Qual</code> for a more intuitive understanding of the columns
* Changed the frequency filter tags from <code>0.05</code>, <code>0.01</code>, <code>0.005</code>, <code>0.0015</code>, <code>0</code> to <code>0.05</code>, <code>0.01</code>, <code>0.001</code>, <code>0.0001</code>, <code>0</code>
* Changed the following criteria:
  * High-quality variants are now defined as variants with a "PASS" FILTER and <code>DP &ge; 2</code> (i.e. <code>F_Qual = 1</code>)
  * Restricted LOF with a splicing mechanism to intronic variants within 2 nucleotides from the splice junction (i.e. <code>"typeseq_priority" %in% c("splicing") & distance_spliceJunction <= 2</code>)
    * Retired <code>F_DamageTier = 1</code> for more suspicious splicing LOF variants, this issue is taken care by the updated splicing LOF definition
    * Retired <code>add_coding_lof_spliceJunction_tag</code>, this issue is taken care by the updated splicing LOF definition
  * Retired <code>sift_score</code>, <code>polyphen_score</code>, and <code>ma_score</code> from the Missense criteria and changed how Missense variants are defined
  * Splicing predictions
    * Removed the SPIDEX-related criteria
    * Added indels with a splicing effect and an unrestricted distance from the splice site (i.e. defaulting to Annovar's cutoff) to the lowest splicing damage tier (i.e. = 1) by adding <code>| (var_type %in% c("del", "ins") & typeseq_priority %in% c("splicing", "exonic;splicing"))</code>
  * Removed the <code>(effect_priority %in% eff_other_sub.chv & (phylopMam_avg >= 1.5 | phylopVert100_avg >= 2.0 | CADD_phred >= 13.0) & is.na (dbsnp) & is.na (dbsnp_region))</code> condition from the criteria for defining damaging variants with type Other Coding
* Changed the <code>phylopMam</code>, <code>phylopVert100</code>, <code>CADD_phred</code> cutoffs for the following sections:
  * 6.2. Missense
  * 6.3. Other coding
  * 6.5. UTR
  * 6.6. Non-coding
* Changed the default type for <code>F_DamageType</code> to "NotDmg" for a more intuitive understanding
* Created a new definition for dominant genes, by grepping "i:AD" from <code>omim_phenotype</code> and "AD" from <code>CGD_inheritance</code>
* Retired <code>add_stats_pheno_filter</code> as the separate OMIM and CGD definitions make less sense
* Changed allele frequency cutoff for dominant, X-linked haploid and predicted dominant to 0.0001
* Fixed zygosity filter for X-linked haploid
* The new script annotates all variants with a maximum allele frequency of 0.05 throughout and only outputs one annotated data set containing rare 0.05 variants. On contrary, the old script annotated the rare 0.05 variants with a "PASS" FILTER throughout and output rare 0.05 and rare 0.05 with "PASS" variants separately. This change was made because only ~4% of variants called by DRAGEN are low-quality (i.e. they do not have a "PASS" filter)
* Summary statistics
  * Removed the extra "0" from the frequencies in the stats names (e.g. <code>050</code> to <code>05</code>, <code>010</code> to <code>01</code>)
  * Changed the frequencies from <code>050</code>, <code>010</code>, <code>005</code>, <code>0015</code>, and <code>000</code> to <code>05</code>, <code>01</code>, <code>001</code>, <code>0001</code>, and <code>00</code>
  * Changed <code>Q1</code> to <code>Q2</code> to <code>LQ</code> and <code>HQ</code>, respectively
  * Changed <code>VarN_Q[12]_Coding_Rare010_LOF</code> to <code>VarN_[LH]Q_Coding_Rare01_LOF_TierLow</code> and added <code>VarN_[LH]Q_Coding_Rare01_LOF_TierHigh</code> which describes the total number of rare 1% Coding LOF variants with a high damage tier (i.e. = 1)
  * Removed <code>VarN_[LH]Q_Coding_Rare010_Dmg_AXD_[HPO|CGD|All]</code>
  * Changed the frequency for Xhap, AXDom, and PDDom (previously HI) from <code>Rare005</code> to <code>Rare0001</code>
  * Changed <code>HI</code> to <code>PDDom</code> in <code>VarN_[LH]Q_Rare0001_DmgT2_HI_PhenoTier[All|Low]</code> and added <code>VarN_[LH]Q_Rare005_DmgT2_PDDom_PhenoTierHigh</code>
  * Added <code>VarN_FM_True_CmpHet</code> and <code>VarN_FS_True_CmpHet</code> to the family pipeline variant stats that indicate how many true compound heterozygotes are there in total using the Main and Secondary Findings, respectively
* Writing session info to file

## :handshake: Contributors
* Dr. Daniele Merico (original creator of the prioritization pipeline)
* [Bhooma Thiruvahindrapuram](https://github.com/bthiruv)
* [Dr. Worrawat Engchuan](https://github.com/naibank)
* [Thomas Nalpathamkalam](https://github.com/TNalpat)
