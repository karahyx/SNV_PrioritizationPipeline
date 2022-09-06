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
    <td rowspan="7">(0) Variables &amp; Cutoffs</td>
    <td>0.1. Input variables</td>
    <td rowspan="7">Modify file locations and cutoffs here</td>
    <td rowspan="7">N/A</td>
  </tr>
  <tr>
    <td>0.2. Output variables</td>
  </tr>
  <tr>
    <td>0.3. Internal variables</td>
  </tr>
  <tr>
    <td>0.4. Cutoffs</td>
  </tr>
  <tr>
    <td>0.4.1. High-quality filter</td>
  </tr>
  <tr>
    <td>0.4.2. Define damage</td>
  </tr>
  <tr>
    <td>0.4.3. Main findings</td>
  </tr>
  <tr>
    <td rowspan="12">(1) Functions<br></td>
    <td>1.1. Frequency filter</td>
    <td rowspan="12">All functions used in the script can be found here</td>
    <td rowspan="12">N/A</td>
  </tr>
  <tr>
    <td>1.2. Quality filter</td>
  </tr>
  <tr>
    <td>1.3. Coding tag</td>
  </tr>
  <tr>
    <td>1.4. Define damage</td>
  </tr>
  <tr>
    <td>1.5. Phenotype filter</td>
  </tr>
  <tr>
    <td>1.6. Main findings</td>
  </tr>
  <tr>
    <td>1.7. Secondary findings</td>
  </tr>
  <tr>
    <td>1.8. Get final results</td>
  </tr>
  <tr>
    <td>1.8.1. Get rare variants</td>
  </tr>
  <tr>
    <td>1.8.2. Get HQ variants (FILTER = "PASS")</td>
  </tr>
  <tr>
    <td>1.8.3. Get HQ rare variants</td>
  </tr>
  <tr>
    <td>1.9. Stats</td>
  </tr>
  <tr>
    <td>(1.5) File Import &amp; Pre-processing</td>
    <td>N/A</td>
    <td>Imports the original variant data and removes those variants with homozygous reference or unknown zygosity</td>
    <td>N/A</td>
  </tr>
  <tr>
    <td rowspan="5">(2) Frequency Filter</td>
    <td>allele frequency cutoff = 0.05</td>
    <td rowspan="5">Add a frequency filter that filters for variants that pass a specific allele frequency cutoff</td>
    <td rowspan="5"><code>F_Rare</code> = <br> the smallest allele frequency cutoff that the variant passes</td>
  </tr>
  <tr>
    <td>allele frequency cutoff = 0.01</td>
  </tr>
  <tr>
    <td>allele frequency cutoff = 0.005</td>
  </tr>
  <tr>
    <td>allele frequency cutoff = 0.0015</td>
  </tr>
  <tr>
    <td>allele frequency cutoff = 0</td>
  </tr>
  <tr>
    <td rowspan="2">(3) Quality Filter</td>
    <td>3.1. Pass tag</td>
    <td>Add a pass tag that indicates whether the variant has FILTER = "PASS"</td>
    <td><code>F_Pass</code> = <br> whether the variant has a "PASS" FILTER</td>
  </tr>
  <tr>
    <td>3.2. Quality tag</td>
    <td>Add a quality tag that indicates whether variants with a "PASS" FILTER pass the DP cutoff</td>
    <td><code>F_Qual_tag</code> = 
      <ul>
        <li>"OK" if the variant has DP &ge; 2</li>
        <li>"LowQuality" if the variant has DP < 2 </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>(4) Coding Tag/td>
    <td>Add coding tags to the variant data</td>
    <td>Add a coding tag that indicates whether the variants' type of sequence overlapped is Coding, ncRNA, or Other</td>
    <td><code>F_Coding</code> = 
      <ul>
        <li>"Coding" if the variant's <code>typeseq_priority</code> is one of <code>exonic</code>, <code>exonic;splicing</code>, or <code>splicing</code></li>
        <li>"ncRNA" if the variant's <code>typeseq_priority</code> is one of <code>ncRNA_exonic</code>, <code>ncRNA_splicing</code>, or <code>ncRNA_exonic;ncRNA_splicing</code></li>
        <li>"Other" otherwise </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td rowspan="7">(5) Define Damage</td>
    <td>5.0. Variable initialization</td>
    <td>Initialize the following columns:<br><code>F_DamageType = "NotDmg"</code><br><code>F_DamageRank = 0</code><br><code>F_S_DamageType = "NotDmg"</code></td>
    <td><code>F_DamageType</code> = a variant's damage type
      <ul>
        <li> one of <code>LOF</code>, <code>Missense</code>, <code>OtherC</code>, <code>Splc</code>, <code>UTR</code>, <code>DmgNcRNA</code>, or <code>NotDmg</code></li>
      </ul>
      <br><code>F_DamageRank</code> = a variant's damage rank &isin; {0, 1, 2}
        <ul>
          <li>Note that the higher the rank (i.e. the bigger the number), the more damaging a variant is</li>
        </ul>
       <br><code>F_S_DamageType</code> = <br> a more stringent Coding LOF damage type tag with the distance from the nearest exon boundary as an additional condition
        <ul>
          <li>Note that F_S_DamageType is specific to the Coding LOF category, thus one of <code>LOF</code> or <code>NotDmg</code></li>
          <li>May be used if more stringent Coding LOF variants are desired</li>
        </ul>
    </td>
  </tr>
  <tr>
    <td>5.1. Coding LOF</td>
    <td><code>add_coding_lof_tag()</code> identify variants of type <code>Coding LOF</code> and change their <code>F_DamageRank</code> tag to 2;<br><br><code>add_coding_lof_spliceJunction_tag()</code> identifies variants of type <code>Coding LOF</code> with <code>distance_spliceJunction</code> &lt; 3; note that no change is made to <code>F_DamageRank</code> here</td>
    <td>The variant is <code>Coding LOF</code> if it
      <ul>
        <li> is coding </li>
        <li> causes frameshift or point mutations; or its type of sequence overlapped is splicing or exonic splicing </li>
      </ul>
      <br>The <code>F_S_DamageType</code> is changed to "LOF" from "NotDmg" when
      <ul>
        <li> the variant is Coding LOF (i.e. satisfies the two conditions above) </li>
        <li> the variant's <code>distance_spliceJunction</code> < 3 </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>5.2. Missense</td>
    <td> Identify variants of type <code>Missense</code> and change its <code>F_DamageRank</code> to either 1 or 2 based on specific conditions. <br><br> The method first compares the variant's SIFT, Polyphen2, MA, phylopMam, phylopVert, and CADD_phred scores to their corresponding cutoffs and documents the results (0 or 1) in a matrix with individual variants on each row. Next, the sum of the scores for each variant are calculated (a variant has a max score of 6) and compared to rank 1 and 2 cutoffs to decide which damage rank it belongs to. 
    </td>
    <td> The variant is <code>Missense</code> and has <code>F_DamageRank</code> = 1 if it
      <ul>
        <li> is coding </li> 
        <li> is a nonsynonymous SNV </li>
        <li> has a sum score &ge; 2 </li>
      </ul>
      <br> The variant is <code>Missense</code> and has <code>F_DamageRank</code> = 2 if it
      <ul>
        <li> is coding </li> 
        <li> is a nonsynonymous SNV </li>
        <li> has a sum score &ge; 4 </li>
      </ul>
      <br> "1" is documented in the sum score matrix if the variant's
      <ul>
        <li> <code>sift_score</code> < 0.05 </li> 
        <li> <code>polyphen_score</code> &ge; 0.90 </li>
        <li> <code>ma_score</code> &ge; 1.90 </li>
        <li> <code>phylopMam_avg</code> &ge; 1.30 </li>
        <li> <code>phylopVert100_avg</code> &ge; 3.90 </li>
        <li> <code>CADD_phred</code> &ge; 21.10 </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>5.3. Other coding</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>5.4. Splicing predictions</td>
    <td>  </td>
    <td></td>
  </tr>
  <tr>
    <td>5.5. UTR</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>5.6. Non-coding</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td rowspan="3">(6) Phenotype Filter<br><br></td>
    <td>6.1. HPO dominant</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>6.2. CGD dominant</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>6.3. Phenotype ranks</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td rowspan="5">(7) Main Findings<br></td>
    <td>7.1. Recessive homozygous</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>7.2. X-linked haploid</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>7.3. Potential compound heterozygous</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>7.4. Dominant</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>7.5. Heterozygous hotzone</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td rowspan="14">(8) Secondary Findings<br></td>
    <td>8.0. Pathogenicity flag</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.1. Rank 1</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.1.1. Dominant, pathogenic</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.1.2. Recessive, homozygous, pathogenic</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.1.3. Recessive, potential compound heterozygous, <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;pathogenic</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.1.4. X-linked, homozygous/haploid, pathogenic</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.1.5. Complex, homozygous, pathogenic</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.1.6. Complex, potential compound heterozygous,<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;pathogenic</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.1.7. Complex, single heterozygous, uncertain</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.1.8. Recessive, single heterozygous</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.1.9. X-linked, heterozygous, carrier</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.2. Rank 2</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.3. Rank 3</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>8.4. ACMG disease</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td rowspan="10">(9) Main<br></td>
    <td>Step 1. File import</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>Step 2. Re-format column names</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>Step 3. Process the original imported variant data</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>Step 4. Free up memory</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>Step 5. Annotate the data</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>Step 6. Get chromosome counts and chromosome-<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;wise zygosity counts</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>Step 7. Get a summary stats list for each data set</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>Step 8. Convert the stats lists to readable data frames</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>Step 9. Get all the summary stats in one data frame</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>Step 10. Output desired results as .txt files</td>
    <td></td>
    <td></td>
  </tr>
</tbody>
</table>

#### Running the Script

### Working with PrioritizationPipeline Version 2 (HPF Version)

#### Script Structure Overview

#### Running the Script

### Stats Summary

## Changes From the Old Script

## Credits

## License
