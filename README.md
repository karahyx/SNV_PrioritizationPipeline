# :dna: PriorizationPipeline

An adaptation of the TCAG Small Variant Prioritization Pipeline originally developed to be used for Broad Institute's [GATK (Genome Analysis Toolkit)](https://gatk.broadinstitute.org/hc/en-us) to the [Illumina DRAGEN (Dynamic Read Analysis for GENomics) Bio-IT Platform](https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/DRAGEN/GPipelineIntro_fDG.htm). 

## :green_book: Description

The PrioritizationPipeline script is used to annotate small variants called by DRAGEN with features including frequency and quality filters, damage type and ranks, as well as main and secondary findings for the purpose of variant prioritization. The script also generates pre-defined summary statistics such as chromosome-wise zygosity type counts and the number of variants in each category of interest. 

PrioritizationPipeline relies on the TCAG Small Variant Annotate Pipeline rev27.7 developed and maintained by Thomas Nalpathamkalam. The documentation can be found here [TCAG_SMALL VARIANT_ANNOTATION_PIPELINE_rev27.7_hg38_JUN2022.pdf](./TCAG_SMALL_VARIANT_ANNOTATION_PIPELINE_rev27.7_hg38_JUN2022.pdf). In comparison, [Pipeline_v16](./R/pipeline_old.R) relies on [TCAG_SMALL VARIANT_ANNOTATION_PIPELINE_rev27.4_hg18_AUG2022.pdf](./TCAG_SMALL_VARIANT_ANNOTATION_PIPELINE_rev27.4_hg18_AUG2022.pdf).

Two versions of PrioritizationPipeline are provided. [Version 1: Pipeline_v17_ILMN_DRAGEN_rev27.7](./R/Pipeline_v17_ILMN_DRAGEN_rev27.7_20220909.R) is a general and more detailed version for users who wish to run the script on one file at a time and on their local machine. The individual sections that add different types of annotations and make up the Main Section are provided. As the file name suggests, [Version 2: Pipeline_v17_ILMN_DRAGEN_rev27.7_HPF](./R/Pipeline_v17_ILMN_DRAGEN_rev27.7_HPF_20220909.R) is designed to be called by a shell script that receives user input and usually used to process multiple files at once on HPF platforms. The shell scripts in the [bash](./bash) folder are used to perform such tasks. Detailed instructions for how to run the file on HPF can be found in the Running the Script sections below.

The main change from the old script is the utilization of functions that reduce the program to smaller, more manageable chunks and allow for reusability and extension. Detailed changes can be found in the Changes From the Old Script section below.

## :desktop_computer: Instructions

### :brain: Working with PrioritizationPipeline v17 Version 1 (General Version)

#### :mag_right: Script Structure Overview

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
    <td rowspan="13">(1) Functions<br></td>
    <td>1.0. Alternate contigs and unlocalized/unplaced sequence filter</td>
    <td rowspan="13">All functions used in the script can be found here</td>
    <td rowspan="13">N/A</td>
  </tr>
  <tr>
    <td>1.1. Frequency filter</td>
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
    <td>(2) Alternate Contigs and Unlocalized/Unplaced Sequence Filter</td>
    <td> N/A </td>
    <td> Remove <a href="https://gatk.broadinstitute.org/hc/en-us/articles/360041155232-Reference-Genome-Components"> alternate contigs and unlocalized/unplaced sequence</a> from the variant data </td>
    <td> N/A
    </td>
  </tr>
  <tr>
    <td rowspan="5">(3) Frequency Filter</td>
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
    <td rowspan="2">(4) Quality Filter</td>
    <td>4.1. Pass tag</td>
    <td>Add a pass tag that indicates whether the variant has FILTER = "PASS"</td>
    <td><code>F_Pass</code> = 
      <ul>
        <li> 1, if the variant has a "PASS" FILTER </li>
        <li> 0, otherwise </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>4.2. Quality tag</td>
    <td>Add a quality tag that indicates whether variants with a "PASS" FILTER pass the DP cutoff</td>
    <td><code>F_Qual_tag</code> = 
      <ul>
        <li>"OK" if the variant has DP &ge; 2</li>
        <li>"LowQuality" if the variant has DP < 2 </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>(5) Coding Tag</td>
    <td>Add coding tags to the variant data</td>
    <td>Add a coding tag that indicates whether the variants' type of sequence overlapped with respect to known genes/transcripts is Coding, ncRNA, or Other</td>
    <td><code>F_Coding</code> = 
      <ul>
        <li>"Coding" if the variant's <code>typeseq_priority</code> is one of <code>exonic</code>, <code>exonic;splicing</code>, or <code>splicing</code></li>
        <li>"ncRNA" if the variant's <code>typeseq_priority</code> is one of <code>ncRNA_exonic</code>, <code>ncRNA_splicing</code>, or <code>ncRNA_exonic;ncRNA_splicing</code></li>
        <li>"Other" otherwise </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td rowspan="7">(6) Define Damage</td>
    <td>6.0. Variable initialization</td>
    <td>Initialize the following columns:<br><code>F_DamageType = "NotDmg"</code><br><code>F_DamageRank = 0</code><br><code>F_S_DamageType = "NotLOF"</code></td>
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
          <li>Note that F_S_DamageType is specific to the Coding LOF category, thus one of <code>LOF</code> or <code>NotLOF</code></li>
          <li>May be used if more stringent Coding LOF variants are desired</li>
        </ul>
    </td>
  </tr>
  <tr>
    <td>6.1. Coding LOF</td>
    <td><code>add_coding_lof_tag()</code> identify variants of type <code>Coding LOF</code> and change their <code>F_DamageRank</code> to 2;<br><br><code>add_coding_lof_spliceJunction_tag()</code> identifies variants of type <code>Coding LOF</code> with <code>distance_spliceJunction &lt; 3</code>; note that no change is made to <code>F_DamageRank</code> here</td>
    <td>The variant is <code>Coding LOF</code> if it
      <ul>
        <li> is coding </li>
        <li> causes frameshift or point mutations in the coding sequence; or its type of sequence overlapped is splicing or exonic splicing </li>
      </ul>
      <br>The <code>F_S_DamageType</code> is changed to "LOF" from "NotLOF" when
      <ul>
        <li> the variant is Coding LOF (i.e. satisfies the two conditions above) </li>
        <li> the variant's <code>distance_spliceJunction < 3</code></li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>6.2. Missense</td>
    <td> Identify variants of type <code>Missense</code> and change its <code>F_DamageRank</code> to either 1 or 2 based on specific conditions. <br><br> The method first compares the variant's SIFT, Polyphen2, MA, phylopMam, phylopVert, and CADD_phred scores to their corresponding cutoffs and documents the results (0 or 1) in a matrix with individual variants on each row. Next, the sum of the scores for each variant are calculated (each variant has a max score of 6) and compared to rank 1 and 2 cutoffs to decide which damage rank it belongs to. 
    </td>
    <td> The variant is <code>Missense</code> and has <code>F_DamageRank = 1</code> if it
      <ul>
        <li> is coding </li> 
        <li> is a nonsynonymous SNV </li>
        <li> has a sum score &ge; 2 </li>
      </ul>
      <br> The variant is <code>Missense</code> and has <code>F_DamageRank = 2</code> if it
      <ul>
        <li> is coding </li> 
        <li> is a nonsynonymous SNV </li>
        <li> has a sum score &ge; 4 </li>
      </ul>
      <br> "1" is documented in the sum score matrix if the variant's
      <ul>
        <li> <code>sift_score < 0.05</code> </li> 
        <li> <code>polyphen_score &ge; 0.90</code> </li>
        <li> <code>ma_score &ge; 1.90</code> </li>
        <li> <code>phylopMam_avg &ge; 1.30</code> </li>
        <li> <code>phylopVert100_avg &ge; 3.90</code> </li>
        <li> <code>CADD_phred &ge; 21.10</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>6.3. Other coding</td>
    <td> Identify variants of type <code>Other coding</code> and change its <code>F_DamageRank</code> to either 1 or 2 based on specific cutoffs </td>
    <td> The variant is <code>Other coding</code> if it
      <ul>
        <li> is coding </li> 
        <li> causes nonframeshift mutations in the coding sequence </li>
        <li> is not an exact match to common dbSNP track UCSC, and <code>phylopMam_avg &ge; phylopMam_cutoff</code> or <code>phylopVert100_avg &ge; phylopVert_cutoff</code> or <code>CADD_phred &ge; CADD_phred_cutoff</code> </li>
      </ul><br>
      <code>F_DamageRank = 1</code> if 
      <ul>
        <li> <code>phylopMam_cutoff = 1.1</code> </li>
        <li> <code>phylopVert_cutoff = 1.6</code> </li>
        <li> <code>CADD_phred_cutoff = 13.7</code> </li>
      </ul><br>
      <code>F_DamageRank</code> = 2 if 
      <ul>
        <li> <code>phylopMam_cutoff = 1.3</code> </li>
        <li> <code>phylopVert_cutoff = 3.9</code> </li>
        <li> <code>CADD_phred_cutoff = 21.1</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>6.4. Splicing predictions</td>
    <td> Predict whether the variant is of type <code>Splicing</code> and change its <code>F_DamageRank</code> to either 1 or 2 based on specific conditions </td>
    <td> The variant is predicted to be "Splicing" and has <code>F_DamageRank = 1</code> if
      <ul>
        <li> its <code>F_DamageType</code> is not "LOF", "Splc" or "Missense" </li>
        <li> its <code>F_DamageRank &ne; 2</code></li>
        <li> satisfies one or more of the following: </li>
          <ul>
            <li> <code>spliceAI_DS_AG > 0.2</code> & <code>|spliceAI_DP_AG| &le; 50</code> </li>
            <li> <code>spliceAI_DS_AL > 0.2</code> & <code>|spliceAI_DP_AL| &le; 50</code> </li>
            <li> <code>spliceAI_DS_DG > 0.2</code> & <code>|spliceAI_DP_DG| &le; 50</code> </li>
            <li> <code>spliceAI_DS_DL > 0.2</code> & <code>|spliceAI_DP_DL| &le; 50</code> </li>
          </ul>
      </ul> <br>
      The variant is predicted to be "Splicing" and has <code>F_DamageRank = 2</code> if
      <ul>
        <li> its <code>F_DamageType</code> is not "LOF" </li>
        <li> satisfies one or more of the following: </li>
          <ul>
            <li> <code>spliceAI_DS_AG > 0.2</code> & <code>|spliceAI_DP_AG| &le; 50</code> </li>
            <li> <code>spliceAI_DS_AL > 0.2</code> & <code>|spliceAI_DP_AL| &le; 50</code> </li>
            <li> <code>spliceAI_DS_DG > 0.2</code> & <code>|spliceAI_DP_DG| &le; 50</code> </li>
            <li> <code>spliceAI_DS_DL > 0.2</code> & <code>|spliceAI_DP_DL| &le; 50</code> </li>
            <li> <code>dbscSNV_ADA_SCORE > 0.6</code> & <code>dbscSNV_RF_SCORE > 0.6</code> </li>
          </ul>
      </ul>
    </td>
  </tr>
  <tr>
    <td>6.5. UTR</td>
    <td> Identify variants of type <code>UTR</code> and change its <code>F_DamageRank</code> to 1 or 2 based on specific cutoffs </td>
    <td> The variant is <code>UTR</code> if 
      <ul>
        <li> the type of sequence overlapped with respect to known genes/transcripts is "UTR3", "UTR5", or both </li>
        <li> has a PhastCons score for the Placental Mammal genome group </li>
        <li> satisfies <code>phylopMam_avg &ge; phylopMam_cutoff</code> or <code>CADD_phred &ge; CADD_phred_cutoff</code> </li>
      </ul> <br>
      <code>F_DamageRank = 1</code> if 
        <ul>
          <li> <code>phylopMam_cutoff = 1.1</code> </li>
          <li> <code>CADD_phred_cutoff = 13.7</code> </li>
        </ul> <br>
      <code>F_DamageRank = 2</code> if 
        <ul>
          <li> <code>phylopMam_cutoff = 1.3</code> </li>
          <li> <code>CADD_phred_cutoff = 21.1</code> </li>
        </ul>
    </td>
  </tr>
  <tr>
    <td>6.6. Non-coding</td>
    <td> Identify variants of type <code>Non-coding</code> and change its <code>F_DamageRank</code> to either 1 or 2 based on specific cutoffs </td>
      <td> The variant is <code>Non-coding</code> if 
        <ul>
          <li> it is identified as "ncRNA" in <code>F_Coding</code> </li>
          <li> its full gene name is not "pseudogene" </li>
          <li> its <code>F_DamageRank = 0</code> (i.e. it's not damaging) </li>
          <li> satisfies one or more of the following: </li>
            <ul>
              <li> has a PhastCons score for the Placental Mammal genome group AND, <code>phylopMam_avg &ge; phylopMam_cutoff</code> OR <code>phylopVert100_avg &ge; phylopVert_cutoff</code> </li>
              <li> <code>CADD_phred &ge; CADD_phred_cutoff</code> </li>
            </ul>
        </ul> <br>
        <code>F_DamageRank = 1</code> if 
          <ul>
            <li> <code>phylopMam_cutoff = 1.1</code> </li>
            <li> <code>phylopVert_cutoff = 1.6</code> </li>
            <li> <code>CADD_phred_cutoff = 13.7</code> </li>
          </ul> <br>
        <code>F_DamageRank = 2</code> if 
        <ul>
            <li> <code>phylopMam_cutoff = 1.3</code> </li>
            <li> <code>phylopVert_cutoff = 3.9</code> </li>
            <li> <code>CADD_phred_cutoff = 21.1</code> </li>
          </ul>
      </td>
  </tr>
  <tr>
    <td rowspan="3">(7) Phenotype Filter</td>
    <td>7.1. HPO dominant</td>
    <td> Add a tag <code>G_AXD_HPO</code> that indicates whether the variant has autosomal dominant (AD) as their mode of inheritance based on the HPO annotations </td>
    <td> <code>G_AXD_HPO</code> = 
      <ul>
        <li> 1, if the pattern "@AD" is found in column <code>HPO</code> </li>
        <li> 0, otherwise </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>7.2. CGD dominant</td>
    <td> Add a tag <code>G_AXD_CGD</code> that indicates whether the variant has autosomal dominant (AD) as their mode of inheritance based on the CGD inheritance annotations </td>
    <td> <code>G_AXD_CGD</code> = 
      <ul>
        <li> 1, if the pattern "AD" is found in column <code>CGD_inheritance</code> </li>
        <li> 0, otherwise </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>7.3. Phenotype ranks</td>
    <td> Add a tag <code>F_PhenoRank</code> that indicates the phenotype rank of a variant. Here, rank 1 and 2 are used to differentiate mouse and human phenotype annotations, respectively. Human phenotype annotations take priority. </td>
    <td> <code>F_PhenoRank</code> =
      <ul>
        <li> 1, if the variant has an MPO annotation imported from MGI (Mouse Genome Informatics) and mapped from an orthologous mouse gene </li>
        <li> 2, if the variant has one or more annotations from OMIM, HPO, or CGD </li>
        <li> 0, otherwise </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td rowspan="5">(8) Main Findings</td>
    <td>8.1. Recessive homozygous</td>
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
    <td>8.2. X-linked haploid</td>
    <td> Add a tag <code>FM_XHAP</code> that indicates whether the variant is an X-linked haploid </td>
    <td> <code>FM_XHAP = 1</code> if the variant
      <ul>
        <li> is a rare variant with a maximum allele frequency of 0.05 </li>
        <li> is damaging, i.e. <code>F_DamageType</code> != "NotDmg" </li>
        <li> has zygosity <code>hom-alt</code> (homozygous alternative) </li>
        <li> is found in chromosome X </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>8.3. Potential compound heterozygous</td>
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
        <li> has <code>F_DamageRank = 2</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>8.4. Dominant</td>
    <td> Add a tag <code>FM_AXDOM</code> that indicates whether the variant is autosomal dominant (AD) </td>
    <td> <code>FM_AXDOM = 1</code> if the variant
      <ul>
        <li> is a rare variant with a maximum allele frequency of 0.005 </li>
        <li> is damaging, i.e. <code>F_DamageType</code> != "NotDmg" </li>
        <li> has either <code>G_AXD_CGD == 1</code> or <code>G_AXD_HPO == 1</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>8.5. Heterozygous hotzone</td>
    <td> Add a tag <code>FM_HZ</code> that indicates whether the variant is part of a heterozygous hotzone </td>
    <td> <code>FM_HZ = 1</code> if the variant
      <ul>
        <li> is a rare variant with a maximum allele frequency of 0.0015 </li>
        <li> has zygosity <code>ref-alt</code> (heterozygous reference) </li>
        <li> has <code>gnomAD_oe_lof_upper < 0.35</code> </li>
        <li> satisfies one of the following: </li>
          <ul>
            <li> <code>F_DamageType</code> is one of "LOF", "Splc", or "OtherC" </li>
            <li> <code>F_DamageType</code> is "Missense" and <code>F_DamageRank = 2</code> </li>
          </ul>
      </ul>
    </td>
  </tr>
  <tr>
    <td rowspan="14">(9) Secondary Findings</td>
    <td>9.0. Pathogenicity flag</td>
    <td> Add two pathogenicity related flags:
      <ul> 
        <li> flag <code>F_Clinvar_Pathg</code> indicates whether the variant has at least one record submitted with pathogenic or likely pathogenic based on ClinVar </li>
        <li> flag <code>F_Clinvar_notPathg</code> indicates whether the variant has no current value of pathogenic based on ClinVar </li>
      </ul>
    </td>
    <td> <code>F_Clinvar_Pathg = 1</code> if <code>Clinvar_SIG_Simple = 1</code>
    <br><br> <code>F_Clinvar_notPathg = 1</code> if <code>Clinvar_SIG_Simple = 0</code>
    </td>
  </tr>
  <tr>
    <td>9.1. Rank 1</td>
    <td> Add a tag <code>FS1_Select</code> that indicates whether a variant belongs to rank 1 for secondary findings </td>
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
    <td>9.1.1. Dominant, pathogenic</td>
    <td> Add a tag <code>FS1_AD_Pathg_Any</code> that indicates whether a variant is dominant and pathogenic </td>
    <td> A variant has <code>FS1_AD_Pathg_Any = 1</code> if
      <ul>
        <li> it belongs to rank 1 for secondary findings, i.e. <code>FS1_Select = 1</code> </li>
        <li> it has autosomal dominant (AD) as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>9.1.2. Recessive, homozygous, pathogenic</td>
    <td> Add a tag <code>FS1_AR_Pathg_Hom</code> that indicates whether a variant is recessive, homozygous and pathogenic </td>
    <td> A variant has <code>FS1_AR_Pathg_Hom</code> if it
      <ul>
        <li> belongs to rank 1 in secondary findings, i.e. <code>FS1_Select = 1</code> </li>
        <li> has autosomal recessive (AR) as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> has zygosity "homozygous alternative" (hom-alt) </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>9.1.3. Recessive, potential compound heterozygous, pathogenic</td>
    <td> Two tags are being added here:
      <ul>
        <li> Add a tag <code>F_CmpHet_S1</code> that indicates whether a variant is a potential compound heterozygote </li>
        <li> Add a tag <code>FS1_AR_Pathg_PotCompHet</code> that indicates whether a variant is recessive, potential compound heterozygous, and pathogenic </li>
      </ul>
    </td>
    <td>
      <code>F_CmpHet_S1</code> =
      <ul>
        <li> 1, if it has a "PASS" FILTER and belongs to rank 1 for secondary findings </li>
        <li> 2, if it has a "PASS" FILTER and DP &ge; 2 and belongs to rank 1 for secondary findings </li>
        <li> 0, otherwise </li>
      </ul><br>
      A variant has <code>FS1_AR_Pathg_PotCompHet = 1</code> if it
      <ul>
        <li> belongs to rank 1 for secondary findings </li>
        <li> has autosomal recessive (AR) as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> has <code>F_CmpHet_S1 &ge; 1</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>9.1.4. X-linked, homozygous/haploid, pathogenic</td>
    <td> Add a tag <code>FS1_XL_Pathg_Hom</code> that indicates whether a variant is X-linked homozygous and pathogenic <br><br>
    Add a tag <code>FS1_XL_Pathg_Hap</code> that indicates whether a variant is X-linked haploid and pathogenic </td>
    <td> <code>FS1_XL_Pathg_Hom = 1</code> or <code>FS1_XL_Pathg_Hap = 1</code> if a variant
      <ul>
        <li> belongs to rank 1 for secondary findings </li>
        <li> has X-linked (XL) as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>9.1.5. Complex, homozygous/haploid, pathogenic</td>
    <td> Add a tag <code>FS1_CX_Pathg_HomHap</code> that indicates whether a variant is complex, homozygous/haploid, and pathogenic </td>
    <td> A variant has <code>FS1_CX_Pathg_HomHap = 1</code> if it
      <ul>
        <li> belongs to rank 1 for secondary findings </li>
        <li> does not have AD, AR, or XL as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> has zygosity <code>hom-alt</code> (homozygous alternative) </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>9.1.6. Complex, potential compound heterozygous, pathogenic</td>
    <td> Add a tag <code>FS1_CX_Pathg_PotCompHet</code> that indicates whether a variant is complex, potential compound heterozygous, and pathogenic </td>
    <td> A variant has <code>FS1_CX_Pathg_HomHap = 1</code> if it
      <ul>
        <li> belongs to rank 1 for secondary findings </li>
        <li> does not have AD, AR, or XL as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> has <code>F_CmpHet_S1 &ge; 1</code> </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>9.1.7. Complex, single heterozygous, uncertain</td>
    <td> Add a tag <code>FS1_CX_Uncertain</code> that indicates whether a variant is complex, single heterozygous, and pathogenicity uncertain </td>
    <td> A variant has <code>FS1_CX_Uncertain = 1</code> if it
      <ul>
        <li> belongs to rank 1 for secondary findings </li>
        <li> does not have AD, AR, or XL as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> is not potential compound heterozygous, i.e. has <code>F_CmpHet_S1 = 0</code> </li>
        <li> has zygosity <code>ref-alt</code> (heterozygous reference) or <code>alt-alt</code> (heterozygous alternate) </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>9.1.8. Recessive, single heterozygous, carrier</td>
    <td> Add a tag <code>FS1_AR_Carrier</code> that indicates whether a variant is recessive, single heterozygous, and a carrier </td>
    <td> A variant has <code>FS1_AR_Carrier = 1</code> if it
      <ul>
        <li> belongs to rank 1 for secondary findings </li>
        <li> has AR as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> is not potential compound heterozygous, i.e. has <code>F_CmpHet_S1 = 0</code> </li>
        <li> has zygosity <code>ref-alt</code> (heterozygous reference) or <code>alt-alt</code> (heterozygous alternate) </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>9.1.9. X-linked, heterozygous, carrier</td>
    <td> Add a tag <code>FS1_XL_Carrier</code> that indicates whether a variant is X-linked, heterozygous, and a carrier </td>
    <td> A variant has <code>FS1_XL_Carrier = 1</code> if it
      <ul>
        <li> belongs to rank 1 for secondary findings </li>
        <li> has XL as one of its modes of inheritance based on <code>CGD_inheritance</code> </li>
        <li> is found in chromosome X </li>
        <li> has zygosity <code>ref-alt</code> (heterozygous reference) or <code>alt-alt</code> (heterozygous alternate) </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>9.2. Rank 2</td>
    <td> Add a tag <code>FS2_Select</code> that indicates whether a variant belongs to rank 2 for secondary findings </td>
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
    <td>9.3. Rank 3</td>
    <td> Add a tag <code>FS3_Select</code> that indicates whether a variant belongs to rank 3 for secondary findings </td>
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
    <td>9.4. ACMG disease</td>
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
  <tr>
    <td rowspan="10">(10) Main</td>
    <td>Step 1. File import</td>
    <td> Imports the original variant data </td>
    <td> N/A </td>
  </tr>
  <tr>
    <td>Step 2. Re-format column names</td>
    <td> Remove <code>{genome_name}:</code> from several columns for easier access </td>
    <td> N/A </td>
  </tr>
  <tr>
    <td>Step 3. Process the original imported variant data</td>
    <td> 
      <ul>
        <li> Remove variants with homozygous reference (hom-ref) or unknown zygosity from the data </li>
        <li> Add a pass tag that indicates whether the variant has FILTER = "PASS" </li>
        <li> Calculate the alternate allele frequency </li>
      </ul>
    </td>
    <td> N/A </td>
  </tr>
  <tr>
    <td>Step 4. Free up memory</td>
  <td> Remove <code>v_full.temp.df</code> from the current workspace </td>
    <td> N/A </td>
  </tr>
  <tr>
    <td>Step 5. Annotate the data</td>
    <td> 
      <ul>
        <li> Obtain variants with a maximum frequency of 0.05 and annotate them with filtering tags </li>
        <li> Obtain high-quality variants with a maximum frequency of 0.05, a "PASS" FILTER and DP &ge; 2 and annotate them with filtering tags </li>
      </ul>
    </td>
    <td> N/A </td>
  </tr>
  <tr>
    <td>Step 6. Get chromosome counts and chromosome-wise zygosity counts</td>
    <td> 
      <ul>
        <li> Obtain the number of variants in each chromosome </li> 
        <li> Obtain the number of alt-alt, hom-alt, ref-alt and the percentage of hom-alt in each chromosome </li>
      </ul>
    </td>
    <td> N/A </td>
  </tr>
  <tr>
    <td>Step 7. Get a summary stats list for each data set</td>
    <td> Obtain a list of pre-defined summary statistics for the data (see how each summary statistic is defined below) </td>
    <td> N/A </td>
  </tr>
  <tr>
    <td>Step 8. Convert the stats lists to readable tables</td>
    <td> See title </td>
    <td> N/A </td>
  </tr>
  <tr>
    <td>Step 9. Get all the summary stats in one table</td>
    <td> Combine the summary statistics for different data sets into one table for easier comparison </td>
    <td> N/A </td>
  </tr>
  <tr>
    <td>Step 10. Output desired results as .txt files</td>
    <td> Output the annotated data sets, chromosome zygosity counts tables, and summary statistics tables to user-defined output directory </td>
    <td> N/A </td>
  </tr>
</tbody>
</table>

#### :arrow_forward: Running the Script
Steps:
1. Change the following parameters in Section (0)Variables & Cutoffs:
  * <code>input_var_genome.file</code> to the path of your input file
  * <code>input_var_genome.name</code> to the name of your sample
  * <code>output_path</code> to the path of your output directory
2. Then run the following sections:
  * Settings
  * (0) Variables and Cutoffs
  * (1) Functions
  * (10) Main (runs sections 2-9)

### :brain: Working with PrioritizationPipeline v17 Version 2 (HPF Version)

#### :mag_right: Script Structure Overview
Version 2 only contains sections <code>(0) Variables & Cutoffs</code>, <code>(1) Functions</code>, and <code>(10) Main</code>, where <code>(10) Main</code> runs Sections (2)-(9).


#### :arrow_forward: Running the Script
If your variant data files are in the same folder:
1. In <code>run_prioritization_tasks.sh</code>, make sure that lines 12-28 are uncommented and lines 32-49 are commented.
2. Change <code>infile_dir</code> on line 12 to the path of your folder that contains the variant data and <code>output_dir</code> on line 13 to your desired output directory.
3. If your files do not end with .tsv, change the <code>'\*.tsv'</code> on line 16 to <code>'\*.{your_file_format}'</code>. 
4. Run <code>qsub ~/run_prioritization_tasks.sh</code> on HPF.
<br>

If your variant data files are in their own folders and the folders are named after their sample name:
1. In <code>run_prioritization_tasks.sh</code>, make sure that lines 32-49 are uncommented and lines 12-28 are commented.
2. Change <code>infile_dir</code> on line 32 to the path of your folder that contains the variant data and <code>output_dir</code> on line 33 to your desired output directory.
3. Change the <code>'\*SUBSET\*'</code> part on line 41 to a part of the file name that's found in all the variant data file names. 
  * For instance, the file names all have the format <code>{sample_name}.hard-filtered.vcf.gz.annovar.out_SUBSET_rev27.7_hg38.tsv</code> in the example.
5. Run <code>qsub ~/run_prioritization_tasks.sh</code> on HPF.

### Stats Summary

## :bulb: Changes From the Old Script
* Added <code>frameshift block substitution</code> and <code>nonframeshift block subsitution</code> to <code>eff_lof.chv</code> and <code>eff_other_sub.chv</code> upon [updates from ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/gene/)
* Replaced tag <code>F_Qual</code> with <code>F_Pass</code> for a more intuitive understanding of the column
* Changed the following criteria:
  * High-quality variants are now defined as variants with a "PASS" FILTER and <code>DP &ge; 2</code> 
  * Removed the <code> (effect_priority %in% eff_other_sub.chv & (phylopMam_avg >= 1.5 | phylopVert100_avg >= 2.0 | CADD_phred >= 13.0) & is.na (dbsnp) & is.na (dbsnp_region))</code> condition from the criteria for defining damaging variants with type Other Coding
* Changed the <code>phylopMam</code>, <code>phylopVert100</code>, <code>CADD_phred</code> cutoffs for the following sections:
  * 6.2. Missense
  * 6.3. Other coding
  * 6.5. UTR
  * 6.6. Non-coding
* Changed the default type for <code>F_DamageType</code> and <code>F_S_DamageType</code> to "NotDmg" and "NotLOF" for a more intuitive understanding
* Used <code>data.table::fread</code> to achieve a faster speed when importing the original variant data
* Added a line to remove column DP immediately after reading in the original variant data - this is because at one point the script removes the <code>{genome_name}.</code> part in columns that start with it, which includes <code>{genome_name}.DP</code>. After removal, there would be two DP columns and the first DP column would be used by default, which is not desired.

## :handshake: Contributors
* Daniele Merico - Original creator of the pipeline
* [Bhooma Thiruvahindrapuram](https://github.com/bthiruv)
* [Dr. Worrawat Engchuan](https://github.com/naibank)
* [Thomas Nalpathamkalam](https://github.com/TNalpat)

