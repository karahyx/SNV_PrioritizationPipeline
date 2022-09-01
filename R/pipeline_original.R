#### RELEASE NOTES
# ** modified for CGI 2.5 on 17 Jan 2014 by Daniele Merico; originally created for CGI 2.2 on 16 Jan 2014 by Daniele Merico
# ** major changes on 02 Apr 2014 by Daniele Merico: incidental finding, Clinvar probably pathogenic, HGMD, CGD, ploidy 1 variants with hom zygosity can be top quality tier
# ** major changes of 08-18 Apr 2014 by Daniele Merico: to restructure the code flow, and also to integrate more annotations added by Thomas, updated phenotypes by Giovanna
# ** major changes on 23 May 2014 by Daniele Merico: to integrate CGI Wellderly and 1000 Genome frequencies, and minor fix in the ClinVar field
# ** major changes on 04 Jun 2014 by Daniele Merico: "cg1KB436_[...]" was replaced with "cg1KG436_", phenotype pre-composed in separate script and imported from RData file
# ** minor changes on 04 Jun 2014 by Daniele Merico: relaxed high quality tier 2 allelic ratio alt / ref to 0.3
# ** major changes on 12 Dec 2014 by Daniele Merico: adptation to annotation pipeline version rev12 
# ** major changes on 07 Jan 2015 by Daniele Merico: adptation to annotation pipeline version rev22.4, added synonymous conserved / CADD at damaging tier-1
# ** major changes on 13 Feb 2015 by Daniele Merico: adptation to annotation pipeline version rev23, including splicing predictions and hot-zone output
# ** minor changes on 16 Mar 2015 by Daniele Merico: fixed bug with missing ethnic groups
# ** minor changes on 06 Nov 2015 by Daniele Merico: revised splicing prediction criteria, inactivate write all rare variants
# ** general updates on 03 May 2016 by Daniele Merico: grouped stoploss with LOF, missense t2 requires 5 preds, added dbscsnv, removed spidex missense, set spidex cutoffs at -2 (t1) and -4 (t2), ...
# ** ..., relaxed non-coding cutoffs, modified gene prioritization for immune sets, replaced HomHet ratio with Hom percentage

# ** major change on 09 Jun 2016 by BTG: adptation to annotation pipeline version rev25 and HAS
# ** major change on 08 Aug 2016 by BTG: adptation to annotation pipeline version rev26 - spidex cutoffs modified based on DM's suggestions, ClinVar entries processed by TCAG
# ** minor change on 08 Aug 2016 by BTG: adptation to annotation pipeline version rev26.2 - ClinVar header changed from Clinvar_CLNDBN to Clinvar_CLNREF 
# ** minor change on 12 Sep 2016 by BTG: adptation to annotation pipeline version rev26.2 for Illumina GATK pipeline
# ** minor change on 4 Nov 2016 by BTG: change # of lines read to 100K instead of 250K
# ** major change on 28 Mar 2017 by BTG: Secondary findings - ClinVar rev26.2
# ** major change on 7 Nov 2017 by BTG: adptation to annotation pipeline version rev26.5 - fixed 0 filter; removed dbSNP filter for rare 0; fixed splidex cutoffs; UTR added to prioritization; chet_dmg
# ** minor change on 11 Nov 2017 by BTG: changing command line parameters
# ** major change on 23 Oct 2019 by BTG: adptation to annotation pipeline version rev27.1 (removed 1000G from freq filter, added oe instead of ExAC_pli, ClinVar_SIG_Simple for filtering likely pathogenic)
# ** major change on 18 Aug 2020 by Thomas: adptation to annotation pipeline version rev27.4 (using max_frq as for frequency filter, HGMD removed)
# ** minor change on 3 Feb 2021 by Bank: improve the script to handle output path better

################# functions
greplist <- function(targetList, pattern){
  return(as.numeric(which(sapply(sapply(sapply(targetList, strsplit, ","), grep, 
                                        pattern = pattern, ignore.case = T), length) != 0)))
} # not really used
#################

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

# (0) INPUT VARIABLES

# 0.1. variables that need to be set for each genome
input_var_genome.file <- "~/Desktop/PrioritizationPipeline/data/starting_data/NA12878.vcf.gz.annovar.out_SUBSET_rev27.6_hg38.tsv"
input_var_genome.name <- "NA12878"
alt_input_var_genome.name <- ifelse(grepl("^[0-9]", input_var_genome.name),paste("X",input_var_genome.name,".",sep=""),paste(input_var_genome.name,".",sep=""))

## 0.2 output related
output_path.ch <- "/hpf/largeprojects/tcagstor/users/kara.han/PrioritizationPipeline/Results"
output_prefix.ch <- paste(output_path.ch,input_var_genome.name,sep="/")

#print input and output
input_var_genome.file
input_var_genome.name
alt_input_var_genome.name
output_path.ch
output_prefix.ch

# (1) IMPORT 

# 1.2. ANNOTATED TSV

v_full.temp.df <- read.table (input_var_genome.file, sep = "\t", header = T, quote = "\"", comment.char = "", stringsAsFactors = F)

names(v_full.temp.df) <- gsub(alt_input_var_genome.name, "", names(v_full.temp.df))
names(v_full.temp.df)

v_full.df <- subset (v_full.temp.df, subset = (Zygosity!="hom-ref" & Zygosity!="unknown"))
v_full.df$alt_fraction <- v_full.df$AD_ALT / (v_full.df$AD_REF + v_full.df$AD_ALT)

rm(v_full.temp.df)
gc (); gc (); gc ()

# (2) INITIALIZE INTERNAL VARIABLES

# 2.1. standard codes

typeseq_coding.chv <- c ("exonic", "exonic;splicing", "splicing")
typeseq_ncrna.chv  <- c ("ncRNA_exonic", "ncRNA_splicing","ncRNA_exonic;ncRNA_splicing")
typeseq_utr.chv    <- c ("UTR3", "UTR5", "UTR3;UTR5", "UTR5;UTR3")

eff_lof.chv    <- c ("frameshift deletion", "frameshift insertion", "frameshift substitution", "stopgain SNV", "stopgain", "stoploss", "stoploss SNV")
eff_missn.chv  <- c ("nonsynonymous SNV")
eff_other_sub.chv  <- c ("nonframeshift deletion", "nonframeshift insertion", "nonframeshift substitution")
eff_syn.chv <- "synonymous SNV"

# (3) STATS: GENE, GENDER, FULL VARIANT TABLE

# 3.1. Stats 

stats.ls <- list ()

# 3.3. Checks: chromosomes, hom-alts, zigosity for autosomal losses

stats.ls$chromosome_counts <- table (v_full.df$X.CHROM)
#      1     10     11     12     13     14     15     16     17     18     19      2     20     21     22      3      4      5      6      7      8      9      M      X 
# 341487 215961 209890 194853 167913 143544 117024 137304 121358 120210  98560 349942  81281  66911  51157 298080 294063 255533 280677 247367 234423 174982     30 148016 

chr_counts.mx <- t (as.data.frame.matrix (table (v_full.df[, c ("Zygosity", "X.CHROM")])))
# good to have
chr_counts.mx <- cbind (chr_counts.mx, (chr_counts.mx[, "hom-alt"]) / (chr_counts.mx[, "ref-alt"] + chr_counts.mx[, "hom-alt"] + chr_counts.mx[, "alt-alt"]) * 100)
colnames (chr_counts.mx)[ncol (chr_counts.mx)] <- "HomPerc"

stats.ls$chr_zigosity_counts_AllQ_AllSeq_AllFreq <- chr_counts.mx; 
rm (chr_counts.mx)

gc (); gc (); gc ()

# 3.4. Variant number stats

stats.ls$VarN_AllQ_AllSeq_AllFreq <- length (unique (v_full.df$Original_VCFKEY))
stats.ls$VarN_AllQ_Coding_AllFreq <- length (unique (subset (v_full.df, subset = typeseq_priority %in% typeseq_coding.chv, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_AllQ_ncRNA_AllFreq  <- length (unique (subset (v_full.df, subset = typeseq_priority %in% typeseq_ncrna.chv,  select = Original_VCFKEY, drop = T)))

gc (); gc ()

# (4) QUALITY AND RARITY FILTER

# 4.2. High quality, tiers

v_full_hq.df <- subset (v_full.df, subset = (FILTER=="PASS"))

v_full_hq.df$F_Qual          <- 2 # what is F_qual? -> look for change in lines below

zhq2.ix <- with (v_full_hq.df, which (DP >= 10 &
                                        ((Zygosity %in% c("ref-alt","alt-alt") & var_type == "snp" & GQ >= 99 & alt_fraction >= 0.3) | 
                                           (Zygosity %in% c("ref-alt","alt-alt") & var_type != "snp" & GQ >= 90 & alt_fraction >= 0.3) | 
                                           (Zygosity %in% c("hom-alt") & GQ >= 25 & alt_fraction >= 0.8)))) 

v_full_hq.df$F_Qual_tag          <- "LowQuality"
v_full_hq.df$F_Qual_tag[zhq2.ix] <- "OK"

rm (zhq2.ix); 
#rm (v_full.df); 

stats.ls$VarN_Q1_AllSeq_AllFreq <- length (unique (v_full_hq.df$Original_VCFKEY))
stats.ls$VarN_Q2_AllSeq_AllFreq <- length (unique (subset (v_full_hq.df, subset = F_Qual_tag != "LowQuality", select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Coding_AllFreq <- length (unique (subset (v_full_hq.df, subset = typeseq_priority %in% typeseq_coding.chv, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_AllFreq <- length (unique (subset (v_full_hq.df, subset = F_Qual_tag != "LowQuality" & typeseq_priority %in% typeseq_coding.chv, select = Original_VCFKEY, drop = T))) 

stats.ls$VarN_Q1_ncRNA_AllFreq <- length (unique (subset (v_full_hq.df, subset = typeseq_priority %in% typeseq_ncrna.chv, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_ncRNA_AllFreq <- length (unique (subset (v_full_hq.df, subset = F_Qual_tag != "LowQuality" & typeseq_priority %in% typeseq_ncrna.chv, select = Original_VCFKEY, drop = T))) 

stats.ls$VarN_Q1_AllSeq_AllFreq_Xhom <- length (unique (subset (v_full_hq.df, subset = X.CHROM == "X" & Zygosity == "hom-alt",                              select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_AllSeq_AllFreq_Xhom <- length (unique (subset (v_full_hq.df, subset = X.CHROM == "X" & Zygosity == "hom-alt" & F_Qual_tag != "LowQuality", select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_AllSeq_AllFreq_Xhet <- length (unique (subset (v_full_hq.df, subset = X.CHROM == "X" & Zygosity == "ref-alt",                              select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_AllSeq_AllFreq_Xhet <- length (unique (subset (v_full_hq.df, subset = X.CHROM == "X" & Zygosity == "ref-alt" & F_Qual_tag != "LowQuality", select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_AllSeq_AllFreq_Hom  <- length (unique (subset (v_full_hq.df, subset = Zygosity == "hom-alt",                                           select = Original_VCFKEY, drop = T)))     
stats.ls$VarN_Q2_AllSeq_AllFreq_Hom  <- length (unique (subset (v_full_hq.df, subset = Zygosity == "hom-alt"              & F_Qual_tag != "LowQuality", select = Original_VCFKEY, drop = T))) 

stats.ls$VarN_Q1_AllSeq_AllFreq_HetR <- length (unique (subset (v_full_hq.df, subset = Zygosity == "ref-alt",                                           select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q2_AllSeq_AllFreq_HetR <- length (unique (subset (v_full_hq.df, subset = Zygosity == "ref-alt"              & F_Qual_tag != "LowQuality", select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q1_AllSeq_AllFreq_HetA <- length (unique (subset (v_full_hq.df, subset = Zygosity == "alt-alt",                                           select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q2_AllSeq_AllFreq_HetA <- length (unique (subset (v_full_hq.df, subset = Zygosity == "alt-alt"              & F_Qual_tag != "LowQuality", select = Original_VCFKEY, drop = T))) 

chr_counts.mx <- t (as.data.frame.matrix (table (v_full_hq.df[, c ("Zygosity", "X.CHROM")])))

chr_counts.mx <- cbind (chr_counts.mx, (chr_counts.mx[, "hom-alt"] ) / (chr_counts.mx[, "ref-alt"] + chr_counts.mx[, "hom-alt"] + chr_counts.mx[, "alt-alt"]) * 100)
colnames (chr_counts.mx)[ncol (chr_counts.mx)] <- "HomPerc"

stats.ls$chr_zigosity_counts_Q1_AllSeq_AllFreq <- chr_counts.mx; rm (chr_counts.mx)

#stats.ls$chr_ploidy_counts_Q1_AllSeq_AllFreq <- t (as.data.frame.matrix (table (v_full_hq.df[, c ("calledPloidy", "chr")])))

# 4.3. Frequency filters
#full set, for secondary findings only
#freq_max
#127     gnomAD_exome_freq_max
#128     gnomAD_genome_freq_max

#obtain variants with at most 0.05 or NA freq_max (rare variants)
v_full_r05.df <- subset (v_full.df, 
                         subset = ((is.na (gnomAD_exome_freq_max)  | gnomAD_exome_freq_max  <= 0.05) &
                                     (is.na (gnomAD_genome_freq_max)  | gnomAD_genome_freq_max  <= 0.05)))

v_full_r05.df$F_Coding <- "Other"
v_full_r05.df$F_Coding[which (v_full_r05.df$typeseq_priority %in% typeseq_ncrna.chv) ] <- "ncRNA"
v_full_r05.df$F_Coding[which (v_full_r05.df$typeseq_priority %in% typeseq_coding.chv)] <- "Coding"

rm (v_full.df);
#obtain high quality variants with at most 0.05 or NA freq_max
v_full_hq_r05.df <- subset (v_full_hq.df, 
                            subset = ((is.na (gnomAD_exome_freq_max)  | gnomAD_exome_freq_max  <= 0.05) &
                                        (is.na (gnomAD_genome_freq_max)  | gnomAD_genome_freq_max  <= 0.05)))
#obtain high quality variants with various freq_max thresholds
zr01.ix <- with (v_full_hq_r05.df, 
                 which ((is.na (gnomAD_exome_freq_max) | gnomAD_exome_freq_max <= 0.01) & 
                          (is.na (gnomAD_genome_freq_max) | gnomAD_genome_freq_max <= 0.01)))

zr005.ix <- with (v_full_hq_r05.df, 
                  which ((is.na (gnomAD_exome_freq_max) | gnomAD_exome_freq_max <= 0.005) & 
                           (is.na (gnomAD_genome_freq_max) | gnomAD_genome_freq_max <= 0.005)))

zr0015.ix <- with (v_full_hq_r05.df, 
                   which ((is.na (gnomAD_exome_freq_max) | gnomAD_exome_freq_max <= 0.0015) & 
                            (is.na (gnomAD_genome_freq_max) | gnomAD_genome_freq_max <= 0.0015)))

zr00.ix <- with (v_full_hq_r05.df, 
                 which ((is.na (gnomAD_exome_freq_max) | gnomAD_exome_freq_max <= 0) & 
                          (is.na (gnomAD_genome_freq_max) | gnomAD_genome_freq_max <= 0)))

v_full_hq_r05.df$F_Rare            <- 0.05
v_full_hq_r05.df$F_Rare[zr01.ix  ] <- 0.01
v_full_hq_r05.df$F_Rare[zr005.ix ] <- 0.005
v_full_hq_r05.df$F_Rare[zr0015.ix] <- 0.0015
v_full_hq_r05.df$F_Rare[zr00.ix  ] <- 0
# of high quality rare variants across different freq_max thresholds
stats.ls$VarN_Q1_AllSeq_Rare050  <- length (unique (v_full_hq_r05.df$Original_VCFKEY))
stats.ls$VarN_Q2_AllSeq_Rare050  <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality",                   select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q1_AllSeq_Rare010  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.01,                               select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_AllSeq_Rare010  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.01 & F_Qual_tag != "LowQuality",  select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q1_AllSeq_Rare005  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.005,                              select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_AllSeq_Rare005  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.005 & F_Qual_tag != "LowQuality", select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q1_AllSeq_Rare0015 <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.0015,                             select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_AllSeq_Rare0015 <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.0015 & F_Qual_tag != "LowQuality",select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q1_AllSeq_Rare000  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare == 0,                                  select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_AllSeq_Rare000  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare == 0 & F_Qual_tag != "LowQuality",     select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Coding_Rare050  <- length (unique (subset (v_full_hq_r05.df, subset = typeseq_priority %in% typeseq_coding.chv,                                                 select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare050  <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual == 2 & typeseq_priority %in% typeseq_coding.chv,                                   select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q1_Coding_Rare010  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.01 & typeseq_priority  %in% typeseq_coding.chv,                               select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare010  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.01 & F_Qual_tag != "LowQuality" & typeseq_priority %in% typeseq_coding.chv,   select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q1_Coding_Rare005  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.005 & typeseq_priority %in% typeseq_coding.chv,                               select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare005  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.005 & F_Qual_tag != "LowQuality" & typeseq_priority %in% typeseq_coding.chv,  select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q1_Coding_Rare0015 <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.0015 & typeseq_priority %in% typeseq_coding.chv,                              select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare0015 <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.0015 & F_Qual_tag != "LowQuality" & typeseq_priority %in% typeseq_coding.chv, select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q1_Coding_Rare000  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare == 0 & typeseq_priority %in% typeseq_coding.chv,                                   select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare000  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare == 0 & F_Qual_tag != "LowQuality" & typeseq_priority %in% typeseq_coding.chv,      select = Original_VCFKEY, drop = T))) 

stats.ls$VarN_Q1_ncRNA_Rare050  <- length (unique (subset (v_full_hq_r05.df, subset = typeseq_priority %in% typeseq_ncrna.chv,                                                   select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_ncRNA_Rare050  <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual == 2 & typeseq_priority %in% typeseq_ncrna.chv,                                     select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q1_ncRNA_Rare010  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.01 & typeseq_priority %in% typeseq_ncrna.chv,                                  select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_ncRNA_Rare010  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.01 & F_Qual_tag != "LowQuality" & typeseq_priority %in% typeseq_ncrna.chv,     select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q1_ncRNA_Rare005  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.005 & typeseq_priority %in% typeseq_ncrna.chv,                                 select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_ncRNA_Rare005  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.005 & F_Qual_tag != "LowQuality" & typeseq_priority %in% typeseq_ncrna.chv,    select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q1_ncRNA_Rare0015 <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.0015 & typeseq_priority %in% typeseq_ncrna.chv,                                select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_ncRNA_Rare0015 <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare <= 0.0015 & F_Qual_tag != "LowQuality" & typeseq_priority %in% typeseq_ncrna.chv,   select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q1_ncRNA_Rare000  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare == 0 & typeseq_priority %in% typeseq_ncrna.chv,                                     select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_ncRNA_Rare000  <- length (unique (subset (v_full_hq_r05.df, subset = F_Rare == 0 & F_Qual_tag != "LowQuality" & typeseq_priority %in% typeseq_ncrna.chv,        select = Original_VCFKEY, drop = T))) 

stats.ls$VarN_Q1_AllSeq_Rare050_Xhom <- length (unique (subset (v_full_hq_r05.df, subset = X.CHROM == "X" & Zygosity == "hom-alt",                                select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_AllSeq_Rare050_Xhom <- length (unique (subset (v_full_hq_r05.df, subset = X.CHROM == "X" & Zygosity == "hom-alt" & F_Qual_tag != "LowQuality",   select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q1_AllSeq_Rare050_Hom  <- length (unique (subset (v_full_hq_r05.df, subset = Zygosity == "hom-alt",                                             select = Original_VCFKEY, drop = T)))     
stats.ls$VarN_Q2_AllSeq_Rare050_Hom  <- length (unique (subset (v_full_hq_r05.df, subset = Zygosity == "hom-alt"     & F_Qual_tag != "LowQuality",            select = Original_VCFKEY, drop = T)))     
stats.ls$VarN_Q1_AllSeq_Rare050_HetR <- length (unique (subset (v_full_hq_r05.df, subset = Zygosity == "ref-alt",                                             select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q2_AllSeq_Rare050_HetR <- length (unique (subset (v_full_hq_r05.df, subset = Zygosity == "ref-alt" & F_Qual_tag != "LowQuality",                select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q1_AllSeq_Rare050_HetA <- length (unique (subset (v_full_hq_r05.df, subset = Zygosity == "alt-alt",                                             select = Original_VCFKEY, drop = T))) 
stats.ls$VarN_Q2_AllSeq_Rare050_HetA <- length (unique (subset (v_full_hq_r05.df, subset = Zygosity == "alt-alt" & F_Qual_tag != "LowQuality",                select = Original_VCFKEY, drop = T))) 

chr_counts.mx <- t (as.data.frame.matrix (table (v_full_hq_r05.df[, c ("Zygosity", "X.CHROM")])))

chr_counts.mx <- cbind (chr_counts.mx, (chr_counts.mx[, "hom-alt"]) / (chr_counts.mx[, "ref-alt"] + chr_counts.mx[, "hom-alt"] + chr_counts.mx[, "alt-alt"]) * 100)
colnames (chr_counts.mx)[ncol (chr_counts.mx)] <- "HomPerc"

stats.ls$chr_zigosity_counts_Q1_AllSeq_Rare050 <- chr_counts.mx; rm (chr_counts.mx)

#stats.ls$chr_ploidy_counts_Q1_AllSeq_Rare050 <- t (as.data.frame.matrix (table (v_full_hq_r05.df[, c ("calledPloidy", "chr")])))

# 4.4. Coding tag

v_full_hq_r05.df$F_Coding <- "Other"
v_full_hq_r05.df$F_Coding[which (v_full_hq_r05.df$typeseq_priority %in% typeseq_ncrna.chv) ] <- "ncRNA"
v_full_hq_r05.df$F_Coding[which (v_full_hq_r05.df$typeseq_priority %in% typeseq_coding.chv)] <- "Coding"

rm (v_full_hq.df); gc (); gc ()

# (5) GENES
#AD = Autosomal Dominant; remove NAs
g_dom.eg <- setdiff (c (v_full_hq_r05.df$entrez_id[grep (v_full_hq_r05.df$HPO, pattern = "@AD")], subset (v_full_hq_r05.df, subset = CGD_inheritance == "AD", select = entrez_id, drop = T)), NA)
#remove NAs from 1) subset of genes with OMIM/HPO/CGD annotations; 2) subset of genes with MPO annotations
g_phm1.eg <- setdiff (c (subset (v_full_hq_r05.df, subset = (! is.na (omim_phenotype) & omim_phenotype != "") | (! is.na (HPO) & HPO != "") | (! is.na (CGD_disease) & CGD_disease != ""), select = entrez_id, drop = T)), NA)
g_phm0.eg <- setdiff (c (subset (v_full_hq_r05.df, subset = (! is.na (MPO) & MPO != ""), select = entrez_id, drop = T)), NA)
# of genes with OMIM/HPO/CGD annotations; # of genes with AD mode of inheritance
stats.ls$GeneN_Pheno_rare050 <- length (intersect (g_phm1.eg, v_full_hq_r05.df$entrez_id))
stats.ls$GeneN_Dominant_rare050 <- length (intersect (g_dom.eg,  v_full_hq_r05.df$entrez_id))

# (6) DEFINE DAMAGE
v_full_r05.df$F_DamageType <- "NO"
v_full_r05.df$F_DamageTier <- 0
v_full_r05.df$F_S_DamageType <- "NO"

v_full_hq_r05.df$F_DamageType <- "NO"
v_full_hq_r05.df$F_DamageTier <- 0
v_full_hq_r05.df$F_S_DamageType <- "NO"

# 6.1. Coding LOF
zLOF.ix <- with (v_full_r05.df, which (F_Coding == "Coding" & (effect_priority %in% eff_lof.chv | typeseq_priority %in% c ("splicing", "exonic;splicing"))))
v_full_r05.df$F_DamageType[zLOF.ix] <- "LOF"
v_full_r05.df$F_DamageTier[zLOF.ix] <- 2
rm (zLOF.ix)

zsLOF.ix <- with (v_full_r05.df, which (F_Coding == "Coding" & (effect_priority %in% eff_lof.chv | (typeseq_priority %in% c ("splicing", "exonic;splicing") & distance_spliceJunction < 3))))
v_full_r05.df$F_S_DamageType[zsLOF.ix] <- "LOF"
rm (zsLOF.ix)
#high quality
zLOF.ix <- with (v_full_hq_r05.df, which (F_Coding == "Coding" & (effect_priority %in% eff_lof.chv | typeseq_priority %in% c ("splicing", "exonic;splicing"))))
v_full_hq_r05.df$F_DamageType[zLOF.ix] <- "LOF"
v_full_hq_r05.df$F_DamageTier[zLOF.ix] <- 2
rm (zLOF.ix)

zsLOF.ix <- with (v_full_hq_r05.df, which (F_Coding == "Coding" & (effect_priority %in% eff_lof.chv | (typeseq_priority %in% c ("splicing", "exonic;splicing") & distance_spliceJunction < 3))))
v_full_hq_r05.df$F_S_DamageType[zsLOF.ix] <- "LOF"
rm (zsLOF.ix)

# 6.2. Missense
#the scores represent predicted protein impact; greater or less than a threshold indicates damaging
zMsDmg.mx <- matrix (data = 0, nrow = nrow (v_full_hq_r05.df), ncol = 6) # 6 is the number of prediction methods
zMsDmg.mx[, 1] <- with (v_full_hq_r05.df, as.numeric (sift_score     <  0.05))
zMsDmg.mx[, 2] <- with (v_full_hq_r05.df, as.numeric (polyphen_score >= 0.90)) # might need to change to 0.95 based on new documentation
zMsDmg.mx[, 3] <- with (v_full_hq_r05.df, as.numeric (ma_score       >= 1.90)) # might need to change to 2
zMsDmg.mx[, 4] <- with (v_full_hq_r05.df, as.numeric (phylopMam_avg >= 2.30)) # does not factor structural variation
zMsDmg.mx[, 5] <- with (v_full_hq_r05.df, as.numeric (phylopVert100_avg >= 4.0)) # does not factor structural variation
zMsDmg.mx[, 6] <- with (v_full_hq_r05.df, as.numeric (CADD_phred     >= 15))
#zMsDmg.mx[, 7] <- with (v_full_hq_r05.df, as.numeric (mt_score       >= 0.5))
zMsDmg.mx[is.na (zMsDmg.mx)] <- 0
zMsDmg.nv <- apply (zMsDmg.mx, 1, sum)
#eff_missn.chv = c("nonsynonymous SNVs") 
zMsDmg1.ix <- with (v_full_hq_r05.df, which (F_Coding == "Coding" & 
                                               (effect_priority %in% eff_missn.chv &  
                                                  (zMsDmg.nv >= 2)))) # Why use 2 and 4 here?
zMsDmg2.ix <- with (v_full_hq_r05.df, which (F_Coding == "Coding" & 
                                               (effect_priority %in% eff_missn.chv & 
                                                  (zMsDmg.nv >= 4))))
v_full_hq_r05.df$F_DamageType[zMsDmg1.ix] <- "Missense"
v_full_hq_r05.df$F_DamageTier[zMsDmg1.ix] <- 1 # Tier 1
v_full_hq_r05.df$F_DamageTier[zMsDmg2.ix] <- 2 # Tier 2
rm (zMsDmg1.ix, zMsDmg2.ix)

# 5.3. Other coding 

zOthDmg1.ix <- with (v_full_hq_r05.df, which (F_Coding == "Coding" & (
  (effect_priority %in% eff_other_sub.chv & ((phylopMam_avg >= 1.2 | phylopVert100_avg >= 2.5 | CADD_phred >= 13.5) & is.na (dbsnp_common) )) | 
    (effect_priority %in% eff_other_sub.chv & ((phylopMam_avg >= 1.5 | phylopVert100_avg >= 2.0 | CADD_phred >= 13.0) & is.na (dbsnp) & is.na (dbsnp_region))) ) ))
v_full_hq_r05.df$F_DamageType[zOthDmg1.ix] <- "OtherC"
v_full_hq_r05.df$F_DamageTier[zOthDmg1.ix] <- 1

zOthDmg2.ix <- with (v_full_hq_r05.df, which (F_Coding == "Coding" & (
  (effect_priority %in% eff_other_sub.chv & ((phylopMam_avg >= 2.0 | phylopVert100_avg >= 3.5 | CADD_phred >= 14  ) & is.na (dbsnp_common) )) | 
    (effect_priority %in% eff_other_sub.chv & ((phylopMam_avg >= 1.5 | phylopVert100_avg >= 2.5 | CADD_phred >= 13.5) & is.na (dbsnp) & is.na (dbsnp_region))) ) ))
v_full_hq_r05.df$F_DamageType[zOthDmg2.ix] <- "OtherC"
v_full_hq_r05.df$F_DamageTier[zOthDmg2.ix] <- 2

rm (zOthDmg2.ix, zOthDmg1.ix)

# 5.5. Splicing predictions
zSpNegT2.ix <- with (v_full_hq_r05.df, which (
  ( (spliceAI_DS_AG > 0.5 & abs (spliceAI_DP_AG) <= 50) | # spliceAI Delta Score: probability of the variant being splice-altering
      (spliceAI_DS_AL > 0.5 & abs (spliceAI_DP_AL) <= 50) |
      (spliceAI_DS_DG > 0.5 & abs (spliceAI_DP_DG) <= 50) |
      (spliceAI_DS_DL > 0.5 & abs (spliceAI_DP_DL) <= 50) |
      (dbscSNV_ADA_SCORE > 0.6 | dbscSNV_RF_SCORE > 0.6) | # splice site prediction scores from dbscSNV
      (spx_sequence_event_type == "splice_region_variant" & spx_dpsi <= -20) | 
      (spx_sequence_event_type %in% c ("synonymous_variant", "missense_variant") & spx_ss_dist <= 3 & spx_dpsi <= -20) |
      (spx_sequence_event_type %in% c ("synonymous_variant", "missense_variant") & spx_ss_dist >  3 & spx_dpsi <= -20) | 
      (spx_sequence_event_type == "coding_transcript_intron_variant" & spx_dpsi <= -20) ) & 
    ! F_DamageType %in% c ("LOF") ) ) # NOTE: spx_sequence_event_type and spx_ss_dist only found in old version
v_full_hq_r05.df$F_DamageType[zSpNegT2.ix] <- "Splc"
v_full_hq_r05.df$F_DamageTier[zSpNegT2.ix] <- 2
# Where can I find info on SPIDEX?
zSpNegT1.ix <- with (v_full_hq_r05.df, which (
  ( (spliceAI_DS_AG > 0.2 & abs (spliceAI_DP_AG) <= 50) |
      (spliceAI_DS_AL > 0.2 & abs (spliceAI_DP_AL) <= 50) |
      (spliceAI_DS_DG > 0.2 & abs (spliceAI_DP_DG) <= 50) |
      (spliceAI_DS_DL > 0.2 & abs (spliceAI_DP_DL) <= 50) |
      (spx_sequence_event_type == "splice_region_variant" & spx_dpsi <= -10) | 
      (spx_sequence_event_type %in% c ("synonymous_variant", "missense_variant") & spx_ss_dist <= 3 & spx_dpsi <= -10) |
      (spx_sequence_event_type %in% c ("synonymous_variant", "missense_variant") & spx_ss_dist >  3 & spx_dpsi <= -10) | 
      (spx_sequence_event_type == "coding_transcript_intron_variant" & spx_dpsi <= -10) ) & 
    ! F_DamageType %in% c ("LOF", "Splc") & ! (F_DamageType %in% "Missense" & F_DamageTier == 2)) )
v_full_hq_r05.df$F_DamageType[zSpNegT1.ix] <- "Splc"
v_full_hq_r05.df$F_DamageTier[zSpNegT1.ix] <- 1

# 5.5. UTR

zUtrDmg2.ix <- with (v_full_hq_r05.df, which (
  typeseq_priority %in% typeseq_utr.chv & 
    (phylopMam_avg >= 1.5 | CADD_phred >= 15) & # stronger conservation based on Mam_avg
    (! is.na (phastCons_placental)) ))
zUtrDmg1.ix <- with (v_full_hq_r05.df, which (
  typeseq_priority %in% typeseq_utr.chv & 
    (phylopMam_avg >  0 | CADD_phred >= 12.5) & 
    (! is.na (phastCons_placental)) ))

v_full_hq_r05.df$F_DamageType[zUtrDmg1.ix] <- "Utr"
v_full_hq_r05.df$F_DamageTier[zUtrDmg1.ix] <- 1
v_full_hq_r05.df$F_DamageTier[zUtrDmg2.ix] <- 2
rm (zUtrDmg1.ix, zUtrDmg2.ix)

# 5.6. Non-coding

zNcrDmg1.ix <- with (v_full_hq_r05.df, which (F_Coding == "ncRNA" & ! gene_desc %in% "pseudogene" & F_DamageTier == 0 & (
  CADD_phred >= 12 | ((phylopMam_avg >= 1.25 | phylopVert100_avg >= 2.00) & (! is.na (phastCons_placental))) )))
zNcrDmg2.ix <- with (v_full_hq_r05.df, which (F_Coding == "ncRNA" & ! gene_desc %in% "pseudogene" & F_DamageTier <= 1 & (
  ((CADD_phred >= 15 | phylopMam_avg >= 2.30 | phylopVert100_avg >= 4.00) & (! is.na (phastCons_placental))) )))
v_full_hq_r05.df$F_DamageType[zNcrDmg1.ix] <- "DmgNcRNA"
v_full_hq_r05.df$F_DamageTier[zNcrDmg1.ix] <- 1
v_full_hq_r05.df$F_DamageTier[zNcrDmg2.ix] <- 2
rm (zNcrDmg1.ix, zNcrDmg2.ix)

# 5.7. Stats

stats.ls$VarN_Q1_Coding_Rare010_LOF <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_DamageType == "LOF", select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare010_LOF <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_DamageType == "LOF", select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Coding_Rare010_MissDmgT1 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_DamageType == "Missense" & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare010_MissDmgT1 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_DamageType == "Missense" & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Coding_Rare010_MissDmgT2 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_DamageType == "Missense" & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare010_MissDmgT2 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_DamageType == "Missense" & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Coding_Rare010_SplcT1 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_DamageType == "Splc" & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare010_SplcT1 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_DamageType == "Splc" & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Coding_Rare010_SplcT2 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_DamageType == "Splc" & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare010_SplcT2 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_DamageType == "Splc" & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

stats.ls$Coding_HQ1_Rare010_OtherDmg <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_DamageType == "OtherC", select = Original_VCFKEY, drop = T)))
stats.ls$Coding_HQ2_Rare010_OtherDmg <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_DamageType == "OtherC", select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_ncRNA_Rare010_DmgT1 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_DamageType == "DmgNcRNA" & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_ncRNA_Rare010_DmgT1 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_DamageType == "DmgNcRNA" & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_ncRNA_Rare010_DmgT2 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_DamageType == "DmgNcRNA" & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_ncRNA_Rare010_DmgT2 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_DamageType == "DmgNcRNA" & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

# (6) MORE GENE ANNOTATIONS

# 6.1. HPO dominant

v_full_hq_r05.df$G_AXD_HPO <- 0
zDomHpo.ix <- grep (v_full_hq_r05.df$HPO, pattern = "@AD")
v_full_hq_r05.df$G_AXD_HPO[zDomHpo.ix] <- 1
rm (zDomHpo.ix)

# 6.2. CGD dominant

# ** NOTE: currently CGD does not have "XD" as a possible mode of inheritance, but if it becomes available it should be included
# still no "XD" in the new version of CGD
v_full_hq_r05.df$G_AXD_CGD <- 0
zDomCgd.ix <- with (v_full_hq_r05.df, which (CGD_inheritance == "AD"))
v_full_hq_r05.df$G_AXD_CGD[zDomCgd.ix] <- 1
rm (zDomCgd.ix)

stats.ls$VarN_Q1_Coding_Rare010_Dmg_AXD_HPO <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_Coding == "Coding" & F_DamageType != "NO" & (G_AXD_HPO == 1), select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare010_Dmg_AXD_HPO <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_Coding == "Coding" & F_DamageType != "NO" & (G_AXD_HPO == 1), select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Coding_Rare010_DmG_AXD_CGD <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_Coding == "Coding" & F_DamageType != "NO" & (G_AXD_CGD == 1), select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare010_DmG_AXD_CGD <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_Coding == "Coding" & F_DamageType != "NO" & (G_AXD_CGD == 1), select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Coding_Rare010_Dmg_AXD_All <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_Coding == "Coding" & F_DamageType != "NO" & (G_AXD_HPO == 1 | G_AXD_CGD == 1), select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare010_Dmg_AXD_All <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_Coding == "Coding" & F_DamageType != "NO" & (G_AXD_HPO == 1 | G_AXD_CGD == 1), select = Original_VCFKEY, drop = T)))

# 6.3. Phenotype tiers

v_full_hq_r05.df$F_PhenoTier <- 0

zPhT1.ix <- with (v_full_hq_r05.df, which (entrez_id %in% g_phm0.eg)) # genes with MPO annotations
zPhT2.ix <- with (v_full_hq_r05.df, which (entrez_id %in% g_phm1.eg)) # genes with OMIM/HPO/CGD annotations

v_full_hq_r05.df$F_PhenoTier[zPhT1.ix] <- 1
v_full_hq_r05.df$F_PhenoTier[zPhT2.ix] <- 2

stats.ls$VarN_Q1_Coding_Rare010_Dmg_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_DamageType != "NO" & F_PhenoTier >= 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare010_Dmg_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_DamageType != "NO" & F_PhenoTier >= 1, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Coding_Rare010_Dmg_PhenoTier2 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual >= 1 & F_Rare <= 0.01 & F_DamageType != "NO" & F_PhenoTier >= 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Coding_Rare010_Dmg_PhenoTier2 <- length (unique (subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & F_Rare <= 0.01 & F_DamageType != "NO" & F_PhenoTier >= 2, select = Original_VCFKEY, drop = T)))

# (7) FINAL GROUPS

# 7.1. RECESSIVE -- HOMOZYGOUS

z.Hom.ix <- with (v_full_hq_r05.df, which (F_Rare <= 0.05 & F_DamageType != "NO" & Zygosity == "hom-alt"))
v_full_hq_r05.df$FM_HOM <- 0
v_full_hq_r05.df$FM_HOM[z.Hom.ix] <- 1
rm (z.Hom.ix)

stats.ls$VarN_Q1_Rare050_DmgT2_Hom_PhenoTier0 <- length (unique (subset (v_full_hq_r05.df, subset = FM_HOM == 1 & F_Qual >= 1 & F_PhenoTier >= 0 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare050_DmgT2_Hom_PhenoTier0 <- length (unique (subset (v_full_hq_r05.df, subset = FM_HOM == 1 & F_Qual_tag != "LowQuality" & F_PhenoTier >= 0 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Rare050_DmgT2_Hom_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = FM_HOM == 1 & F_Qual >= 1 & F_PhenoTier >= 1 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare050_DmgT2_Hom_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = FM_HOM == 1 & F_Qual_tag != "LowQuality" & F_PhenoTier >= 1 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Rare050_DmgT2_Hom_PhenoTier2 <- length (unique (subset (v_full_hq_r05.df, subset = FM_HOM == 1 & F_Qual >= 1 & F_PhenoTier >= 2 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare050_DmgT2_Hom_PhenoTier2 <- length (unique (subset (v_full_hq_r05.df, subset = FM_HOM == 1 & F_Qual_tag != "LowQuality" & F_PhenoTier >= 2 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

# 7.2. X-LINKED HAPLOID

z.Xhom.ix <- with (v_full_hq_r05.df, which (F_Rare <= 0.05 & F_DamageType != "NO" & Zygosity == "hom-alt" & X.CHROM == "X"))
v_full_hq_r05.df$FM_XHAP <- 0
v_full_hq_r05.df$FM_XHAP[z.Xhom.ix] <- 1
rm (z.Xhom.ix)

stats.ls$VarN_Q1_Rare050_DmgT2_Xhom_PhenoTier0 <- length (unique (subset (v_full_hq_r05.df, subset = FM_XHAP == 1 & F_Qual >= 1 & F_PhenoTier >= 0 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare050_DmgT2_Xhom_PhenoTier0 <- length (unique (subset (v_full_hq_r05.df, subset = FM_XHAP == 1 & F_Qual_tag != "LowQuality" & F_PhenoTier >= 0 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Rare050_DmgT2_Xhom_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = FM_XHAP == 1 & F_Qual >= 1 & F_PhenoTier >= 1 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare050_DmgT2_Xhom_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = FM_XHAP == 1 & F_Qual_tag != "LowQuality" & F_PhenoTier >= 1 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Rare050_DmgT2_Xhom_PhenoTier2 <- length (unique (subset (v_full_hq_r05.df, subset = FM_XHAP == 1 & F_Qual >= 1 & F_PhenoTier >= 2 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare050_DmgT2_Xhom_PhenoTier2 <- length (unique (subset (v_full_hq_r05.df, subset = FM_XHAP == 1 & F_Qual_tag != "LowQuality" & F_PhenoTier >= 2 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

# 7.3. POTENTIAL COMPOUND HETS

v_full_hq_r05.df$FM_PCHET <- 0

zchet_q1.df <- subset (v_full_hq_r05.df, subset = F_Rare <= 0.05 & F_DamageType != "NO" & F_Qual >= 1, select = c (gene_symbol, Original_VCFKEY))
zchet_q1.df <- zchet_q1.df[! duplicated (zchet_q1.df), ]
zchet_q1_f.df <- subset (zchet_q1.df, subset = gene_symbol %in% zchet_q1.df$gene_symbol[duplicated (zchet_q1.df$gene_symbol)], select = c (Original_VCFKEY, gene_symbol))
v_full_hq_r05.df$FM_PCHET[which (v_full_hq_r05.df$Original_VCFKEY %in% zchet_q1_f.df$Original_VCFKEY & v_full_hq_r05.df$gene_symbol %in% zchet_q1_f.df$gene_symbol)] <- 1
# high quality
zchet_q2.df <- subset (v_full_hq_r05.df, subset = F_Rare <= 0.05 & F_DamageType != "NO" & F_Qual_tag != "LowQuality", select = c (gene_symbol, Original_VCFKEY))
zchet_q2.df <- zchet_q2.df[! duplicated (zchet_q2.df), ]
zchet_q2_f.df <- subset (zchet_q2.df, subset = gene_symbol %in% zchet_q2.df$gene_symbol[duplicated (zchet_q2.df$gene_symbol)], select = c (Original_VCFKEY, gene_symbol))
v_full_hq_r05.df$FM_PCHET[which (v_full_hq_r05.df$Original_VCFKEY %in% zchet_q2_f.df$Original_VCFKEY & v_full_hq_r05.df$gene_symbol %in% zchet_q2_f.df$gene_symbol)] <- 2

v_full_hq_r05.df$FM_PCHET_DMG <- 0 

zchet_q2d2.sy <- subset (v_full_hq_r05.df, subset = FM_PCHET == 2 & F_DamageTier == 2, select = gene_symbol, drop = T)
zchet_qd2f.sy <- zchet_q2d2.sy[which (duplicated (zchet_q2d2.sy))]

v_full_hq_r05.df$FM_PCHET_DMG[with (v_full_hq_r05.df, which (gene_symbol %in% zchet_qd2f.sy & FM_PCHET == 2))] <- 1

rm (zchet_q1.df, zchet_q1_f.df, zchet_q2.df, zchet_q2_f.df)

stats.ls$VarN_Q1_Rare050_DmgT1_CmpHet_PhenoTier0 <- length (unique (subset (v_full_hq_r05.df, subset = FM_PCHET >= 1 & F_PhenoTier >= 0 & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare050_DmgT1_CmpHet_PhenoTier0 <- length (unique (subset (v_full_hq_r05.df, subset = FM_PCHET == 2 & F_PhenoTier >= 0 & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Rare050_DmgT1_CmpHet_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = FM_PCHET >= 1 & F_PhenoTier >= 1 & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare050_DmgT1_CmpHet_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = FM_PCHET == 2 & F_PhenoTier >= 1 & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Rare050_DmgT1_CmpHet_PhenoTier2 <- length (unique (subset (v_full_hq_r05.df, subset = FM_PCHET >= 1 & F_PhenoTier >= 2 & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare050_DmgT1_CmpHet_PhenoTier2 <- length (unique (subset (v_full_hq_r05.df, subset = FM_PCHET == 2 & F_PhenoTier >= 2 & F_DamageTier >= 1, select = Original_VCFKEY, drop = T)))

# 7.4. DOMINANT

v_full_hq_r05.df$FM_AXDOM <- 0

z.Dom.ix <- with (v_full_hq_r05.df, which (F_Rare <= 0.005 & F_DamageType != "NO" & (G_AXD_CGD == 1 | G_AXD_HPO == 1)))
v_full_hq_r05.df$FM_AXDOM <- 0
v_full_hq_r05.df$FM_AXDOM[z.Dom.ix] <- 1
rm (z.Dom.ix)

stats.ls$VarN_Q1_Rare005_DmgT2_AXDom_PhenoTier0 <- length (unique (subset (v_full_hq_r05.df, subset = FM_AXDOM == 1 & F_Qual >= 1 & F_PhenoTier >= 0 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare005_DmgT2_AXDom_PhenoTier0 <- length (unique (subset (v_full_hq_r05.df, subset = FM_AXDOM == 1 & F_Qual_tag != "LowQuality" & F_PhenoTier >= 0 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Rare005_DmgT2_AXDom_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = FM_AXDOM == 1 & F_Qual >= 1 & F_PhenoTier >= 1 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare005_DmgT2_AXDom_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = FM_AXDOM == 1 & F_Qual_tag != "LowQuality" & F_PhenoTier >= 1 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Rare005_DmgT2_AXDom_PhenoTier2 <- length (unique (subset (v_full_hq_r05.df, subset = FM_AXDOM == 1 & F_Qual >= 1 & F_PhenoTier >= 2 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare005_DmgT2_AXDom_PhenoTier2 <- length (unique (subset (v_full_hq_r05.df, subset = FM_AXDOM == 1 & F_Qual_tag != "LowQuality" & F_PhenoTier >= 2 & F_DamageTier == 2, select = Original_VCFKEY, drop = T)))

# 7.5. HETEROZYGOUS HOTZONE
# What is heterozygous hotzone? A high-interest region with lots of heterozygous?
v_full_hq_r05.df$FM_HZ <- 0

z.HZ.ix <- with (v_full_hq_r05.df, which (Zygosity == "ref-alt" & F_Rare <= 0.0015 & (
  (gnomAD_oe_lof_upper < 0.35 & (F_DamageType %in% c ("LOF", "Splc", "OtherC") | (F_DamageType == "Missense" & F_DamageTier == 2) )) )))
v_full_hq_r05.df$FM_HZ <- 0
v_full_hq_r05.df$FM_HZ[z.HZ.ix] <- 1
rm (z.HZ.ix)

stats.ls$VarN_Q1_Rare005_DmgT2_HI_PhenoTier0 <- length (unique (subset (v_full_hq_r05.df, subset = FM_HZ == 1 & F_Qual >= 1 & F_PhenoTier >= 0, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare005_DmgT2_HI_PhenoTier0 <- length (unique (subset (v_full_hq_r05.df, subset = FM_HZ == 1 & F_Qual_tag != "LowQuality" & F_PhenoTier >= 0, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_Q1_Rare005_DmgT2_HI_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = FM_HZ == 1 & F_Qual >= 1 & F_PhenoTier >= 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_Q2_Rare005_DmgT2_HI_PhenoTier1 <- length (unique (subset (v_full_hq_r05.df, subset = FM_HZ == 1 & F_Qual_tag != "LowQuality" & F_PhenoTier >= 1, select = Original_VCFKEY, drop = T)))

# (8) SECONDARY FINDINGS

# 8.0. Flag ClinVar field reported at least once as not pathogenic or probable not pathogenic
#Full set
#non-pathogenic
v_full_r05.df$F_Clinvar_notPathg <- 0
#zClinvNotPathg.ix <- with (v_full_r05.df, which (Clinvar_SIG %in% c ("Benign", "Likely benign", "Uncertain significance")))
zClinvNotPathg.ix <- with (v_full_r05.df, which (Clinvar_SIG_Simple == 0))
v_full_r05.df$F_Clinvar_notPathg[zClinvNotPathg.ix] <- 1

#pathogenic
v_full_r05.df$F_Clinvar_Pathg <- 0
#zClinvPathg.ix <- greplist(v_full_r05.df$Clinvar_SIG, pattern = paste(c("^Likely pathogenic$", "^Pathogenic$","^Pathologic$","^probable-pathogenic$","^probably pathogenic$","^risk factor$"), collapse="|"))
zClinvPathg.ix <- with (v_full_r05.df, which (Clinvar_SIG_Simple == 1))
v_full_r05.df$F_Clinvar_Pathg[zClinvPathg.ix] <- 1

# 8.1. Tier 1
v_full_r05.df$FS1_Select <- 0
zSF1dmg.ix <- with (v_full_r05.df, 
                    which ( 
                      CGD_disease != "" & 
                        (F_S_DamageType == "LOF") |
                        (F_Clinvar_Pathg == 1) 
                    ) ) # Variants with 1) CGD annotations; 2) LOF; 3) pathogenic
v_full_r05.df$FS1_Select[zSF1dmg.ix] <- 1 
rm (zSF1dmg.ix)
###########################################
#HQ set
# non-pathogenic
v_full_hq_r05.df$F_Clinvar_notPathg <- 0
#zClinvNotPathg.ix <- with (v_full_hq_r05.df, which (Clinvar_SIG %in% c ("Benign", "Likely benign", "Uncertain significance")))
zClinvNotPathg.ix <- with (v_full_hq_r05.df, which (Clinvar_SIG_Simple == 0))
v_full_hq_r05.df$F_Clinvar_notPathg[zClinvNotPathg.ix] <- 1

# pathogenic
v_full_hq_r05.df$F_Clinvar_Pathg <- 0
#zClinvPathg.ix <- greplist(v_full_hq_r05.df$Clinvar_SIG, pattern = paste(c("^Likely pathogenic$", "^Pathogenic$","^Pathologic$","^probable-pathogenic$","^probably pathogenic$","^risk factor$"), collapse="|"))
zClinvPathg.ix <- with (v_full_hq_r05.df, which (Clinvar_SIG_Simple == 1))
v_full_hq_r05.df$F_Clinvar_Pathg[zClinvPathg.ix] <- 1

# 8.1. Tier 1 
v_full_hq_r05.df$FS1_Select <- 0 # high quality
zSF1dmg.ix <- with (v_full_hq_r05.df, 
                    which ( 
                      CGD_disease != "" & 
                        (F_S_DamageType == "LOF") |
                        (F_Clinvar_Pathg == 1) 
                    ) )
v_full_hq_r05.df$FS1_Select[zSF1dmg.ix] <- 1
rm (zSF1dmg.ix)
############################################
v_full_hq_r05.df$F_CmpHet_S1 <- 0 # compound heterozygous again?

zchet_q1.df <- subset (v_full_hq_r05.df, subset = F_Qual >= 1 & FS1_Select == 1, select = c (gene_symbol, Original_VCFKEY))
zchet_q1.df <- zchet_q1.df[! duplicated (zchet_q1.df), ]
zchet_q1.Original_VCFKEY <- subset (zchet_q1.df, subset = gene_symbol %in% zchet_q1.df$gene_symbol[duplicated (zchet_q1.df$gene_symbol)], select = Original_VCFKEY, drop = T)
v_full_hq_r05.df$F_CmpHet_S1[which (v_full_hq_r05.df$Original_VCFKEY %in% zchet_q1.Original_VCFKEY)] <- 1

zchet_q2.df <- subset (v_full_hq_r05.df, subset = F_Qual_tag != "LowQuality" & FS1_Select == 1, select = c (gene_symbol, Original_VCFKEY))
zchet_q2.df <- zchet_q2.df[! duplicated (zchet_q2.df), ]
zchet_q2.Original_VCFKEY <- subset (zchet_q2.df, subset = gene_symbol %in% zchet_q2.df$gene_symbol[duplicated (zchet_q2.df$gene_symbol)], select = Original_VCFKEY, drop = T)
v_full_hq_r05.df$F_CmpHet_S1[which (v_full_hq_r05.df$Original_VCFKEY %in% zchet_q2.Original_VCFKEY)] <- 2

rm (zchet_q1.df, zchet_q1.Original_VCFKEY, zchet_q2.df, zchet_q2.Original_VCFKEY)
#############################################
# dominant: pathogenic
z_d_any.ix <- with (v_full_hq_r05.df, which (FS1_Select == 1 & CGD_inheritance == "AD"))

# recessive + hom: pathogenic
z_r_hom.ix <- with (v_full_hq_r05.df, which (FS1_Select == 1 & CGD_inheritance == "AR" & Zygosity == "hom-alt"))

# recessive + potential compound het: pathogenic
z_r_cht.ix <- with (v_full_hq_r05.df, which (FS1_Select == 1 & CGD_inheritance == "AR" & F_CmpHet_S1 >= 1))

# recessive + single het: carrier
z_r_car.ix <- with (v_full_hq_r05.df, which (FS1_Select == 1 & CGD_inheritance == "AR" & Zygosity %in% c ("ref-alt", "alt-alt") & F_CmpHet_S1 == 0))

# X-linked + hom: pathogenic
z_x_hom.ix <- with (v_full_hq_r05.df, which (FS1_Select == 1 & CGD_inheritance == "XL" & Zygosity == "hom-alt" & X.CHROM == "X"))

# X-linked + het: carrier
z_x_car.ix <- with (v_full_hq_r05.df, which (FS1_Select == 1 & CGD_inheritance == "XL" & Zygosity %in% c ("ref-alt", "alt-alt") & X.CHROM == "X"))

# complex + hom / hap: pathogenic
z_cx_homhap.ix <- with (v_full_hq_r05.df, which (FS1_Select == 1 & (! CGD_inheritance %in% c ("AD", "AR", "XL")) & (Zygosity == "hom-alt")))

# complex + potential compound het: pathogenic
z_cx_cht.ix <- with (v_full_hq_r05.df, which (FS1_Select == 1 & (! CGD_inheritance %in% c ("AD", "AR", "XL")) & F_CmpHet_S1 >= 1))

# complex + single het: uncertain
z_cx_unc.ix <- with (v_full_hq_r05.df, which (FS1_Select == 1 & (! CGD_inheritance %in% c ("AD", "AR", "XL")) & Zygosity %in% c ("ref-alt", "alt-alt") & F_CmpHet_S1 == 0))

# ACMG_disease

v_full_hq_r05.df$FS1_AD_Pathg_Any <- 0
v_full_hq_r05.df$FS1_AR_Pathg_Hom <- 0
v_full_hq_r05.df$FS1_AR_Pathg_PotCompHet <- 0
v_full_hq_r05.df$FS1_XL_Pathg_Hap <- 0
v_full_hq_r05.df$FS1_XL_Pathg_Hom <- 0
v_full_hq_r05.df$FS1_CX_Pathg_HomHap     <- 0 
v_full_hq_r05.df$FS1_CX_Pathg_PotCompHet <- 0 
v_full_hq_r05.df$FS1_CX_Uncertain <- 0 
v_full_hq_r05.df$FS1_AR_Carrier   <- 0
v_full_hq_r05.df$FS1_XL_Carrier   <- 0

v_full_hq_r05.df$FS1_AD_Pathg_Any[z_d_any.ix] <- 1
v_full_hq_r05.df$FS1_AR_Pathg_Hom[z_r_hom.ix] <- 1
v_full_hq_r05.df$FS1_AR_Pathg_PotCompHet[z_r_cht.ix] <- 1
v_full_hq_r05.df$FS1_XL_Pathg_Hap[z_x_hom.ix] <- 1 # haploid the same as homozygous in this dataset? 
v_full_hq_r05.df$FS1_XL_Pathg_Hom[z_x_hom.ix] <- 1 
v_full_hq_r05.df$FS1_CX_Pathg_HomHap[z_cx_homhap.ix]  <- 1 
v_full_hq_r05.df$FS1_CX_Pathg_PotCompHet[z_cx_cht.ix] <- 1 
v_full_hq_r05.df$FS1_CX_Uncertain[z_cx_unc.ix] <- 1 
v_full_hq_r05.df$FS1_AR_Carrier[z_r_car.ix]   <- 1
v_full_hq_r05.df$FS1_XL_Carrier[z_x_car.ix]   <- 1

rm (z_d_any.ix, z_r_hom.ix, z_r_cht.ix, z_x_hom.ix, z_cx_homhap.ix, z_cx_cht.ix, z_cx_unc.ix, z_r_car.ix, z_x_car.ix)

# 8.1. Tier 2 and 3

v_full_hq_r05.df$FS2_Select <- 0
zSF2dmg.ix <- with (v_full_hq_r05.df, 
                    which ( 
                      CGD_disease != "" & F_Qual_tag != "LowQuality"  & F_Rare <= 0.01 & (
                        (F_DamageType == "LOF") |
                          (F_Clinvar_Pathg == 1) | 
                          ( F_Clinvar_notPathg == 0) )
                    ) ) # genes with 1) CGD annotations; 2) hq; 3) rare; 4) not non-pathogenic
v_full_hq_r05.df$FS2_Select[zSF2dmg.ix] <- 1
rm (zSF2dmg.ix)

v_full_hq_r05.df$FS3_Select <- 0
zSF3dmg.ix <- with (v_full_hq_r05.df, 
                    which ( 
                      CGD_disease != "" & F_Qual_tag != "LowQuality"  & F_Rare <= 0.01 & F_DamageType != "NO" & (
                        (F_DamageType == "LOF") |
                          (F_Clinvar_Pathg == 1) | 
                          ( F_Clinvar_notPathg == 0) )
                    ) ) # genes with 1) CGD annotations; 2) hq; 3) rare; 4) not non-pathogenic; 5) damaging
v_full_hq_r05.df$FS3_Select[zSF3dmg.ix] <- 1
rm (zSF3dmg.ix)

stats.ls$VarN_FS1_Q1_Rare050_Tot     <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS2_Q2_Rare010_Tot     <- length (unique (subset (v_full_hq_r05.df, subset = FS2_Select == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS3_Q2_Rare010_Dmg_Tot <- length (unique (subset (v_full_hq_r05.df, subset = FS3_Select == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS1_Q1_Rare050_Tot     <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_FS1_Q1_Rare050_AD_Pathg_Any        <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1 & FS1_AD_Pathg_Any        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS1_Q1_Rare050_AR_Pathg_Hom        <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1 & FS1_AR_Pathg_Hom        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS1_Q1_Rare050_AR_Pathg_PotCompHet <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1 & FS1_AR_Pathg_PotCompHet == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS1_Q1_Rare050_XL_Pathg_Hap        <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1 & FS1_XL_Pathg_Hap        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS1_Q1_Rare050_XL_Pathg_Hom        <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1 & FS1_XL_Pathg_Hom        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS1_Q1_Rare050_CX_Pathg_HomHap     <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1 & FS1_CX_Pathg_HomHap     == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS1_Q1_Rare050_CX_Pathg_PotCompHet <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1 & FS1_CX_Pathg_PotCompHet == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS1_Q1_Rare050_CX_Uncertain        <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1 & FS1_CX_Uncertain        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS1_Q1_Rare050_AR_Carrier          <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1 & FS1_AR_Carrier          == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS1_Q1_Rare050_XL_Carrier          <- length (unique (subset (v_full_hq_r05.df, subset = FS1_Select == 1 & FS1_XL_Carrier          == 1, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_FS2_Q2_Rare010_AD_Pathg_Any        <- length (unique (subset (v_full_hq_r05.df, subset = FS2_Select == 1 & FS1_AD_Pathg_Any        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS2_Q2_Rare010_AR_Pathg_Hom        <- length (unique (subset (v_full_hq_r05.df, subset = FS2_Select == 1 & FS1_AR_Pathg_Hom        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS2_Q2_Rare010_AR_Pathg_PotCompHet <- length (unique (subset (v_full_hq_r05.df, subset = FS2_Select == 1 & FS1_AR_Pathg_PotCompHet == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS2_Q2_Rare010_XL_Pathg_Hap        <- length (unique (subset (v_full_hq_r05.df, subset = FS2_Select == 1 & FS1_XL_Pathg_Hap        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS2_Q2_Rare010_XL_Pathg_Hom        <- length (unique (subset (v_full_hq_r05.df, subset = FS2_Select == 1 & FS1_XL_Pathg_Hom        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS2_Q2_Rare010_CX_Pathg_HomHap     <- length (unique (subset (v_full_hq_r05.df, subset = FS2_Select == 1 & FS1_CX_Pathg_HomHap     == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS2_Q2_Rare010_CX_Pathg_PotCompHet <- length (unique (subset (v_full_hq_r05.df, subset = FS2_Select == 1 & FS1_CX_Pathg_PotCompHet == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS2_Q2_Rare010_CX_Uncertain        <- length (unique (subset (v_full_hq_r05.df, subset = FS2_Select == 1 & FS1_CX_Uncertain        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS2_Q2_Rare010_AR_Carrier          <- length (unique (subset (v_full_hq_r05.df, subset = FS2_Select == 1 & FS1_AR_Carrier          == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS2_Q2_Rare010_XL_Carrier          <- length (unique (subset (v_full_hq_r05.df, subset = FS2_Select == 1 & FS1_XL_Carrier          == 1, select = Original_VCFKEY, drop = T)))

stats.ls$VarN_FS3_Q2_Rare010_Dmg_AD_Pathg_Any        <- length (unique (subset (v_full_hq_r05.df, subset = FS3_Select == 1 & FS1_AD_Pathg_Any        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS3_Q2_Rare010_Dmg_AR_Pathg_Hom        <- length (unique (subset (v_full_hq_r05.df, subset = FS3_Select == 1 & FS1_AR_Pathg_Hom        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS3_Q2_Rare010_Dmg_AR_Pathg_PotCompHet <- length (unique (subset (v_full_hq_r05.df, subset = FS3_Select == 1 & FS1_AR_Pathg_PotCompHet == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS3_Q2_Rare010_Dmg_XL_Pathg_Hap        <- length (unique (subset (v_full_hq_r05.df, subset = FS3_Select == 1 & FS1_XL_Pathg_Hap        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS3_Q2_Rare010_Dmg_XL_Pathg_Hom        <- length (unique (subset (v_full_hq_r05.df, subset = FS3_Select == 1 & FS1_XL_Pathg_Hom        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS3_Q2_Rare010_Dmg_CX_Pathg_HomHap     <- length (unique (subset (v_full_hq_r05.df, subset = FS3_Select == 1 & FS1_CX_Pathg_HomHap     == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS3_Q2_Rare010_Dmg_CX_Pathg_PotCompHet <- length (unique (subset (v_full_hq_r05.df, subset = FS3_Select == 1 & FS1_CX_Pathg_PotCompHet == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS3_Q2_Rare010_Dmg_CX_Uncertain        <- length (unique (subset (v_full_hq_r05.df, subset = FS3_Select == 1 & FS1_CX_Uncertain        == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS3_Q2_Rare010_Dmg_AR_Carrier          <- length (unique (subset (v_full_hq_r05.df, subset = FS3_Select == 1 & FS1_AR_Carrier          == 1, select = Original_VCFKEY, drop = T)))
stats.ls$VarN_FS3_Q2_Rare010_Dmg_XL_Carrier          <- length (unique (subset (v_full_hq_r05.df, subset = FS3_Select == 1 & FS1_XL_Carrier          == 1, select = Original_VCFKEY, drop = T)))

gc (); gc ()

# (9) REFORMAT STATS.LS

# 9.1. split structured and atomic stats for printing

## stats1.chv <- c ("check_typeseq_cc", "check_typeseq_nc", "check_effect_lof", "check_effect_mis", "check_effect_ot1", "check_effect_ot2",
##                           "chromosome_counts")

stats1.chv <- "chromosome_counts"

stats1.ls <- stats.ls[stats1.chv]
stats2.ls <- stats.ls[setdiff (names (stats.ls), 
                               c (stats1.chv, "chromosome_counts", "chr_zigosity_counts_AllQ_AllSeq_AllFreq", "chr_ploidy_counts_Allq_AllSeq_AllFreq", 
                                  "chr_zigosity_counts_Q1_AllSeq_AllFreq",   "chr_ploidy_counts_Q1_AllSeq_AllFreq", "chr_zigosity_counts_Q1_AllSeq_Rare050", "chr_ploidy_counts_Q1_AllSeq_Rare050"))]

# 9.2. reformat structured into character vector

# this works for one-level lists with any type of atomic or vectorial-1D object in the slots
f.listCollapse <- function (slot.obj, slot.name)
{
  if (! is.null (names (slot.obj)))
  {out.chv <- paste (slot.name, paste (paste (names (slot.obj), slot.obj, sep = ": "), collapse = "; "), sep = " = ")}
  else
  {out.chv <- paste (slot.name, paste (slot.obj, collapse = "; "), sep = " = ")}
  return (out.chv)
}

stat1.chv <- mapply (FUN = f.listCollapse, stats1.ls, names (stats1.ls))

# 9.3. reformat structured into data frame

stats2.df <- data.frame (Stat = names (stats2.ls), Value = unlist (stats2.ls), stringsAsFactors = F)

# (10) WRITE OUT

# 10.0. Setwd
if(!exists(output_path.ch))
  dir.create(output_path.ch, recursive = T)
setwd (output_path.ch)

# 10.1. Chromosome counts

write.table (cbind (data.frame (Chr_X_Zygosity = rownames (stats.ls$chr_zigosity_counts_AllQ_AllSeq_AllFreq)), as.data.frame (stats.ls$chr_zigosity_counts_AllQ_AllSeq_AllFreq)), 
             col.names = T, row.names = F, quote = F, sep = "\t", 
             file = paste (output_prefix.ch, "_Stats_chr_zigosity_counts_AllQ_AllSeq_AllFreq", ".txt", sep = ""))

write.table (cbind (data.frame (Chr_X_Zygosity = rownames (stats.ls$chr_zigosity_counts_Q1_AllSeq_AllFreq)), as.data.frame (stats.ls$chr_zigosity_counts_Q1_AllSeq_AllFreq)), 
             col.names = T, row.names = F, quote = F, sep = "\t", 
             file = paste (output_prefix.ch, "_Stats_chr_zigosity_counts_Q1_AllSeq_AllFreq", ".txt", sep = ""))

write.table (cbind (data.frame (Chr_X_Zygosity = rownames (stats.ls$chr_zigosity_counts_Q1_AllSeq_Rare050)), as.data.frame (stats.ls$chr_zigosity_counts_Q1_AllSeq_Rare050)), 
             col.names = T, row.names = F, quote = F, sep = "\t", 
             file = paste (output_prefix.ch, "_Stats_chr_zigosity_counts_Q1_AllSeq_Rare05", ".txt", sep = ""))

# 10.2. Variant Counts and Other Stats

cat (stat1.chv, sep = "\n", file = paste (output_prefix.ch, "_Stats_checks", ".txt", sep = ""))

write.table (stats2.df, col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.ch, "_Stats_var_counts", ".txt", sep = ""))

# 10.3. Variants

write.table (v_full_hq_r05.df, col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.ch, "_Full_Rare05", ".txt", sep = ""))
write.table (subset (v_full_hq_r05.df, F_DamageType != "NO"), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.ch, "_Full_Rare05_Dmg", ".txt", sep = ""))

write.table (subset (v_full_r05.df, FS1_Select == 1), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.ch, "_Full_SecondaryFindings_Rare05", ".txt", sep = ""))
write.table (subset (v_full_hq_r05.df, FS1_Select == 1), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.ch, "_SecondaryFindings_Rare05", ".txt", sep = ""))

write.table (subset (v_full_hq_r05.df, FM_HOM     == 1), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.ch, "_Main_RecessiveHom_Rare05",  ".txt", sep = ""))
write.table (subset (v_full_hq_r05.df, FM_XHAP    == 1), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.ch, "_Main_RecessiveXHap_Rare05", ".txt", sep = ""))
write.table (subset (v_full_hq_r05.df, FM_PCHET   >= 1), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.ch, "_Main_RecessiveCHet_Rare05", ".txt", sep = ""))
write.table (subset (v_full_hq_r05.df, FM_AXDOM   == 1), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.ch, "_Main_AXDominant_Rare005",    ".txt", sep = ""))
write.table (subset (v_full_hq_r05.df, FM_HZ == 1), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.ch, "_Main_HetHotzone_Rare0015",   ".txt", sep = ""))