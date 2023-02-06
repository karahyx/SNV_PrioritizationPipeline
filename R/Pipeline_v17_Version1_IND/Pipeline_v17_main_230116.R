# ---
# Title: Variant Prioritization
# Purpose: An adaption of the original small variant prioritization pipeline
#          applicable to Broad Institute's GATK to the Illumina DRAGEN Bio-IT platform
# Author: Kara Han <kara.han@sickkids.ca>
# Adapted from: Prioritizaiton pipeline developed by Dr. Daniele M., Bhooma T., 
#               Thomas N., and Dr. Worrawat E. at TCAG
# Date script created: 2022-07-20 14:32:21 EDT
# Date last modified: 2023-01-16 13:35:54 EST
# Version: v17
# TCAG annotation pipeline version: rev27.7 hg38
# Depends: 
#     R (>= 3.4.0)
# Imports:
#     data.table_1.14.4
#     dplyr_1.0.10
#     rmarkdown_2.17
#     tibble_3.1.8
# ---

# (0) SETTINGS ----------------------------------------------------------------

# 0.1. Libraries

if (!require(data.table)) {
  install.packages("data.table")
  library(data.table)
}

if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}

if (!require(rmarkdown)) {
  install.packages("rmarkdown")
  library(rmarkdown)
}

if (!require(tibble)) {
  install.packages("tibble")
  library(tibble)
}

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)

# 0.2. Input variables

# Arguments
# args[1]: root path
# args[2]: sub-path and filename of second R code file defining functions
# args[3]: variant input file path
# args[4]: genome name used at columns
# args[5]: sub-path and filename of stats visualization Rmd file
# args[6]: output sub-path

# args <- character(6)
# args[1] <- "/Users/karahan/Desktop/PrioritizationPipeline/"
# args[2] <- "script/R/pipeline/v230116/Pipeline_v17_funcs_230116.R"
# 
# # individual:
# args[3] <- "/Users/karahan/Desktop/PrioritizationPipeline/data/testing_data/dragen/NA12878.hard-filtered.vcf.gz.annovar.out_SUBSET_rev27.7_hg38.tsv"
# args[4] <- "NA12878"
# # multisample:
# # args[3] <- "/Users/karahan/Desktop/PrioritizationPipeline/data/family/AZ/AZ.hard-filtered.vcf.gz.annovar.out_SUBSET_rev27.7_hg38.tsv"
# # args[4] <- "NA24385_b"
# #
# args[5] <- "script/R/pipeline/stats/Pipeline_stats.Rmd"
# args[6] <- "tests/NA12878"

root.path <- args[1]
fun.path_filename     <- paste (root.path, args[2], sep = "")
input_var.path        <- args[3]
input_var_genome.name <- args[4]
input_stats_vis.path_filename <- paste (root.path, args[5], sep = "")
output_prefix.path_filename <- paste (root.path, args[6], "/", input_var_genome.name, sep = "")

alt_input_var_genome.name <- paste0 (input_var_genome.name, ":")

cat ("\n")
cat ("R functions path and filename ");         cat (fun.path_filename);             cat ("\n")
cat ("Variant data path and filename ");        cat (input_var.path);                cat ("\n")
cat ("Genome name ");                           cat (input_var_genome.name);         cat ("\n")
cat ("Stats visualization path and filename "); cat (input_stats_vis.path_filename); cat ("\n")
cat ("Output path and filename prefix ");       cat (output_prefix.path_filename);   cat ("\n")
cat ("\n")

# 0.3. Functions

source(fun.path_filename)

# (1) VARIABLES & CUTOFFS -----------------------------------------------------

# 1.1. Internal Variables

typeseq_coding.chv <- c ("exonic", "exonic;splicing", "splicing")
typeseq_ncrna.chv  <- c ("ncRNA_exonic", "ncRNA_splicing","ncRNA_exonic;ncRNA_splicing")
typeseq_utr.chv    <- c ("UTR3", "UTR5", "UTR3;UTR5", "UTR5;UTR3")

eff_lof.chv    <- c ("frameshift deletion", "frameshift insertion", "frameshift substitution", "frameshift block substitution", "stopgain", "stoploss", "stopgain SNV", "stoploss SNV")
eff_missn.chv  <- c ("nonsynonymous SNV")
eff_other_sub.chv  <- c ("nonframeshift deletion", "nonframeshift insertion", "nonframeshift substitution", "nonframeshift block substitution")
eff_syn.chv <- "synonymous SNV"

# 1.2. Cutoffs

# 1.2.1. High-quality Filter

DP_cutoff <- 2

# 1.2.2. Define Damage

# Missense
REVEL_t1_cutoff <- 0.25
REVEL_t2_cutoff <- 0.45 
phylopMam_missense_cutoff <- 1.3
phylopVert_missense_cutoff <- 3.9
CADD_phred_missense_cutoff <- 30
MPC_cutoff <- 2

# Other coding - in order: phylopMam_avg, phylopVert100_avg, CADD_phred
otherc_t1_cutoffs <- c(1.1, 1.6, 13.7)
otherc_t2_cutoffs <- c(1.3, 3.9, 21.1)

# Splicing predictions
spliceAI_DS_AG_t1_cutoff <- 0.2
spliceAI_DP_AG_t1_cutoff <- 50
spliceAI_DS_AL_t1_cutoff <- 0.2
spliceAI_DP_AL_t1_cutoff <- 50
spliceAI_DS_DG_t1_cutoff <- 0.2
spliceAI_DP_DG_t1_cutoff <- 50
spliceAI_DS_DL_t1_cutoff <- 0.2
spliceAI_DP_DL_t1_cutoff <- 50

spliceAI_DS_AG_t2_cutoff <- 0.5
spliceAI_DP_AG_t2_cutoff <- 50
spliceAI_DS_AL_t2_cutoff <- 0.5
spliceAI_DP_AL_t2_cutoff <- 50
spliceAI_DS_DG_t2_cutoff <- 0.5
spliceAI_DP_DG_t2_cutoff <- 50
spliceAI_DS_DL_t2_cutoff <- 0.5
spliceAI_DP_DL_t2_cutoff <- 50

dbscSNV_ADA_SCORE_cutoff <- 0.6
dbscSNV_RF_SCORE_cutoff <- 0.6

# UTR
phylopMam_utr_t1_cutoff <- 1.0
CADD_phred_utr_t1_cutoff <- 12.5
phylopMam_utr_t2_cutoff <- 1.5
CADD_phred_utr_t2_cutoff <- 15

# Non-coding (nc) - in order: phylopMam_avg, phylopVert100_avg, CADD_phred
nc_t1_cutoffs <- c(1.1, 1.6, 13.7)
nc_t2_cutoffs <- c(1.3, 3.9, 21.1)

# 1.2.3 Main Findings

# Heterozygous Hotzone
gnomAD_oe_lof_cutoff <- 0.15
gnomAD_oe_mis_cutoff <- 0.7

# (2) Main ----------------------------------------------------------------

# 2.1. File import
v_full.temp.df <- data.table::fread(input_var.path, data.table = F, na.strings = c ("NA", "."))
v_full.temp.df <- subset(v_full.temp.df, select = -c(DP, cg_freq_max)) # DP is removed to avoid conflict after removing "genome_name:" for processing

# 2.1.1. Check if multi-sample
multisample <- FALSE
sample_cols <- names(v_full.temp.df)[grepl(":", names(v_full.temp.df))]
sample_names <- vapply(strsplit(sample_cols, ":"), `[`, 1, FUN.VALUE=character(1))
if (length(unique(sample_names)) > 1) multisample <- TRUE

# 2.2. Re-format column names
input_var_cols <- grep(input_var_genome.name, names(v_full.temp.df))
names(v_full.temp.df) <- gsub(alt_input_var_genome.name, "", names(v_full.temp.df)) # remove "genome_name:" from columns
names(v_full.temp.df)[1] <- "CHROM"
names(v_full.temp.df)

# 2.3. Process the full data
v_full.df <- subset(v_full.temp.df, subset = (Zygosity != "hom-ref" & Zygosity != "unknown"))
v_full.df$alt_fraction <- as.numeric(v_full.df$AD_ALT) / ( as.numeric(v_full.df$AD_REF) + as.numeric(v_full.df$AD_ALT) )

# 2.3.1. Add some fields for allele frequency
v_full_AF.df <- subset (v_full.df, select = c (gnomAD_exome211_AF_afr, gnomAD_exome211_AF_amr, gnomAD_exome211_AF_eas, gnomAD_exome211_AF_nfe, gnomAD_exome211_AF_oth, gnomAD_exome211_AF_sas,
                                               gnomAD_genome31_AF_afr, gnomAD_genome31_AF_amr, gnomAD_genome31_AF_eas, gnomAD_genome31_AF_nfe, gnomAD_genome31_AF_oth, gnomAD_genome31_AF_sas) )
v_full.df$FreqMaxSimple_AfrAmrEasNfeSasOth <- apply (v_full_AF.df, 1, max, na.rm = T)
rm (v_full_AF.df)
v_full_HomCount.df <- subset (v_full.df, select = c (gnomAD_exome211_nhomalt, gnomAD_genome31_nhomalt) )
v_full.df$FreqHomCount_AfrAmrEasNfeSasOth <- apply (v_full_HomCount.df, 1, max, na.rm = T)
rm (v_full_HomCount.df)
v_full.df$dbsnp_region_count <- sapply (lapply (strsplit (v_full.df$dbsnp_region, split = ","), setdiff, NA), length)

# 2.3.2. Add some fields for gene constraint
v_full.df$F_GeneConstr <- with (v_full.df, as.numeric (gnomAD_oe_lof <= gnomAD_oe_lof_cutoff | gnomAD_oe_mis <= gnomAD_oe_mis_cutoff))

# 2.4. Free up memory
rm(v_full.temp.df)
gc()

# 2.5. Add filters
v_full.df <- remove_alt_contigs (v_full.df)
v_full_r05.df <- get_rare05_variants (v_full.df, multisample)
v_full_hq.df <- get_hq_variants (v_full.df, multisample) # used for stats.ls

# 2.6. Get chromosome counts + chromosome-wise zygosity counts
chr_zygosity_stats_full <- get_chr_zygosity_stats(v_full.df[, c("CHROM", "Zygosity")])
chr_zygosity_stats_full_hq <- get_chr_zygosity_stats(v_full_hq.df[, c("CHROM", "Zygosity")])
chr_zygosity_stats_full_r05 <- get_chr_zygosity_stats(v_full_r05.df[, c("CHROM", "Zygosity")])

# 2.7. Get all variant stats
stats.df.all <- get_all_stats(v_full.df, v_full_hq.df, v_full_r05.df)

# 2.8. Process the results before outputting
v_full_r05.df <- change_tier_levels(v_full_r05.df) # change tier levels from 0, 1 and 2 to "-", "Low" and "High"
names(v_full_r05.df)[input_var_cols] <- paste0(alt_input_var_genome.name, names(v_full_r05.df)[input_var_cols]) # add genome_name back to columns

# 2.9. Output

rmarkdown::render(input = input_stats_vis.path_filename,
                  output_format = "html_document",
                  output_file = paste(output_prefix.path_filename, "Stats_visualization", sep = "_"))

v_sub_hom.df  <- subset (v_full_r05.df, FM_HOM   == 1 & F_Qual >= 1)
v_sub_chet.df <- subset (v_full_r05.df, FM_PCHET == 2)
v_sub_xhap.df <- subset (v_full_r05.df, FM_XHAP  == 1 & F_Qual >= 1)
v_sub_axd.df  <- subset (v_full_r05.df, FM_AXDOM == 1 & F_Qual >= 1)
v_sub_pdom.df <- subset (v_full_r05.df, FM_PDDOM == 1 & F_Qual >= 1 & FM_AXDOM == 0)

write.table (v_full_r05.df, col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Full_Rare05", ".txt", sep = ""))
write.table (subset(v_full_r05.df, F_DamageType != "NotDmg"), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Full_Rare05_Dmg", ".txt", sep = ""))
write.table (subset(v_full_r05.df, FS1_Select == 1 | FS2_Select == 1 | FS3_Select == 1), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Full_Rare05_SecondaryFindings", ".txt", sep = ""))

write.table (v_sub_hom.df,  col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Full_Rare05_Set1_HOM",     ".txt", sep = ""))
write.table (v_sub_xhap.df, col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Full_Rare05_Set2_XHAP",    ".txt", sep = ""))
write.table (v_sub_chet.df, col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Full_Rare05_Set3_CHET",    ".txt", sep = ""))
write.table (v_sub_axd.df,  col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Full_Rare05_Set4_AXDOM",   ".txt", sep = ""))
write.table (v_sub_pdom.df, col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Full_Rare05_Set5_PredDom", ".txt", sep = ""))

write.table (chr_zygosity_stats_full,     col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Stats_chr_zygosity_counts_AllQ_AllSeq_AllFreq", ".txt", sep = ""))
write.table (chr_zygosity_stats_full_hq,  col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Stats_chr_zygosity_counts_HQ_AllSeq_AllFreq",   ".txt", sep = ""))
write.table (chr_zygosity_stats_full_r05, col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Stats_chr_zygosity_counts_AllQ_AllSeq_Rare05",  ".txt", sep = ""))

write.table (stats.df.all, col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix.path_filename, "_Stats_var_counts", ".txt", sep = ""))

writeLines (capture.output (sessionInfo ()), con = paste (output_prefix.path_filename, "_RSessionInfo", ".txt", sep = ""))