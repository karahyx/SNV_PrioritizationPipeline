

# (0) SETTINGS ----------------------------------------------------------------

if (!require(data.table)) {
  install.packages("data.table")
  library(data.table)
}

if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}

if (!require(tibble)) {
  install.packages("tibble")
  library(tibble)
}

if (!require(stringr)) {
  install.packages("stringr")
  library(stringr)
}

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)

source("~/Desktop/PrioritizationPipeline/Rscript/pipeline/pipeline_fam_funcs.R")
# source(args[1])

# (1) VARIABLES & CUTOFFS -----------------------------------------------------

# 1.1. Input Variables
input_var_fam.file <- "~/Desktop/PrioritizationPipeline/data/dragen/AZ/AZ.hard-filtered.vcf.gz.annovar.out_SUBSET_rev27.7_hg38.tsv"
input_child.id <- "NA24385_b"

# input_var_fam.file <- args[2]
# input_child.id <- args[3]
child_column.name <- paste0(input_child.id, ":")

# 1.2. Output Variables
output_path <- args[4]
output_prefix <- paste(output_path, input_child.id, sep = "/")

# print input and output
input_var_fam.file
input_child.id
child_column.name
# input_pedigree.file
output_path
output_prefix

# 1.3. Internal Variables

typeseq_coding.chv <- c("exonic", "exonic;splicing", "splicing")
typeseq_ncrna.chv  <- c("ncRNA_exonic", "ncRNA_splicing","ncRNA_exonic;ncRNA_splicing")
typeseq_utr.chv    <- c("UTR3", "UTR5", "UTR3;UTR5", "UTR5;UTR3")

eff_lof.chv    <- c("frameshift deletion", "frameshift insertion", "frameshift substitution", "frameshift block substitution", "stopgain", "stoploss", "stopgain SNV", "stoploss SNV")
eff_missn.chv  <- c("nonsynonymous SNV")
eff_other_sub.chv  <- c("nonframeshift deletion", "nonframeshift insertion", "nonframeshift substitution", "nonframeshift block substitution")
eff_syn.chv <- "synonymous SNV"

# 1.4. Cutoffs

# 1.4.1. High-quality Filter

DP_cutoff <- 2

# 1.4.2. Define Damage

# Missense
sift_cutoff <- 0.05
polyphen_cutoff <- 0.9
ma_cutoff <- 1.9
phylopMam_missense_cutoff <- 1.3 
phylopVert_missense_cutoff <- 3.9
CADD_phred_missense_cutoff <- 21.1
REVEL_cutoff <- 0.75 
MPC_cutoff <- 2

missense_t1_cutoff <- 2
missense_t2_cutoff <- 4

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
phylopMam_utr_t1_cutoff <- 1.1
CADD_phred_utr_t1_cutoff <- 13.7
phylopMam_utr_t2_cutoff <- 1.3
CADD_phred_utr_t2_cutoff <- 21.1

# Non-coding (nc) - in order: phylopMam_avg, phylopVert100_avg, CADD_phred
nc_t1_cutoffs <- c(1.1, 1.6, 13.7)
nc_t2_cutoffs <- c(1.3, 3.9, 21.1)

# 1.4.3 Main Findings

# Heterozygous Hotzone
gnomAD_oe_lof_upper_cutoff <- 0.35

# (2) Main ----------------------------------------------------------------

# 2.1. Import pedigree file
# pedigree <- data.table::fread(input_pedigree.file)
pedigree <- data.frame(rep("AZ", 3),
                       c("NA24149_b", "NA24143_b", "NA24385_b"),
                       c("unknown", "unknown", "NA24149_b"),
                       c("unknown", "unknown", "NA24143_b"),
                       c("male", "female", "male"),
                       c("0", "0", "0"))
colnames(pedigree) <- c("FAMILY", "PERSON", "FATHER", "MOTHER", "SEX", "AFFECTION")

# 2.2. Find the parents
child <- pedigree %>% subset(PERSON == input_child.id)
parents <- pedigree %>% subset(PERSON == child$FATHER | PERSON == child$MOTHER)

father.id <- parents %>% subset(SEX == "male") %>% pull(PERSON)
mother.id <- parents %>% subset(SEX == "female") %>% pull(PERSON)
father_column.name <- paste0(father.id, ":")
mother_column.name <- paste0(mother.id, ":")

fam_dat <- data.table(name = c("child", "father", "mother"), 
                      id = c(input_child.id, father.id, mother.id),
                      sex = c(child$SEX, "male", "female"),
                      gt = character(3),
                      zygosity = character(3))

# 2.3. Import variant data
v_full.temp.df <- data.table::fread(input_var_fam.file, data.table = F)
check_samples_exist(v_full.temp.df, child_column.name, father_column.name, mother_column.name)
v_full.temp.df <- subset(v_full.temp.df, select = -c(DP, cg_freq_max)) # DP is removed because there would be two DP columns after removing "genomeName:"

# 2.4. Re-format column names
names(v_full.temp.df) <- gsub(child_column.name, "", names(v_full.temp.df)) # remove "genomeName:" from the child's columns
names(v_full.temp.df)[1] <- "CHROM"
names(v_full.temp.df)

# 2.5. Process the full data
v_full.df <- subset(v_full.temp.df, subset = (Zygosity != "hom-ref" & Zygosity != "unknown"))
v_full.df$alt_fraction <- as.numeric(v_full.df$AD_ALT) / ( as.numeric(v_full.df$AD_REF) + as.numeric(v_full.df$AD_ALT) )

# 2.6. Free up memory
rm(v_full.temp.df)
gc()

# 2.7. Add filters
v_full.df <- remove_alt_contigs(v_full.df)
v_full_r05.df <- get_rare05_variants(v_full.df)
v_full_hq.df <- get_hq_variants(v_full.df) # used for stats.ls

# 2.8. Get chromosome counts + chromosome-wise zygosity counts
chr_zygosity_stats_full <- get_chr_zygosity_stats(v_full.df[, c("CHROM", "Zygosity")])
chr_zygosity_stats_full_hq <- get_chr_zygosity_stats(v_full_hq.df[, c("CHROM", "Zygosity")])
chr_zygosity_stats_full_r05 <- get_chr_zygosity_stats(v_full_r05.df[, c("CHROM", "Zygosity")])

# 2.9. Get all variant stats
stats.df.all <- get_all_stats(v_full.df, v_full_hq.df, v_full_r05.df)

# 2.10. Change tier levels from 1 and 2 to "Low" and "High"
v_full_r05.df <- change_tier_levels(v_full_r05.df)

# 2.11. Output
write.table(v_full_r05.df, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Full_Rare05", ".txt", sep = ""))
write.table(subset(v_full_r05.df, F_DamageType != "NotDmg"), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix, "_Full_Rare05_Dmg", ".txt", sep = ""))
write.table(subset(v_full_r05.df, FS1_Select == 1 | FS2_Select == 1 | FS3_Select == 1), col.names = T, row.names = F, quote = F, sep = "\t", file = paste (output_prefix, "_Full_Rare05_SecondaryFindings", ".txt", sep = ""))

write.table(chr_zygosity_stats_full, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Stats_chr_zygosity_counts_AllQ_AllSeq_AllFreq", ".txt", sep = ""))
write.table(chr_zygosity_stats_full_hq, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Stats_chr_zygosity_counts_HQ_AllSeq_AllFreq", ".txt", sep = ""))
write.table(chr_zygosity_stats_full_r05, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Stats_chr_zygosity_counts_AllQ_AllSeq_Rare05", ".txt", sep = ""))

write.table(stats.df.all, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Stats_var_counts", ".txt", sep = ""))

sessionInfo()