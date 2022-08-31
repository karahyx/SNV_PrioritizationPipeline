# ---
# Title: Variant Prioritization
# Purpose: An adaption of the original variant prioritization pipeline to DRAGEN
# Author: Kara Han <kara.han@sickkids.ca>
# Date script created: 2022-07-20 14:32:21 EDT
# Version: 0.1.0
# Depends: 
#     R (>= 3.4.0)
# Imports:
#     data.table_1.14.2
#     dplyr_1.0.9
#     tibble_3.1.8
# ---

# SETTINGS ----------------------------------------------------------------

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

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)

# (0) VARIABLES & CUTOFFS -----------------------------------------------------

# 0.1. Input Variables

input_var_genome.file <- args[1]
input_var_genome.name <- args[2]
alt_input_var_genome.name <- paste(input_var_genome.name, ":", sep = "")

# 0.2. Output Variables

output_path <- args[3]
output_prefix <- paste(output_path, input_var_genome.name, sep = "/")

# print input and output
input_var_genome.file
input_var_genome.name
alt_input_var_genome.name
output_path
output_prefix

# 0.3. Internal Variables

typeseq_coding.chv <- c("exonic", "exonic;splicing", "splicing")
typeseq_ncrna.chv  <- c("ncRNA_exonic", "ncRNA_splicing","ncRNA_exonic;ncRNA_splicing")
typeseq_utr.chv    <- c("UTR3", "UTR5", "UTR3;UTR5", "UTR5;UTR3")

eff_lof.chv    <- c("frameshift deletion", "frameshift insertion", "frameshift substitution", "frameshift block substitution", "stopgain", "stoploss", "stopgain SNV", "stoploss SNV")
eff_missn.chv  <- c("nonsynonymous SNV")
eff_other_sub.chv  <- c("nonframeshift deletion", "nonframeshift insertion", "nonframeshift substitution", "nonframeshift block substitution")
eff_syn.chv <- "synonymous SNV"

# 0.4. Cutoffs

# 0.4.1. High-quality Filter
DP_cutoff <- 2

# 0.4.2. Define Damage

# Missense
sift_cutoff <- 0.05
polyphen_cutoff <- 0.90
ma_cutoff <- 1.90
phylopMam_missense_cutoff <- 1.3 
phylopVert_missense_cutoff <- 3.9
CADD_phred_missense_cutoff <- 21.1

missense_rank1_cutoff <- 2
missense_rank2_cutoff <- 4

# Other coding - in order: phylopMam_avg, phylopVert100_avg, CADD_phred
otherc_rk1_cr1_cutoffs <- c(1.2, 2.5, 13.5)
otherc_rk1_cr2_cutoffs <- c(1.5, 2.0, 13.0)
otherc_rk2_cr1_cutoffs <- c(2.0, 3.5, 14)
otherc_rk2_cr2_cutoffs <- c(1.5, 2.5, 13.5)

# Splicing predictions
spliceAI_DS_AG_r1_cutoff <- 0.5
spliceAI_DP_AG_r1_cutoff <- 50
spliceAI_DS_AL_r1_cutoff <- 0.5
spliceAI_DP_AL_r1_cutoff <- 50
spliceAI_DS_DG_r1_cutoff <- 0.5
spliceAI_DP_DG_r1_cutoff <- 50
spliceAI_DS_DL_r1_cutoff <- 0.5
spliceAI_DP_DL_r1_cutoff <- 50

spliceAI_DS_AG_r2_cutoff <- 0.5
spliceAI_DP_AG_r2_cutoff <- 50
spliceAI_DS_AL_r2_cutoff <- 0.5
spliceAI_DP_AL_r2_cutoff <- 50
spliceAI_DS_DG_r2_cutoff <- 0.5
spliceAI_DP_DG_r2_cutoff <- 50
spliceAI_DS_DL_r2_cutoff <- 0.5
spliceAI_DP_DL_r2_cutoff <- 50
dbscSNV_ADA_SCORE_cutoff <- 0.6
dbscSNV_RF_SCORE_cutoff <- 0.6

# UTR
phylopMam_utr_rk1_cutoff <- 1.1
CADD_phred_utr_rk1_cutoff <- 13.7
phylopMam_utr_rk2_cutoff <- 1.3
CADD_phred_utr_rk2_cutoff <- 21.1

# Non-coding (nc) - in order: phylopMam_avg, phylopVert100_avg, CADD_phred
nc_rk1_cutoffs <- c(1.1, 1.6, 13.7)
nc_rk2_cutoffs <- c(1.3, 3.9, 21.1)

# 0.4.3 Main Findings

# Heterozygous Hotzone
gnomAD_oe_lof_upper_cutoff <- 0.35


# (1) FUNCTIONS -----------------------------------------------------------

# 1.1. Frequency filter

add_freq_tag <- function(data, freq_max_cutoff) {
  
  freq_max_var <- with(data,
                       which((is.na (gnomAD_exome_freq_max)  | gnomAD_exome_freq_max  <= freq_max_cutoff) &
                               (is.na (gnomAD_genome_freq_max)  | gnomAD_genome_freq_max  <= freq_max_cutoff)))
  
  data$F_Rare[freq_max_var] <- freq_max_cutoff
  
  return(data)
}

# 1.2. Quality filter

add_pass_tag <- function(data) {
  
  passed_var <- with(data, which(data$FILTER == "PASS"))
  
  data$F_Pass <- 0
  data$F_Pass[passed_var] <- 1
  
  return(data)
}

#' Assumption: data is high-quality (FILTER == "PASS") variants 
add_qual_tag <- function(data, DP_cutoff) {
  
  ok_qual_var <- with(data, which(DP >= DP_cutoff))
  
  data$F_Qual_tag <- "LowQuality"
  data$F_Qual_tag[ok_qual_var] <- "OK"
  
  return(data)
}


# 1.3. Coding tag

add_coding_tag <- function(data) {
  
  data$F_Coding <- "Other"
  data$F_Coding[which(data$typeseq_priority %in% typeseq_ncrna.chv)] <- "ncRNA"
  data$F_Coding[which(data$typeseq_priority %in% typeseq_coding.chv)] <- "Coding"
  
  return(data)
}


# 1.4. Define damage

# Coding LOF
add_coding_lof_tag <- function(data, eff_lof.chv) {
  
  lof_var <- with(data, which(F_Coding == "Coding" & 
                                (effect_priority %in% eff_lof.chv | typeseq_priority %in% c ("splicing", "exonic;splicing"))))
  
  data$F_DamageType[lof_var] <- "LOF"
  data$F_DamageRank[lof_var] <- 2
  
  return(data)
}

add_coding_lof_spliceJunction_tag <- function(data, eff_lof.chv) {
  
  lof_spliceJunction_var <- with(data, which(F_Coding == "Coding" 
                                             & (effect_priority %in% eff_lof.chv 
                                                | (typeseq_priority %in% c ("splicing", "exonic;splicing") 
                                                   & distance_spliceJunction < 3))))
  
  data$F_S_DamageType[lof_spliceJunction_var] <- "LOF"
  
  return(data)
}

# Missense
add_missense_tag <- function(data, sift_cutoff, polyphen_cutoff, ma_cutoff, 
                             phylopMam_missense_cutoff, phylopVert_missense_cutoff, 
                             CADD_phred_missense_cutoff, missense_rank1_cutoff,
                             missense_rank2_cutoff) {
  
  if (missense_rank1_cutoff > missense_rank2_cutoff) {
    stop("missense_rank1_cutoff should be <= missense_rank2_cutoff, otherwise 
         F_DamageType would not be correctly labelled as 'Missense'")
  }
  
  missense_mx <- matrix(data = 0, nrow = nrow(data), ncol = 6)
  missense_mx[, 1] <- with(data, as.numeric(sift_score     <= sift_cutoff))
  missense_mx[, 2] <- with(data, as.numeric(polyphen_score >= polyphen_cutoff)) 
  missense_mx[, 3] <- with(data, as.numeric(ma_score       >= ma_cutoff)) 
  missense_mx[, 4] <- with(data, as.numeric(phylopMam_avg  >= phylopMam_missense_cutoff)) 
  missense_mx[, 5] <- with(data, as.numeric(phylopVert100_avg >= phylopVert_missense_cutoff))
  missense_mx[, 6] <- with(data, as.numeric(CADD_phred     >= CADD_phred_missense_cutoff))
  
  missense_mx[is.na(missense_mx)] <- 0
  missense_score <- apply(missense_mx, 1, sum)
  
  missense_rank1_var <- with(data, which(F_Coding == "Coding" & 
                                           (effect_priority %in% eff_missn.chv &  
                                              (missense_score >= missense_rank1_cutoff))))
  
  missense_rank2_var <- with(data, which(F_Coding == "Coding" & 
                                           (effect_priority %in% eff_missn.chv & 
                                              (missense_score >= missense_rank2_cutoff))))
  
  data$F_DamageType[missense_rank1_var] <- "Missense"
  data$F_DamageRank[missense_rank1_var] <- 1 
  data$F_DamageRank[missense_rank2_var] <- 2 
  
  return(data)
}


# Other coding
add_otherc_tag <- function(data, cutoffs1, cutoffs2, rank) {
  
  otherc_dmg_var <- with(data, which(F_Coding == "Coding" & (
    (effect_priority %in% eff_other_sub.chv & 
       ((phylopMam_avg >= cutoffs1[1] | phylopVert100_avg >= cutoffs1[2] | CADD_phred >= cutoffs1[3]) & is.na (dbsnp_common) )) | 
      (effect_priority %in% eff_other_sub.chv & 
         ((phylopMam_avg >= cutoffs2[1] | phylopVert100_avg >= cutoffs2[2] | CADD_phred >= cutoffs2[3]) & is.na (dbsnp) & is.na (dbsnp_region))) ) ))
  
  data$F_DamageType[otherc_dmg_var] <- "OtherC"
  data$F_DamageRank[otherc_dmg_var] <- rank
  
  return(data)
} 

# Splicing predictions
add_splicing_tag <- function(data, spliceAI_DS_AG_cutoff, spliceAI_DP_AG_cutoff, 
                             spliceAI_DS_AL_cutoff, spliceAI_DP_AL_cutoff, 
                             spliceAI_DS_DG_cutoff, spliceAI_DP_DG_cutoff,
                             spliceAI_DS_DL_cutoff, spliceAI_DP_DL_cutoff, 
                             dbsc_SNV_ADA_SCORE_cutoff = 0, 
                             dbscSNV_RF_SCORE_cutoff = 0, rank) {
  
  if (rank == 1) {
    splicing_var <- with(data, which (
      ((spliceAI_DS_AG > spliceAI_DS_AG_cutoff & abs (spliceAI_DP_AG) <= spliceAI_DP_AG_cutoff) |
         (spliceAI_DS_AL > spliceAI_DS_AL_cutoff & abs (spliceAI_DP_AL) <= spliceAI_DP_AL_cutoff) |
         (spliceAI_DS_DG > spliceAI_DS_DG_cutoff & abs (spliceAI_DP_DG) <= spliceAI_DP_DG_cutoff) |
         (spliceAI_DS_DL > spliceAI_DS_DL_cutoff & abs (spliceAI_DP_DL) <= spliceAI_DP_DL_cutoff) ) & 
        ! F_DamageType %in% c ("LOF", "Splc") & ! (F_DamageType %in% "Missense" & F_DamageRank == 2)))
  }
  else if (rank == 2) {
    splicing_var <- with(data, which (
      ((spliceAI_DS_AG > spliceAI_DS_AG_cutoff & abs (spliceAI_DP_AG) <= spliceAI_DP_AG_cutoff) | 
         (spliceAI_DS_AL > spliceAI_DS_AL_cutoff & abs (spliceAI_DP_AL) <= spliceAI_DP_AL_cutoff) |
         (spliceAI_DS_DG > spliceAI_DS_DG_cutoff & abs (spliceAI_DP_DG) <= spliceAI_DP_DG_cutoff) |
         (spliceAI_DS_DL > spliceAI_DS_DL_cutoff & abs (spliceAI_DP_DL) <= spliceAI_DP_DL_cutoff) |
         (dbscSNV_ADA_SCORE > dbscSNV_ADA_SCORE_cutoff | dbscSNV_RF_SCORE > dbscSNV_RF_SCORE_cutoff)) & 
        ! F_DamageType %in% c("LOF") ))
  }
  
  data$F_DamageType[splicing_var] <- "Splc"
  data$F_DamageRank[splicing_var] <- rank
  
  return(data)
}


# UTR
add_utr_dmg_tag <- function(data, phylopMam_cutoff, CADD_phred_cutoff, rank) {
  
  utr_dmg_var <- with(data, which(
    typeseq_priority %in% typeseq_utr.chv & 
      (phylopMam_avg >= phylopMam_cutoff | CADD_phred >= CADD_phred_cutoff) &
      (! is.na (phastCons_placental)) ))
  
  data$F_DamageType[utr_dmg_var] <- "UTR"
  data$F_DamageRank[utr_dmg_var] <- rank
  
  return(data)
}


# Non-coding
add_nc_dmg_tag <- function(data, cutoffs, rank) {
  
  nc_dmg_var <- with(data, which(F_Coding == "ncRNA" & ! gene_desc %in% "pseudogene" & F_DamageRank == 0 & (
    ((phylopMam_avg >= cutoffs[1] | phylopVert100_avg >= cutoffs[2]) & (! is.na (phastCons_placental))) | CADD_phred >= cutoffs[3])))
  
  data$F_DamageType[nc_dmg_var] <- "DmgNcRNA"
  data$F_DamageRank[nc_dmg_var] <- rank
  
  return(data)
}


# 1.4. Phenotype filter 

add_HPO_CGD_dominant_tag <- function(data, database, pattern) {
  
  data <- data.frame(data)
  col_name <- paste("G_AXD_", database, sep = "")
  data[, col_name] <- 0
  
  if (database == "CGD") {
    database <- "CGD_inheritance"
  }
  
  database_dom_var <- grep(data[, database], pattern = pattern)
  data[, col_name][database_dom_var] <- 1
  
  return(data)
}

add_phenotype_rank_tag <- function(data) {
  
  data$F_PhenoRank <- 0
  
  # Remove NAs from:
  # 1) subset of genes with MPO annotations
  pheno_mouse_anno <- setdiff(c(subset(data, subset = (! is.na (MPO) & MPO != ""), 
                                       select = entrez_id, drop = T)), NA)
  # 2) subset of genes with OMIM/HPO/CGD annotations 
  pheno_human_anno <- setdiff(c(subset(data, subset = (! is.na (omim_phenotype) & omim_phenotype != "") | 
                                         (! is.na (HPO) & HPO != "") | 
                                         (! is.na (CGD_disease) & CGD_disease != ""), 
                                       select = entrez_id, drop = T)), NA)
  
  pheno_mouse_var <- with(data, which(entrez_id %in% pheno_mouse_anno))
  pheno_human_var <- with(data, which(entrez_id %in% pheno_human_anno)) 
  
  data$F_PhenoRank[pheno_mouse_var] <- 1
  data$F_PhenoRank[pheno_human_var] <- 2
  
  return(data)
}


# 1.5. Main Findings

# RECESSIVE -- HOMOZYGOUS
add_recessive_homozygous_tag <- function(data) {
  
  rec_hom_var <- with(data, which(F_Rare <= 0.05 & F_DamageType != "NotDmg" & Zygosity == "hom-alt"))
  data$FM_HOM <- 0
  data$FM_HOM[rec_hom_var] <- 1
  
  return(data)
}

# X-LINKED HAPLOID
add_X_haploid_tag <- function(data) {
  
  X_haploid_var <- with(data, which(F_Rare <= 0.05 & F_DamageType != "NotDmg" & Zygosity == "hom-alt" & CHROM == "chrX"))
  data$FM_XHAP <- 0
  data$FM_XHAP[X_haploid_var] <- 1
  
  return(data)
}

# POTENTIAL COMPOUND HETS

get_potential_cmphet_df <- function(cmp_hets_df) { # help function for add_potential_compound_heterozygous_tag
  
  cmp_hets_df <- unique(cmp_hets_df)
  cmp_hets_df_final <- subset(cmp_hets_df, 
                              subset = gene_symbol %in% na.omit(cmp_hets_df$gene_symbol[duplicated(cmp_hets_df$gene_symbol)]), 
                              select = c(Original_VCFKEY, gene_symbol))
  
  return(cmp_hets_df_final)
}

#' Assumption: data is v_full_hq_r05.df for rank == 2 because of F_Qual_tag
#' This function is used for both main and secondary findings
add_potential_compound_heterozygous_tag <- function(data, secondary_findings = FALSE) {
  
  if (secondary_findings == TRUE) {
    col_name <- "F_CmpHet_S1"
    # used F_Pass == 1 instead of F_Qual > 1 in the original script for variants with PASS FILTER
    cmp_hets_r1_df <- subset(data, subset = F_Pass == 1 & FS1_Select == 1, select = c(gene_symbol, Original_VCFKEY))
    cmp_hets_r2_df <- subset(data, subset = F_Qual_tag != "LowQuality" & FS1_Select == 1, select = c(gene_symbol, Original_VCFKEY))
  }
  else {
    col_name <- "FM_PCHET"
    cmp_hets_r1_df <- subset(data, 
                             subset = F_Rare <= 0.05 & F_DamageType != "NotDmg", 
                             select = c(gene_symbol, Original_VCFKEY))
    cmp_hets_r2_df <- subset(data, 
                             subset = F_Rare <= 0.05 & F_DamageType != "NotDmg" & F_Qual_tag != "LowQuality", 
                             select = c(gene_symbol, Original_VCFKEY))
  }
  
  cmp_hets_r1_df_final <- get_potential_cmphet_df(cmp_hets_r1_df)
  cmp_hets_r2_df_final <- get_potential_cmphet_df(cmp_hets_r2_df)
  
  data[, col_name] <- 0
  data[, col_name][which(data$Original_VCFKEY %in% cmp_hets_r1_df_final$Original_VCFKEY & data$gene_symbol %in% cmp_hets_r1_df_final$gene_symbol)] <- 1
  data[, col_name][which(data$Original_VCFKEY %in% cmp_hets_r2_df_final$Original_VCFKEY & data$gene_symbol %in% cmp_hets_r2_df_final$gene_symbol)] <- 2
  
  return(data)
}

# POTENTIAL COMPOUND HETS (DAMAGING)
add_potential_dmg_compound_heterozygous_tag <- function(data) {
  
  dmg_cmp_hets_df <- subset(data, subset = FM_PCHET == 2 & F_DamageRank == 2, select = gene_symbol, drop = T)
  dmg_cmp_hets_df_final <- dmg_cmp_hets_df[duplicated(dmg_cmp_hets_df)]
  
  data$FM_PCHET_DMG <- 0 
  data$FM_PCHET_DMG[with(data, which(gene_symbol %in% dmg_cmp_hets_df_final & FM_PCHET == 2))] <- 1 
  
  return(data)
}

# DOMINANT
add_dom_tag <- function(data) {
  
  dom_var <- with(data, which(F_Rare <= 0.005 & F_DamageType != "NotDmg" & (G_AXD_CGD == 1 | G_AXD_HPO == 1)))
  
  data$FM_AXDOM <- 0
  data$FM_AXDOM[dom_var] <- 1
  
  return(data)
}

# HETEROZYGOUS HOTZONE
add_het_hotzone_tag <- function(data, gnomAD_oe_lof_upper_cutoff) {
  
  het_hotzone_var <- with(data, which(Zygosity == "ref-alt" & F_Rare <= 0.0015 & (
    (gnomAD_oe_lof_upper < gnomAD_oe_lof_upper_cutoff & 
       (F_DamageType %in% c ("LOF", "Splc", "OtherC") | (F_DamageType == "Missense" & F_DamageRank == 2) )) )))
  
  data$FM_HZ <- 0
  data$FM_HZ[het_hotzone_var] <- 1
  
  return(data)
}


# 1.6. Secondary Findings

# Pathogenicity flag
add_pathogenic_tag <- function(data) {
  
  clinvar_not_pathg_var <- with(data, which(Clinvar_SIG_Simple == 0))
  clinvar_pathg_var <- with(data, which(Clinvar_SIG_Simple == 1))
  
  data$F_Clinvar_notPathg <- 0
  data$F_Clinvar_Pathg <- 0
  data$F_Clinvar_notPathg[clinvar_not_pathg_var] <- 1
  data$F_Clinvar_Pathg[clinvar_pathg_var] <- 1
  
  return(data)
}

# Rank 1
add_secondary_findings_rank1_tag <- function(data) {
  
  sf1_var <- with(data, which(CGD_disease != "" & (F_S_DamageType == "LOF") | (F_Clinvar_Pathg == 1)))
  
  data$FS1_Select <- 0
  data$FS1_Select[sf1_var] <- 1
  
  return(data)
}

# Rank 2
add_secondary_findings_rank2_tag <- function(data) {
  
  # genes with 1) CGD annotations; 2) hq; 3) rare 1%; 4) not non-pathogenic
  sf2_var <- with(data, which(CGD_disease != "" & F_Qual_tag != "LowQuality"  
                              & F_Rare <= 0.01 & ((F_DamageType == "LOF") | 
                                                    (F_Clinvar_Pathg == 1) | 
                                                    (F_Clinvar_notPathg == 0)) ))
  data$FS2_Select <- 0
  data$FS2_Select[sf2_var] <- 1
  
  return(data)
}

# Rank 3
add_secondary_findings_rank3_tag <- function(data) {
  
  # genes with 1) CGD annotations; 2) hq; 3) rare 1%; 4) not non-pathogenic; 5) damaging
  sf2_var <- with(data, which(CGD_disease != "" & F_Qual_tag != "LowQuality" & F_Rare <= 0.01 
                              & F_DamageType != "NotDmg" & 
                                ((F_DamageType == "LOF") | 
                                   (F_Clinvar_Pathg == 1) | 
                                   (F_Clinvar_notPathg == 0)) )) 
  data$FS3_Select <- 0
  data$FS3_Select[sf2_var] <- 1
  
  return(data)
}

# dominant: pathogenic
add_dominant_pathogenic_tag <- function(data) {
  
  dom_pathg_var <- with(data, which(FS1_Select == 1 & grepl("AD", CGD_inheritance)))
  
  data$FS1_AD_Pathg_Any <- 0
  data$FS1_AD_Pathg_Any[dom_pathg_var] <- 1
  
  return(data)
}

# recessive + hom: pathogenic
add_recessive_hom_pathg_tag <- function(data) {
  
  rec_hom_pathg_var <- with(data, which(FS1_Select == 1 & grepl("AR", CGD_inheritance) & Zygosity == "hom-alt"))
  
  data$FS1_AR_Pathg_Hom <- 0
  data$FS1_AR_Pathg_Hom[rec_hom_pathg_var] <- 1
  
  return(data)
}

# recessive + potential compound het: pathogenic
add_recessive_cmphet_pathg_tag <- function(data) {
  
  rec_cmphet_pathg_var <- with(data, which(FS1_Select == 1 & grepl("AR", CGD_inheritance) & F_CmpHet_S1 >= 1))
  
  data$FS1_AR_Pathg_PotCompHet <- 0
  data$FS1_AR_Pathg_PotCompHet[rec_cmphet_pathg_var] <- 1
  
  return(data)
}

# X-linked + hom/hap: pathogenic
add_X_hom_or_hap_pathg_tag <- function(data, type) {
  
  X_hom_or_hap_pathg_var <- with(data, which(FS1_Select == 1 & grepl("XL", CGD_inheritance) & Zygosity == "hom-alt" & CHROM == "chrX"))
  
  if (type == "hom") {
    col_name <- "FS1_XL_Pathg_Hom"
  }
  else if (type == "hap") {
    col_name <- "FS1_XL_Pathg_Hap"
  }
  
  data[, col_name] <- 0
  data[, col_name][X_hom_or_hap_pathg_var] <- 1
  
  return(data)
}

# complex + hom / hap: pathogenic
add_complex_hom_hap_pathg_tag <- function(data) {
  
  complex_hom_hap_pathg_tag <- with(data, which(FS1_Select == 1 & (! grepl("AD|AR|XL", CGD_inheritance)) & (Zygosity == "hom-alt")))
  
  data$FS1_CX_Pathg_HomHap <- 0
  data$FS1_CX_Pathg_HomHap[complex_hom_hap_pathg_tag] <- 1
  
  return(data)
}

# complex + potential compound het: pathogenic
add_complex_cmphet_pathg_tag <- function(data) {
  
  complex_cmphet_pathg_var <- with(data, which(FS1_Select == 1 & (! grepl("AD|AR|XL", CGD_inheritance)) & F_CmpHet_S1 >= 1))
  
  data$FS1_CX_Pathg_PotCompHet <- 0
  data$FS1_CX_Pathg_PotCompHet[complex_cmphet_pathg_var] <- 1
  
  return(data)
}

# complex + single het: uncertain
add_complex_single_het_uncertain_tag <- function(data) {
  
  complex_single_het_uncertain_var <- with(data, 
                                           which(FS1_Select == 1 & (! grepl("AD|AR|XL", CGD_inheritance)) & 
                                                   Zygosity %in% c ("ref-alt", "alt-alt") & F_CmpHet_S1 == 0))
  
  data$FS1_CX_Uncertain <- 0
  data$FS1_CX_Uncertain[complex_single_het_uncertain_var] <- 1
  
  return(data)
}

# recessive + single het: carrier
add_recessive_single_het_carrier_tag <- function(data) {
  
  rec_single_het_var <- with(data, which(FS1_Select == 1 & grepl("AR", CGD_inheritance) & 
                                           Zygosity %in% c ("ref-alt", "alt-alt") & F_CmpHet_S1 == 0))
  
  data$FS1_AR_Carrier <- 0
  data$FS1_AR_Carrier[rec_single_het_var] <- 1
  
  return(data)
}

# X-linked + het: carrier
add_X_het_carrier_tag <- function(data){
  
  X_het_carrier_var <- with(data, which(FS1_Select == 1 & grepl("XL", CGD_inheritance) & 
                                          Zygosity %in% c ("ref-alt", "alt-alt") & CHROM == "chrX"))
  
  data$FS1_XL_Carrier <- 0
  data$FS1_XL_Carrier[X_het_carrier_var] <- 1
  
  return(data)
}

# ACMG Disease
add_ACMG_tag <- function(data) {
  
  ACMG_var <- with(data, which(!is.na(ACMG_disease)))
  
  data$F_ACMG <- 0
  data$F_ACMG[ACMG_var] <- 1
  
  return(data)
}

add_ACMG_coding_tag <- function(data, typeseq_coding.chv) {
  
  ACMG_coding_var <- with(data, which(!is.na(ACMG_disease) & typeseq_priority %in% typeseq_coding.chv))
  
  data$F_ACMG_Coding <- 0
  data$F_ACMG_Coding[ACMG_coding_var] <- 1
  
  return(data)
}

# 1.7. Get final results

# 1.7.1. Rare Variants

# input: data = v_full.df
get_rare05_variants <- function(data) {
  
  data <- add_freq_tag(data, freq_max_cutoff = 0.05)
  data <- add_freq_tag(data, freq_max_cutoff = 0.01)
  data <- add_freq_tag(data, freq_max_cutoff = 0.005)
  data <- add_freq_tag(data, freq_max_cutoff = 0.0015)
  data <- add_freq_tag(data, freq_max_cutoff = 0)
  
  data <- add_pass_tag(data)
  data <- add_coding_tag(data)
  
  data$F_DamageType <- "NotDmg" 
  data$F_S_DamageType <- "NotDmg"
  data$F_DamageRank <- 0
  
  data <- add_coding_lof_tag(data, eff_lof.chv)
  data <- add_coding_lof_spliceJunction_tag(data, eff_lof.chv)
  
  data <- add_pathogenic_tag(data)
  data <- add_secondary_findings_rank1_tag(data)
  
  return(data)
}

# 1.7.2. HQ Variants (FILTER == "PASS")
get_hq_variants <- function(data) {
  
  data <- subset(data, subset = (FILTER == "PASS"))
  data <- add_qual_tag(data, DP_cutoff) 
  
  return(data)
}

# 1.7.3. HQ Rare Variants

# input: data = v_full_r05.df
get_hq_rare05_variants <- function(data) {
  
  data <- subset(data, subset = (FILTER == "PASS"))
  data <- add_qual_tag(data, DP_cutoff)
  data <- add_coding_tag(data)
  
  data$F_DamageType <- "NotDmg"
  data$F_S_DamageType <- "NotDmg"
  data$F_DamageRank <- 0
  
  data <- add_coding_lof_tag(data, eff_lof.chv)
  data <- add_coding_lof_spliceJunction_tag(data, eff_lof.chv)
  
  data <- add_missense_tag(data, sift_cutoff, polyphen_cutoff, ma_cutoff, phylopMam_missense_cutoff, phylopVert_missense_cutoff, 
                           CADD_phred_missense_cutoff, missense_rank1_cutoff, missense_rank2_cutoff)
  
  data <- add_otherc_tag(data, otherc_rk1_cr1_cutoffs, otherc_rk1_cr2_cutoffs, rank = 1)
  data <- add_otherc_tag(data, otherc_rk2_cr1_cutoffs, otherc_rk2_cr2_cutoffs, rank = 2)
  
  data <- add_splicing_tag(data, spliceAI_DS_AG_r1_cutoff, spliceAI_DP_AG_r1_cutoff, 
                           spliceAI_DS_AL_r1_cutoff, spliceAI_DP_AL_r1_cutoff,
                           spliceAI_DS_DG_r1_cutoff, spliceAI_DP_DG_r1_cutoff,
                           spliceAI_DS_DL_r1_cutoff, spliceAI_DP_DL_r1_cutoff, rank = 1)
  
  data <- add_splicing_tag(data, spliceAI_DS_AG_r2_cutoff, spliceAI_DP_AG_r2_cutoff,
                           spliceAI_DS_AL_r2_cutoff, spliceAI_DP_AL_r2_cutoff,
                           spliceAI_DS_DG_r2_cutoff, spliceAI_DP_DG_r2_cutoff,
                           spliceAI_DS_DL_r2_cutoff, spliceAI_DP_DL_r2_cutoff, 
                           dbscSNV_ADA_SCORE_cutoff, dbscSNV_RF_SCORE_cutoff, rank = 2)
  
  data <- add_utr_dmg_tag(data, phylopMam_utr_rk1_cutoff,
                          CADD_phred_utr_rk1_cutoff, rank = 1)
  data <- add_utr_dmg_tag(data, phylopMam_utr_rk2_cutoff,
                          CADD_phred_utr_rk2_cutoff, rank = 2)
  
  data <- add_nc_dmg_tag(data, nc_rk1_cutoffs, rank = 1)
  data <- add_nc_dmg_tag(data, nc_rk2_cutoffs, rank = 2)
  
  data <- add_HPO_CGD_dominant_tag(data, database = "HPO", pattern = "@AD")
  data <- add_HPO_CGD_dominant_tag(data, database = "CGD", pattern = "AD")
  data <- add_phenotype_rank_tag(data)
  
  data <- add_recessive_homozygous_tag(data)
  data <- add_X_haploid_tag(data)
  data <- add_potential_compound_heterozygous_tag(data, secondary_findings = F)
  data <- add_potential_dmg_compound_heterozygous_tag(data)
  data <- add_dom_tag(data)
  data <- add_het_hotzone_tag(data, gnomAD_oe_lof_upper_cutoff)
  
  data <- add_pathogenic_tag(data)
  data <- add_secondary_findings_rank1_tag(data)
  data <- add_secondary_findings_rank2_tag(data)
  data <- add_secondary_findings_rank3_tag(data)
  
  data <- add_dominant_pathogenic_tag(data)
  data <- add_recessive_hom_pathg_tag(data)
  data <- add_potential_compound_heterozygous_tag(data, secondary_findings = T)
  data <- add_recessive_cmphet_pathg_tag(data)
  data <- add_X_hom_or_hap_pathg_tag(data, type = "hom")
  data <- add_X_hom_or_hap_pathg_tag(data, type = "hap")
  data <- add_complex_hom_hap_pathg_tag(data)
  data <- add_complex_cmphet_pathg_tag(data)
  data <- add_complex_single_het_uncertain_tag(data)
  data <- add_recessive_single_het_carrier_tag(data)
  data <- add_X_het_carrier_tag(data)
  
  data <- add_ACMG_tag(data)
  data <- add_ACMG_coding_tag(data, typeseq_coding.chv)
  
  return(data)
}

# 1.8. Stats

get_chr_zygosity_stats <- function(data) {
  
  chrom_list <- paste("chr", c(1:22, "X", "Y", "M"), sep = "")
  data <- as.data.table(data)
  data <- data[CHROM %in% chrom_list]
  
  chr_counts.mx <- table(data$CHROM) %>% as.matrix()
  zygosity_counts.mx <- table(data[, c("CHROM", "Zygosity")]) %>% as.data.frame.matrix()
  
  chr_zygosity_counts.mx <- cbind(chr_counts.mx, zygosity_counts.mx)
  chr_zygosity_counts.mx <- chr_zygosity_counts.mx[order(factor(rownames(chr_zygosity_counts.mx), levels = chrom_list)), ]
  
  chr_zygosity_counts.mx$HomPerc <- (zygosity_counts.mx[, "hom-alt"]) / 
    (zygosity_counts.mx[, "ref-alt"] + zygosity_counts.mx[, "hom-alt"] + zygosity_counts.mx[, "alt-alt"]) * 100
  
  colnames(chr_zygosity_counts.mx)[1] <- "Freq"
  chr_zygosity_counts.mx <- tibble::rownames_to_column(chr_zygosity_counts.mx, "Chrom")
  
  return(chr_zygosity_counts.mx)
}

add_stats_unique_vcfkey <- function(data, high_quality = FALSE) {
  
  if (high_quality == TRUE) {
    num_var <- length(unique(subset(data, subset = F_Qual_tag != "LowQuality", 
                                    select = Original_VCFKEY, drop = T))) # Didn't use a wrapper here for better readability
  } else {
    num_var <- length(unique(data$Original_VCFKEY))
  }
  
  return(num_var)
}

get_num_var <- function(data, cond) {
  
  envir_use <- parent.frame()
  num_var <- data %>% subset(eval(cond, envir = data, enclos = envir_use)) %>% select(Original_VCFKEY) %>% unique() %>% nrow()
  
  return(num_var)
}

add_stats_freq <- function(data, freq_cutoff, high_quality = FALSE) {
  
  if (high_quality == TRUE) {
    cond <- expression(F_Rare <= freq_cutoff & F_Qual_tag != "LowQuality")
  } else {
    cond <- expression(F_Rare <= freq_cutoff)
  }
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_typeseq_tag <- function(data, typeseq.chv, freq_cutoff = -1, high_quality = FALSE) {
  
  if (high_quality == TRUE) {
    if (freq_cutoff != -1) {
      cond <- expression(typeseq_priority %in% typeseq.chv & F_Qual_tag != "LowQuality" & F_Rare <= freq_cutoff)
    } else {
      cond <- expression(typeseq_priority %in% typeseq.chv & F_Qual_tag != "LowQuality")
    }
  } else {
    if (freq_cutoff != -1) {
      cond <- expression(typeseq_priority %in% typeseq.chv & F_Rare <= freq_cutoff)
    } else {
      cond <- expression(typeseq_priority %in% typeseq.chv)
    }
  }
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_zygosity_tag <- function(data, zygosity, chrX = FALSE, high_quality = FALSE) {
  
  if (high_quality == TRUE) {
    if (chrX == TRUE) {
      cond <- expression(Zygosity == zygosity & F_Qual_tag != "LowQuality" & CHROM == "chrX")
    } else {
      cond <- expression(Zygosity == zygosity & F_Qual_tag != "LowQuality")
    }
  } else {
    if (chrX == TRUE) {
      cond <- expression(Zygosity == zygosity & CHROM == "chrX")
    } else {
      cond <- expression(Zygosity == zygosity)
    }
  }
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_dmg_tag <- function(data, freq_cutoff, dmg_type, dmg_rank = 0, high_quality) {
  
  if (dmg_rank > 0) {
    if (high_quality == TRUE) {
      cond <- expression(F_Rare <= freq_cutoff & F_DamageType == dmg_type & F_DamageRank >= dmg_rank & F_Qual_tag != "LowQuality")
    } else {
      # Note that F_Pass is not necessarily needed if this function is only ever used for v_full_hq.df and v_full_hq_r05.df
      # F_Pass is added to ensure that the function works for v_full.df and v_full_r05.df (which contain FAIL filters)
      cond <- expression(F_Rare <= freq_cutoff & F_DamageType == dmg_type & F_DamageRank >= dmg_rank & F_Pass == 1)
    }
  } else {
    if (high_quality == TRUE) {
      cond <- expression(F_Rare <= freq_cutoff & F_DamageType == dmg_type & F_Qual_tag != "LowQuality")
    } else {
      cond <- expression(F_Rare <= freq_cutoff & F_DamageType == dmg_type & F_Pass == 1)
    }
  }
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_pheno_filter <- function(data, freq_cutoff, database, high_quality = FALSE) {
  
  switch(database,
         HPO = db_cond <- expression(G_AXD_HPO == 1),
         CGD = db_cond <- expression(G_AXD_CGD == 1),
         all = db_cond <- expression(G_AXD_HPO == 1 | G_AXD_CGD == 1))
  
  if (high_quality == TRUE) {
    cond <- expression(F_Rare <= freq_cutoff & F_Coding == "Coding" & F_DamageType != "NotDmg" & eval(db_cond) & F_Qual_tag != "LowQuality")
  } else {
    cond <- expression(F_Rare <= freq_cutoff & F_Coding == "Coding" & F_DamageType != "NotDmg" & eval(db_cond) & F_Pass == 1)
  }
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_pheno_rank <- function(data, freq_cutoff, pheno_rank, high_quality = FALSE) {
  
  if (high_quality == TRUE) {
    cond <- expression(F_Rare <= freq_cutoff & F_DamageType != "NotDmg" & F_PhenoRank >= pheno_rank & F_Qual_tag != "LowQuality")
  } else {
    cond <- expression(F_Rare <= freq_cutoff & F_DamageType != "NotDmg" & F_PhenoRank >= pheno_rank & F_Pass == 1)
  }
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_main_findings <- function(data, type, pheno_rank, dmg_rank = 0, high_quality = FALSE) {
  
  switch(type,
         hom = type_cond <- expression(FM_HOM == 1),
         Xhom = type_cond <- expression(FM_XHAP == 1),
         cmp_het = ifelse(high_quality == TRUE, type_cond <- expression(FM_PCHET >= 2), type_cond <- expression(FM_PCHET >= 1)),
         dominant = type_cond <- expression(FM_AXDOM == 1),
         het_hotzone = type_cond <- expression(FM_HZ == 1))
  
  if (type == "cmp_het") {
    cond <- expression(eval(type_cond) & F_PhenoRank >= pheno_rank & F_DamageRank >= dmg_rank)
  } else {
    if (high_quality == TRUE) {
      cond <- expression(F_PhenoRank >= pheno_rank & F_DamageRank >= dmg_rank & eval(type_cond) & F_Qual_tag != "LowQuality")
    } else {
      cond <- expression(F_PhenoRank >= pheno_rank & F_DamageRank >= dmg_rank & eval(type_cond) & F_Pass == 1)
    }
  }
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_secondary_findings <- function(data, rank, type = NA) {
  
  switch(rank,
         "1" = rank_cond <- expression(FS1_Select == 1),
         "2" = rank_cond <- expression(FS2_Select == 1),
         "3" = rank_cond <- expression(FS3_Select == 1))
  
  if (is.na(type)) {
    cond <- rank_cond
    num_var <- get_num_var(data, cond)
    
    return(num_var)
  }
  
  switch(type,
         dominant = type_cond <- expression(FS1_AD_Pathg_Any == 1),
         recessive_hom = type_cond <- expression(FS1_AR_Pathg_Hom == 1),
         recessive_cmp_het = type_cond <- expression(FS1_AR_Pathg_PotCompHet == 1),
         Xhap = type_cond <- expression(FS1_XL_Pathg_Hap == 1),
         Xhom = type_cond <- expression(FS1_XL_Pathg_Hom == 1),
         complex_hom_hap = type_cond <- expression(FS1_CX_Pathg_HomHap == 1),
         complex_cmp_het = type_cond <- expression(FS1_CX_Pathg_PotCompHet == 1),
         complex_single_het = type_cond <- expression(FS1_CX_Uncertain == 1),
         recessive_single_het = type_cond <- expression(FS1_AR_Carrier == 1),
         Xhet = type_cond <- expression(FS1_XL_Carrier == 1))
  
  cond <- expression(eval(rank_cond) & eval(type_cond))
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

#' input: data = v_full.df
add_stats_full_df <- function(data) {
  
  stats.ls <- list()
  
  stats.ls$VarN_AllQ_AllSeq_AllFreq <- add_stats_unique_vcfkey(data, high_quality = F)
  stats.ls$VarN_AllQ_Coding_AllFreq <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, high_quality = F)
  stats.ls$VarN_AllQ_ncRNA_AllFreq  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, high_quality = F)
  
  return(stats.ls)
}

#' input: data = v_full_hq.df
add_stats_full_hq_df <- function(data) {
  
  stats.ls <- list()
  
  stats.ls$VarN_Q1_AllSeq_AllFreq <- add_stats_unique_vcfkey(data, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_AllFreq <- add_stats_unique_vcfkey(data, high_quality = T)
  
  stats.ls$VarN_Q1_Coding_AllFreq <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, high_quality = F)
  stats.ls$VarN_Q2_Coding_AllFreq <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, high_quality = T)
  
  stats.ls$VarN_Q1_ncRNA_AllFreq <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, high_quality = F)
  stats.ls$VarN_Q2_ncRNA_AllFreq <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, high_quality = T)
  
  stats.ls$VarN_Q1_AllSeq_AllFreq_Xhom <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = T, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_AllFreq_Xhom <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = T, high_quality = T)
  
  stats.ls$VarN_Q1_AllSeq_AllFreq_Xhet <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = T, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_AllFreq_Xhet <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = T, high_quality = T)
  
  stats.ls$VarN_Q1_AllSeq_AllFreq_Hom  <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = F, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_AllFreq_Hom  <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = F, high_quality = T)
  
  stats.ls$VarN_Q1_AllSeq_AllFreq_HetR <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = F, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_AllFreq_HetR <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = F, high_quality = T)
  stats.ls$VarN_Q1_AllSeq_AllFreq_HetA <- add_stats_zygosity_tag(data, zygosity = "alt-alt", chrX = F, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_AllFreq_HetA <- add_stats_zygosity_tag(data, zygosity = "alt-alt", chrX = F, high_quality = T)
  
  return(stats.ls)
}

#' input: data = v_full_hq_r05.df
add_stats_full_hq_r05_df <- function(data) {
  
  stats.ls <- list()
  
  stats.ls$VarN_Q1_AllSeq_Rare050 <- add_stats_freq(data, freq_cutoff = 0.05, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_Rare050 <- add_stats_freq(data, freq_cutoff = 0.05, high_quality = T)
  stats.ls$VarN_Q1_AllSeq_Rare010 <- add_stats_freq(data, freq_cutoff = 0.01, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_Rare010 <- add_stats_freq(data, freq_cutoff = 0.01, high_quality = T)
  stats.ls$VarN_Q1_AllSeq_Rare005 <- add_stats_freq(data, freq_cutoff = 0.005, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_Rare005 <- add_stats_freq(data, freq_cutoff = 0.005, high_quality = T)
  stats.ls$VarN_Q1_AllSeq_Rare0015 <- add_stats_freq(data, freq_cutoff = 0.0015, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_Rare0015 <- add_stats_freq(data, freq_cutoff = 0.0015, high_quality = T)
  stats.ls$VarN_Q1_AllSeq_Rare000 <- add_stats_freq(data, freq_cutoff = 0, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_Rare000 <- add_stats_freq(data, freq_cutoff = 0, high_quality = T)
  
  stats.ls$VarN_Q1_Coding_Rare050  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.05, high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare050  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.05, high_quality = T)
  stats.ls$VarN_Q1_Coding_Rare010  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.01, high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare010  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.01, high_quality = T)
  stats.ls$VarN_Q1_Coding_Rare005  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.005, high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare005  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.005, high_quality = T)
  stats.ls$VarN_Q1_Coding_Rare0015 <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.0015, high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare0015 <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.0015, high_quality = T)
  stats.ls$VarN_Q1_Coding_Rare000  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0, high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare000  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0, high_quality = T)
  
  stats.ls$VarN_Q1_ncRNA_Rare050  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.05, high_quality = F)
  stats.ls$VarN_Q2_ncRNA_Rare050  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.05, high_quality = T)
  stats.ls$VarN_Q1_ncRNA_Rare010  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.01, high_quality = F)
  stats.ls$VarN_Q2_ncRNA_Rare010  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.01, high_quality = T)
  stats.ls$VarN_Q1_ncRNA_Rare005  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.005, high_quality = F)
  stats.ls$VarN_Q2_ncRNA_Rare005  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.005, high_quality = T)
  stats.ls$VarN_Q1_ncRNA_Rare0015 <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.0015, high_quality = F)
  stats.ls$VarN_Q2_ncRNA_Rare0015 <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.0015, high_quality = T)
  stats.ls$VarN_Q1_ncRNA_Rare000  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0, high_quality = F)
  stats.ls$VarN_Q2_ncRNA_Rare000  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0, high_quality = T)
  
  stats.ls$VarN_Q1_AllSeq_Rare050_Xhom <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = T, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_Rare050_Xhom <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = T, high_quality = T)
  stats.ls$VarN_Q1_AllSeq_Rare050_Hom  <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = F, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_Rare050_Hom  <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = F, high_quality = T)
  stats.ls$VarN_Q1_AllSeq_Rare050_HetR <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = F, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_Rare050_HetR <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = F, high_quality = T)
  stats.ls$VarN_Q1_AllSeq_Rare050_HetA <- add_stats_zygosity_tag(data, zygosity = "alt-alt", chrX = F, high_quality = F)
  stats.ls$VarN_Q2_AllSeq_Rare050_HetA <- add_stats_zygosity_tag(data, zygosity = "alt-alt", chrX = F, high_quality = T)
  
  stats.ls$GeneN_Pheno_rare050 <- length (unique (data[data$F_PhenoRank == 2, ]$entrez_id))
  stats.ls$GeneN_Dominant_rare050 <- length (setdiff (c (data$entrez_id[grep ("@AD", data$HPO)], data$entrez_id[grep ("AD", data$CGD_inheritance)]), NA))
  
  stats.ls$VarN_Q1_Coding_Rare010_LOF <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "LOF", high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare010_LOF <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "LOF", high_quality = T)
  
  stats.ls$VarN_Q1_Coding_Rare010_MissDmgR1 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Missense", dmg_rank = 1, high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare010_MissDmgR1 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Missense", dmg_rank = 1, high_quality = T)
  
  stats.ls$VarN_Q1_Coding_Rare010_MissDmgR2 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Missense", dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare010_MissDmgR2 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Missense", dmg_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Coding_Rare010_SplcR1 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Splc", dmg_rank = 1, high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare010_SplcR1 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Splc", dmg_rank = 1, high_quality = T)
  
  stats.ls$VarN_Q1_Coding_Rare010_SplcR2 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Splc", dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare010_SplcR2 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Splc", dmg_rank = 2, high_quality = T)
  
  stats.ls$Coding_HQ1_Rare010_OtherDmg <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "OtherC", high_quality = F)
  stats.ls$Coding_HQ2_Rare010_OtherDmg <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "OtherC", high_quality = T)
  
  stats.ls$VarN_Q1_ncRNA_Rare010_DmgR1 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "DmgNcRNA", dmg_rank = 1, high_quality = F)
  stats.ls$VarN_Q2_ncRNA_Rare010_DmgR1 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "DmgNcRNA", dmg_rank = 1, high_quality = T)
  
  stats.ls$VarN_Q1_ncRNA_Rare010_DmgR2 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "DmgNcRNA", dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_ncRNA_Rare010_DmgR2 <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "DmgNcRNA", dmg_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Coding_Rare010_Dmg_AXD_HPO <- add_stats_pheno_filter(data, freq_cutoff = 0.01, database = "HPO", high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare010_Dmg_AXD_HPO <- add_stats_pheno_filter(data, freq_cutoff = 0.01, database = "HPO", high_quality = T)
  
  stats.ls$VarN_Q1_Coding_Rare010_DmG_AXD_CGD <- add_stats_pheno_filter(data, freq_cutoff = 0.01, database = "CGD", high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare010_DmG_AXD_CGD <- add_stats_pheno_filter(data, freq_cutoff = 0.01, database = "CGD", high_quality = T)
  
  stats.ls$VarN_Q1_Coding_Rare010_Dmg_AXD_All <- add_stats_pheno_filter(data, freq_cutoff = 0.01, database = "all", high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare010_Dmg_AXD_All <- add_stats_pheno_filter(data, freq_cutoff = 0.01, database = "all", high_quality = T)
  
  stats.ls$VarN_Q1_Coding_Rare010_Dmg_PhenoRank1 <- add_stats_pheno_rank(data, freq_cutoff = 0.01, pheno_rank = 1, high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare010_Dmg_PhenoRank1 <- add_stats_pheno_rank(data, freq_cutoff = 0.01, pheno_rank = 1, high_quality = T)
  
  stats.ls$VarN_Q1_Coding_Rare010_Dmg_PhenoRank2 <- add_stats_pheno_rank(data, freq_cutoff = 0.01, pheno_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Coding_Rare010_Dmg_PhenoRank2 <- add_stats_pheno_rank(data, freq_cutoff = 0.01, pheno_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Rare050_DmgR2_Hom_PhenoRank0 <- add_stats_main_findings(data, type = "hom", pheno_rank = 0, dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Rare050_DmgR2_Hom_PhenoRank0 <- add_stats_main_findings(data, type = "hom", pheno_rank = 0, dmg_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Rare050_DmgR2_Hom_PhenoRank1 <- add_stats_main_findings(data, type = "hom", pheno_rank = 1, dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Rare050_DmgR2_Hom_PhenoRank1 <- add_stats_main_findings(data, type = "hom", pheno_rank = 1, dmg_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Rare050_DmgR2_Hom_PhenoRank2 <- add_stats_main_findings(data, type = "hom", pheno_rank = 2, dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Rare050_DmgR2_Hom_PhenoRank2 <- add_stats_main_findings(data, type = "hom", pheno_rank = 2, dmg_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Rare050_DmgR2_Xhom_PhenoRank0 <- add_stats_main_findings(data, type = "Xhom", pheno_rank = 0, dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Rare050_DmgR2_Xhom_PhenoRank0 <- add_stats_main_findings(data, type = "Xhom", pheno_rank = 0, dmg_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Rare050_DmgR2_Xhom_PhenoRank1 <- add_stats_main_findings(data, type = "Xhom", pheno_rank = 1, dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Rare050_DmgR2_Xhom_PhenoRank1 <- add_stats_main_findings(data, type = "Xhom", pheno_rank = 1, dmg_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Rare050_DmgR2_Xhom_PhenoRank2 <- add_stats_main_findings(data, type = "Xhom", pheno_rank = 2, dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Rare050_DmgR2_Xhom_PhenoRank2 <- add_stats_main_findings(data, type = "Xhom", pheno_rank = 2, dmg_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Rare050_DmgR1_CmpHet_PhenoRank0 <- add_stats_main_findings(data, type = "cmp_het", pheno_rank = 0, dmg_rank = 1, high_quality = F)
  stats.ls$VarN_Q2_Rare050_DmgR1_CmpHet_PhenoRank0 <- add_stats_main_findings(data, type = "cmp_het", pheno_rank = 0, dmg_rank = 1, high_quality = T)
  
  stats.ls$VarN_Q1_Rare050_DmgR1_CmpHet_PhenoRank1 <- add_stats_main_findings(data, type = "cmp_het", pheno_rank = 1, dmg_rank = 1, high_quality = F)
  stats.ls$VarN_Q2_Rare050_DmgR1_CmpHet_PhenoRank1 <- add_stats_main_findings(data, type = "cmp_het", pheno_rank = 1, dmg_rank = 1, high_quality = T)
  
  stats.ls$VarN_Q1_Rare050_DmgR1_CmpHet_PhenoRank2 <- add_stats_main_findings(data, type = "cmp_het", pheno_rank = 2, dmg_rank = 1, high_quality = F)
  stats.ls$VarN_Q2_Rare050_DmgR1_CmpHet_PhenoRank2 <- add_stats_main_findings(data, type = "cmp_het", pheno_rank = 2, dmg_rank = 1, high_quality = T)
  
  stats.ls$VarN_Q1_Rare005_DmgR2_AXDom_PhenoRank0 <- add_stats_main_findings(data, type = "dominant", pheno_rank = 0, dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Rare005_DmgR2_AXDom_PhenoRank0 <- add_stats_main_findings(data, type = "dominant", pheno_rank = 0, dmg_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Rare005_DmgR2_AXDom_PhenoRank1 <- add_stats_main_findings(data, type = "dominant", pheno_rank = 1, dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Rare005_DmgR2_AXDom_PhenoRank1 <- add_stats_main_findings(data, type = "dominant", pheno_rank = 1, dmg_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Rare005_DmgR2_AXDom_PhenoRank2 <- add_stats_main_findings(data, type = "dominant", pheno_rank = 2, dmg_rank = 2, high_quality = F)
  stats.ls$VarN_Q2_Rare005_DmgR2_AXDom_PhenoRank2 <- add_stats_main_findings(data, type = "dominant", pheno_rank = 2, dmg_rank = 2, high_quality = T)
  
  stats.ls$VarN_Q1_Rare005_DmgR2_HI_PhenoRank0 <- add_stats_main_findings(data, type = "het_hotzone", pheno_rank = 0, high_quality = F)
  stats.ls$VarN_Q2_Rare005_DmgR2_HI_PhenoRank0 <- add_stats_main_findings(data, type = "het_hotzone", pheno_rank = 0, high_quality = T)
  
  stats.ls$VarN_Q1_Rare005_DmgR2_HI_PhenoRank1 <- add_stats_main_findings(data, type = "het_hotzone", pheno_rank = 1, high_quality = F)
  stats.ls$VarN_Q2_Rare005_DmgR2_HI_PhenoRank1 <- add_stats_main_findings(data, type = "het_hotzone", pheno_rank = 1, high_quality = T)
  
  stats.ls$VarN_FS1_Q1_Rare050_Tot     <- add_stats_secondary_findings(data, rank = 1)
  stats.ls$VarN_FS2_Q2_Rare010_Tot     <- add_stats_secondary_findings(data, rank = 2)
  stats.ls$VarN_FS3_Q2_Rare010_Dmg_Tot <- add_stats_secondary_findings(data, rank = 3)
  
  stats.ls$VarN_FS1_Q1_Rare050_AD_Pathg_Any        <- add_stats_secondary_findings(data, rank = 1, type = "dominant")
  stats.ls$VarN_FS1_Q1_Rare050_AR_Pathg_Hom        <- add_stats_secondary_findings(data, rank = 1, type = "recessive_hom")
  stats.ls$VarN_FS1_Q1_Rare050_AR_Pathg_PotCompHet <- add_stats_secondary_findings(data, rank = 1, type = "recessive_cmp_het")
  stats.ls$VarN_FS1_Q1_Rare050_XL_Pathg_Hap        <- add_stats_secondary_findings(data, rank = 1, type = "Xhap")
  stats.ls$VarN_FS1_Q1_Rare050_XL_Pathg_Hom        <- add_stats_secondary_findings(data, rank = 1, type = "Xhom")
  stats.ls$VarN_FS1_Q1_Rare050_CX_Pathg_HomHap     <- add_stats_secondary_findings(data, rank = 1, type = "complex_hom_hap")
  stats.ls$VarN_FS1_Q1_Rare050_CX_Pathg_PotCompHet <- add_stats_secondary_findings(data, rank = 1, type = "complex_cmp_het")
  stats.ls$VarN_FS1_Q1_Rare050_CX_Uncertain        <- add_stats_secondary_findings(data, rank = 1, type = "complex_single_het")
  stats.ls$VarN_FS1_Q1_Rare050_AR_Carrier          <- add_stats_secondary_findings(data, rank = 1, type = "recessive_single_het")
  stats.ls$VarN_FS1_Q1_Rare050_XL_Carrier          <- add_stats_secondary_findings(data, rank = 1, type = "Xhet")
  
  stats.ls$VarN_FS2_Q2_Rare010_AD_Pathg_Any        <- add_stats_secondary_findings(data, rank = 2, type = "dominant")
  stats.ls$VarN_FS2_Q2_Rare010_AR_Pathg_Hom        <- add_stats_secondary_findings(data, rank = 2, type = "recessive_hom")
  stats.ls$VarN_FS2_Q2_Rare010_AR_Pathg_PotCompHet <- add_stats_secondary_findings(data, rank = 2, type = "recessive_cmp_het")
  stats.ls$VarN_FS2_Q2_Rare010_XL_Pathg_Hap        <- add_stats_secondary_findings(data, rank = 2, type = "Xhap")
  stats.ls$VarN_FS2_Q2_Rare010_XL_Pathg_Hom        <- add_stats_secondary_findings(data, rank = 2, type = "Xhom")
  stats.ls$VarN_FS2_Q2_Rare010_CX_Pathg_HomHap     <- add_stats_secondary_findings(data, rank = 2, type = "complex_hom_hap")
  stats.ls$VarN_FS2_Q2_Rare010_CX_Pathg_PotCompHet <- add_stats_secondary_findings(data, rank = 2, type = "complex_cmp_het")
  stats.ls$VarN_FS2_Q2_Rare010_CX_Uncertain        <- add_stats_secondary_findings(data, rank = 2, type = "complex_single_het")
  stats.ls$VarN_FS2_Q2_Rare010_AR_Carrier          <- add_stats_secondary_findings(data, rank = 2, type = "recessive_single_het")
  stats.ls$VarN_FS2_Q2_Rare010_XL_Carrier          <- add_stats_secondary_findings(data, rank = 2, type = "Xhet")
  
  stats.ls$VarN_FS3_Q2_Rare010_Dmg_AD_Pathg_Any        <- add_stats_secondary_findings(data, rank = 3, type = "dominant")
  stats.ls$VarN_FS3_Q2_Rare010_Dmg_AR_Pathg_Hom        <- add_stats_secondary_findings(data, rank = 3, type = "recessive_hom")
  stats.ls$VarN_FS3_Q2_Rare010_Dmg_AR_Pathg_PotCompHet <- add_stats_secondary_findings(data, rank = 3, type = "recessive_cmp_het")
  stats.ls$VarN_FS3_Q2_Rare010_Dmg_XL_Pathg_Hap        <- add_stats_secondary_findings(data, rank = 3, type = "Xhap")
  stats.ls$VarN_FS3_Q2_Rare010_Dmg_XL_Pathg_Hom        <- add_stats_secondary_findings(data, rank = 3, type = "Xhom")
  stats.ls$VarN_FS3_Q2_Rare010_Dmg_CX_Pathg_HomHap     <- add_stats_secondary_findings(data, rank = 3, type = "complex_hom_hap")
  stats.ls$VarN_FS3_Q2_Rare010_Dmg_CX_Pathg_PotCompHet <- add_stats_secondary_findings(data, rank = 3, type = "complex_cmp_het")
  stats.ls$VarN_FS3_Q2_Rare010_Dmg_CX_Uncertain        <- add_stats_secondary_findings(data, rank = 3, type = "complex_single_het")
  stats.ls$VarN_FS3_Q2_Rare010_Dmg_AR_Carrier          <- add_stats_secondary_findings(data, rank = 3, type = "recessive_single_het")
  stats.ls$VarN_FS3_Q2_Rare010_Dmg_XL_Carrier          <- add_stats_secondary_findings(data, rank = 3, type = "Xhet")
  
  return(stats.ls)
}

convert_stats_list_to_df <- function(stats.ls, df_name = NA) {
  
  vals <- sapply(seq_along(stats.ls), function(i) {stats.ls[i][[1]]})
  
  if (!is.na(df_name)) {
    stats.df <- cbind(names(stats.ls), vals, rep(df_name, length(stats.ls))) %>% as.data.frame() %>% setNames(c("Name", "Count", "Source"))
    return(stats.df)
  }
  stats.df <- cbind(names(stats.ls), vals) %>% as.data.frame() %>% setNames(c("Name", "Count"))
  
  return(stats.df)
}

get_all_stats <- function() {
  
  stats.df.all <- convert_stats_list_to_df(stats.ls.full, "v_full.df")
  stats.df.all <- rbind(stats.df.all, convert_stats_list_to_df(stats.ls.hq, "v_full_hq.df"))
  stats.df.all <- rbind(stats.df.all, convert_stats_list_to_df(stats.ls.hq.r05, "v_full_hq_r05.df"))
  
  return(stats.df.all)
}

# (2) Main ----------------------------------------------------------------

# File import
v_full.temp.df <- data.table::fread(input_var_genome.file, data.table = F)
v_full.temp.df <- subset(v_full.temp.df, select = -c(DP, cg_freq_max)) # DP is removed because there would be two DP columns after removing "genomeName:"

# Re-format column names
names(v_full.temp.df) <- gsub(alt_input_var_genome.name, "", names(v_full.temp.df)) # remove "genomeName:" from columns
names(v_full.temp.df)[1] <- "CHROM"
names(v_full.temp.df)

# Process the full data
v_full.df <- subset(v_full.temp.df, subset = (Zygosity != "hom-ref" & Zygosity != "unknown"))
v_full.df <- add_pass_tag(v_full.df)
v_full.df$alt_fraction <- v_full.df$AD_ALT/(v_full.df$AD_REF + v_full.df$AD_ALT)

# Free up memory
rm(v_full.temp.df)
gc()

# Annotate data
v_full_r05.df <- get_rare05_variants(v_full.df)
v_full_hq_r05.df <- get_hq_rare05_variants(v_full_r05.df)
v_full_hq.df <- get_hq_variants(v_full.df) # used for stats.ls

# Get chromosome counts + chromosome-wise zygosity counts
chr_zygosity_stats_full <- get_chr_zygosity_stats(v_full.df[, c("CHROM", "Zygosity")])
chr_zygosity_stats_full_hq <- get_chr_zygosity_stats(v_full_hq.df[, c("CHROM", "Zygosity")])
chr_zygosity_stats_full_hq_r05 <- get_chr_zygosity_stats(v_full_hq_r05.df[, c("CHROM", "Zygosity")])

# Get stats list for each data frame
stats.ls.full <- add_stats_full_df(v_full.df)
stats.ls.hq <- add_stats_full_hq_df(v_full_hq.df)
stats.ls.hq.r05 <- add_stats_full_hq_r05_df(v_full_hq_r05.df)

# Convert stats lists to readable data frames
stats.df.full <- convert_stats_list_to_df(stats.ls.full)
stats.df.hq <- convert_stats_list_to_df(stats.ls.hq)
stats.df.hq.r05 <- convert_stats_list_to_df(stats.ls.hq.r05)

# Get all stats
stats.df.all <- get_all_stats()

# Output
write.table(v_full_r05.df, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Full_Rare05", ".txt", sep = ""))
write.table(v_full_hq_r05.df, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Full_HQ_Rare05", ".txt", sep = ""))

write.table(chr_zygosity_stats_full, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Stats_chr_zygosity_counts_AllQ_AllSeq_AllFreq", ".txt", sep = ""))
write.table(chr_zygosity_stats_full_hq, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Stats_chr_zygosity_counts_HQ_AllSeq_AllFreq", ".txt", sep = ""))
write.table(chr_zygosity_stats_full_hq_r05, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Stats_chr_zygosity_counts_HQ_AllSeq_Rare05", ".txt", sep = ""))

write.table(stats.df.full, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Stats_Full", ".txt", sep = ""))
write.table(stats.df.hq, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Stats_Full_HQ", ".txt", sep = ""))
write.table(stats.df.hq.r05, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Stats_Full_HQ_Rare05", ".txt", sep = ""))
write.table(stats.df.all, col.names = T, row.names = F, quote = F, sep = "\t", file = paste(output_prefix, "_Stats_All", ".txt", sep = ""))

sessionInfo()