# ---
# Title: Variant Prioritization Functions
# Purpose: Contains functions used in the small variant prioritization pipeline v17
# Author: Kara Han <kara.han@sickkids.ca>
# Adapted from: Prioritizaiton pipeline developed by Dr. Daniele M., Bhooma T., 
#               Thomas N., and Dr. Worrawat E. at TCAG
# Date script created: 2022-07-20 14:32:21 EDT
# Date last modified: 2023-01-16 14:04:50 EST
# Version: v17
# TCAG annotation pipeline version: rev27.7 hg38
# ---

# (3) FUNCTIONS -----------------------------------------------------------

# 3.0. Alternate contig and unlocalized/unplaced sequence filter

remove_alt_contigs <- function(data) {
  
  filtered_var <- data[!grepl("alt|random|chrUn", data$CHROM), ]
  
  return(filtered_var)
}

# 3.1. Frequency filter

add_freq_tag <- function(data, freq_max_cutoff) {
  
  freq_max_var <- with(data,
                       which(freq_max <= freq_max_cutoff))
  
  data$F_Rare[freq_max_var] <- freq_max_cutoff
  
  return(data)
}

# 3.2. Quality filter

# 3.2.1. Pass tag
add_pass_tag <- function(data) {
  
  passed_var <- with(data, which(data$FILTER == "PASS"))
  
  data$F_Pass <- 0
  data$F_Pass[passed_var] <- 1
  
  return(data)
}

# 3.2.2. Quality tag
add_qual_tag <- function(data, DP_cutoff) {
  
  high_qual_var <- with(data, which(FILTER == "PASS" & DP >= DP_cutoff))
  
  data$F_Qual <- 0
  data$F_Qual[high_qual_var] <- 1
  
  return(data)
}

add_qual_tag_multisample <- function(data, DP_cutoff) {

  qual_var_1 <- with (data, which (FT == "PASS" & DP >= DP_cutoff))
  qual_var_2 <- with (data, which (FT == "PASS" & DP >= DP_cutoff & is.na (SegDup) & (
    (Zygosity == "hap"     &  alt_fraction >= 0.8) |
      (Zygosity == "hom-alt" &  alt_fraction >= 0.8) |
      (Zygosity == "ref-alt" & (alt_fraction > 0.3 & alt_fraction < 0.7)) ) ))

  data$F_Qual <- 0
  data$F_Qual[qual_var_1] <- 1
  data$F_Qual[qual_var_2] <- 2

  return(data)
}

# 3.3. Coding tag

add_coding_tag <- function(data) {
  
  data$F_Coding <- "Other"
  data$F_Coding[which(data$typeseq_priority %in% typeseq_ncrna.chv)] <- "ncRNA"
  data$F_Coding[which(data$typeseq_priority %in% typeseq_coding.chv)] <- "Coding"
  
  return(data)
}

# 3.4. Define damage

# 3.4.1. Coding LOF
add_coding_lof_tag <- function(data, eff_lof.chv) {
  
  lof_var <- with (data, which (F_Coding == "Coding" & (
    (effect_priority %in% eff_lof.chv) | 
      (typeseq_priority %in% c ("splicing") & distance_spliceJunction <= 2) )))
  data$F_DamageType[lof_var] <- "LOF"
  
  data$F_DamageTier[lof_var] <- 2
  
  return(data)
}

# 3.4.2. Missense
add_missense_tag <- function(data, REVEL_t1_cutoff, REVEL_t2_cutoff, phylopMam_missense_cutoff, 
                             phylopVert_missense_cutoff, CADD_phred_missense_cutoff, MPC_cutoff) {
  
  missense_tier1_var <- with(data, which( F_Coding == "Coding" &
                                            (effect_priority %in% eff_missn.chv &
                                               (REVEL_score >= REVEL_t1_cutoff | 
                                                  ((phylopMam_avg >= phylopMam_missense_cutoff & phylopVert100_avg >= phylopVert_missense_cutoff) | 
                                                     CADD_phred >= CADD_phred_missense_cutoff | MPC_score >= MPC_cutoff) 
                                               )) ) ) 
  
  missense_tier2_var <- with(data, which( F_Coding == "Coding" &
                                            (effect_priority %in% eff_missn.chv &
                                               REVEL_score >= REVEL_t2_cutoff) ))
  
  data$F_DamageType[missense_tier1_var] <- "Missense"
  data$F_DamageTier[missense_tier1_var] <- 1
  data$F_DamageTier[missense_tier2_var] <- 2
  
  return(data)
}

# 3.4.3. Other coding
add_otherc_tag <- function(data, cutoffs, tier) {
  
  otherc_dmg_var <- with(data, which(F_Coding == "Coding" & 
                                       effect_priority %in% eff_other_sub.chv & 
                                       (phylopMam_avg >= cutoffs[1] | phylopVert100_avg >= cutoffs[2] | CADD_phred >= cutoffs[3]) & 
                                       is.na (dbsnp_common)   ))
  
  data$F_DamageType[otherc_dmg_var] <- "OtherC"
  data$F_DamageTier[otherc_dmg_var] <- tier
  
  return(data)
} 

# 3.4.4. Splicing predictions
add_splicing_tag <- function(data, spliceAI_DS_AG_cutoff, spliceAI_DP_AG_cutoff, 
                             spliceAI_DS_AL_cutoff, spliceAI_DP_AL_cutoff, 
                             spliceAI_DS_DG_cutoff, spliceAI_DP_DG_cutoff,
                             spliceAI_DS_DL_cutoff, spliceAI_DP_DL_cutoff, 
                             dbsc_SNV_ADA_SCORE_cutoff = 0, 
                             dbscSNV_RF_SCORE_cutoff = 0, tier) {
  
  if (tier == 1) {
    splicing_var <- with(data, which (
      ((spliceAI_DS_AG > spliceAI_DS_AG_cutoff & abs (spliceAI_DP_AG) <= spliceAI_DP_AG_cutoff) |
         (spliceAI_DS_AL > spliceAI_DS_AL_cutoff & abs (spliceAI_DP_AL) <= spliceAI_DP_AL_cutoff) |
         (spliceAI_DS_DG > spliceAI_DS_DG_cutoff & abs (spliceAI_DP_DG) <= spliceAI_DP_DG_cutoff) |
         (spliceAI_DS_DL > spliceAI_DS_DL_cutoff & abs (spliceAI_DP_DL) <= spliceAI_DP_DL_cutoff) |
         (var_type %in% c ("del", "ins") & typeseq_priority %in% c ("splicing", "exonic;splicing")) ) & 
        ! F_DamageType %in% c ("LOF", "Splc") & ! (F_DamageType %in% "Missense" & F_DamageTier == 2) ))
  }
  else if (tier == 2) {
    splicing_var <- with(data, which (
      ((spliceAI_DS_AG > spliceAI_DS_AG_cutoff & abs (spliceAI_DP_AG) <= spliceAI_DP_AG_cutoff) | 
         (spliceAI_DS_AL > spliceAI_DS_AL_cutoff & abs (spliceAI_DP_AL) <= spliceAI_DP_AL_cutoff) |
         (spliceAI_DS_DG > spliceAI_DS_DG_cutoff & abs (spliceAI_DP_DG) <= spliceAI_DP_DG_cutoff) |
         (spliceAI_DS_DL > spliceAI_DS_DL_cutoff & abs (spliceAI_DP_DL) <= spliceAI_DP_DL_cutoff) |
         (dbscSNV_ADA_SCORE > dbscSNV_ADA_SCORE_cutoff | dbscSNV_RF_SCORE > dbscSNV_RF_SCORE_cutoff)) & 
        ! F_DamageType %in% c ("LOF") ))
  }
  
  data$F_DamageType[splicing_var] <- "Splc"
  data$F_DamageTier[splicing_var] <- tier
  
  return(data)
}

# 3.4.5. UTR
add_utr_dmg_tag <- function(data, phylopMam_cutoff, CADD_phred_cutoff, tier) {
  
  utr_dmg_var <- with(data, which(
    typeseq_priority %in% typeseq_utr.chv & 
      (phylopMam_avg >= phylopMam_cutoff | CADD_phred >= CADD_phred_cutoff) &
      (! is.na (phastCons_placental)) ))
  
  data$F_DamageType[utr_dmg_var] <- "UTR"
  data$F_DamageTier[utr_dmg_var] <- tier
  
  return(data)
}

# 3.4.6. Non-coding
add_nc_dmg_tag <- function(data, cutoffs, tier) {
  
  nc_dmg_var <- with(data, which(F_Coding == "ncRNA" & ! gene_desc %in% "pseudogene" & F_DamageTier == 0 & (
    ((phylopMam_avg >= cutoffs[1] | phylopVert100_avg >= cutoffs[2]) & (! is.na (phastCons_placental))) | CADD_phred >= cutoffs[3])))
  
  data$F_DamageType[nc_dmg_var] <- "DmgNcRNA"
  data$F_DamageTier[nc_dmg_var] <- tier
  
  return(data)
}

# 3.5. Phenotype filter 

# 3.5.1. and 3.5.2. OMIM / CGD dominant
add_dominant_gene_tag <- function(data, database, pattern) {
  
  data$F_AXD <- 0
  
  dom_var <- intersect (
    grep (data$omim_phenotype,  pattern = "i:AD", fixed = T),
    grep (data$CGD_inheritance, pattern = "AD",   fixed = T) )
  data$F_AXD[dom_var] <- 1
  
  return(data)
}

# 3.5.3. Phenotype tiers
add_phenotype_tier_tag <- function(data) {

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
  
   data$F_PhenoTier <- 0
   data$F_PhenoTier[pheno_mouse_var] <- 1
   data$F_PhenoTier[pheno_human_var] <- 2
   
   return(data)
}

# 3.6. Main Findings

# 3.6.1. RECESSIVE HOMOZYGOUS
add_recessive_homozygous_tag <- function(data) {
  
  rec_hom_var <- with(data, which(F_Rare <= 0.05 & F_DamageType != "NotDmg" & Zygosity == "hom-alt"))
  data$FM_HOM <- 0
  data$FM_HOM[rec_hom_var] <- 1
  
  return(data)
}

# 3.6.2. X-LINKED HAPLOID
add_X_haploid_tag <- function(data) {
  
  X_haploid_var <- with(data, which(F_Rare <= 1e-4 & F_DamageType != "NotDmg" & Zygosity == "hap" & CHROM == "chrX"))
  data$FM_XHAP <- 0
  data$FM_XHAP[X_haploid_var] <- 1
  
  return(data)
}

# 3.6.3. POTENTIAL COMPOUND HETS

get_potential_cmphet_df <- function(cmp_hets_df) { # help function for add_potential_compound_heterozygous_tag
  
  cmp_hets_df <- unique(cmp_hets_df)
  cmp_hets_df_final <- subset(cmp_hets_df, 
                              subset = gene_symbol %in% na.omit(cmp_hets_df$gene_symbol[duplicated(cmp_hets_df$gene_symbol)]), 
                              select = c(Original_VCFKEY, gene_symbol))
  
  return(cmp_hets_df_final)
}

#' This function is used for both main and secondary findings
add_potential_compound_heterozygous_tag <- function(data, secondary_findings = FALSE) {
  
  if (secondary_findings == TRUE) {
    col_name <- "F_CmpHet_S1"
    # used F_Pass == 1 instead of F_Qual > 1 in the original script for variants with PASS FILTER
    cmp_hets_r1_df <- subset(data, subset = F_Pass == 1 & FS1_Select == 1, select = c(gene_symbol, Original_VCFKEY))
    cmp_hets_r2_df <- subset(data, subset = F_Qual == 1 & FS1_Select == 1, select = c(gene_symbol, Original_VCFKEY))
  }
  else {
    col_name <- "FM_PCHET"
    cmp_hets_r1_df <- subset(data, 
                             subset = F_Rare <= 0.05 & F_DamageType != "NotDmg", 
                             select = c(gene_symbol, Original_VCFKEY))
    cmp_hets_r2_df <- subset(data, 
                             subset = F_Rare <= 0.05 & F_DamageType != "NotDmg" & F_Qual == 1, 
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
  
  dmg_cmp_hets_df <- subset(data, subset = FM_PCHET == 2 & F_DamageTier == 2, select = gene_symbol, drop = T)
  dmg_cmp_hets_df_final <- dmg_cmp_hets_df[duplicated(dmg_cmp_hets_df)]
  
  data$FM_PCHET_DMG <- 0 
  data$FM_PCHET_DMG[with(data, which(gene_symbol %in% dmg_cmp_hets_df_final & FM_PCHET == 2))] <- 1 
  
  return(data)
}

# 3.6.4. DOMINANT
add_dominant_tag <- function(data) {
  
  dom_var <- with(data, which(F_Rare <= 1e-4 & F_DamageType != "NotDmg" & F_AXD == 1))
  
  data$FM_AXDOM <- 0
  data$FM_AXDOM[dom_var] <- 1
  
  return(data)
}

# 3.6.5. PREDICTED DOMINANT
add_predicted_dominant_tag <- function(data, gnomAD_oe_lof_upper_cutoff) {
  
  pred_dom_var <- with(data, which(Zygosity == "ref-alt" & F_Rare <= 1e-4 & (
    (F_GeneConstr == 1 & 
       (F_DamageType %in% c ("LOF", "Splc", "OtherC") | (F_DamageType == "Missense" & F_DamageTier == 2) )) )))
  
  data$FM_PDDOM <- 0
  data$FM_PDDOM[pred_dom_var] <- 1
  
  return(data)
}

# 3.7. Secondary Findings

# 3.7.0. Pathogenicity flag
add_pathogenic_tag <- function(data) {
  
  clinvar_not_pathg_var <- with(data, which(Clinvar_SIG_Simple == 0))
  clinvar_pathg_var <- with(data, which(Clinvar_SIG_Simple == 1))
  
  data$F_Clinvar_notPathg <- 0
  data$F_Clinvar_Pathg <- 0
  data$F_Clinvar_notPathg[clinvar_not_pathg_var] <- 1
  data$F_Clinvar_Pathg[clinvar_pathg_var] <- 1
  
  return(data)
}

# 3.7.1. Tier 1
add_secondary_findings_tier1_tag <- function(data) {
  
  sf1_var <- with(data, which(CGD_disease != "" & ((F_DamageType == "LOF") | (F_Clinvar_Pathg == 1))))
  
  data$FS1_Select <- 0
  data$FS1_Select[sf1_var] <- 1
  
  return(data)
}

# 3.7.2. Tier 2
add_secondary_findings_tier2_tag <- function(data) {
  
  # genes with 1) CGD annotations; 2) hq; 3) rare 1%; 4) not non-pathogenic
  sf2_var <- with(data, which(CGD_disease != "" & F_Qual == 1  
                              & F_Rare <= 0.01 & ((F_DamageType == "LOF") | 
                                                    (F_Clinvar_Pathg == 1) | 
                                                    (F_Clinvar_notPathg == 0)) ))
  data$FS2_Select <- 0
  data$FS2_Select[sf2_var] <- 1
  
  return(data)
}

# 3.7.3. Tier 3
add_secondary_findings_tier3_tag <- function(data) {
  
  # genes with 1) CGD annotations; 2) hq; 3) rare 1%; 4) not non-pathogenic; 5) damaging
  sf3_var <- with(data, which(CGD_disease != "" & F_Qual == 1 & F_Rare <= 0.01 
                              & F_DamageType != "NotDmg" & 
                                ((F_DamageType == "LOF") | 
                                   (F_Clinvar_Pathg == 1) | 
                                   (F_Clinvar_notPathg == 0)) )) 
  data$FS3_Select <- 0
  data$FS3_Select[sf3_var] <- 1
  
  return(data)
}

# 3.7.4. dominant: pathogenic
add_dominant_pathogenic_tag <- function(data) {
  
  dom_pathg_var <- with(data, which(FS1_Select == 1 & grepl("AD", CGD_inheritance)))
  
  data$FS1_AD_Pathg_Any <- 0
  data$FS1_AD_Pathg_Any[dom_pathg_var] <- 1
  
  return(data)
}

# 3.7.5. recessive + hom: pathogenic
add_recessive_hom_pathg_tag <- function(data) {
  
  rec_hom_pathg_var <- with(data, which(FS1_Select == 1 & grepl("AR", CGD_inheritance) & Zygosity == "hom-alt"))
  
  data$FS1_AR_Pathg_Hom <- 0
  data$FS1_AR_Pathg_Hom[rec_hom_pathg_var] <- 1
  
  return(data)
}

# 3.7.6. recessive + potential compound het: pathogenic
add_recessive_cmphet_pathg_tag <- function(data) {
  
  rec_cmphet_pathg_var <- with(data, which(FS1_Select == 1 & grepl("AR", CGD_inheritance) & F_CmpHet_S1 >= 1))
  
  data$FS1_AR_Pathg_PotCompHet <- 0
  data$FS1_AR_Pathg_PotCompHet[rec_cmphet_pathg_var] <- 1
  
  return(data)
}

# 3.7.7. X-linked + hom/hap: pathogenic
add_X_hom_or_hap_pathg_tag <- function(data, type) {
  
  X_hom_or_hap_pathg_var <- with(data, which(FS1_Select == 1 & grepl("XL", CGD_inheritance) & Zygosity == "hom-alt" & CHROM == "chrX"))
  ifelse(type == "hom", col_name <- "FS1_XL_Pathg_Hom", col_name <- "FS1_XL_Pathg_Hap")
  
  data[, col_name] <- 0
  data[, col_name][X_hom_or_hap_pathg_var] <- 1
  
  return(data)
}

# 3.7.8. complex + hom / hap: pathogenic
add_complex_hom_hap_pathg_tag <- function(data) {
  
  complex_hom_hap_pathg_tag <- with(data, which(FS1_Select == 1 & (! grepl("AD|AR|XL", CGD_inheritance)) & (Zygosity == "hom-alt")))
  
  data$FS1_CX_Pathg_HomHap <- 0
  data$FS1_CX_Pathg_HomHap[complex_hom_hap_pathg_tag] <- 1
  
  return(data)
}

# 3.7.9. complex + potential compound het: pathogenic
add_complex_cmphet_pathg_tag <- function(data) {
  
  complex_cmphet_pathg_var <- with(data, which(FS1_Select == 1 & (! grepl("AD|AR|XL", CGD_inheritance)) & F_CmpHet_S1 >= 1))
  
  data$FS1_CX_Pathg_PotCompHet <- 0
  data$FS1_CX_Pathg_PotCompHet[complex_cmphet_pathg_var] <- 1
  
  return(data)
}

# 3.7.10. complex + single het: uncertain
add_complex_single_het_uncertain_tag <- function(data) {
  
  complex_single_het_uncertain_var <- with(data, 
                                           which(FS1_Select == 1 & (! grepl("AD|AR|XL", CGD_inheritance)) & 
                                                   Zygosity %in% c ("ref-alt", "alt-alt") & F_CmpHet_S1 == 0))
  
  data$FS1_CX_Uncertain <- 0
  data$FS1_CX_Uncertain[complex_single_het_uncertain_var] <- 1
  
  return(data)
}

# 3.7.11. recessive + single het: carrier
add_recessive_single_het_carrier_tag <- function(data) {
  
  rec_single_het_var <- with(data, which(FS1_Select == 1 & grepl("AR", CGD_inheritance) & 
                                           Zygosity %in% c ("ref-alt", "alt-alt") & F_CmpHet_S1 == 0))
  
  data$FS1_AR_Carrier <- 0
  data$FS1_AR_Carrier[rec_single_het_var] <- 1
  
  return(data)
}

# 3.7.12. X-linked + het: carrier
add_X_het_carrier_tag <- function(data){
  
  X_het_carrier_var <- with(data, which(FS1_Select == 1 & grepl("XL", CGD_inheritance) & 
                                          Zygosity %in% c ("ref-alt", "alt-alt") & CHROM == "chrX"))
  
  data$FS1_XL_Carrier <- 0
  data$FS1_XL_Carrier[X_het_carrier_var] <- 1
  
  return(data)
}

# 3.7.13. ACMG Disease
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

# 3.8. Get final results

# 3.8.1. Rare05 Variants

get_rare05_variants <- function(data, multisample) {
  
  data$F_Rare <- 2
  data <- add_freq_tag(data, freq_max_cutoff = 0.05)
  data <- add_freq_tag(data, freq_max_cutoff = 0.01)
  data <- add_freq_tag(data, freq_max_cutoff = 0.001)
  data <- add_freq_tag(data, freq_max_cutoff = 0.0001)
  data <- add_freq_tag(data, freq_max_cutoff = 0)
  
  data <- subset(data, subset = F_Rare <= 0.05)
  data <- add_pass_tag(data)
  
  ifelse(multisample,
         data <- add_qual_tag_multisample(data, DP_cutoff),
         data <- add_qual_tag(data, DP_cutoff))
  
  data <- add_coding_tag(data)
  
  data$F_DamageType <- "NotDmg"
  data$F_DamageTier <- 0
  
  data <- add_coding_lof_tag(data, eff_lof.chv)
  
  data <- add_missense_tag(data, REVEL_t1_cutoff, REVEL_t2_cutoff, phylopMam_missense_cutoff, 
                           phylopVert_missense_cutoff, CADD_phred_missense_cutoff, MPC_cutoff)
  
  data <- add_otherc_tag(data, otherc_t1_cutoffs, tier = 1)
  data <- add_otherc_tag(data, otherc_t2_cutoffs, tier = 2)
  
  data <- add_splicing_tag(data, spliceAI_DS_AG_t1_cutoff, spliceAI_DP_AG_t1_cutoff, 
                           spliceAI_DS_AL_t1_cutoff, spliceAI_DP_AL_t1_cutoff,
                           spliceAI_DS_DG_t1_cutoff, spliceAI_DP_DG_t1_cutoff,
                           spliceAI_DS_DL_t1_cutoff, spliceAI_DP_DL_t1_cutoff, tier = 1)
  
  data <- add_splicing_tag(data, spliceAI_DS_AG_t2_cutoff, spliceAI_DP_AG_t2_cutoff,
                           spliceAI_DS_AL_t2_cutoff, spliceAI_DP_AL_t2_cutoff,
                           spliceAI_DS_DG_t2_cutoff, spliceAI_DP_DG_t2_cutoff,
                           spliceAI_DS_DL_t2_cutoff, spliceAI_DP_DL_t2_cutoff, 
                           dbscSNV_ADA_SCORE_cutoff, dbscSNV_RF_SCORE_cutoff, tier = 2)
  
  data <- add_utr_dmg_tag(data, phylopMam_utr_t1_cutoff,
                          CADD_phred_utr_t1_cutoff, tier = 1)
  data <- add_utr_dmg_tag(data, phylopMam_utr_t2_cutoff,
                          CADD_phred_utr_t2_cutoff, tier = 2)
  
  data <- add_nc_dmg_tag(data, nc_t1_cutoffs, tier = 1)
  data <- add_nc_dmg_tag(data, nc_t2_cutoffs, tier = 2)
  
  data <- add_dominant_gene_tag(data)
  
  data <- add_phenotype_tier_tag(data)
  
  data <- add_recessive_homozygous_tag(data)
  data <- add_X_haploid_tag(data)
  data <- add_potential_compound_heterozygous_tag(data, secondary_findings = F)
  data <- add_potential_dmg_compound_heterozygous_tag(data)
  data <- add_dominant_tag(data)
  data <- add_predicted_dominant_tag(data, gnomAD_oe_lof_upper_cutoff)
  
  data <- add_pathogenic_tag(data)
  data <- add_secondary_findings_tier1_tag(data)
  data <- add_secondary_findings_tier2_tag(data)
  data <- add_secondary_findings_tier3_tag(data)
  
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

# 3.8.2. HQ Variants (FILTER == "PASS")

get_hq_variants <- function(data, multisample) {
  
  if (multisample) {
    data <- subset(data, subset = (FT == "PASS"))
    data <- add_qual_tag_multisample(data, DP_cutoff)
  } else {
    data <- subset(data, subset = (FILTER == "PASS"))
    data <- add_qual_tag(data, DP_cutoff) 
  }
  
  return(data)
}

# 3.9. Stats

get_chr_zygosity_stats <- function(data) {
  
  chrom_list <- paste("chr", c(1:22, "X", "Y", "M"), sep = "")
  data <- as.data.table(data)
  data <- data[CHROM %in% chrom_list]
  
  chr_counts.mx <- table(data$CHROM) %>% as.matrix()
  zygosity_counts.mx <- table(data[, c("CHROM", "Zygosity")]) %>% as.data.frame.matrix()
  
  chr_zygosity_counts.mx <- cbind(chr_counts.mx, zygosity_counts.mx)
  chr_zygosity_counts.mx <- chr_zygosity_counts.mx[order(factor(rownames(chr_zygosity_counts.mx), levels = chrom_list)), ]
  
  chr_zygosity_counts.mx$HomPerc <- (chr_zygosity_counts.mx[, "hom-alt"]) / 
    (chr_zygosity_counts.mx[, "ref-alt"] + chr_zygosity_counts.mx[, "hom-alt"] + chr_zygosity_counts.mx[, "alt-alt"]) * 100
  
  colnames(chr_zygosity_counts.mx)[1] <- "Freq"
  chr_zygosity_counts.mx <- tibble::rownames_to_column(chr_zygosity_counts.mx, "Chrom")
  
  return(chr_zygosity_counts.mx)
}

add_stats_unique_vcfkey <- function(data, high_quality = FALSE) {
  
  if (high_quality == TRUE) {
    num_var <- length(unique(subset(data, subset = F_Qual >= 1, 
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
    cond <- expression(F_Rare <= freq_cutoff & F_Qual >= 1)
  } else {
    cond <- expression(F_Rare <= freq_cutoff)
  }
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_typeseq_tag <- function(data, typeseq.chv, freq_cutoff = -1, high_quality = FALSE) {
  
  if (high_quality == TRUE) {
    if (freq_cutoff != -1) {
      cond <- expression(typeseq_priority %in% typeseq.chv & F_Qual >= 1 & F_Rare <= freq_cutoff)
    } else {
      cond <- expression(typeseq_priority %in% typeseq.chv & F_Qual >= 1)
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
      cond <- expression(Zygosity == zygosity & F_Qual >= 1 & CHROM == "chrX")
    } else {
      cond <- expression(Zygosity == zygosity & F_Qual >= 1)
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

add_stats_dmg_tag <- function(data, freq_cutoff, dmg_type, dmg_tier = 0, high_quality) {
  
  if (dmg_tier > 0) {
    if (high_quality == TRUE) {
      cond <- expression(F_Rare <= freq_cutoff & F_DamageType == dmg_type & F_DamageTier >= dmg_tier & F_Qual >= 1)
    } else {
      cond <- expression(F_Rare <= freq_cutoff & F_DamageType == dmg_type & F_DamageTier >= dmg_tier & F_Pass == 1)
    }
  } else {
    if (high_quality == TRUE) {
      cond <- expression(F_Rare <= freq_cutoff & F_DamageType == dmg_type & F_Qual >= 1)
    } else {
      cond <- expression(F_Rare <= freq_cutoff & F_DamageType == dmg_type & F_Pass == 1)
    }
  }
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_pheno_tier <- function(data, freq_cutoff, pheno_tier, high_quality = FALSE) {
  
  if (high_quality == TRUE) {
    cond <- expression(F_Rare <= freq_cutoff & F_DamageType != "NotDmg" & F_PhenoTier >= pheno_tier & F_Qual >= 1)
  } else {
    cond <- expression(F_Rare <= freq_cutoff & F_DamageType != "NotDmg" & F_PhenoTier >= pheno_tier & F_Pass == 1)
  }
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_main_findings <- function(data, type, pheno_tier, dmg_tier = 0, high_quality = FALSE) {
  
  switch(type,
         hom = type_cond <- expression(FM_HOM == 1),
         Xhap = type_cond <- expression(FM_XHAP == 1),
         cmp_het = ifelse(high_quality == TRUE, type_cond <- expression(FM_PCHET >= 2), type_cond <- expression(FM_PCHET >= 1)),
         dominant = type_cond <- expression(FM_AXDOM == 1),
         pred_dom = type_cond <- expression(FM_PDDOM == 1))
  
  if (type == "cmp_het") {
    cond <- expression(eval(type_cond) & F_PhenoTier >= pheno_tier & F_DamageTier >= dmg_tier)
  } else {
    if (high_quality == TRUE) {
      cond <- expression(F_PhenoTier >= pheno_tier & F_DamageTier >= dmg_tier & eval(type_cond) & F_Qual >= 1)
    } else {
      cond <- expression(F_PhenoTier >= pheno_tier & F_DamageTier >= dmg_tier & eval(type_cond) & F_Pass == 1)
    }
  }
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_secondary_findings <- function(data, tier, type = NA) {
  
  switch(tier,
         "1" = tier_cond <- expression(FS1_Select == 1),
         "2" = tier_cond <- expression(FS2_Select == 1),
         "3" = tier_cond <- expression(FS3_Select == 1))
  
  if (is.na(type)) {
    cond <- tier_cond
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
  
  cond <- expression(eval(tier_cond) & eval(type_cond))
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

add_stats_acmg <- function(data, coding = FALSE) {
  
  ifelse(coding, 
         cond <- expression(F_ACMG_Coding == 1),
         cond <- expression(F_ACMG == 1))
  
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
  
  stats.ls$VarN_LQ_AllSeq_AllFreq <- add_stats_unique_vcfkey(data, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_AllFreq <- add_stats_unique_vcfkey(data, high_quality = T)
  
  stats.ls$VarN_LQ_Coding_AllFreq <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, high_quality = F)
  stats.ls$VarN_HQ_Coding_AllFreq <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, high_quality = T)
  
  stats.ls$VarN_LQ_ncRNA_AllFreq <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, high_quality = F)
  stats.ls$VarN_HQ_ncRNA_AllFreq <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, high_quality = T)
  
  stats.ls$VarN_LQ_AllSeq_AllFreq_Xhom <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = T, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_AllFreq_Xhom <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = T, high_quality = T)
  
  stats.ls$VarN_LQ_AllSeq_AllFreq_Xhet <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = T, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_AllFreq_Xhet <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = T, high_quality = T)
  
  stats.ls$VarN_LQ_AllSeq_AllFreq_Hom  <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = F, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_AllFreq_Hom  <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = F, high_quality = T)
  
  stats.ls$VarN_LQ_AllSeq_AllFreq_HetR <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = F, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_AllFreq_HetR <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = F, high_quality = T)
  stats.ls$VarN_LQ_AllSeq_AllFreq_HetA <- add_stats_zygosity_tag(data, zygosity = "alt-alt", chrX = F, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_AllFreq_HetA <- add_stats_zygosity_tag(data, zygosity = "alt-alt", chrX = F, high_quality = T)
  
  return(stats.ls)
}

#' input: data = v_full_r05.df
add_stats_full_r05_df <- function(data) {
  
  stats.ls <- list()
  
  stats.ls$VarN_LQ_AllSeq_Rare05 <- add_stats_freq(data, freq_cutoff = 0.05, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_Rare05 <- add_stats_freq(data, freq_cutoff = 0.05, high_quality = T)
  stats.ls$VarN_LQ_AllSeq_Rare01 <- add_stats_freq(data, freq_cutoff = 0.01, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_Rare01 <- add_stats_freq(data, freq_cutoff = 0.01, high_quality = T)
  stats.ls$VarN_LQ_AllSeq_Rare001 <- add_stats_freq(data, freq_cutoff = 0.001, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_Rare001 <- add_stats_freq(data, freq_cutoff = 0.001, high_quality = T)
  stats.ls$VarN_LQ_AllSeq_Rare0001 <- add_stats_freq(data, freq_cutoff = 0.0001, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_Rare0001 <- add_stats_freq(data, freq_cutoff = 0.0001, high_quality = T)
  stats.ls$VarN_LQ_AllSeq_Rare00 <- add_stats_freq(data, freq_cutoff = 0, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_Rare00 <- add_stats_freq(data, freq_cutoff = 0, high_quality = T)
  
  stats.ls$VarN_LQ_Coding_Rare05  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.05, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare05  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.05, high_quality = T)
  stats.ls$VarN_LQ_Coding_Rare01  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.01, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare01  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.01, high_quality = T)
  stats.ls$VarN_LQ_Coding_Rare001  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.001, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare001  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.001, high_quality = T)
  stats.ls$VarN_LQ_Coding_Rare0001 <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.0001, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare0001 <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0.0001, high_quality = T)
  stats.ls$VarN_LQ_Coding_Rare00  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare00  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_coding.chv, freq_cutoff = 0, high_quality = T)
  
  stats.ls$VarN_LQ_ncRNA_Rare05  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.05, high_quality = F)
  stats.ls$VarN_HQ_ncRNA_Rare05  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.05, high_quality = T)
  stats.ls$VarN_LQ_ncRNA_Rare01  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.01, high_quality = F)
  stats.ls$VarN_HQ_ncRNA_Rare01  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.01, high_quality = T)
  stats.ls$VarN_LQ_ncRNA_Rare001  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.001, high_quality = F)
  stats.ls$VarN_HQ_ncRNA_Rare001  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.001, high_quality = T)
  stats.ls$VarN_LQ_ncRNA_Rare0001 <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.0001, high_quality = F)
  stats.ls$VarN_HQ_ncRNA_Rare0001 <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0.0001, high_quality = T)
  stats.ls$VarN_LQ_ncRNA_Rare00  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0, high_quality = F)
  stats.ls$VarN_HQ_ncRNA_Rare00  <- add_stats_typeseq_tag(data, typeseq.chv = typeseq_ncrna.chv, freq_cutoff = 0, high_quality = T)
  
  stats.ls$VarN_LQ_AllSeq_Rare05_Xhom <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = T, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_Rare05_Xhom <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = T, high_quality = T)
  stats.ls$VarN_LQ_AllSeq_Rare05_Hom  <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = F, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_Rare05_Hom  <- add_stats_zygosity_tag(data, zygosity = "hom-alt", chrX = F, high_quality = T)
  stats.ls$VarN_LQ_AllSeq_Rare05_HetR <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = F, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_Rare05_HetR <- add_stats_zygosity_tag(data, zygosity = "ref-alt", chrX = F, high_quality = T)
  stats.ls$VarN_LQ_AllSeq_Rare05_HetA <- add_stats_zygosity_tag(data, zygosity = "alt-alt", chrX = F, high_quality = F)
  stats.ls$VarN_HQ_AllSeq_Rare05_HetA <- add_stats_zygosity_tag(data, zygosity = "alt-alt", chrX = F, high_quality = T)
  
  stats.ls$GeneN_Pheno_Rare05 <- length (unique (data[data$F_PhenoTier == 2, ]$entrez_id))
  stats.ls$GeneN_Dominant_Rare05 <- length (setdiff (c (data$entrez_id[grep ("i:AD", data$omim_phenotype)], 
                                                        data$entrez_id[grep ("AD", data$CGD_inheritance)]), NA))
  
  stats.ls$VarN_LQ_Coding_Rare01_LOF_TierLow <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "LOF", dmg_tier = 1, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare01_LOF_TierLow <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "LOF", dmg_tier = 1, high_quality = T)
  
  stats.ls$VarN_LQ_Coding_Rare01_LOF_TierHigh <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "LOF", dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare01_LOF_TierHigh <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "LOF", dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Coding_Rare01_MissDmg_TierLow <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Missense", dmg_tier = 1, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare01_MissDmg_TierLow <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Missense", dmg_tier = 1, high_quality = T)
  
  stats.ls$VarN_LQ_Coding_Rare01_MissDmg_TierHigh <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Missense", dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare01_MissDmg_TierHigh <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Missense", dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Coding_Rare01_Splc_TierLow <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Splc", dmg_tier = 1, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare01_Splc_TierLow <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Splc", dmg_tier = 1, high_quality = T)
  
  stats.ls$VarN_LQ_Coding_Rare01_Splc_TierHigh <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Splc", dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare01_Splc_TierHigh <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "Splc", dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Coding_Rare01_OtherDmg <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "OtherC", high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare01_OtherDmg <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "OtherC", high_quality = T)
  
  stats.ls$VarN_LQ_ncRNA_Rare01_Dmg_TierLow <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "DmgNcRNA", dmg_tier = 1, high_quality = F)
  stats.ls$VarN_HQ_ncRNA_Rare01_Dmg_TierLow <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "DmgNcRNA", dmg_tier = 1, high_quality = T)
  
  stats.ls$VarN_LQ_ncRNA_Rare01_Dmg_TierHigh <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "DmgNcRNA", dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_ncRNA_Rare01_Dmg_TierHigh <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "DmgNcRNA", dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Coding_Rare01_Dmg_PhenoTierLow <- add_stats_pheno_tier(data, freq_cutoff = 0.01, pheno_tier = 1, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare01_Dmg_PhenoTierLow <- add_stats_pheno_tier(data, freq_cutoff = 0.01, pheno_tier = 1, high_quality = T)
  
  stats.ls$VarN_LQ_Coding_Rare01_Dmg_PhenoTierHigh <- add_stats_pheno_tier(data, freq_cutoff = 0.01, pheno_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Coding_Rare01_Dmg_PhenoTierHigh <- add_stats_pheno_tier(data, freq_cutoff = 0.01, pheno_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Rare05_DmgT2_Hom_PhenoTierAll <- add_stats_main_findings(data, type = "hom", pheno_tier = 0, dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Rare05_DmgT2_Hom_PhenoTierAll <- add_stats_main_findings(data, type = "hom", pheno_tier = 0, dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Rare05_DmgT2_Hom_PhenoTierLow <- add_stats_main_findings(data, type = "hom", pheno_tier = 1, dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Rare05_DmgT2_Hom_PhenoTierLow <- add_stats_main_findings(data, type = "hom", pheno_tier = 1, dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Rare05_DmgT2_Hom_PhenoTierHigh <- add_stats_main_findings(data, type = "hom", pheno_tier = 2, dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Rare05_DmgT2_Hom_PhenoTierHigh <- add_stats_main_findings(data, type = "hom", pheno_tier = 2, dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Rare0001_DmgT2_Xhap_PhenoTierAll <- add_stats_main_findings(data, type = "Xhap", pheno_tier = 0, dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Rare0001_DmgT2_Xhap_PhenoTierAll <- add_stats_main_findings(data, type = "Xhap", pheno_tier = 0, dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Rare0001_DmgT2_Xhap_PhenoTierLow <- add_stats_main_findings(data, type = "Xhap", pheno_tier = 1, dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Rare0001_DmgT2_Xhap_PhenoTierLow <- add_stats_main_findings(data, type = "Xhap", pheno_tier = 1, dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Rare0001_DmgT2_Xhap_PhenoTierHigh <- add_stats_main_findings(data, type = "Xhap", pheno_tier = 2, dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Rare0001_DmgT2_Xhap_PhenoTierHigh <- add_stats_main_findings(data, type = "Xhap", pheno_tier = 2, dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Rare05_DmgT1_CmpHet_PhenoTierAll <- add_stats_main_findings(data, type = "cmp_het", pheno_tier = 0, dmg_tier = 1, high_quality = F)
  stats.ls$VarN_HQ_Rare05_DmgT1_CmpHet_PhenoTierAll <- add_stats_main_findings(data, type = "cmp_het", pheno_tier = 0, dmg_tier = 1, high_quality = T)
  
  stats.ls$VarN_LQ_Rare05_DmgT1_CmpHet_PhenoTierLow <- add_stats_main_findings(data, type = "cmp_het", pheno_tier = 1, dmg_tier = 1, high_quality = F)
  stats.ls$VarN_HQ_Rare05_DmgT1_CmpHet_PhenoTierLow <- add_stats_main_findings(data, type = "cmp_het", pheno_tier = 1, dmg_tier = 1, high_quality = T)
  
  stats.ls$VarN_LQ_Rare05_DmgT1_CmpHet_PhenoTierHigh <- add_stats_main_findings(data, type = "cmp_het", pheno_tier = 2, dmg_tier = 1, high_quality = F)
  stats.ls$VarN_HQ_Rare05_DmgT1_CmpHet_PhenoTierHigh <- add_stats_main_findings(data, type = "cmp_het", pheno_tier = 2, dmg_tier = 1, high_quality = T)
  
  stats.ls$VarN_LQ_Rare0001_DmgT2_AXDom_PhenoTierAll <- add_stats_main_findings(data, type = "dominant", pheno_tier = 0, dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Rare0001_DmgT2_AXDom_PhenoTierAll <- add_stats_main_findings(data, type = "dominant", pheno_tier = 0, dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Rare0001_DmgT2_AXDom_PhenoTierLow <- add_stats_main_findings(data, type = "dominant", pheno_tier = 1, dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Rare0001_DmgT2_AXDom_PhenoTierLow <- add_stats_main_findings(data, type = "dominant", pheno_tier = 1, dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Rare0001_DmgT2_AXDom_PhenoTierHigh <- add_stats_main_findings(data, type = "dominant", pheno_tier = 2, dmg_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Rare0001_DmgT2_AXDom_PhenoTierHigh <- add_stats_main_findings(data, type = "dominant", pheno_tier = 2, dmg_tier = 2, high_quality = T)
  
  stats.ls$VarN_LQ_Rare0001_DmgT2_PDDom_PhenoTierAll <- add_stats_main_findings(data, type = "pred_dom", pheno_tier = 0, high_quality = F)
  stats.ls$VarN_HQ_Rare0001_DmgT2_PDDom_PhenoTierAll <- add_stats_main_findings(data, type = "pred_dom", pheno_tier = 0, high_quality = T)
  
  stats.ls$VarN_LQ_Rare0001_DmgT2_PDDom_PhenoTierLow <- add_stats_main_findings(data, type = "pred_dom", pheno_tier = 1, high_quality = F)
  stats.ls$VarN_HQ_Rare0001_DmgT2_PDDom_PhenoTierLow <- add_stats_main_findings(data, type = "pred_dom", pheno_tier = 1, high_quality = T)
  
  stats.ls$VarN_LQ_Rare0001_DmgT2_PDDom_PhenoTierHigh <- add_stats_main_findings(data, type = "pred_dom", pheno_tier = 2, high_quality = F)
  stats.ls$VarN_HQ_Rare0001_DmgT2_PDDom_PhenoTierHigh <- add_stats_main_findings(data, type = "pred_dom", pheno_tier = 2, high_quality = T)
  
  stats.ls$VarN_FS1_LQ_Rare05_Tot     <- add_stats_secondary_findings(data, tier = 1)
  stats.ls$VarN_FS2_HQ_Rare01_Tot     <- add_stats_secondary_findings(data, tier = 2)
  stats.ls$VarN_FS3_HQ_Rare01_Dmg_Tot <- add_stats_secondary_findings(data, tier = 3)
  
  stats.ls$VarN_FS1_LQ_Rare05_AD_Pathg_Any        <- add_stats_secondary_findings(data, tier = 1, type = "dominant")
  stats.ls$VarN_FS1_LQ_Rare05_AR_Pathg_Hom        <- add_stats_secondary_findings(data, tier = 1, type = "recessive_hom")
  stats.ls$VarN_FS1_LQ_Rare05_AR_Pathg_PotCompHet <- add_stats_secondary_findings(data, tier = 1, type = "recessive_cmp_het")
  stats.ls$VarN_FS1_LQ_Rare05_XL_Pathg_Hap        <- add_stats_secondary_findings(data, tier = 1, type = "Xhap")
  stats.ls$VarN_FS1_LQ_Rare05_XL_Pathg_Hom        <- add_stats_secondary_findings(data, tier = 1, type = "Xhom")
  stats.ls$VarN_FS1_LQ_Rare05_CX_Pathg_HomHap     <- add_stats_secondary_findings(data, tier = 1, type = "complex_hom_hap")
  stats.ls$VarN_FS1_LQ_Rare05_CX_Pathg_PotCompHet <- add_stats_secondary_findings(data, tier = 1, type = "complex_cmp_het")
  stats.ls$VarN_FS1_LQ_Rare05_CX_Uncertain        <- add_stats_secondary_findings(data, tier = 1, type = "complex_single_het")
  stats.ls$VarN_FS1_LQ_Rare05_AR_Carrier          <- add_stats_secondary_findings(data, tier = 1, type = "recessive_single_het")
  stats.ls$VarN_FS1_LQ_Rare05_XL_Carrier          <- add_stats_secondary_findings(data, tier = 1, type = "Xhet")
  
  stats.ls$VarN_FS2_HQ_Rare01_AD_Pathg_Any        <- add_stats_secondary_findings(data, tier = 2, type = "dominant")
  stats.ls$VarN_FS2_HQ_Rare01_AR_Pathg_Hom        <- add_stats_secondary_findings(data, tier = 2, type = "recessive_hom")
  stats.ls$VarN_FS2_HQ_Rare01_AR_Pathg_PotCompHet <- add_stats_secondary_findings(data, tier = 2, type = "recessive_cmp_het")
  stats.ls$VarN_FS2_HQ_Rare01_XL_Pathg_Hap        <- add_stats_secondary_findings(data, tier = 2, type = "Xhap")
  stats.ls$VarN_FS2_HQ_Rare01_XL_Pathg_Hom        <- add_stats_secondary_findings(data, tier = 2, type = "Xhom")
  stats.ls$VarN_FS2_HQ_Rare01_CX_Pathg_HomHap     <- add_stats_secondary_findings(data, tier = 2, type = "complex_hom_hap")
  stats.ls$VarN_FS2_HQ_Rare01_CX_Pathg_PotCompHet <- add_stats_secondary_findings(data, tier = 2, type = "complex_cmp_het")
  stats.ls$VarN_FS2_HQ_Rare01_CX_Uncertain        <- add_stats_secondary_findings(data, tier = 2, type = "complex_single_het")
  stats.ls$VarN_FS2_HQ_Rare01_AR_Carrier          <- add_stats_secondary_findings(data, tier = 2, type = "recessive_single_het")
  stats.ls$VarN_FS2_HQ_Rare01_XL_Carrier          <- add_stats_secondary_findings(data, tier = 2, type = "Xhet")
  
  stats.ls$VarN_FS3_HQ_Rare01_Dmg_AD_Pathg_Any        <- add_stats_secondary_findings(data, tier = 3, type = "dominant")
  stats.ls$VarN_FS3_HQ_Rare01_Dmg_AR_Pathg_Hom        <- add_stats_secondary_findings(data, tier = 3, type = "recessive_hom")
  stats.ls$VarN_FS3_HQ_Rare01_Dmg_AR_Pathg_PotCompHet <- add_stats_secondary_findings(data, tier = 3, type = "recessive_cmp_het")
  stats.ls$VarN_FS3_HQ_Rare01_Dmg_XL_Pathg_Hap        <- add_stats_secondary_findings(data, tier = 3, type = "Xhap")
  stats.ls$VarN_FS3_HQ_Rare01_Dmg_XL_Pathg_Hom        <- add_stats_secondary_findings(data, tier = 3, type = "Xhom")
  stats.ls$VarN_FS3_HQ_Rare01_Dmg_CX_Pathg_HomHap     <- add_stats_secondary_findings(data, tier = 3, type = "complex_hom_hap")
  stats.ls$VarN_FS3_HQ_Rare01_Dmg_CX_Pathg_PotCompHet <- add_stats_secondary_findings(data, tier = 3, type = "complex_cmp_het")
  stats.ls$VarN_FS3_HQ_Rare01_Dmg_CX_Uncertain        <- add_stats_secondary_findings(data, tier = 3, type = "complex_single_het")
  stats.ls$VarN_FS3_HQ_Rare01_Dmg_AR_Carrier          <- add_stats_secondary_findings(data, tier = 3, type = "recessive_single_het")
  stats.ls$VarN_FS3_HQ_Rare01_Dmg_XL_Carrier          <- add_stats_secondary_findings(data, tier = 3, type = "Xhet")
  
  stats.ls$VarN_Rare05_ACMG <- add_stats_acmg(data, coding = F)
  stats.ls$VarN_Rare05_ACMG_Coding <- add_stats_acmg(data, coding = T)
  
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

get_all_stats <- function(v_full.df, v_full_hq.df, v_full_r05.df) {
  
  stats.ls.full <- add_stats_full_df(v_full.df)
  stats.ls.hq <- add_stats_full_hq_df(v_full_hq.df)
  stats.ls.r05 <- add_stats_full_r05_df(v_full_r05.df)
  
  stats.df.all <- convert_stats_list_to_df(stats.ls.full, "v_full.df")
  stats.df.all <- rbind(stats.df.all, convert_stats_list_to_df(stats.ls.hq, "v_full_hq.df"))
  stats.df.all <- rbind(stats.df.all, convert_stats_list_to_df(stats.ls.r05, "v_full_r05.df"))
  
  return(stats.df.all)
}

# 3.10. Change tier names to "Low" and "High"

translate_tier_levels <- function(data, col) {
  
  data[, col] <- gsub("1", "Low", as.character(data[, col]))
  data[, col] <- gsub("2", "High", as.character(data[, col]))
  data[, col] <- gsub("0", "-", as.character(data[, col]))
  
  return(data)
}

change_tier_levels <- function(data) {
  
  data <- translate_tier_levels(data, "F_DamageTier")
  data <- translate_tier_levels(data, "F_PhenoTier")
  
  return(data)
}