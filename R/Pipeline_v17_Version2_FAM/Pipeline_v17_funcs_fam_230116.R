# ---
# Title: Family-based Variant Prioritization Functions
# Purpose: Contains functions used in the small variant prioritization family-based 
#          pipeline v17
# Author: Kara Han <kara.han@sickkids.ca>
# Adapted from: Prioritizaiton pipeline developed by Dr. Daniele M., Bhooma T., 
#               Thomas N., and Dr. Worrawat E. at TCAG
# Date script created: 2022-09-22 14:32:21 EDT
# Date last modified: 2023-01-16 13:35:54 EST
# Version: v17
# TCAG annotation pipeline version: rev27.7 hg38
# ---

# source("/Users/karahan/Desktop/PrioritizationPipeline/script/R/pipeline/v230116/Pipeline_v17_funcs_230116.R")
source("/hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/new_script/script/R/pipeline/v230116/Pipeline_v17_funcs_230116.R")

# (4) FAMILY-BASED FUNCTIONS ----------------------------------------------

get_rare05_variants <- function(data, output_prefix.path_filename) {
  
  data$F_Rare <- 2
  data <- add_freq_tag(data, freq_max_cutoff = 0.05)
  data <- add_freq_tag(data, freq_max_cutoff = 0.01)
  data <- add_freq_tag(data, freq_max_cutoff = 0.001)
  data <- add_freq_tag(data, freq_max_cutoff = 0.0001)
  data <- add_freq_tag(data, freq_max_cutoff = 0)
  
  data <- subset(data, subset = F_Rare <= 0.05)
  data <- add_pass_tag(data)
  data <- add_qual_tag_multisample(data, DP_cutoff)
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
  
  # For family-wise prioritization only
  data <- add_compound_heterozygous_tag(data, father_column.name, mother_column.name, 
                                        secondary_findings = F, output_prefix.path_filename)
  data <- add_compound_heterozygous_tag(data, father_column.name, mother_column.name, 
                                        secondary_findings = T, output_prefix.path_filename)
  
  return(data)
}

get_hq_variants <- function(data) {
  
  data <- subset(data, subset = (FT == "PASS"))
  data <- add_qual_tag_multisample(data, DP_cutoff)
  
  return(data)
}

add_stats_true_cmphet <- function(data, secondary_findings = FALSE) {
  
  ifelse(secondary_findings, 
         cond <- expression(FS_Fam_CmpHet == "True"),
         cond <- expression(FM_Fam_CmpHet == "True"))
  
  num_var <- get_num_var(data, cond)
  
  return(num_var)
}

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
  
  stats.ls$Coding_LQ_Rare01_OtherDmg <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "OtherC", high_quality = F)
  stats.ls$Coding_HQ_Rare01_OtherDmg <- add_stats_dmg_tag(data, freq_cutoff = 0.01, dmg_type = "OtherC", high_quality = T)
  
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
  
  stats.ls$VarN_FM_True_CmpHet <- add_stats_true_cmphet(data, secondary_findings = F)
  stats.ls$VarN_FS_True_CmpHet <- add_stats_true_cmphet(data, secondary_findings = T)
  
  return(stats.ls)
}

check_samples_exist <- function(data, child_column.name, father_column.name, mother_column.name) {
  
  col_names <- colnames(data)
  
  if (! ( TRUE %in% grepl(child_column.name, col_names) & TRUE %in% grepl(child_column.name, col_names) & TRUE %in% grepl(child_column.name, col_names) )) {
    stop("At least one of the samples in the family trio is missing from the variant data")
  } else {
    message("All members of the family trio are found in the variant data, proceed to next step")
  }
}

sep_gt <- function(gt) {
  
  gt_sep <- unlist(strsplit(gt, "/|\\Q|\\E"))
  
  return(gt_sep)
}

get_var_origin <- function(child_gt_sep) {
  
  if (sum(child_gt_sep == 0) == 1) {
    origin <- child_gt_sep[child_gt_sep != 0] %>% names()
    return(origin)
  } else if (sum(child_gt_sep == 0) == 0) {
    origin <- "both"
    return(origin)
  }
}

get_each_var_origin <- function(i, origins, fam_dat, CHROM, child_gts, father_gts, mother_gts, father_gt.col, mother_gt.col) {
  
  if (fam_dat$sex[1] == "male" & (CHROM[i] == "chrX" | CHROM[i] == "chrY")) {
    if (child_gts[i] != 0) {
      ifelse(CHROM[i] == "chrX", origins[i] <- "mother", origins[i] <- "father")
      return(origins[i])
    }
    return(origins[i])
  }
  
  child_gt <- child_gts[i]
  father_gt <- father_gts[i]
  mother_gt <- mother_gts[i]
  fam_dat$gt <- c(child_gt, father_gt, mother_gt)
  
  child_gt_sep <- sep_gt(child_gt)
  father_gt_sep <- sep_gt(father_gt)
  mother_gt_sep <- sep_gt(mother_gt)
  
  if (setequal(father_gt_sep, mother_gt_sep) & setequal(father_gt_sep, child_gt_sep) & setequal(mother_gt_sep, child_gt_sep)) {
    
    return(origins[i]) # if all family members have the same gt (order doesn't matter), origin stays as "none"
    
  } else {
    
    cond1 <- child_gt_sep[1] %in% father_gt_sep & child_gt_sep[2] %in% mother_gt_sep
    cond2 <- child_gt_sep[1] %in% mother_gt_sep & child_gt_sep[2] %in% father_gt_sep
    
    if (cond1) { # if the first copy comes from the father and the second the mother
      
      names(child_gt_sep) <- c("father", "mother")
      origins[i] <- get_var_origin(child_gt_sep)
      return(origins[i])
      
    } else if (cond2) { # if the first copy comes from the mother and the second the father
      
      names(child_gt_sep) <- c("mother", "father")
      origins[i] <- get_var_origin(child_gt_sep)
      return(origins[i])
    }
  }
  return(origins[i]) # in all other cases the origin stays as "none"
}

get_cmp_het_subset <- function(data, father_gt.col, mother_gt.col, secondary_findings) {
  
  ifelse(secondary_findings,
         cmp_het <- data %>% subset(F_CmpHet_S1 == 1 | F_CmpHet_S1 == 2),
         cmp_het <- data %>% subset(FM_PCHET == 1 | FM_PCHET == 2))
  
  cmp_het <- unique(cmp_het[, c("CHROM", "start", "Original_VCFKEY", "var_type", 
                                father_gt.col, mother_gt.col, "GT_PreNorm", 
                                "DN", "gene_symbol")])
  cmp_het <- cmp_het %>% subset(!grepl("\\.(/\\.){0,}", cmp_het[, father_gt.col]) & 
                                  !grepl("\\./\\.", cmp_het[, mother_gt.col]) &
                                  !grepl("\\.(/\\.){0,}", GT_PreNorm))
  cmp_het <- cmp_het %>% group_by(gene_symbol) %>% filter(n() > 1) # remove genes with less than 2 variants
  
  return(cmp_het)
}

get_var_origins <- function(cmp_het, father_gt.col, mother_gt.col) {
  
  CHROM <- cmp_het$CHROM
  child_gts <- cmp_het$GT_PreNorm
  father_gts <- cmp_het[, father_gt.col] %>% deframe()
  mother_gts <- cmp_het[, mother_gt.col] %>% deframe()
  
  origins <- rep("none", nrow(cmp_het))
  origins <- sapply(seq_along(1:nrow(cmp_het)), get_each_var_origin, origin = origins, 
                    fam_dat = fam_dat, CHROM = CHROM, child_gts = child_gts, 
                    father_gts = father_gts, mother_gts = mother_gts, 
                    father_gt.col = father_gt.col, mother_gt.col = mother_gt.col)
  
  return(origins)
}

find_true_cmp_het <- function(cmp_het, output_prefix, secondary_findings) {
  
  same_pos <- cmp_het %>% group_by(gene_symbol) %>% filter(n() == 2 & "father" %in% origin & "mother" %in% origin) %>% mutate(unique_pos = n_distinct(start)) %>% filter(unique_pos == 1)
  
  cmp_het <- cmp_het %>% filter(! gene_symbol %in% same_pos$gene_symbol) %>% group_by(gene_symbol) %>% mutate(is_cmp_het = 
                                                                                                                ("father" %in% origin & "mother" %in% origin) | 
                                                                                                                ("both" %in% origin)) %>% arrange(gene_symbol)
  true_cmp_het <- cmp_het %>% subset(is_cmp_het == TRUE)
  
  if (secondary_findings) {
    write.table(true_cmp_het, col.names = T, row.names = F, quote = F, sep = "\t", 
                file = paste(output_prefix, "_True_CmpHet_SecondaryFindings", ".txt", sep = ""))
  } else {
    write.table(true_cmp_het, col.names = T, row.names = F, quote = F, sep = "\t", 
                file = paste(output_prefix, "_True_CmpHet_MainFindings", ".txt", sep = ""))
  }
  
  return(true_cmp_het)
}

initialize_cmp_het_col <- function(data, col_name, secondary_findings) {
  
  ifelse(secondary_findings,
         cond <- expression(data$F_CmpHet_S1 == 1 | data$F_CmpHet_S1 == 2),
         cond <- expression(data$FM_PCHET == 1 | data$FM_PCHET == 2))
  
  vals <- data.table::fifelse(eval(cond), 
                              data[, col_name] <- "Potential",
                              fifelse(data$DN == "LowDQ" | data$DN == "DeNovo", 
                                      data[, col_name] <- data$DN,
                                      data[, col_name] <- "-"))
  return(vals)
}

add_compound_heterozygous_tag <- function(data, father_column.name, mother_column.name,
                                          secondary_findings = FALSE, output_prefix) {
  
  father_gt.col <- paste0(father_column.name, "GT_PreNorm")
  mother_gt.col <- paste0(mother_column.name, "GT_PreNorm")
  ifelse(secondary_findings, col_name <- "FS_Fam_CmpHet", col_name <- "FM_Fam_CmpHet")
  
  cmp_het <- get_cmp_het_subset(data, father_gt.col, mother_gt.col, secondary_findings)
  if (nrow(cmp_het) == 0) {
    data[, col_name] <- "-"
    return(data)
  }
  
  cmp_het$origin <- get_var_origins(cmp_het, father_gt.col, mother_gt.col)
  true_cmp_het <- find_true_cmp_het(cmp_het, output_prefix, secondary_findings)
  data[, col_name] <- initialize_cmp_het_col(data, col_name, secondary_findings)
  data <- merge(data, true_cmp_het[, c("Original_VCFKEY", "gene_symbol", "is_cmp_het")], 
                by = c("Original_VCFKEY", "gene_symbol"), all.x = TRUE)
  data[, col_name][data$is_cmp_het == TRUE] <- "True"
  data <- subset(data, select = -c(is_cmp_het))
  
  return(data)
}
