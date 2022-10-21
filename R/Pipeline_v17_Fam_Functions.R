source("~/Desktop/PrioritizationPipeline/Rscript/pipeline/pipeline_new_funcs.R")


# (4) FAMILY-BASED FUNCTIONS ----------------------------------------------

get_rare05_variants <- function(data) {
  
  data$F_Rare <- 2
  data <- add_freq_tag(data, freq_max_cutoff = 0.05)
  data <- add_freq_tag(data, freq_max_cutoff = 0.01)
  data <- add_freq_tag(data, freq_max_cutoff = 0.005)
  data <- add_freq_tag(data, freq_max_cutoff = 0.0015)
  data <- add_freq_tag(data, freq_max_cutoff = 0)
  
  data <- subset(data, subset = F_Rare <= 0.05)
  data <- add_pass_tag(data)
  data <- add_qual_tag(data, DP_cutoff)
  data <- add_coding_tag(data)
  
  data$F_DamageType <- "NotDmg"
  data$F_S_DamageType <- "NotLOF"
  data$F_DamageTier <- 0
  
  data <- add_coding_lof_tag(data, eff_lof.chv)
  data <- add_coding_lof_spliceJunction_tag(data, eff_lof.chv)
  
  data <- add_missense_tag(data, sift_cutoff, polyphen_cutoff, ma_cutoff,
                           phylopMam_missense_cutoff, phylopVert_missense_cutoff,
                           CADD_phred_missense_cutoff, REVEL_cutoff, MPC_cutoff, 
                           missense_t1_cutoff, missense_t2_cutoff)
  
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
  
  data <- add_HPO_CGD_dominant_tag(data, database = "HPO", pattern = "@AD")
  data <- add_HPO_CGD_dominant_tag(data, database = "CGD", pattern = "AD")
  data <- add_phenotype_tier_tag(data)
  
  data <- add_recessive_homozygous_tag(data)
  data <- add_X_haploid_tag(data)
  data <- add_potential_compound_heterozygous_tag(data, secondary_findings = F)
  data <- add_potential_dmg_compound_heterozygous_tag(data)
  data <- add_dom_tag(data)
  data <- add_het_hotzone_tag(data, gnomAD_oe_lof_upper_cutoff)
  
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
  data <- add_compound_heterozygous_tag(data, father_column.name, 
                                        mother_column.name, secondary_findings = F)
  data <- add_compound_heterozygous_tag(data, father_column.name, 
                                        mother_column.name, secondary_findings = T)
  
  return(data)
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

which_zygosity <- function(gt) {
  
  gt_sep <- sep_gt(gt)
  
  if (gt_sep[1] == gt_sep[2]) {
    return("hom")
  } else {
    return("het")
  }
}

get_zygosity <- function(child_gt, father_gt, mother_gt) {
  
  child_zg <- which_zygosity(child_gt)
  father_zg <- which_zygosity(father_gt)
  mother_zg <- which_zygosity(mother_gt)
  
  return(c(child_zg, father_zg, mother_zg))
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
  fam_dat$zygosity <- get_zygosity(child_gt, father_gt, mother_gt)
  
  child_gt_sep <- sep_gt(child_gt)
  father_gt_sep <- sep_gt(father_gt)
  mother_gt_sep <- sep_gt(mother_gt)
  
  if (fam_dat$zygosity[1] == "hom" ) {
    
    if (child_gt_sep[1] %in% father_gt_sep & child_gt_sep[1] %in% mother_gt_sep) { # if the copy can be found in both parents
      origins[i] <- "both"
      return(origins[i])
    } 
    
  } else if (setequal(father_gt_sep, mother_gt_sep) & setequal(father_gt_sep, child_gt_sep) & setequal(mother_gt_sep, child_gt_sep)) {
    
    return(origins[i]) # if all family members have the same gt (order doesn't matter), origin stays as "none"
    
  } else if (fam_dat$zygosity[1] == "het") {
    
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
  cmp_het <- cmp_het %>% subset(!grepl("\\./\\.", cmp_het[, father_gt.col]) & 
                                  !grepl("\\./\\.", cmp_het[, mother_gt.col]) &
                                  !grepl("\\./\\.", GT_PreNorm))
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

find_true_cmp_het <- function(cmp_het) {
  
  same_pos <- cmp_het %>% group_by(gene_symbol) %>% filter(n() == 2 & "father" %in% origin & "mother" %in% origin) %>% mutate(unique_pos = n_distinct(start)) %>% filter(unique_pos == 1)
  
  cmp_het <- cmp_het %>% filter(! gene_symbol %in% same_pos$gene_symbol) %>% group_by(gene_symbol) %>% mutate(is_cmp_het = 
                                                                                                                ("father" %in% origin & "mother" %in% origin) | 
                                                                                                                ("both" %in% origin)) %>% arrange(gene_symbol)
  true_cmp_het <- cmp_het %>% subset(is_cmp_het == TRUE)
  
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
                                          secondary_findings = FALSE) {
  
  father_gt.col <- paste0(father_column.name, "GT_PreNorm")
  mother_gt.col <- paste0(mother_column.name, "GT_PreNorm")
  
  cmp_het <- get_cmp_het_subset(data, father_gt.col, mother_gt.col, secondary_findings)
  cmp_het$origin <- get_var_origins(cmp_het, father_gt.col, mother_gt.col)
  true_cmp_het <- find_true_cmp_het(cmp_het)
  
  ifelse(secondary_findings, col_name <- "FS_Fam_CmpHet", col_name <- "FM_Fam_CmpHet")
  data[, col_name] <- initialize_cmp_het_col(data, col_name, secondary_findings)
  data <- merge(data, true_cmp_het[, c("Original_VCFKEY", "gene_symbol", "is_cmp_het")], 
                by = c("Original_VCFKEY", "gene_symbol"), all.x = TRUE)
  data[, col_name][data$is_cmp_het == TRUE] <- "True"
  data <- subset(data, select = -c(is_cmp_het))
  
  return(data)
}
