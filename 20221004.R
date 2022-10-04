library("TwoSampleMR")
library("MRPRESSO")
library("MendelianRandomization")
library("data.table")
library("dplyr")
setwd("E:/MR of Graves")
gut.data<-read.table("MBG.allHits.p1e4.txt",header = TRUE)
gut.data_1<-format_data(
  gut.data,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "rsID",
  beta_col = "beta",
  se_col = "SE",
  eaf_col = "eaf",
  effect_allele_col = "eff.allele",
  other_allele_col = "ref.allele",
  pval_col = "P.weightedSumZ",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "N",
  gene_col = "bac",
  id_col = "id",
  min_pval = 1e-200,
  z_col = "Z.weightedSumZ",
  info_col = "info",
  chr_col = "chr",
  pos_col = "bp",
  log_pval = FALSE
)
#when bac is exposure
gut.data_7<-gut.data_1[which(gut.data_1$pval.exposure < 1e-5), ]
gut.data_10 <- clump_data(
  gut.data_7,
  clump_kb = 500,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EAS"
)
gut.data_11<-gut.data_10
gut.data_11$id.exposure<-gut.data_11$gene.exposure
outcome_dat_3 <- extract_outcome_data(
  snps = gut.data_11$SNP,
  outcomes = 'bbj-a-123'
)
data_guttograves_3 <- harmonise_data(
  exposure_dat = gut.data_11, 
  outcome_dat = outcome_dat_3
)
res_guttograves_5 <- mr(data_guttograves_3)
res_guttograves_6<-res_guttograves_5[which(res_guttograves_5$pval < 0.05), ]
or_res_guttograves_5 <- generate_odds_ratios(res_guttograves_6)
data_guttograves_3_mr_presso<-mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = data_guttograves_3, NbDistribution = 1500,  SignifThreshold = 0.05)
mr_heterogeneity<-mr_heterogeneity(data_guttograves_3)
mr_leaveoneout<-mr_leaveoneout(data_guttograves_3)

#when defaut the p<5e-8 and bbj-a-123 is exposure
iv_graves <- extract_instruments(outcomes='bbj-a-123')
iv_graves_1 <- clump_data(iv_graves,
                          pop = "EAS")
iv_outcome_dat <- extract_outcome_data(
  snps = iv_graves_1$SNP,
  outcomes = gut.data_1
)
#We use the family.Oxalobacteraceae.id.2966 as example
#The causal effect between other types of bac with Graves' disease could refer to the following codes
family.Oxalobacteraceae.id.2966.summary<-fread("family.Oxalobacteraceae.id.2966.summary.txt")
family.Oxalobacteraceae.id.2966.summary<-rename(family.Oxalobacteraceae.id.2966.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Oxalobacteraceae.id.2966.summary,file = "family.Oxalobacteraceae.id.2966.summary.csv")
family.Oxalobacteraceae.id.2966.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Oxalobacteraceae.id.2966.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Oxalobacteraceae.id.2966.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Oxalobacteraceae.id.2966.summary
)
family.Oxalobacteraceae.id.2966.summary_res_1 <- mr(family.Oxalobacteraceae.id.2966.summary_data_1)
or_family.Oxalobacteraceae.id.2966.summary <- generate_odds_ratios(family.Oxalobacteraceae.id.2966.summary_res_1)
fwrite(family.Oxalobacteraceae.id.2966.summary_res_1,file = "family.Oxalobacteraceae.id.2966.summary_res_1.csv")
fwrite(family.Oxalobacteraceae.id.2966.summary_data_1,file = "family.Oxalobacteraceae.id.2966.summary_data_1.csv")
p_family.Oxalobacteraceae.id.2966.summary_data_1 <- mr_scatter_plot(family.Oxalobacteraceae.id.2966.summary_res_1, family.Oxalobacteraceae.id.2966.summary_data_1)
p_family.Oxalobacteraceae.id.2966.summary_data_1[[1]]
heterogeneity_family.Oxalobacteraceae.id.2966.summary_data_1<-mr_heterogeneity(family.Oxalobacteraceae.id.2966.summary_data_1)
fwrite(heterogeneity_family.Oxalobacteraceae.id.2966.summary_data_1,file = "heterogeneity_family.Oxalobacteraceae.id.2966.summary_data_1.csv")
family.Oxalobacteraceae.id.2966.summary_data_1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                                   SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                                   data = family.Oxalobacteraceae.id.2966.summary_data_1, 
                                                                   NbDistribution = 1000,  SignifThreshold = 0.05)
leaveoneout_family.Oxalobacteraceae.id.2966.summary_data_1 <- mr_leaveoneout(family.Oxalobacteraceae.id.2966.summary_data_1)
p_family.Oxalobacteraceae.id.2966.summary_data_leaveoneout_plot_1 <- mr_leaveoneout_plot(leaveoneout_family.Oxalobacteraceae.id.2966.summary_data_1)
p_family.Oxalobacteraceae.id.2966.summary_data_leaveoneout_plot_1[[1]]
family.Oxalobacteraceae.id.2966.summary_data_1_presso