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
#when defaut the p<1e-5 and clump_data clump_r2 = 0.01
gut.data_7<-gut.data_1[which(gut.data_1$pval.exposure < 1e-5), ]
gut.data_10 <- clump_data(
  gut.data_7,
  clump_kb = 500,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EAS"
)
fwrite(gut.data_10,file = "gut.data_10.csv")
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
fwrite(or_res_guttograves_5,file = "or_res_guttograves_5.csv")
data_guttograves_3_mr_presso<-mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = data_guttograves_3, NbDistribution = 1500,  SignifThreshold = 0.05)
mr_heterogeneity<-mr_heterogeneity(data_guttograves_3)
fwrite(mr_heterogeneity,file = "mr_heterogeneity.csv")
fwrite(data_guttograves_3,file = "data_guttograves_3.csv")
figure_1<-read.csv("figure_1.csv",header = TRUE)
res_figure_1 <- mr(figure_1)
figure_2 <- mr_scatter_plot(res_figure_1, figure_1)
figure_2[1]
class.Deltaproteobacteria.id.3087_presso<-read.csv("class.Deltaproteobacteria.id.3087_presso.csv",header = TRUE)
res_class.Deltaproteobacteria.id.3087_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                          data = class.Deltaproteobacteria.id.3087_presso, 
                                                          NbDistribution = 1000,  SignifThreshold = 0.05)
mr_leaveoneout<-mr_leaveoneout(data_guttograves_3)
leaveoneout_class.Deltaproteobacteria.id.3087<-read.csv("leaveoneout_class.Deltaproteobacteria.id.3087.csv",header = TRUE)
p1_leaveoneout <- mr_leaveoneout_plot(leaveoneout_class.Deltaproteobacteria.id.3087)
p1_leaveoneout[[1]]
#when defaut the p<5e-8 and bbj-a-123 is exposure
iv_graves <- extract_instruments(outcomes='bbj-a-123')
iv_graves_1 <- clump_data(iv_graves,
                          pop = "EAS")
iv_outcome_dat <- extract_outcome_data(
  snps = iv_graves_1$SNP,
  outcomes = gut.data_1
)
write.csv(iv_graves_1,file = "iv_graves_1.csv")
head(gut.data_1)
gut.data_rv<-gut.data_1
gut.data_rv1<-rename(gut.data_rv, id.outcome = id.exposure,
                     outcome = exposure,
                     beta.outcome = beta.exposure,
                     se.outcome = se.exposure,
                     effect_allele.outcome = effect_allele.exposure,
                     other_allele.outcome = other_allele.exposure)
grave_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = gut.data_rv1
)
#for example
phylum.Actinobacteria.id.400.summary<-read.table("phylum.Actinobacteria.id.400.summary.txt",header = TRUE)
Actinobacteria<-phylum.Actinobacteria.id.400.summary
Actinobacteria<-rename(Actinobacteria, 
                       id.outcome = bac,
                       outcome = bac,
                       beta.outcome = beta,
                       se.outcome = SE,
                       effect_allele.outcome = eff.allele,
                       other_allele.outcome = ref.allele,
                       SNP = rsID)
Actinobacteria_1<-Actinobacteria[,1]
Actinobacteria_2<-cbind(Actinobacteria_1,Actinobacteria)
Actinobacteria_2<-rename(Actinobacteria_2,
                         id.outcome = Actinobacteria_1)
write.csv(Actinobacteria_2,file = "Actinobacteria_2.csv")
Actinobacteria_3 <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Actinobacteria_2.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Actinobacteria_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Actinobacteria_3
)
Actinobacteria_res_1 <- mr(Actinobacteria_data_1)