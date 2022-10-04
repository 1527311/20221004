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
#when defaut the p<5e-8
gut.data_2<-gut.data_1[which(gut.data_1$pval.exposure < 5e-8), ]
gut.data_3 <- clump_data(gut.data_2,
                         pop = "EAS")
outcome_dat <- extract_outcome_data(
  snps = gut.data_3$SNP,
  outcomes = 'bbj-a-123'
)
data_guttograves <- harmonise_data(
  exposure_dat = gut.data_3, 
  outcome_dat = outcome_dat
)
res_guttograves <- mr(data_guttograves)
res_single_guttograves <- mr_singlesnp(data_guttograves, all_method="mr_two_sample_ml")
#when defaut the p<1e-5
gut.data_4<-gut.data_1[which(gut.data_1$pval.exposure < 1e-5), ]
gut.data_5 <- clump_data(gut.data_4,
                         pop = "EAS")
gut.data_6<-gut.data_5
gut.data_6$id.exposure<-gut.data_6$gene.exposure
outcome_dat_1 <- extract_outcome_data(
  snps = gut.data_6$SNP,
  outcomes = 'bbj-a-123'
)
data_guttograves_1 <- harmonise_data(
  exposure_dat = gut.data_6, 
  outcome_dat = outcome_dat_1
)
res_guttograves_1 <- mr(data_guttograves_1)
res_guttograves_2<-res_guttograves_1[which(res_guttograves_1$pval < 0.05), ]
or_res_guttograves_2 <- generate_odds_ratios(res_guttograves_2)
#when defaut the p<1e-5 and clump_data clump_r2 = 0.1
gut.data_7<-gut.data_1[which(gut.data_1$pval.exposure < 1e-5), ]
gut.data_8 <- clump_data(
  gut.data_7,
  clump_kb = 500,
  clump_r2 = 0.1,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EAS"
)
gut.data_9<-gut.data_8
gut.data_9$id.exposure<-gut.data_9$gene.exposure
outcome_dat_2 <- extract_outcome_data(
  snps = gut.data_9$SNP,
  outcomes = 'bbj-a-123'
)
data_guttograves_2 <- harmonise_data(
  exposure_dat = gut.data_9, 
  outcome_dat = outcome_dat_2
)
res_guttograves_3 <- mr(data_guttograves_2)
res_guttograves_4<-res_guttograves_3[which(res_guttograves_3$pval < 0.05), ]
or_res_guttograves_4 <- generate_odds_ratios(res_guttograves_4)
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
??mrpresso
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
??mr_scatter_plot
figure_2[8]
class.Deltaproteobacteria.id.3087_presso<-read.csv("class.Deltaproteobacteria.id.3087_presso.csv",header = TRUE)
res_class.Deltaproteobacteria.id.3087_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                           SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                           OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                           data = class.Deltaproteobacteria.id.3087_presso, 
                           NbDistribution = 1000,  SignifThreshold = 0.05)
res_class.Deltaproteobacteria.id.3087_presso
class.Mollicutes.id.3920_presso<-read.csv("class.Mollicutes.id.3920_presso.csv",header = TRUE)
res_class.Mollicutes.id.3920_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                          data = class.Mollicutes.id.3920_presso, 
                                                          NbDistribution = 1000,  SignifThreshold = 0.05)
family.Peptococcaceae.id.2024_presso<-read.csv("family.Peptococcaceae.id.2024_presso.csv",header = TRUE)
res_family.Peptococcaceae.id.2024_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                          data = family.Peptococcaceae.id.2024_presso, 
                                                          NbDistribution = 1000,  SignifThreshold = 0.05)
family.unknownfamily.id.1000006161_presso<-read.csv("family.unknownfamily.id.1000006161_presso.csv",header = TRUE)
res_family.unknownfamily.id.1000006161_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                          data = family.unknownfamily.id.1000006161_presso, 
                                                          NbDistribution = 1000,  SignifThreshold = 0.05)
genus..Ruminococcustorquesgroup.id.14377_presso<-read.csv("genus..Ruminococcustorquesgroup.id.14377_presso.csv",header = TRUE)
res_genus..Ruminococcustorquesgroup.id.14377_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                          data = genus..Ruminococcustorquesgroup.id.14377_presso, 
                                                          NbDistribution = 1000,  SignifThreshold = 0.05)
genus.Anaerostipes.id.1991_presso<-read.csv("genus.Anaerostipes.id.1991_presso.csv",header = TRUE)
res_genus.Anaerostipes.id.1991_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                          data = genus.Anaerostipes.id.1991_presso, 
                                                          NbDistribution = 1000,  SignifThreshold = 0.05)
genus.Oxalobacter.id.2978_presso<-read.csv("genus.Oxalobacter.id.2978_presso.csv",header = TRUE)
res_genus.Oxalobacter.id.2978_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                          data = genus.Oxalobacter.id.2978_presso, 
                                                          NbDistribution = 1000,  SignifThreshold = 0.05)
genus.RuminococcaceaeUCG011.id.11368_presso<-read.csv("genus.RuminococcaceaeUCG011.id.11368_presso.csv",header = TRUE)
res_genus.RuminococcaceaeUCG011.id.11368_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                          data = genus.RuminococcaceaeUCG011.id.11368_presso, 
                                                          NbDistribution = 1000,  SignifThreshold = 0.05)
mr_leaveoneout<-mr_leaveoneout(data_guttograves_3)
fwrite(mr_leaveoneout,file = "mr_leaveoneout.csv")
leaveoneout_class.Deltaproteobacteria.id.3087<-read.csv("leaveoneout_class.Deltaproteobacteria.id.3087.csv",header = TRUE)
p1_leaveoneout <- mr_leaveoneout_plot(leaveoneout_class.Deltaproteobacteria.id.3087)
p1_leaveoneout[[1]]
leaveoneout_class.Mollicutes.id.3920<-read.csv("leaveoneout_class.Mollicutes.id.3920.csv",header = TRUE)
p2_leaveoneout <- mr_leaveoneout_plot(leaveoneout_class.Mollicutes.id.3920)
p2_leaveoneout[[1]]
leaveoneout_family.Peptococcaceae.id.2024<-read.csv("leaveoneout_family.Peptococcaceae.id.2024.csv",header = TRUE)
p3_leaveoneout <- mr_leaveoneout_plot(leaveoneout_family.Peptococcaceae.id.2024)
p3_leaveoneout[[1]]
leaveoneout_family.unknownfamily.id.1000006161<-read.csv("leaveoneout_family.unknownfamily.id.1000006161.csv",header = TRUE)
p4_leaveoneout <- mr_leaveoneout_plot(leaveoneout_family.unknownfamily.id.1000006161)
p4_leaveoneout[[1]]
leaveoneout_genus..Ruminococcustorquesgroup.id.14377<-read.csv("leaveoneout_genus..Ruminococcustorquesgroup.id.14377.csv",header = TRUE)
p5_leaveoneout <- mr_leaveoneout_plot(leaveoneout_genus..Ruminococcustorquesgroup.id.14377)
p5_leaveoneout[[1]]
leaveoneout_genus.Anaerostipes.id.1991<-read.csv("leaveoneout_genus.Anaerostipes.id.1991.csv",header = TRUE)
p6_leaveoneout <- mr_leaveoneout_plot(leaveoneout_genus.Anaerostipes.id.1991)
p6_leaveoneout[[1]]
leaveoneout_genus.RuminococcaceaeUCG011.id.11368<-read.csv("leaveoneout_genus.RuminococcaceaeUCG011.id.11368.csv",header = TRUE)
p7_leaveoneout <- mr_leaveoneout_plot(leaveoneout_genus.RuminococcaceaeUCG011.id.11368)
p7_leaveoneout[[1]]
#when defaut the p<5e-8 and bbj-a-123 is exposure
iv_graves <- extract_instruments(outcomes='bbj-a-123')
iv_graves_1 <- clump_data(iv_graves,
                         pop = "EAS")
iv_outcome_dat <- extract_outcome_data(
  snps = iv_graves_1$SNP,
  outcomes = gut.data_1
)
write.csv(iv_graves_1,file = "iv_graves_1.csv")
#The following required columns are missing from outcome: 
#id.outcome, outcome, beta.outcome, se.outcome, effect_allele.outcome, other_allele.outcome
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
#Actinobacteria
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
#Firmicutesphylum.Firmicutes.id.1672.summary.txt<-read.table("phylum.Firmicutes.id.1672.summary.txt",header = TRUE)
Firmicutes<-phylum.Firmicutes.id.1672.summary.txt
Firmicutes<-rename(Firmicutes, 
                       id.outcome = bac,
                       outcome = bac,
                       beta.outcome = beta,
                       se.outcome = SE,
                       effect_allele.outcome = eff.allele,
                       other_allele.outcome = ref.allele,
                       SNP = rsID)
Firmicutes_1<-Firmicutes[,1]
Firmicutes_2<-cbind(Firmicutes_1,Firmicutes)
Firmicutes_2<-rename(Firmicutes_2,
                         id.outcome = Firmicutes_1)
write.csv(Firmicutes_2,file = "Firmicutes_2.csv")
Firmicutes_3 <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Firmicutes_2.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Firmicutes_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Firmicutes_3
)
Firmicutes_res_1 <- mr(Firmicutes_data_1)
#phylum.Bacteroidetes.id.905.summary.txt
Bacteroidetes<-read.table("phylum.Bacteroidetes.id.905.summary.txt",header = TRUE)
Bacteroidetes<-rename(Bacteroidetes, 
                   id.outcome = bac,
                   outcome = bac,
                   beta.outcome = beta,
                   se.outcome = SE,
                   effect_allele.outcome = eff.allele,
                   other_allele.outcome = ref.allele,
                   SNP = rsID)
Bacteroidetes_1<-Bacteroidetes[,1]
Bacteroidetes_2<-cbind(Bacteroidetes_1,Bacteroidetes)
Bacteroidetes_2<-rename(Bacteroidetes_2,
                     id.outcome = Bacteroidetes_1)
write.csv(Bacteroidetes_2,file = "Bacteroidetes_2.csv")
Bacteroidetes_3 <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Bacteroidetes_2.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Bacteroidetes_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Bacteroidetes_3
)
Bacteroidetes_res_1 <- mr(Bacteroidetes_data_1)
remove(Bacteroidetes)
remove(Bacteroidetes_1)
remove(Bacteroidetes_2)
#phylum.Cyanobacteria.id.1500.summary.txt
Cyanobacteria<-read.table("phylum.Cyanobacteria.id.1500.summary.txt",header = TRUE)
Cyanobacteria<-rename(Cyanobacteria, 
                      id.outcome = bac,
                      outcome = bac,
                      beta.outcome = beta,
                      se.outcome = SE,
                      effect_allele.outcome = eff.allele,
                      other_allele.outcome = ref.allele,
                      SNP = rsID)
Cyanobacteria_1<-Cyanobacteria[,1]
Cyanobacteria_2<-cbind(Cyanobacteria_1,Cyanobacteria)
Cyanobacteria_2<-rename(Cyanobacteria_2,
                        id.outcome = Cyanobacteria_1)
write.csv(Cyanobacteria_2,file = "Cyanobacteria_2.csv")
Cyanobacteria_3 <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Cyanobacteria_2.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Cyanobacteria_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Cyanobacteria_3
)
Cyanobacteria_res_1 <- mr(Cyanobacteria_data_1)
remove(Cyanobacteria)
remove(Cyanobacteria_1)
remove(Cyanobacteria_2)
#phylum.Euryarchaeota.id.55.summary.txt
Euryarchaeota<-read.table("phylum.Euryarchaeota.id.55.summary.txt",header = TRUE)
Euryarchaeota<-rename(Euryarchaeota, 
                      id.outcome = bac,
                      outcome = bac,
                      beta.outcome = beta,
                      se.outcome = SE,
                      effect_allele.outcome = eff.allele,
                      other_allele.outcome = ref.allele,
                      SNP = rsID)
Euryarchaeota_1<-Euryarchaeota[,1]
Euryarchaeota_2<-cbind(Euryarchaeota_1,Euryarchaeota)
Euryarchaeota_2<-rename(Euryarchaeota_2,
                        id.outcome = Euryarchaeota_1)
write.csv(Euryarchaeota_2,file = "Euryarchaeota_2.csv")
Euryarchaeota_3 <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Euryarchaeota_2.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Euryarchaeota_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Euryarchaeota_3
)
Euryarchaeota_res_1 <- mr(Euryarchaeota_data_1)
remove(Euryarchaeota)
remove(Euryarchaeota_1)
remove(Euryarchaeota_2)
#phylum.Lentisphaerae.id.2238.summary.txt
Lentisphaerae<-read.table("phylum.Lentisphaerae.id.2238.summary.txt",header = TRUE)
Lentisphaerae<-rename(Lentisphaerae, 
                      id.outcome = bac,
                      outcome = bac,
                      beta.outcome = beta,
                      se.outcome = SE,
                      effect_allele.outcome = eff.allele,
                      other_allele.outcome = ref.allele,
                      SNP = rsID)
Lentisphaerae_1<-Lentisphaerae[,1]
Lentisphaerae_2<-cbind(Lentisphaerae_1,Lentisphaerae)
Lentisphaerae_2<-rename(Lentisphaerae_2,
                        id.outcome = Lentisphaerae_1)
write.csv(Lentisphaerae_2,file = "Lentisphaerae_2.csv")
Lentisphaerae_3 <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Lentisphaerae_2.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Lentisphaerae_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Lentisphaerae_3
)
Lentisphaerae_res_1 <- mr(Lentisphaerae_data_1)
remove(Lentisphaerae)
remove(Lentisphaerae_1)
remove(Lentisphaerae_2)
#phylum.Proteobacteria.id.2375.summary.txt
Proteobacteria<-read.table("phylum.Proteobacteria.id.2375.summary.txt",header = TRUE)
Proteobacteria<-rename(Proteobacteria, 
                      id.outcome = bac,
                      outcome = bac,
                      beta.outcome = beta,
                      se.outcome = SE,
                      effect_allele.outcome = eff.allele,
                      other_allele.outcome = ref.allele,
                      SNP = rsID)
Proteobacteria_1<-Proteobacteria[,1]
Proteobacteria_2<-cbind(Proteobacteria_1,Proteobacteria)
Proteobacteria_2<-rename(Proteobacteria_2,
                        id.outcome = Proteobacteria_1)
write.csv(Proteobacteria_2,file = "Proteobacteria_2.csv")
Proteobacteria_3 <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Proteobacteria_2.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Proteobacteria_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Proteobacteria_3
)
Proteobacteria_res_1 <- mr(Proteobacteria_data_1)
remove(Proteobacteria)
remove(Proteobacteria_1)
remove(Proteobacteria_2)
#phylum.Tenericutes.id.3919.summary.txt
Tenericutes<-read.table("phylum.Tenericutes.id.3919.summary.txt",header = TRUE)
Tenericutes<-rename(Tenericutes, 
                       id.outcome = bac,
                       outcome = bac,
                       beta.outcome = beta,
                       se.outcome = SE,
                       effect_allele.outcome = eff.allele,
                       other_allele.outcome = ref.allele,
                       SNP = rsID)
Tenericutes_1<-Tenericutes[,1]
Tenericutes_2<-cbind(Tenericutes_1,Tenericutes)
Tenericutes_2<-rename(Tenericutes_2,
                         id.outcome = Tenericutes_1)
write.csv(Tenericutes_2,file = "Tenericutes_2.csv")
Tenericutes_3 <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Tenericutes_2.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Tenericutes_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Tenericutes_3
)
Tenericutes_res_1 <- mr(Tenericutes_data_1)
remove(Tenericutes)
remove(Tenericutes_1)
remove(Tenericutes_2)
#phylum.Verrucomicrobia.id.3982.summary.txt
Verrucomicrobia<-read.table("phylum.Verrucomicrobia.id.3982.summary.txt",header = TRUE)
Verrucomicrobia<-rename(Verrucomicrobia, 
                    id.outcome = bac,
                    outcome = bac,
                    beta.outcome = beta,
                    se.outcome = SE,
                    effect_allele.outcome = eff.allele,
                    other_allele.outcome = ref.allele,
                    SNP = rsID)
Verrucomicrobia_1<-Verrucomicrobia[,1]
Verrucomicrobia_2<-cbind(Verrucomicrobia_1,Verrucomicrobia)
Verrucomicrobia_2<-rename(Verrucomicrobia_2,
                      id.outcome = Verrucomicrobia_1)
write.csv(Verrucomicrobia_2,file = "Verrucomicrobia_2.csv")
Verrucomicrobia_3 <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Verrucomicrobia_2.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Verrucomicrobia_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Verrucomicrobia_3
)
Verrucomicrobia_res_1 <- mr(Verrucomicrobia_data_1)
remove(Verrucomicrobia)
remove(Verrucomicrobia_1)
remove(Verrucomicrobia_2)
??read_outcome_data
setwd("G:/gut/MiBioGen_QmbQTL_summary_class/MiBioGen_QmbQTL_summary_class")
#class.Actinobacteria.id.419.summary.txt
remove(Alphaproteobacteria_2)
#class.Bacilli.id.1673.summary.txt
setwd("G:/gut/MiBioGen_QmbQTL_summary_class/MiBioGen_QmbQTL_summary_class")
Bacilli<-read.table("class.Bacilli.id.1673.summary.txt",header = TRUE)
Bacilli<-rename(Bacilli, 
                        id.outcome = bac,
                        outcome = bac,
                        beta.outcome = beta,
                        se.outcome = SE,
                        effect_allele.outcome = eff.allele,
                        other_allele.outcome = ref.allele,
                        SNP = rsID)
Bacilli_1<-Bacilli[,1]
Bacilli<-cbind(Bacilli_1,Bacilli)
Bacilli<-rename(Bacilli,
                          id.outcome = Bacilli_1)
write.csv(Bacilli,file = "Bacilli.csv")
Bacilli <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Bacilli.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Bacilli_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Bacilli
)
Bacilli_res_1 <- mr(Bacilli_data_1)
remove(Bacilli_1)
#class.Betaproteobacteria.id.912.summary.txt
Betaproteobacteria<-read.table("class.Betaproteobacteria.id.912.summary.txt",header = TRUE)
Betaproteobacteria<-rename(Betaproteobacteria, 
                id.outcome = bac,
                outcome = bac,
                beta.outcome = beta,
                se.outcome = SE,
                effect_allele.outcome = eff.allele,
                other_allele.outcome = ref.allele,
                SNP = rsID)
Betaproteobacteria_1<-Betaproteobacteria[,1]
Betaproteobacteria<-cbind(Betaproteobacteria_1,Betaproteobacteria)
Betaproteobacteria<-rename(Betaproteobacteria,
                id.outcome = Bacteroidia_1)
write.csv(Bacteroidia,file = "Bacteroidia.csv")
Bacteroidia <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Bacteroidia.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Bacteroidia_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Bacteroidia
)
Bacteroidia_res_1 <- mr(Bacteroidia_data_1)
remove(Bacteroidia_1)
#class.Betaproteobacteria.id.2867.summary.txt
Betaproteobacteria<-read.table("class.Betaproteobacteria.id.2867.summary.txt",header = TRUE)
Betaproteobacteria<-rename(Betaproteobacteria, 
                    id.outcome = bac,
                    outcome = bac,
                    beta.outcome = beta,
                    se.outcome = SE,
                    effect_allele.outcome = eff.allele,
                    other_allele.outcome = ref.allele,
                    SNP = rsID)
Betaproteobacteria_1<-Betaproteobacteria[,1]
Betaproteobacteria<-cbind(Betaproteobacteria_1,Betaproteobacteria)
Betaproteobacteria<-rename(Betaproteobacteria,
                    id.outcome = Betaproteobacteria_1)
write.csv(Betaproteobacteria,file = "Betaproteobacteria.csv")
Betaproteobacteria <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Betaproteobacteria.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Betaproteobacteria_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Betaproteobacteria
)
Betaproteobacteria_res_1 <- mr(Betaproteobacteria_data_1)
remove(Betaproteobacteria_1)
#class.Clostridia.id.1859.summary.txt
Clostridia<-read.table("class.Clostridia.id.1859.summary.txt",header = TRUE)
Clostridia<-rename(Clostridia, 
                           id.outcome = bac,
                           outcome = bac,
                           beta.outcome = beta,
                           se.outcome = SE,
                           effect_allele.outcome = eff.allele,
                           other_allele.outcome = ref.allele,
                           SNP = rsID)
Clostridia_1<-Clostridia[,1]
Clostridia<-cbind(Clostridia_1,Clostridia)
Clostridia<-rename(Clostridia,
                           id.outcome = Clostridia_1)
write.csv(Clostridia,file = "Clostridia.csv")
Clostridia <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Clostridia.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Clostridia_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Clostridia
)
Clostridia_res_1 <- mr(Clostridia_data_1)
remove(Clostridia_1)
setwd("G:/gut/MiBioGen_QmbQTL_summary_class/MiBioGen_QmbQTL_summary_class")
#class.Coriobacteriia.id.809.summary.txt
Coriobacteriia<-fread("class.Coriobacteriia.id.809.summary.txt",header = TRUE)
Coriobacteriia<-rename(Coriobacteriia, 
                           id.outcome = bac,
                           outcome = bac,
                           beta.outcome = beta,
                           se.outcome = SE,
                           effect_allele.outcome = eff.allele,
                           other_allele.outcome = ref.allele,
                           SNP = rsID)
Coriobacteriia_1<-Coriobacteriia[,1]
Coriobacteriia_1<-rename(Coriobacteriia_1,
                           id.outcome = outcome)
Coriobacteriia<-cbind(Coriobacteriia_1,Coriobacteriia)
fwrite(Coriobacteriia,file = "Coriobacteriia.csv")
Coriobacteriia <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Coriobacteriia.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Coriobacteriia_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Coriobacteriia
)
Coriobacteriia_res_1 <- mr(Coriobacteriia_data_1)
remove(Coriobacteriia_1)
#class.Deltaproteobacteria.id.3087.summary.txt
Deltaproteobacteria<-fread("class.Deltaproteobacteria.id.3087.summary.txt",header = TRUE)
Deltaproteobacteria<-rename(Deltaproteobacteria, 
                           id.outcome = bac,
                           outcome = bac,
                           beta.outcome = beta,
                           se.outcome = SE,
                           effect_allele.outcome = eff.allele,
                           other_allele.outcome = ref.allele,
                           SNP = rsID)
Deltaproteobacteria_1<-Deltaproteobacteria[,1]
Deltaproteobacteria_1<-rename(Deltaproteobacteria_1,
                           id.outcome = bac)
Deltaproteobacteria<-cbind(Deltaproteobacteria_1,Deltaproteobacteria)
??fwrite
fwrite(Deltaproteobacteria,file = "Deltaproteobacteria.csv")
Coriobacteriia <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Coriobacteriia.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Deltaproteobacteria_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Deltaproteobacteria
)
Deltaproteobacteria_res_1 <- mr(Deltaproteobacteria_data_1)
remove(Deltaproteobacteria)
#class.Erysipelotrichia.id.2147.summary.txt
Erysipelotrichia<-fread("class.Erysipelotrichia.id.2147.summary.txt")
Erysipelotrichia<-rename(Erysipelotrichia, 
                   id.outcome = bac,
                   beta.outcome = beta,
                   se.outcome = SE,
                   effect_allele.outcome = eff.allele,
                   other_allele.outcome = ref.allele,
                   SNP = rsID)
fwrite(Erysipelotrichia,file = "Erysipelotrichia.csv")
Erysipelotrichia <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Erysipelotrichia.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Erysipelotrichia_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Erysipelotrichia
)
Erysipelotrichia_res_1 <- mr(Erysipelotrichia_data_1)
remove(Erysipelotrichia_1)
#class.Gammaproteobacteria.id.3303.summary.txt
Gammaproteobacteria<-fread("class.Gammaproteobacteria.id.3303.summary.txt")
Gammaproteobacteria<-rename(Gammaproteobacteria, 
                         id.outcome = bac,
                         beta.outcome = beta,
                         se.outcome = SE,
                         effect_allele.outcome = eff.allele,
                         other_allele.outcome = ref.allele,
                         SNP = rsID)
fwrite(Gammaproteobacteria,file = "Gammaproteobacteria.csv")
Gammaproteobacteria <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "Gammaproteobacteria.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
Gammaproteobacteria_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = Gammaproteobacteria
)
Gammaproteobacteria_res_1 <- mr(Gammaproteobacteria_data_1)
#class.Lentisphaeria.id.2250.summary.txt
class.Lentisphaeria.id.2250.summary<-fread("class.Lentisphaeria.id.2250.summary.txt")
class.Lentisphaeria.id.2250.summary<-rename(class.Lentisphaeria.id.2250.summary, 
                            id.outcome = bac,
                            beta.outcome = beta,
                            se.outcome = SE,
                            effect_allele.outcome = eff.allele,
                            other_allele.outcome = ref.allele,
                            SNP = rsID)
fwrite(class.Lentisphaeria.id.2250.summary,file = "class.Lentisphaeria.id.2250.summary.csv")
class.Lentisphaeria.id.2250.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "class.Lentisphaeria.id.2250.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
class.Lentisphaeria.id.2250.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = class.Lentisphaeria.id.2250.summary
)
class.Lentisphaeria.id.2250.summary_res_1 <- mr(class.Lentisphaeria.id.2250.summary_data_1)
#class.Melainabacteria.id.1589.summary.txt
class.Melainabacteria.id.1589.summary<-fread("class.Melainabacteria.id.1589.summary.txt")
class.Melainabacteria.id.1589.summary<-rename(class.Melainabacteria.id.1589.summary, 
                                            id.outcome = bac,
                                            beta.outcome = beta,
                                            se.outcome = SE,
                                            effect_allele.outcome = eff.allele,
                                            other_allele.outcome = ref.allele,
                                            SNP = rsID)
fwrite(class.Melainabacteria.id.1589.summary,file = "class.Melainabacteria.id.1589.summary.csv")
class.Melainabacteria.id.1589.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "class.Melainabacteria.id.1589.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
class.Melainabacteria.id.1589.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = class.Melainabacteria.id.1589.summary
)
class.Melainabacteria.id.1589.summary_res_1 <- mr(class.Melainabacteria.id.1589.summary_data_1)
#class.Methanobacteria.id.119.summary.txt
class.Methanobacteria.id.119.summary<-fread("class.Methanobacteria.id.119.summary.txt")
class.Methanobacteria.id.119.summary<-rename(class.Methanobacteria.id.119.summary, 
                                            id.outcome = bac,
                                            beta.outcome = beta,
                                            se.outcome = SE,
                                            effect_allele.outcome = eff.allele,
                                            other_allele.outcome = ref.allele,
                                            SNP = rsID)
fwrite(class.Methanobacteria.id.119.summary,file = "class.Methanobacteria.id.119.summary.csv")
class.Methanobacteria.id.119.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "class.Methanobacteria.id.119.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
class.Methanobacteria.id.119.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = class.Methanobacteria.id.119.summary
)
class.Methanobacteria.id.119.summary_res_1 <- mr(class.Methanobacteria.id.119.summary_data_1)
#class.Mollicutes.id.3920.summary.txt
class.Mollicutes.id.3920.summary<-fread("class.Mollicutes.id.3920.summary.txt")
class.Mollicutes.id.3920.summary<-rename(class.Mollicutes.id.3920.summary, 
                                            id.outcome = bac,
                                            beta.outcome = beta,
                                            se.outcome = SE,
                                            effect_allele.outcome = eff.allele,
                                            other_allele.outcome = ref.allele,
                                            SNP = rsID)
fwrite(class.Mollicutes.id.3920.summary,file = "class.Mollicutes.id.3920.summary.csv")
class.Mollicutes.id.3920.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "class.Mollicutes.id.3920.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
class.Mollicutes.id.3920.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = class.Mollicutes.id.3920.summary
)
class.Mollicutes.id.3920.summary_res_1 <- mr(class.Mollicutes.id.3920.summary_data_1)
#class.Negativicutes.id.2164.summary.txt
class.Negativicutes.id.2164.summary<-fread("class.Negativicutes.id.2164.summary.txt")
class.Negativicutes.id.2164.summary<-rename(class.Negativicutes.id.2164.summary, 
                                            id.outcome = bac,
                                            beta.outcome = beta,
                                            se.outcome = SE,
                                            effect_allele.outcome = eff.allele,
                                            other_allele.outcome = ref.allele,
                                            SNP = rsID)
fwrite(class.Negativicutes.id.2164.summary,file = "class.Negativicutes.id.2164.summary.csv")
class.Negativicutes.id.2164.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "class.Negativicutes.id.2164.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
class.Negativicutes.id.2164.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = class.Negativicutes.id.2164.summary
)
class.Negativicutes.id.2164.summary_res_1 <- mr(class.Negativicutes.id.2164.summary_data_1)
#class.Verrucomicrobiae.id.4029.summary.txt
class.Verrucomicrobiae.id.4029.summary<-fread("class.Verrucomicrobiae.id.4029.summary.txt")
class.Verrucomicrobiae.id.4029.summary<-rename(class.Verrucomicrobiae.id.4029.summary, 
                                            id.outcome = bac,
                                            beta.outcome = beta,
                                            se.outcome = SE,
                                            effect_allele.outcome = eff.allele,
                                            other_allele.outcome = ref.allele,
                                            SNP = rsID)
fwrite(class.Verrucomicrobiae.id.4029.summary,file = "class.Verrucomicrobiae.id.4029.summary.csv")
class.Verrucomicrobiae.id.4029.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "class.Verrucomicrobiae.id.4029.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
class.Verrucomicrobiae.id.4029.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = class.Verrucomicrobiae.id.4029.summary
)
class.Verrucomicrobiae.id.4029.summary_res_1 <- mr(class.Verrucomicrobiae.id.4029.summary_data_1)
#setwd("G:/gut/MiBioGen_QmbQTL_summary_order/order")
setwd("G:/gut/MiBioGen_QmbQTL_summary_order/order")
#order.Actinomycetales.id.420.summary.txt
order.Actinomycetales.id.420.summary<-fread("order.Actinomycetales.id.420.summary.txt")
order.Actinomycetales.id.420.summary<-rename(order.Actinomycetales.id.420.summary, 
                                               id.outcome = bac,
                                               beta.outcome = beta,
                                               se.outcome = SE,
                                               effect_allele.outcome = eff.allele,
                                               other_allele.outcome = ref.allele,
                                               SNP = rsID)
fwrite(order.Actinomycetales.id.420.summary,file = "order.Actinomycetales.id.420.summary.csv")
order.Actinomycetales.id.420.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Actinomycetales.id.420.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Actinomycetales.id.420.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Actinomycetales.id.420.summary
)
order.Actinomycetales.id.420.summary_res_1 <- mr(order.Actinomycetales.id.420.summary_data_1)
fwrite(order.Actinomycetales.id.420.summary_res_1,file = "order.Actinomycetales.id.420.summary_res_1.csv")
fwrite(order.Actinomycetales.id.420.summary_data_1,file = "order.Actinomycetales.id.420.summary_data_1.csv")
#order.Bacillales.id.1674.summary.txt
order.Bacillales.id.1674.summary<-fread("order.Bacillales.id.1674.summary.txt")
order.Bacillales.id.1674.summary<-rename(order.Bacillales.id.1674.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Bacillales.id.1674.summary,file = "order.Bacillales.id.1674.summary.csv")
order.Bacillales.id.1674.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Bacillales.id.1674.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Bacillales.id.1674.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Bacillales.id.1674.summary
)
order.Bacillales.id.1674.summary_res_1 <- mr(order.Bacillales.id.1674.summary_data_1)
fwrite(order.Bacillales.id.1674.summary_res_1,file = "order.Bacillales.id.1674.summary_res_1.csv")
fwrite(order.Bacillales.id.1674.summary_data_1,file = "order.Bacillales.id.1674.summary_data_1.csv")
#order.Bacteroidales.id.913.summary.txt
order.Bacteroidales.id.913.summary<-fread("order.Bacteroidales.id.913.summary.txt")
order.Bacteroidales.id.913.summary<-rename(order.Bacteroidales.id.913.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Bacteroidales.id.913.summary,file = "order.Bacteroidales.id.913.summary.csv")
order.Bacteroidales.id.913.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Bacteroidales.id.913.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Bacteroidales.id.913.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Bacteroidales.id.913.summary
)
order.Bacteroidales.id.913.summary_res_1 <- mr(order.Bacteroidales.id.913.summary_data_1)
fwrite(order.Bacteroidales.id.913.summary_res_1,file = "order.Bacteroidales.id.913.summary_res_1.csv")
fwrite(order.Bacteroidales.id.913.summary_data_1,file = "order.Bacteroidales.id.913.summary_data_1.csv")
#order.Bifidobacteriales.id.432.summary.txt
order.Bifidobacteriales.id.432.summary<-fread("order.Bifidobacteriales.id.432.summary.txt")
order.Bifidobacteriales.id.432.summary<-rename(order.Bifidobacteriales.id.432.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Bifidobacteriales.id.432.summary,file = "order.Bifidobacteriales.id.432.summary.csv")
order.Bifidobacteriales.id.432.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Bifidobacteriales.id.432.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Bifidobacteriales.id.432.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Bifidobacteriales.id.432.summary
)
order.Bifidobacteriales.id.432.summary_res_1 <- mr(order.Bifidobacteriales.id.432.summary_data_1)
fwrite(order.Bifidobacteriales.id.432.summary_res_1,file = "order.Bifidobacteriales.id.432.summary_res_1.csv")
fwrite(order.Bifidobacteriales.id.432.summary_data_1,file = "order.Bifidobacteriales.id.432.summary_data_1.csv")
#order.Burkholderiales.id.2874.summary.txt
order.Burkholderiales.id.2874.summary<-fread("order.Burkholderiales.id.2874.summary.txt")
order.Burkholderiales.id.2874.summary<-rename(order.Burkholderiales.id.2874.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Burkholderiales.id.2874.summary,file = "order.Burkholderiales.id.2874.summary.csv")
order.Burkholderiales.id.2874.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Burkholderiales.id.2874.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Burkholderiales.id.2874.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Burkholderiales.id.2874.summary
)
order.Burkholderiales.id.2874.summary_res_1 <- mr(order.Burkholderiales.id.2874.summary_data_1)
fwrite(order.Burkholderiales.id.2874.summary_res_1,file = "order.Burkholderiales.id.2874.summary_res_1.csv")
fwrite(order.Burkholderiales.id.2874.summary_data_1,file = "order.Burkholderiales.id.2874.summary_data_1.csv")
#order.Clostridiales.id.1863.summary.txt
order.Clostridiales.id.1863.summary<-fread("order.Clostridiales.id.1863.summary.txt")
order.Clostridiales.id.1863.summary<-rename(order.Clostridiales.id.1863.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Clostridiales.id.1863.summary,file = "order.Clostridiales.id.1863.summary.csv")
order.Clostridiales.id.1863.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Clostridiales.id.1863.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Clostridiales.id.1863.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Clostridiales.id.1863.summary
)
order.Clostridiales.id.1863.summary_res_1 <- mr(order.Clostridiales.id.1863.summary_data_1)
fwrite(order.Clostridiales.id.1863.summary_res_1,file = "order.Clostridiales.id.1863.summary_res_1.csv")
fwrite(order.Clostridiales.id.1863.summary_data_1,file = "order.Clostridiales.id.1863.summary_data_1.csv")
#order.Coriobacteriales.id.810.summary.txt
order.Coriobacteriales.id.810.summary<-fread("order.Coriobacteriales.id.810.summary.txt")
order.Coriobacteriales.id.810.summary<-rename(order.Coriobacteriales.id.810.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Coriobacteriales.id.810.summary,file = "order.Coriobacteriales.id.810.summary.csv")
order.Coriobacteriales.id.810.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Coriobacteriales.id.810.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Coriobacteriales.id.810.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Coriobacteriales.id.810.summary
)
order.Coriobacteriales.id.810.summary_res_1 <- mr(order.Coriobacteriales.id.810.summary_data_1)
fwrite(order.Coriobacteriales.id.810.summary_res_1,file = "order.Coriobacteriales.id.810.summary_res_1.csv")
fwrite(order.Coriobacteriales.id.810.summary_data_1,file = "order.Coriobacteriales.id.810.summary_data_1.csv")
#order.Desulfovibrionales.id.3156.summary.txt
order.Desulfovibrionales.id.3156.summary<-fread("order.Desulfovibrionales.id.3156.summary.txt")
order.Desulfovibrionales.id.3156.summary<-rename(order.Desulfovibrionales.id.3156.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Desulfovibrionales.id.3156.summary,file = "order.Desulfovibrionales.id.3156.summary.csv")
order.Desulfovibrionales.id.3156.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Desulfovibrionales.id.3156.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Desulfovibrionales.id.3156.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Desulfovibrionales.id.3156.summary
)
order.Desulfovibrionales.id.3156.summary_res_1 <- mr(order.Desulfovibrionales.id.3156.summary_data_1)
fwrite(order.Desulfovibrionales.id.3156.summary_res_1,file = "order.Desulfovibrionales.id.3156.summary_res_1.csv")
fwrite(order.Desulfovibrionales.id.3156.summary_data_1,file = "order.Desulfovibrionales.id.3156.summary_data_1.csv")
#order.Enterobacteriales.id.3468.summary.txt
order.Enterobacteriales.id.3468.summary<-fread("order.Enterobacteriales.id.3468.summary.txt")
order.Enterobacteriales.id.3468.summary<-rename(order.Enterobacteriales.id.3468.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Enterobacteriales.id.3468.summary,file = "order.Enterobacteriales.id.3468.summary.csv")
order.Enterobacteriales.id.3468.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Enterobacteriales.id.3468.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Enterobacteriales.id.3468.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Enterobacteriales.id.3468.summary
)
order.Enterobacteriales.id.3468.summary_res_1 <- mr(order.Enterobacteriales.id.3468.summary_data_1)
fwrite(order.Enterobacteriales.id.3468.summary_res_1,file = "order.Enterobacteriales.id.3468.summary_res_1.csv")
fwrite(order.Enterobacteriales.id.3468.summary_data_1,file = "order.Enterobacteriales.id.3468.summary_data_1.csv")
#order.Erysipelotrichales.id.2148.summary.txt
order.Erysipelotrichales.id.2148.summary<-fread("order.Erysipelotrichales.id.2148.summary.txt")
order.Erysipelotrichales.id.2148.summary<-rename(order.Erysipelotrichales.id.2148.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Erysipelotrichales.id.2148.summary,file = "order.Erysipelotrichales.id.2148.summary.csv")
order.Erysipelotrichales.id.2148.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Erysipelotrichales.id.2148.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Erysipelotrichales.id.2148.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Erysipelotrichales.id.2148.summary
)
order.Erysipelotrichales.id.2148.summary_res_1 <- mr(order.Erysipelotrichales.id.2148.summary_data_1)
fwrite(order.Erysipelotrichales.id.2148.summary_res_1,file = "order.Erysipelotrichales.id.2148.summary_res_1.csv")
fwrite(order.Erysipelotrichales.id.2148.summary_data_1,file = "order.Erysipelotrichales.id.2148.summary_data_1.csv")
#order.Gastranaerophilales.id.1591.summary.txt
order.Gastranaerophilales.id.1591.summary<-fread("order.Gastranaerophilales.id.1591.summary.txt")
order.Gastranaerophilales.id.1591.summary<-rename(order.Gastranaerophilales.id.1591.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Gastranaerophilales.id.1591.summary,file = "order.Gastranaerophilales.id.1591.summary.csv")
order.Gastranaerophilales.id.1591.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Gastranaerophilales.id.1591.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Gastranaerophilales.id.1591.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Gastranaerophilales.id.1591.summary
)
order.Gastranaerophilales.id.1591.summary_res_1 <- mr(order.Gastranaerophilales.id.1591.summary_data_1)
fwrite(order.Gastranaerophilales.id.1591.summary_res_1,file = "order.Gastranaerophilales.id.1591.summary_res_1.csv")
fwrite(order.Gastranaerophilales.id.1591.summary_data_1,file = "order.Gastranaerophilales.id.1591.summary_data_1.csv")
#order.Lactobacillales.id.1800.summary.txt
order.Lactobacillales.id.1800.summary<-fread("order.Lactobacillales.id.1800.summary.txt")
order.Lactobacillales.id.1800.summary<-rename(order.Lactobacillales.id.1800.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Lactobacillales.id.1800.summary,file = "order.Lactobacillales.id.1800.summary.csv")
order.Lactobacillales.id.1800.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Lactobacillales.id.1800.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Lactobacillales.id.1800.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Lactobacillales.id.1800.summary
)
order.Lactobacillales.id.1800.summary_res_1 <- mr(order.Lactobacillales.id.1800.summary_data_1)
fwrite(order.Lactobacillales.id.1800.summary_res_1,file = "order.Lactobacillales.id.1800.summary_res_1.csv")
fwrite(order.Lactobacillales.id.1800.summary_data_1,file = "order.Lactobacillales.id.1800.summary_data_1.csv")
#order.Methanobacteriales.id.120.summary.txt
order.Methanobacteriales.id.120.summary<-fread("order.Methanobacteriales.id.120.summary.txt")
order.Methanobacteriales.id.120.summary<-rename(order.Methanobacteriales.id.120.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Methanobacteriales.id.120.summary,file = "order.Methanobacteriales.id.120.summary.csv")
order.Methanobacteriales.id.120.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Methanobacteriales.id.120.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Methanobacteriales.id.120.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Methanobacteriales.id.120.summary
)
order.Methanobacteriales.id.120.summary_res_1 <- mr(order.Methanobacteriales.id.120.summary_data_1)
fwrite(order.Methanobacteriales.id.120.summary_res_1,file = "order.Methanobacteriales.id.120.summary_res_1.csv")
fwrite(order.Methanobacteriales.id.120.summary_data_1,file = "order.Methanobacteriales.id.120.summary_data_1.csv")
#order.MollicutesRF9.id.11579.summary.txt
order.MollicutesRF9.id.11579.summary<-fread("order.MollicutesRF9.id.11579.summary.txt")
order.MollicutesRF9.id.11579.summary<-rename(order.MollicutesRF9.id.11579.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.MollicutesRF9.id.11579.summary,file = "order.MollicutesRF9.id.11579.summary.csv")
order.MollicutesRF9.id.11579.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.MollicutesRF9.id.11579.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.MollicutesRF9.id.11579.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.MollicutesRF9.id.11579.summary
)
order.MollicutesRF9.id.11579.summary_res_1 <- mr(order.MollicutesRF9.id.11579.summary_data_1)
fwrite(order.MollicutesRF9.id.11579.summary_res_1,file = "order.MollicutesRF9.id.11579.summary_res_1.csv")
fwrite(order.MollicutesRF9.id.11579.summary_data_1,file = "order.MollicutesRF9.id.11579.summary_data_1.csv")
#order.NB1n.id.3953.summary.txt
order.NB1n.id.3953.summary<-fread("order.NB1n.id.3953.summary.txt")
order.NB1n.id.3953.summary<-rename(order.NB1n.id.3953.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.NB1n.id.3953.summary,file = "order.NB1n.id.3953.summary.csv")
order.NB1n.id.3953.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.NB1n.id.3953.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.NB1n.id.3953.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.NB1n.id.3953.summary
)
order.NB1n.id.3953.summary_res_1 <- mr(order.NB1n.id.3953.summary_data_1)
fwrite(order.NB1n.id.3953.summary_res_1,file = "order.NB1n.id.3953.summary_res_1.csv")
fwrite(order.NB1n.id.3953.summary_data_1,file = "order.NB1n.id.3953.summary_data_1.csv")
#order.Pasteurellales.id.3688.summary.txt
order.Pasteurellales.id.3688.summary<-fread("order.Pasteurellales.id.3688.summary.txt")
order.Pasteurellales.id.3688.summary<-rename(order.Pasteurellales.id.3688.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Pasteurellales.id.3688.summary,file = "order.Pasteurellales.id.3688.summary.csv")
order.Pasteurellales.id.3688.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Pasteurellales.id.3688.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Pasteurellales.id.3688.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Pasteurellales.id.3688.summary
)
order.Pasteurellales.id.3688.summary_res_1 <- mr(order.Pasteurellales.id.3688.summary_data_1)
fwrite(order.Pasteurellales.id.3688.summary_res_1,file = "order.Pasteurellales.id.3688.summary_res_1.csv")
fwrite(order.Pasteurellales.id.3688.summary_data_1,file = "order.Pasteurellales.id.3688.summary_data_1.csv")
#order.Rhodospirillales.id.2667.summary.txt
order.Rhodospirillales.id.2667.summary<-fread("order.Rhodospirillales.id.2667.summary.txt")
order.Rhodospirillales.id.2667.summary<-rename(order.Rhodospirillales.id.2667.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Rhodospirillales.id.2667.summary,file = "order.Rhodospirillales.id.2667.summary.csv")
order.Rhodospirillales.id.2667.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Rhodospirillales.id.2667.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Rhodospirillales.id.2667.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Rhodospirillales.id.2667.summary
)
order.Rhodospirillales.id.2667.summary_res_1 <- mr(order.Rhodospirillales.id.2667.summary_data_1)
fwrite(order.Rhodospirillales.id.2667.summary_res_1,file = "order.Rhodospirillales.id.2667.summary_res_1.csv")
fwrite(order.Rhodospirillales.id.2667.summary_data_1,file = "order.Rhodospirillales.id.2667.summary_data_1.csv")
#order.Selenomonadales.id.2165.summary.txt
order.Selenomonadales.id.2165.summary<-fread("order.Selenomonadales.id.2165.summary.txt")
order.Selenomonadales.id.2165.summary<-rename(order.Selenomonadales.id.2165.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Selenomonadales.id.2165.summary,file = "order.Selenomonadales.id.2165.summary.csv")
order.Selenomonadales.id.2165.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Selenomonadales.id.2165.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Selenomonadales.id.2165.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Selenomonadales.id.2165.summary
)
order.Selenomonadales.id.2165.summary_res_1 <- mr(order.Selenomonadales.id.2165.summary_data_1)
fwrite(order.Selenomonadales.id.2165.summary_res_1,file = "order.Selenomonadales.id.2165.summary_res_1.csv")
fwrite(order.Selenomonadales.id.2165.summary_data_1,file = "order.Selenomonadales.id.2165.summary_data_1.csv")
#order.Verrucomicrobiales.id.4030.summary.txt
order.Verrucomicrobiales.id.4030.summary<-fread("order.Verrucomicrobiales.id.4030.summary.txt")
order.Verrucomicrobiales.id.4030.summary<-rename(order.Verrucomicrobiales.id.4030.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Verrucomicrobiales.id.4030.summary,file = "order.Verrucomicrobiales.id.4030.summary.csv")
order.Verrucomicrobiales.id.4030.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Verrucomicrobiales.id.4030.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Verrucomicrobiales.id.4030.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Verrucomicrobiales.id.4030.summary
)
order.Verrucomicrobiales.id.4030.summary_res_1 <- mr(order.Verrucomicrobiales.id.4030.summary_data_1)
fwrite(order.Verrucomicrobiales.id.4030.summary_res_1,file = "order.Verrucomicrobiales.id.4030.summary_res_1.csv")
fwrite(order.Verrucomicrobiales.id.4030.summary_data_1,file = "order.Verrucomicrobiales.id.4030.summary_data_1.csv")
#order.Victivallales.id.2254.summary.txt
order.Victivallales.id.2254.summary<-fread("order.Victivallales.id.2254.summary.txt")
order.Victivallales.id.2254.summary<-rename(order.Victivallales.id.2254.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(order.Victivallales.id.2254.summary,file = "order.Victivallales.id.2254.summary.csv")
order.Victivallales.id.2254.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "order.Victivallales.id.2254.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
order.Victivallales.id.2254.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = order.Victivallales.id.2254.summary
)
order.Victivallales.id.2254.summary_res_1 <- mr(order.Victivallales.id.2254.summary_data_1)
fwrite(order.Victivallales.id.2254.summary_res_1,file = "order.Victivallales.id.2254.summary_res_1.csv")
fwrite(order.Victivallales.id.2254.summary_data_1,file = "order.Victivallales.id.2254.summary_data_1.csv")
#setwd("G:/gut/MiBioGen_QmbQTL_summary_family/MiBioGen_QmbQTL_summary_family")
setwd("G:/gut/MiBioGen_QmbQTL_summary_family/MiBioGen_QmbQTL_summary_family")
#family.Acidaminococcaceae.id.2166.summary.txt
family.Acidaminococcaceae.id.2166.summary<-fread("family.Acidaminococcaceae.id.2166.summary.txt")
family.Acidaminococcaceae.id.2166.summary<-rename(family.Acidaminococcaceae.id.2166.summary, 
                                                  id.outcome = bac,
                                                  beta.outcome = beta,
                                                  se.outcome = SE,
                                                  effect_allele.outcome = eff.allele,
                                                  other_allele.outcome = ref.allele,
                                                  SNP = rsID)
fwrite(family.Acidaminococcaceae.id.2166.summary,file = "family.Acidaminococcaceae.id.2166.summary.csv")
family.Acidaminococcaceae.id.2166.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Acidaminococcaceae.id.2166.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Acidaminococcaceae.id.2166.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Acidaminococcaceae.id.2166.summary
)
family.Acidaminococcaceae.id.2166.summary_res_1 <- mr(family.Acidaminococcaceae.id.2166.summary_data_1)
fwrite(family.Acidaminococcaceae.id.2166.summary_res_1,file = "family.Acidaminococcaceae.id.2166.summary_res_1.csv")
fwrite(family.Acidaminococcaceae.id.2166.summary_data_1,file = "family.Acidaminococcaceae.id.2166.summary_data_1.csv")
#family.Actinomycetaceae.id.421.summary.txt
family.Actinomycetaceae.id.421.summary<-fread("family.Actinomycetaceae.id.421.summary.txt")
family.Actinomycetaceae.id.421.summary<-rename(family.Actinomycetaceae.id.421.summary, 
                                               id.outcome = bac,
                                               beta.outcome = beta,
                                               se.outcome = SE,
                                               effect_allele.outcome = eff.allele,
                                               other_allele.outcome = ref.allele,
                                               SNP = rsID)
fwrite(family.Actinomycetaceae.id.421.summary,file = "family.Actinomycetaceae.id.421.summary.csv")
family.Actinomycetaceae.id.421.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Actinomycetaceae.id.421.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Actinomycetaceae.id.421.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Actinomycetaceae.id.421.summary
)
family.Actinomycetaceae.id.421.summary_res_1 <- mr(family.Actinomycetaceae.id.421.summary_data_1)
fwrite(family.Actinomycetaceae.id.421.summary_res_1,file = "family.Actinomycetaceae.id.421.summary_res_1.csv")
fwrite(family.Actinomycetaceae.id.421.summary_data_1,file = "family.Actinomycetaceae.id.421.summary_data_1.csv")
#family.Alcaligenaceae.id.2875.summary.txt
family.Alcaligenaceae.id.2875.summary<-fread("family.Alcaligenaceae.id.2875.summary.txt")
family.Alcaligenaceae.id.2875.summary<-rename(family.Alcaligenaceae.id.2875.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(family.Alcaligenaceae.id.2875.summary,file = "family.Alcaligenaceae.id.2875.summary.csv")
family.Alcaligenaceae.id.2875.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Alcaligenaceae.id.2875.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Alcaligenaceae.id.2875.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Alcaligenaceae.id.2875.summary
)
family.Alcaligenaceae.id.2875.summary_res_1 <- mr(family.Alcaligenaceae.id.2875.summary_data_1)
fwrite(family.Alcaligenaceae.id.2875.summary_res_1,file = "family.Alcaligenaceae.id.2875.summary_res_1.csv")
fwrite(family.Alcaligenaceae.id.2875.summary_data_1,file = "family.Alcaligenaceae.id.2875.summary_data_1.csv")
#family.Bacteroidaceae.id.917.summary.txt
family.Bacteroidaceae.id.917.summary<-fread("family.Bacteroidaceae.id.917.summary.txt")
family.Bacteroidaceae.id.917.summary<-rename(family.Bacteroidaceae.id.917.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(family.Bacteroidaceae.id.917.summary,file = "family.Bacteroidaceae.id.917.summary.csv")
family.Bacteroidaceae.id.917.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Bacteroidaceae.id.917.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Bacteroidaceae.id.917.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Bacteroidaceae.id.917.summary
)
family.Bacteroidaceae.id.917.summary_res_1 <- mr(family.Bacteroidaceae.id.917.summary_data_1)
fwrite(family.Bacteroidaceae.id.917.summary_res_1,file = "family.Bacteroidaceae.id.917.summary_res_1.csv")
fwrite(family.Bacteroidaceae.id.917.summary_data_1,file = "family.Bacteroidaceae.id.917.summary_data_1.csv")
#family.BacteroidalesS24.7group.id.11173.summary.txt
family.BacteroidalesS24.7group.id.11173.summary<-fread("family.BacteroidalesS24.7group.id.11173.summary.txt")
family.BacteroidalesS24.7group.id.11173.summary<-rename(family.BacteroidalesS24.7group.id.11173.summary, 
                                                        id.outcome = bac,
                                                        beta.outcome = beta,
                                                        se.outcome = SE,
                                                        effect_allele.outcome = eff.allele,
                                                        other_allele.outcome = ref.allele,
                                                        SNP = rsID)
fwrite(family.BacteroidalesS24.7group.id.11173.summary,file = "family.BacteroidalesS24.7group.id.11173.summary.csv")
family.BacteroidalesS24.7group.id.11173.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.BacteroidalesS24.7group.id.11173.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.BacteroidalesS24.7group.id.11173.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.BacteroidalesS24.7group.id.11173.summary
)
family.BacteroidalesS24.7group.id.11173.summary_res_1 <- mr(family.BacteroidalesS24.7group.id.11173.summary_data_1)
fwrite(family.BacteroidalesS24.7group.id.11173.summary_res_1,file = "family.BacteroidalesS24.7group.id.11173.summary_res_1.csv")
fwrite(family.BacteroidalesS24.7group.id.11173.summary_data_1,file = "family.BacteroidalesS24.7group.id.11173.summary_data_1.csv")
#family.Bifidobacteriaceae.id.433.summary.txt
family.Bifidobacteriaceae.id.433.summary<-fread("family.Bifidobacteriaceae.id.433.summary.txt")
family.Bifidobacteriaceae.id.433.summary<-rename(family.Bifidobacteriaceae.id.433.summary, 
                                                 id.outcome = bac,
                                                 beta.outcome = beta,
                                                 se.outcome = SE,
                                                 effect_allele.outcome = eff.allele,
                                                 other_allele.outcome = ref.allele,
                                                 SNP = rsID)
fwrite(family.Bifidobacteriaceae.id.433.summary,file = "family.Bifidobacteriaceae.id.433.summary.csv")
family.Bifidobacteriaceae.id.433.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Bifidobacteriaceae.id.433.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Bifidobacteriaceae.id.433.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Bifidobacteriaceae.id.433.summary
)
family.Bifidobacteriaceae.id.433.summary_res_1 <- mr(family.Bifidobacteriaceae.id.433.summary_data_1)
fwrite(family.Bifidobacteriaceae.id.433.summary_res_1,file = "family.Bifidobacteriaceae.id.433.summary_res_1.csv")
fwrite(family.Bifidobacteriaceae.id.433.summary_data_1,file = "family.Bifidobacteriaceae.id.433.summary_data_1.csv")
#family.Christensenellaceae.id.1866.summary.txt
family.Christensenellaceae.id.1866.summary<-fread("family.Christensenellaceae.id.1866.summary.txt")
family.Christensenellaceae.id.1866.summary<-rename(family.Christensenellaceae.id.1866.summary, 
                                                 id.outcome = bac,
                                                 beta.outcome = beta,
                                                 se.outcome = SE,
                                                 effect_allele.outcome = eff.allele,
                                                 other_allele.outcome = ref.allele,
                                                 SNP = rsID)
fwrite(family.Christensenellaceae.id.1866.summary,file = "family.Christensenellaceae.id.1866.summary.csv")
family.Christensenellaceae.id.1866.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Christensenellaceae.id.1866.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Christensenellaceae.id.1866.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Christensenellaceae.id.1866.summary
)
family.Christensenellaceae.id.1866.summary_res_1 <- mr(family.Christensenellaceae.id.1866.summary_data_1)
fwrite(family.Christensenellaceae.id.1866.summary_res_1,file = "family.Christensenellaceae.id.1866.summary_res_1.csv")
fwrite(family.Christensenellaceae.id.1866.summary_data_1,file = "family.Christensenellaceae.id.1866.summary_data_1.csv")
#family.Clostridiaceae1.id.1869.summary.txt
family.Clostridiaceae1.id.1869.summary<-fread("family.Clostridiaceae1.id.1869.summary.txt")
family.Clostridiaceae1.id.1869.summary<-rename(family.Clostridiaceae1.id.1869.summary, 
                                               id.outcome = bac,
                                               beta.outcome = beta,
                                               se.outcome = SE,
                                               effect_allele.outcome = eff.allele,
                                               other_allele.outcome = ref.allele,
                                               SNP = rsID)
fwrite(family.Clostridiaceae1.id.1869.summary,file = "family.Clostridiaceae1.id.1869.summary.csv")
family.Clostridiaceae1.id.1869.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Clostridiaceae1.id.1869.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Clostridiaceae1.id.1869.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Clostridiaceae1.id.1869.summary
)
family.Clostridiaceae1.id.1869.summary_res_1 <- mr(family.Clostridiaceae1.id.1869.summary_data_1)
fwrite(family.Clostridiaceae1.id.1869.summary_res_1,file = "family.Clostridiaceae1.id.1869.summary_res_1.csv")
fwrite(family.Clostridiaceae1.id.1869.summary_data_1,file = "family.Clostridiaceae1.id.1869.summary_data_1.csv")#family.Christensenellaceae.id.1866.summary.txt
#family.ClostridialesvadinBB60group.id.11286.summary.txt
family.ClostridialesvadinBB60group.id.11286.summary<-fread("family.ClostridialesvadinBB60group.id.11286.summary.txt")
family.ClostridialesvadinBB60group.id.11286.summary<-rename(family.ClostridialesvadinBB60group.id.11286.summary, 
                                                            id.outcome = bac,
                                                            beta.outcome = beta,
                                                            se.outcome = SE,
                                                            effect_allele.outcome = eff.allele,
                                                            other_allele.outcome = ref.allele,
                                                            SNP = rsID)
fwrite(family.ClostridialesvadinBB60group.id.11286.summary,file = "family.ClostridialesvadinBB60group.id.11286.summary.csv")
family.ClostridialesvadinBB60group.id.11286.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.ClostridialesvadinBB60group.id.11286.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.ClostridialesvadinBB60group.id.11286.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.ClostridialesvadinBB60group.id.11286.summary
)
family.ClostridialesvadinBB60group.id.11286.summary_res_1 <- mr(family.ClostridialesvadinBB60group.id.11286.summary_data_1)
fwrite(family.ClostridialesvadinBB60group.id.11286.summary_res_1,file = "family.ClostridialesvadinBB60group.id.11286.summary_res_1.csv")
fwrite(family.ClostridialesvadinBB60group.id.11286.summary_data_1,file = "family.ClostridialesvadinBB60group.id.11286.summary_data_1.csv")#family.Christensenellaceae.id.1866.summary.txt
#family.Coriobacteriaceae.id.811.summary.txt
family.Coriobacteriaceae.id.811.summary<-fread("family.Coriobacteriaceae.id.811.summary.txt")
family.Coriobacteriaceae.id.811.summary<-rename(family.Coriobacteriaceae.id.811.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Coriobacteriaceae.id.811.summary,file = "family.Coriobacteriaceae.id.811.summary.csv")
family.Coriobacteriaceae.id.811.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Coriobacteriaceae.id.811.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Coriobacteriaceae.id.811.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Coriobacteriaceae.id.811.summary
)
family.Coriobacteriaceae.id.811.summary_res_1 <- mr(family.Coriobacteriaceae.id.811.summary_data_1)
fwrite(family.Coriobacteriaceae.id.811.summary_res_1,file = "family.Coriobacteriaceae.id.811.summary_res_1.csv")
fwrite(family.Coriobacteriaceae.id.811.summary_data_1,file = "family.Coriobacteriaceae.id.811.summary_data_1.csv")#family.Christensenellaceae.id.1866.summary.txt
#family.Defluviitaleaceae.id.1924.summary.txt
family.Defluviitaleaceae.id.1924.summary<-fread("family.Defluviitaleaceae.id.1924.summary.txt")
family.Defluviitaleaceae.id.1924.summary<-rename(family.Defluviitaleaceae.id.1924.summary, 
                                                 id.outcome = bac,
                                                 beta.outcome = beta,
                                                 se.outcome = SE,
                                                 effect_allele.outcome = eff.allele,
                                                 other_allele.outcome = ref.allele,
                                                 SNP = rsID)
fwrite(family.Defluviitaleaceae.id.1924.summary,file = "family.Defluviitaleaceae.id.1924.summary.csv")
family.Defluviitaleaceae.id.1924.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Defluviitaleaceae.id.1924.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Defluviitaleaceae.id.1924.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Defluviitaleaceae.id.1924.summary
)
family.Defluviitaleaceae.id.1924.summary_res_1 <- mr(family.Defluviitaleaceae.id.1924.summary_data_1)
fwrite(family.Defluviitaleaceae.id.1924.summary_res_1,file = "family.Defluviitaleaceae.id.1924.summary_res_1.csv")
fwrite(family.Defluviitaleaceae.id.1924.summary_data_1,file = "family.Defluviitaleaceae.id.1924.summary_data_1.csv")#family.Christensenellaceae.id.1866.summary.txt
#family.Desulfovibrionaceae.id.3169.summary.txt
family.Desulfovibrionaceae.id.3169.summary<-fread("family.Desulfovibrionaceae.id.3169.summary.txt")
family.Desulfovibrionaceae.id.3169.summary<-rename(family.Desulfovibrionaceae.id.3169.summary, 
                                                   id.outcome = bac,
                                                   beta.outcome = beta,
                                                   se.outcome = SE,
                                                   effect_allele.outcome = eff.allele,
                                                   other_allele.outcome = ref.allele,
                                                   SNP = rsID)
fwrite(family.Desulfovibrionaceae.id.3169.summary,file = "family.Desulfovibrionaceae.id.3169.summary.csv")
family.Desulfovibrionaceae.id.3169.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Desulfovibrionaceae.id.3169.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Desulfovibrionaceae.id.3169.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Desulfovibrionaceae.id.3169.summary
)
family.Desulfovibrionaceae.id.3169.summary_res_1 <- mr(family.Desulfovibrionaceae.id.3169.summary_data_1)
fwrite(family.Desulfovibrionaceae.id.3169.summary_res_1,file = "family.Desulfovibrionaceae.id.3169.summary_res_1.csv")
fwrite(family.Desulfovibrionaceae.id.3169.summary_data_1,file = "family.Desulfovibrionaceae.id.3169.summary_data_1.csv")#family.Christensenellaceae.id.1866.summary.txt
#family.Enterobacteriaceae.id.3469.summary.txt
family.Enterobacteriaceae.id.3469.summary<-fread("family.Enterobacteriaceae.id.3469.summary.txt")
family.Enterobacteriaceae.id.3469.summary<-rename(family.Enterobacteriaceae.id.3469.summary, 
                                                  id.outcome = bac,
                                                  beta.outcome = beta,
                                                  se.outcome = SE,
                                                  effect_allele.outcome = eff.allele,
                                                  other_allele.outcome = ref.allele,
                                                  SNP = rsID)
fwrite(family.Enterobacteriaceae.id.3469.summary,file = "family.Enterobacteriaceae.id.3469.summary.csv")
family.Enterobacteriaceae.id.3469.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Enterobacteriaceae.id.3469.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Enterobacteriaceae.id.3469.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Enterobacteriaceae.id.3469.summary
)
family.Enterobacteriaceae.id.3469.summary_res_1 <- mr(family.Enterobacteriaceae.id.3469.summary_data_1)
fwrite(family.Enterobacteriaceae.id.3469.summary_res_1,file = "family.Enterobacteriaceae.id.3469.summary_res_1.csv")
fwrite(family.Enterobacteriaceae.id.3469.summary_data_1,file = "family.Enterobacteriaceae.id.3469.summary_data_1.csv")#family.Christensenellaceae.id.1866.summary.txt
#family.Erysipelotrichaceae.id.2149.summary.txt
family.Erysipelotrichaceae.id.2149.summary<-fread("family.Erysipelotrichaceae.id.2149.summary.txt")
family.Erysipelotrichaceae.id.2149.summary<-rename(family.Erysipelotrichaceae.id.2149.summary, 
                                                   id.outcome = bac,
                                                   beta.outcome = beta,
                                                   se.outcome = SE,
                                                   effect_allele.outcome = eff.allele,
                                                   other_allele.outcome = ref.allele,
                                                   SNP = rsID)
fwrite(family.Erysipelotrichaceae.id.2149.summary,file = "family.Erysipelotrichaceae.id.2149.summary.csv")
family.Erysipelotrichaceae.id.2149.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Erysipelotrichaceae.id.2149.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Erysipelotrichaceae.id.2149.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Erysipelotrichaceae.id.2149.summary
)
family.Erysipelotrichaceae.id.2149.summary_res_1 <- mr(family.Erysipelotrichaceae.id.2149.summary_data_1)
fwrite(family.Erysipelotrichaceae.id.2149.summary_res_1,file = "family.Erysipelotrichaceae.id.2149.summary_res_1.csv")
fwrite(family.Erysipelotrichaceae.id.2149.summary_data_1,file = "family.Erysipelotrichaceae.id.2149.summary_data_1.csv")#family.Christensenellaceae.id.1866.summary.txt
#family.FamilyXI.id.1936.summary.txt
family.FamilyXI.id.1936.summary<-fread("family.FamilyXI.id.1936.summary.txt")
family.FamilyXI.id.1936.summary<-rename(family.FamilyXI.id.1936.summary, 
                                                   id.outcome = bac,
                                                   beta.outcome = beta,
                                                   se.outcome = SE,
                                                   effect_allele.outcome = eff.allele,
                                                   other_allele.outcome = ref.allele,
                                                   SNP = rsID)
fwrite(family.FamilyXI.id.1936.summary,file = "family.FamilyXI.id.1936.summary.csv")
family.FamilyXI.id.1936.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.FamilyXI.id.1936.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.FamilyXI.id.1936.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.FamilyXI.id.1936.summary
)
family.FamilyXI.id.1936.summary_res_1 <- mr(family.FamilyXI.id.1936.summary_data_1)
fwrite(family.FamilyXI.id.1936.summary_res_1,file = "family.FamilyXI.id.1936.summary_res_1.csv")
fwrite(family.FamilyXI.id.1936.summary_data_1,file = "family.FamilyXI.id.1936.summary_data_1.csv")
#family.familyXIII.id.1957.summary.txt
family.familyXIII.id.1957.summary<-fread("family.familyXIII.id.1957.summary.txt")
family.familyXIII.id.1957.summary<-rename(family.familyXIII.id.1957.summary, 
                                        id.outcome = bac,
                                        beta.outcome = beta,
                                        se.outcome = SE,
                                        effect_allele.outcome = eff.allele,
                                        other_allele.outcome = ref.allele,
                                        SNP = rsID)
fwrite(family.familyXIII.id.1957.summary,file = "family.familyXIII.id.1957.summary.csv")
family.familyXIII.id.1957.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.familyXIII.id.1957.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.familyXIII.id.1957.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.familyXIII.id.1957.summary
)
family.familyXIII.id.1957.summary_res_1 <- mr(family.familyXIII.id.1957.summary_data_1)
fwrite(family.familyXIII.id.1957.summary_res_1,file = "family.familyXIII.id.1957.summary_res_1.csv")
fwrite(family.familyXIII.id.1957.summary_data_1,file = "family.familyXIII.id.1957.summary_data_1.csv")
#family.Lachnospiraceae.id.1987.summary.txt
family.Lachnospiraceae.id.1987.summary<-fread("family.Lachnospiraceae.id.1987.summary.txt")
family.Lachnospiraceae.id.1987.summary<-rename(family.Lachnospiraceae.id.1987.summary, 
                                        id.outcome = bac,
                                        beta.outcome = beta,
                                        se.outcome = SE,
                                        effect_allele.outcome = eff.allele,
                                        other_allele.outcome = ref.allele,
                                        SNP = rsID)
fwrite(family.Lachnospiraceae.id.1987.summary,file = "family.Lachnospiraceae.id.1987.summary.csv")
family.Lachnospiraceae.id.1987.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Lachnospiraceae.id.1987.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Lachnospiraceae.id.1987.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Lachnospiraceae.id.1987.summary
)
family.Lachnospiraceae.id.1987.summary_res_1 <- mr(family.Lachnospiraceae.id.1987.summary_data_1)
fwrite(family.Lachnospiraceae.id.1987.summary_res_1,file = "family.Lachnospiraceae.id.1987.summary_res_1.csv")
fwrite(family.Lachnospiraceae.id.1987.summary_data_1,file = "family.Lachnospiraceae.id.1987.summary_data_1.csv")
#family.Lactobacillaceae.id.1836.summary.txt
family.Lactobacillaceae.id.1836.summary<-fread("family.Lactobacillaceae.id.1836.summary.txt")
family.Lactobacillaceae.id.1836.summary<-rename(family.Lactobacillaceae.id.1836.summary, 
                                               id.outcome = bac,
                                               beta.outcome = beta,
                                               se.outcome = SE,
                                               effect_allele.outcome = eff.allele,
                                               other_allele.outcome = ref.allele,
                                               SNP = rsID)
fwrite(family.Lactobacillaceae.id.1836.summary,file = "family.Lactobacillaceae.id.1836.summary.csv")
family.Lactobacillaceae.id.1836.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Lactobacillaceae.id.1836.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Lactobacillaceae.id.1836.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Lactobacillaceae.id.1836.summary
)
family.Lactobacillaceae.id.1836.summary_res_1 <- mr(family.Lactobacillaceae.id.1836.summary_data_1)
fwrite(family.Lactobacillaceae.id.1836.summary_res_1,file = "family.Lactobacillaceae.id.1836.summary_res_1.csv")
fwrite(family.Lactobacillaceae.id.1836.summary_data_1,file = "family.Lactobacillaceae.id.1836.summary_data_1.csv")
#family.Methanobacteriaceae.id.121.summary.txt
family.Methanobacteriaceae.id.121.summary<-fread("family.Methanobacteriaceae.id.121.summary.txt")
family.Methanobacteriaceae.id.121.summary<-rename(family.Methanobacteriaceae.id.121.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Methanobacteriaceae.id.121.summary,file = "family.Methanobacteriaceae.id.121.summary.csv")
family.Methanobacteriaceae.id.121.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Methanobacteriaceae.id.121.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Methanobacteriaceae.id.121.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Methanobacteriaceae.id.121.summary
)
family.Methanobacteriaceae.id.121.summary_res_1 <- mr(family.Methanobacteriaceae.id.121.summary_data_1)
fwrite(family.Methanobacteriaceae.id.121.summary_res_1,file = "family.Methanobacteriaceae.id.121.summary_res_1.csv")
fwrite(family.Methanobacteriaceae.id.121.summary_data_1,file = "family.Methanobacteriaceae.id.121.summary_data_1.csv")
#family.Oxalobacteraceae.id.2966.summary.txt
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
#family.Pasteurellaceae.id.3689.summary.txt
family.Pasteurellaceae.id.3689.summary<-fread("family.Pasteurellaceae.id.3689.summary.txt")
family.Pasteurellaceae.id.3689.summary<-rename(family.Pasteurellaceae.id.3689.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Pasteurellaceae.id.3689.summary,file = "family.Pasteurellaceae.id.3689.summary.csv")
family.Pasteurellaceae.id.3689.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Pasteurellaceae.id.3689.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Pasteurellaceae.id.3689.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Pasteurellaceae.id.3689.summary
)
family.Pasteurellaceae.id.3689.summary_res_1 <- mr(family.Pasteurellaceae.id.3689.summary_data_1)
fwrite(family.Pasteurellaceae.id.3689.summary_res_1,file = "family.Pasteurellaceae.id.3689.summary_res_1.csv")
fwrite(family.Pasteurellaceae.id.3689.summary_data_1,file = "family.Pasteurellaceae.id.3689.summary_data_1.csv")
#family.Peptococcaceae.id.2024.summary.txt
family.Peptococcaceae.id.2024.summary<-fread("family.Peptococcaceae.id.2024.summary.txt")
family.Peptococcaceae.id.2024.summary<-rename(family.Peptococcaceae.id.2024.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Peptococcaceae.id.2024.summary,file = "family.Peptococcaceae.id.2024.summary.csv")
family.Peptococcaceae.id.2024.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Peptococcaceae.id.2024.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Peptococcaceae.id.2024.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Peptococcaceae.id.2024.summary
)
family.Peptococcaceae.id.2024.summary_res_1 <- mr(family.Peptococcaceae.id.2024.summary_data_1)
fwrite(family.Peptococcaceae.id.2024.summary_res_1,file = "family.Peptococcaceae.id.2024.summary_res_1.csv")
fwrite(family.Peptococcaceae.id.2024.summary_data_1,file = "family.Peptococcaceae.id.2024.summary_data_1.csv")
#family.Peptostreptococcaceae.id.2042.summary.txt
family.Peptostreptococcaceae.id.2042.summary<-fread("family.Peptostreptococcaceae.id.2042.summary.txt")
family.Peptostreptococcaceae.id.2042.summary<-rename(family.Peptostreptococcaceae.id.2042.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Peptostreptococcaceae.id.2042.summary,file = "family.Peptostreptococcaceae.id.2042.summary.csv")
family.Peptostreptococcaceae.id.2042.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Peptostreptococcaceae.id.2042.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Peptostreptococcaceae.id.2042.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Peptostreptococcaceae.id.2042.summary
)
family.Peptostreptococcaceae.id.2042.summary_res_1 <- mr(family.Peptostreptococcaceae.id.2042.summary_data_1)
fwrite(family.Peptostreptococcaceae.id.2042.summary_res_1,file = "family.Peptostreptococcaceae.id.2042.summary_res_1.csv")
fwrite(family.Peptostreptococcaceae.id.2042.summary_data_1,file = "family.Peptostreptococcaceae.id.2042.summary_data_1.csv")
#family.Porphyromonadaceae.id.943.summary.txt
family.Porphyromonadaceae.id.943.summary<-fread("family.Porphyromonadaceae.id.943.summary.txt")
family.Porphyromonadaceae.id.943.summary<-rename(family.Porphyromonadaceae.id.943.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Porphyromonadaceae.id.943.summary,file = "family.Porphyromonadaceae.id.943.summary.csv")
family.Porphyromonadaceae.id.943.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Porphyromonadaceae.id.943.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Porphyromonadaceae.id.943.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Porphyromonadaceae.id.943.summary
)
family.Porphyromonadaceae.id.943.summary_res_1 <- mr(family.Porphyromonadaceae.id.943.summary_data_1)
fwrite(family.Porphyromonadaceae.id.943.summary_res_1,file = "family.Porphyromonadaceae.id.943.summary_res_1.csv")
fwrite(family.Porphyromonadaceae.id.943.summary_data_1,file = "family.Porphyromonadaceae.id.943.summary_data_1.csv")
#family.Prevotellaceae.id.960.summary.txt
family.Prevotellaceae.id.960.summary<-fread("family.Prevotellaceae.id.960.summary.txt")
family.Prevotellaceae.id.960.summary<-rename(family.Prevotellaceae.id.960.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Prevotellaceae.id.960.summary,file = "family.Prevotellaceae.id.960.summary.csv")
family.Prevotellaceae.id.960.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Prevotellaceae.id.960.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Prevotellaceae.id.960.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Prevotellaceae.id.960.summary
)
family.Prevotellaceae.id.960.summary_res_1 <- mr(family.Prevotellaceae.id.960.summary_data_1)
fwrite(family.Prevotellaceae.id.960.summary_res_1,file = "family.Prevotellaceae.id.960.summary_res_1.csv")
fwrite(family.Prevotellaceae.id.960.summary_data_1,file = "family.Prevotellaceae.id.960.summary_data_1.csv")
#family.Rhodospirillaceae.id.2717.summary.txt
family.Rhodospirillaceae.id.2717.summary<-fread("family.Rhodospirillaceae.id.2717.summary.txt")
family.Rhodospirillaceae.id.2717.summary<-rename(family.Rhodospirillaceae.id.2717.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Rhodospirillaceae.id.2717.summary,file = "family.Rhodospirillaceae.id.2717.summary.csv")
family.Rhodospirillaceae.id.2717.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Rhodospirillaceae.id.2717.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Rhodospirillaceae.id.2717.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Rhodospirillaceae.id.2717.summary
)
family.Rhodospirillaceae.id.2717.summary_res_1 <- mr(family.Rhodospirillaceae.id.2717.summary_data_1)
fwrite(family.Rhodospirillaceae.id.2717.summary_res_1,file = "family.Rhodospirillaceae.id.2717.summary_res_1.csv")
fwrite(family.Rhodospirillaceae.id.2717.summary_data_1,file = "family.Rhodospirillaceae.id.2717.summary_data_1.csv")
#family.Rikenellaceae.id.967.summary.txt
family.Rikenellaceae.id.967.summary<-fread("family.Rikenellaceae.id.967.summary.txt")
family.Rikenellaceae.id.967.summary<-rename(family.Rikenellaceae.id.967.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Rikenellaceae.id.967.summary,file = "family.Rikenellaceae.id.967.summary.csv")
family.Rikenellaceae.id.967.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Rikenellaceae.id.967.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Rikenellaceae.id.967.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Rikenellaceae.id.967.summary
)
family.Rikenellaceae.id.967.summary_res_1 <- mr(family.Rikenellaceae.id.967.summary_data_1)
fwrite(family.Rikenellaceae.id.967.summary_res_1,file = "family.Rikenellaceae.id.967.summary_res_1.csv")
fwrite(family.Rikenellaceae.id.967.summary_data_1,file = "family.Rikenellaceae.id.967.summary_data_1.csv")
#family.Ruminococcaceae.id.2050.summary.txt
family.Ruminococcaceae.id.2050.summary<-fread("family.Ruminococcaceae.id.2050.summary.txt")
family.Ruminococcaceae.id.2050.summary<-rename(family.Ruminococcaceae.id.2050.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Ruminococcaceae.id.2050.summary,file = "family.Ruminococcaceae.id.2050.summary.csv")
family.Ruminococcaceae.id.2050.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Ruminococcaceae.id.2050.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Ruminococcaceae.id.2050.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Ruminococcaceae.id.2050.summary
)
family.Ruminococcaceae.id.2050.summary_res_1 <- mr(family.Ruminococcaceae.id.2050.summary_data_1)
fwrite(family.Ruminococcaceae.id.2050.summary_res_1,file = "family.Ruminococcaceae.id.2050.summary_res_1.csv")
fwrite(family.Ruminococcaceae.id.2050.summary_data_1,file = "family.Ruminococcaceae.id.2050.summary_data_1.csv")
#family.Streptococcaceae.id.1850.summary.txt
family.Streptococcaceae.id.1850.summary<-fread("family.Streptococcaceae.id.1850.summary.txt")
family.Streptococcaceae.id.1850.summary<-rename(family.Streptococcaceae.id.1850.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Streptococcaceae.id.1850.summary,file = "family.Streptococcaceae.id.1850.summary.csv")
family.Streptococcaceae.id.1850.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Streptococcaceae.id.1850.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Streptococcaceae.id.1850.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Streptococcaceae.id.1850.summary
)
family.Streptococcaceae.id.1850.summary_res_1 <- mr(family.Streptococcaceae.id.1850.summary_data_1)
fwrite(family.Streptococcaceae.id.1850.summary_res_1,file = "family.Streptococcaceae.id.1850.summary_res_1.csv")
fwrite(family.Streptococcaceae.id.1850.summary_data_1,file = "family.Streptococcaceae.id.1850.summary_data_1.csv")
#family.unknownfamily.id.1000001214.summary.txt
family.unknownfamily.id.1000001214.summary<-fread("family.unknownfamily.id.1000001214.summary.txt")
family.unknownfamily.id.1000001214.summary<-rename(family.unknownfamily.id.1000001214.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.unknownfamily.id.1000001214.summary,file = "family.unknownfamily.id.1000001214.summary.csv")
family.unknownfamily.id.1000001214.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.unknownfamily.id.1000001214.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.unknownfamily.id.1000001214.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.unknownfamily.id.1000001214.summary
)
family.unknownfamily.id.1000001214.summary_res_1 <- mr(family.unknownfamily.id.1000001214.summary_data_1)
fwrite(family.unknownfamily.id.1000001214.summary_res_1,file = "family.unknownfamily.id.1000001214.summary_res_1.csv")
fwrite(family.unknownfamily.id.1000001214.summary_data_1,file = "family.unknownfamily.id.1000001214.summary_data_1.csv")
#family.unknownfamily.id.1000005471.summary.txt
family.unknownfamily.id.1000005471.summary<-fread("family.unknownfamily.id.1000005471.summary.txt")
family.unknownfamily.id.1000005471.summary<-rename(family.unknownfamily.id.1000005471.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.unknownfamily.id.1000005471.summary,file = "family.unknownfamily.id.1000005471.summary.csv")
family.unknownfamily.id.1000005471.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.unknownfamily.id.1000005471.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.unknownfamily.id.1000005471.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.unknownfamily.id.1000005471.summary
)
family.unknownfamily.id.1000005471.summary_res_1 <- mr(family.unknownfamily.id.1000005471.summary_data_1)
fwrite(family.unknownfamily.id.1000005471.summary_res_1,file = "family.unknownfamily.id.1000005471.summary_res_1.csv")
fwrite(family.unknownfamily.id.1000005471.summary_data_1,file = "family.unknownfamily.id.1000005471.summary_data_1.csv")
#family.unknownfamily.id.1000006161.summary.txt
family.unknownfamily.id.1000006161.summary<-fread("family.unknownfamily.id.1000006161.summary.txt")
family.unknownfamily.id.1000006161.summary<-rename(family.unknownfamily.id.1000006161.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.unknownfamily.id.1000006161.summary,file = "family.unknownfamily.id.1000006161.summary.csv")
family.unknownfamily.id.1000006161.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.unknownfamily.id.1000006161.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.unknownfamily.id.1000006161.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.unknownfamily.id.1000006161.summary
)
family.unknownfamily.id.1000006161.summary_res_1 <- mr(family.unknownfamily.id.1000006161.summary_data_1)
fwrite(family.unknownfamily.id.1000006161.summary_res_1,file = "family.unknownfamily.id.1000006161.summary_res_1.csv")
fwrite(family.unknownfamily.id.1000006161.summary_data_1,file = "family.unknownfamily.id.1000006161.summary_data_1.csv")
#family.Veillonellaceae.id.2172.summary.txt
family.Veillonellaceae.id.2172.summary<-fread("family.Veillonellaceae.id.2172.summary.txt")
family.Veillonellaceae.id.2172.summary<-rename(family.Veillonellaceae.id.2172.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Veillonellaceae.id.2172.summary,file = "family.Veillonellaceae.id.2172.summary.csv")
family.Veillonellaceae.id.2172.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Veillonellaceae.id.2172.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Veillonellaceae.id.2172.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Veillonellaceae.id.2172.summary
)
family.Veillonellaceae.id.2172.summary_res_1 <- mr(family.Veillonellaceae.id.2172.summary_data_1)
fwrite(family.Veillonellaceae.id.2172.summary_res_1,file = "family.Veillonellaceae.id.2172.summary_res_1.csv")
fwrite(family.Veillonellaceae.id.2172.summary_data_1,file = "family.Veillonellaceae.id.2172.summary_data_1.csv")
#family.Verrucomicrobiaceae.id.4036.summary.txt
family.Verrucomicrobiaceae.id.4036.summary<-fread("family.Verrucomicrobiaceae.id.4036.summary.txt")
family.Verrucomicrobiaceae.id.4036.summary<-rename(family.Verrucomicrobiaceae.id.4036.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Verrucomicrobiaceae.id.4036.summary,file = "family.Verrucomicrobiaceae.id.4036.summary.csv")
family.Verrucomicrobiaceae.id.4036.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Verrucomicrobiaceae.id.4036.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Verrucomicrobiaceae.id.4036.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Verrucomicrobiaceae.id.4036.summary
)
family.Verrucomicrobiaceae.id.4036.summary_res_1 <- mr(family.Verrucomicrobiaceae.id.4036.summary_data_1)
fwrite(family.Verrucomicrobiaceae.id.4036.summary_res_1,file = "family.Verrucomicrobiaceae.id.4036.summary_res_1.csv")
fwrite(family.Verrucomicrobiaceae.id.4036.summary_data_1,file = "family.Verrucomicrobiaceae.id.4036.summary_data_1.csv")
#family.Victivallaceae.id.2255.summary.txt
family.Victivallaceae.id.2255.summary<-fread("family.Victivallaceae.id.2255.summary.txt")
family.Victivallaceae.id.2255.summary<-rename(family.Victivallaceae.id.2255.summary, 
                                                id.outcome = bac,
                                                beta.outcome = beta,
                                                se.outcome = SE,
                                                effect_allele.outcome = eff.allele,
                                                other_allele.outcome = ref.allele,
                                                SNP = rsID)
fwrite(family.Victivallaceae.id.2255.summary,file = "family.Victivallaceae.id.2255.summary.csv")
family.Victivallaceae.id.2255.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "family.Victivallaceae.id.2255.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
family.Victivallaceae.id.2255.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = family.Victivallaceae.id.2255.summary
)
family.Victivallaceae.id.2255.summary_res_1 <- mr(family.Victivallaceae.id.2255.summary_data_1)
fwrite(family.Victivallaceae.id.2255.summary_res_1,file = "family.Victivallaceae.id.2255.summary_res_1.csv")
fwrite(family.Victivallaceae.id.2255.summary_data_1,file = "family.Victivallaceae.id.2255.summary_data_1.csv")
#G:/gut/MiBioGen_QmbQTL_summary_genus/MiBioGen_QmbQTL_summary_genus
setwd("G:/gut/MiBioGen_QmbQTL_summary_genus/MiBioGen_QmbQTL_summary_genus")
#genus..Clostridiuminnocuumgroup.id.14397.summary.txt
genus..Clostridiuminnocuumgroup.id.14397.summary<-fread("genus..Clostridiuminnocuumgroup.id.14397.summary.txt")
genus..Clostridiuminnocuumgroup.id.14397.summary<-rename(genus..Clostridiuminnocuumgroup.id.14397.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Clostridiuminnocuumgroup.id.14397.summary,file = "genus..Clostridiuminnocuumgroup.id.14397.summary.csv")
genus..Clostridiuminnocuumgroup.id.14397.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Clostridiuminnocuumgroup.id.14397.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Clostridiuminnocuumgroup.id.14397.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Clostridiuminnocuumgroup.id.14397.summary
)
genus..Clostridiuminnocuumgroup.id.14397.summary_res_1 <- mr(genus..Clostridiuminnocuumgroup.id.14397.summary_data_1)
or_genus..Clostridiuminnocuumgroup.id.14397.summary_res_1 <- generate_odds_ratios(genus..Clostridiuminnocuumgroup.id.14397.summary_res_1)
fwrite(genus..Clostridiuminnocuumgroup.id.14397.summary_res_1,file = "genus..Clostridiuminnocuumgroup.id.14397.summary_res_1.csv")
fwrite(genus..Clostridiuminnocuumgroup.id.14397.summary_data_1,file = "genus..Clostridiuminnocuumgroup.id.14397.summary_data_1.csv")
ps_genus..Clostridiuminnocuumgroup.id.14397.summary_data_1 <- mr_scatter_plot(genus..Clostridiuminnocuumgroup.id.14397.summary_res_1, genus..Clostridiuminnocuumgroup.id.14397.summary_data_1)
ps_genus..Clostridiuminnocuumgroup.id.14397.summary_data_1[[1]]
heterogeneity_genus..Clostridiuminnocuumgroup.id.14397.summary_data_1<-mr_heterogeneity(genus..Clostridiuminnocuumgroup.id.14397.summary_data_1)
fwrite(heterogeneity_genus..Clostridiuminnocuumgroup.id.14397.summary_data_1,file = "heterogeneity_genus..Clostridiuminnocuumgroup.id.14397.summary_data_1.csv")
genus..Clostridiuminnocuumgroup.id.14397.summary_data_1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                                   SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                                   data = genus..Clostridiuminnocuumgroup.id.14397.summary_data_1, 
                                                                   NbDistribution = 1000,  SignifThreshold = 0.05)
leaveoneout_genus..Clostridiuminnocuumgroup.id.14397.summary_data_1 <- mr_leaveoneout(genus..Clostridiuminnocuumgroup.id.14397.summary_data_1)
p_genus..Clostridiuminnocuumgroup.id.14397.summary_data_1 <- mr_leaveoneout_plot(leaveoneout_genus..Clostridiuminnocuumgroup.id.14397.summary_data_1)
p_genus..Clostridiuminnocuumgroup.id.14397.summary_data_1[[1]]
genus..Clostridiuminnocuumgroup.id.14397.summary_data_1_presso
#genes..Eubacteriumbrachygroup.id.11296.summary.txt
genus..Eubacteriumbrachygroup.id.11296.summary<-fread("genus..Eubacteriumbrachygroup.id.11296.summary.txt")
genus..Eubacteriumbrachygroup.id.11296.summary<-rename(genus..Eubacteriumbrachygroup.id.11296.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Eubacteriumbrachygroup.id.11296.summary,file = "genus..Eubacteriumbrachygroup.id.11296.summary.csv")
genus..Eubacteriumbrachygroup.id.11296.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Eubacteriumbrachygroup.id.11296.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Eubacteriumbrachygroup.id.11296.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Eubacteriumbrachygroup.id.11296.summary
)
genus..Eubacteriumbrachygroup.id.11296.summary_res_1 <- mr(genus..Eubacteriumbrachygroup.id.11296.summary_data_1)
fwrite(genus..Eubacteriumbrachygroup.id.11296.summary_res_1,file = "genus..Eubacteriumbrachygroup.id.11296.summary_res_1.csv")
fwrite(genus..Eubacteriumbrachygroup.id.11296.summary_data_1,file = "genus..Eubacteriumbrachygroup.id.11296.summary_data_1.csv")
#genes..Eubacteriumcoprostanoligenesgroup.id.11375.summary.txt
genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary<-fread("genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary.txt")
genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary<-rename(genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary,file = "genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary.csv")
genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary
)
genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary_res_1 <- mr(genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary_data_1)
fwrite(genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary_res_1,file = "genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary_res_1.csv")
fwrite(genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary_data_1,file = "genus..Eubacteriumcoprostanoligenesgroup.id.11375.summary_data_1.csv")
#genes..Eubacteriumeligensgroup.id.14372.summary.txt
genus..Eubacteriumeligensgroup.id.14372.summary<-fread("genus..Eubacteriumeligensgroup.id.14372.summary.txt")
genus..Eubacteriumeligensgroup.id.14372.summary<-rename(genus..Eubacteriumeligensgroup.id.14372.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Eubacteriumeligensgroup.id.14372.summary,file = "genus..Eubacteriumeligensgroup.id.14372.summary.csv")
genus..Eubacteriumeligensgroup.id.14372.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Eubacteriumeligensgroup.id.14372.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Eubacteriumeligensgroup.id.14372.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Eubacteriumeligensgroup.id.14372.summary
)
genus..Eubacteriumeligensgroup.id.14372.summary_res_1 <- mr(genus..Eubacteriumeligensgroup.id.14372.summary_data_1)
fwrite(genus..Eubacteriumeligensgroup.id.14372.summary_res_1,file = "genus..Eubacteriumeligensgroup.id.14372.summary_res_1.csv")
fwrite(genus..Eubacteriumeligensgroup.id.14372.summary_data_1,file = "genus..Eubacteriumeligensgroup.id.14372.summary_data_1.csv")
#genes..Eubacteriumfissicatenagroup.id.14373.summary.txt
genus..Eubacteriumfissicatenagroup.id.14373.summary<-fread("genus..Eubacteriumfissicatenagroup.id.14373.summary.txt")
genus..Eubacteriumfissicatenagroup.id.14373.summary<-rename(genus..Eubacteriumfissicatenagroup.id.14373.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Eubacteriumfissicatenagroup.id.14373.summary,file = "genus..Eubacteriumfissicatenagroup.id.14373.summary.csv")
genus..Eubacteriumfissicatenagroup.id.14373.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Eubacteriumfissicatenagroup.id.14373.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Eubacteriumfissicatenagroup.id.14373.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Eubacteriumfissicatenagroup.id.14373.summary
)
genus..Eubacteriumfissicatenagroup.id.14373.summary_res_1 <- mr(genus..Eubacteriumfissicatenagroup.id.14373.summary_data_1)
fwrite(genus..Eubacteriumfissicatenagroup.id.14373.summary_res_1,file = "genus..Eubacteriumfissicatenagroup.id.14373.summary_res_1.csv")
fwrite(genus..Eubacteriumfissicatenagroup.id.14373.summary_data_1,file = "genus..Eubacteriumfissicatenagroup.id.14373.summary_data_1.csv")
#genes..Eubacteriumhalliigroup.id.11338.summary.txt
genus..Eubacteriumhalliigroup.id.11338.summary<-fread("genus..Eubacteriumhalliigroup.id.11338.summary.txt")
genus..Eubacteriumhalliigroup.id.11338.summary<-rename(genus..Eubacteriumhalliigroup.id.11338.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Eubacteriumhalliigroup.id.11338.summary,file = "genus..Eubacteriumhalliigroup.id.11338.summary.csv")
genus..Eubacteriumhalliigroup.id.11338.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Eubacteriumhalliigroup.id.11338.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Eubacteriumhalliigroup.id.11338.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Eubacteriumhalliigroup.id.11338.summary
)
genus..Eubacteriumhalliigroup.id.11338.summary_res_1 <- mr(genus..Eubacteriumhalliigroup.id.11338.summary_data_1)
fwrite(genus..Eubacteriumhalliigroup.id.11338.summary_res_1,file = "genus..Eubacteriumhalliigroup.id.11338.summary_res_1.csv")
fwrite(genus..Eubacteriumhalliigroup.id.11338.summary_data_1,file = "genus..Eubacteriumhalliigroup.id.11338.summary_data_1.csv")
#genes..Eubacteriumnodatumgroup.id.11297.summary.txt
genus..Eubacteriumnodatumgroup.id.11297.summary<-fread("genus..Eubacteriumnodatumgroup.id.11297.summary.txt")
genus..Eubacteriumnodatumgroup.id.11297.summary<-rename(genus..Eubacteriumnodatumgroup.id.11297.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Eubacteriumnodatumgroup.id.11297.summary,file = "genus..Eubacteriumnodatumgroup.id.11297.summary.csv")
genus..Eubacteriumnodatumgroup.id.11297.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Eubacteriumnodatumgroup.id.11297.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Eubacteriumnodatumgroup.id.11297.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Eubacteriumnodatumgroup.id.11297.summary
)
genus..Eubacteriumnodatumgroup.id.11297.summary_res_1 <- mr(genus..Eubacteriumnodatumgroup.id.11297.summary_data_1)
fwrite(genus..Eubacteriumnodatumgroup.id.11297.summary_res_1,file = "genus..Eubacteriumnodatumgroup.id.11297.summary_res_1.csv")
fwrite(genus..Eubacteriumnodatumgroup.id.11297.summary_data_1,file = "genus..Eubacteriumnodatumgroup.id.11297.summary_data_1.csv")
#genes..Eubacteriumoxidoreducensgroup.id.11339.summary.txt
genus..Eubacteriumoxidoreducensgroup.id.11339.summary<-fread("genus..Eubacteriumoxidoreducensgroup.id.11339.summary.txt")
genus..Eubacteriumoxidoreducensgroup.id.11339.summary<-rename(genus..Eubacteriumoxidoreducensgroup.id.11339.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Eubacteriumoxidoreducensgroup.id.11339.summary,file = "genus..Eubacteriumoxidoreducensgroup.id.11339.summary.csv")
genus..Eubacteriumoxidoreducensgroup.id.11339.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Eubacteriumoxidoreducensgroup.id.11339.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Eubacteriumoxidoreducensgroup.id.11339.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Eubacteriumoxidoreducensgroup.id.11339.summary
)
genus..Eubacteriumoxidoreducensgroup.id.11339.summary_res_1 <- mr(genus..Eubacteriumoxidoreducensgroup.id.11339.summary_data_1)
fwrite(genus..Eubacteriumoxidoreducensgroup.id.11339.summary_res_1,file = "genus..Eubacteriumoxidoreducensgroup.id.11339.summary_res_1.csv")
fwrite(genus..Eubacteriumoxidoreducensgroup.id.11339.summary_data_1,file = "genus..Eubacteriumoxidoreducensgroup.id.11339.summary_data_1.csv")
#genes..Eubacteriumrectalegroup.id.14374.summary.txt
genus..Eubacteriumrectalegroup.id.14374.summary<-fread("genus..Eubacteriumrectalegroup.id.14374.summary.txt")
genus..Eubacteriumrectalegroup.id.14374.summary<-rename(genus..Eubacteriumrectalegroup.id.14374.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Eubacteriumrectalegroup.id.14374.summary,file = "genus..Eubacteriumrectalegroup.id.14374.summary.csv")
genus..Eubacteriumrectalegroup.id.14374.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Eubacteriumrectalegroup.id.14374.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Eubacteriumrectalegroup.id.14374.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Eubacteriumrectalegroup.id.14374.summary
)
genus..Eubacteriumrectalegroup.id.14374.summary_res_1 <- mr(genus..Eubacteriumrectalegroup.id.14374.summary_data_1)
fwrite(genus..Eubacteriumrectalegroup.id.14374.summary_res_1,file = "genus..Eubacteriumrectalegroup.id.14374.summary_res_1.csv")
fwrite(genus..Eubacteriumrectalegroup.id.14374.summary_data_1,file = "genus..Eubacteriumrectalegroup.id.14374.summary_data_1.csv")
#genes..Eubacteriumruminantiumgroup.id.11340.summary.txt
genus..Eubacteriumruminantiumgroup.id.11340.summary<-fread("genus..Eubacteriumruminantiumgroup.id.11340.summary.txt")
genus..Eubacteriumruminantiumgroup.id.11340.summary<-rename(genus..Eubacteriumruminantiumgroup.id.11340.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Eubacteriumruminantiumgroup.id.11340.summary,file = "genus..Eubacteriumruminantiumgroup.id.11340.summary.csv")
genus..Eubacteriumruminantiumgroup.id.11340.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Eubacteriumruminantiumgroup.id.11340.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Eubacteriumruminantiumgroup.id.11340.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Eubacteriumruminantiumgroup.id.11340.summary
)
genus..Eubacteriumruminantiumgroup.id.11340.summary_res_1 <- mr(genus..Eubacteriumruminantiumgroup.id.11340.summary_data_1)
fwrite(genus..Eubacteriumruminantiumgroup.id.11340.summary_res_1,file = "genus..Eubacteriumruminantiumgroup.id.11340.summary_res_1.csv")
fwrite(genus..Eubacteriumruminantiumgroup.id.11340.summary_data_1,file = "genus..Eubacteriumruminantiumgroup.id.11340.summary_data_1.csv")
#genes..Eubacteriumventriosumgroup.id.11341.summary.txt
genus..Eubacteriumventriosumgroup.id.11341.summary<-fread("genus..Eubacteriumventriosumgroup.id.11341.summary.txt")
genus..Eubacteriumventriosumgroup.id.11341.summary<-rename(genus..Eubacteriumventriosumgroup.id.11341.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Eubacteriumventriosumgroup.id.11341.summary,file = "genus..Eubacteriumventriosumgroup.id.11341.summary.csv")
genus..Eubacteriumventriosumgroup.id.11341.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Eubacteriumventriosumgroup.id.11341.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Eubacteriumventriosumgroup.id.11341.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Eubacteriumventriosumgroup.id.11341.summary
)
genus..Eubacteriumventriosumgroup.id.11341.summary_res_1 <- mr(genus..Eubacteriumventriosumgroup.id.11341.summary_data_1)
fwrite(genus..Eubacteriumventriosumgroup.id.11341.summary_res_1,file = "genus..Eubacteriumventriosumgroup.id.11341.summary_res_1.csv")
fwrite(genus..Eubacteriumventriosumgroup.id.11341.summary_data_1,file = "genus..Eubacteriumventriosumgroup.id.11341.summary_data_1.csv")
#genes..Eubacteriumxylanophilumgroup.id.14375.summary.txt
genus..Eubacteriumxylanophilumgroup.id.14375.summary<-fread("genus..Eubacteriumxylanophilumgroup.id.14375.summary.txt")
genus..Eubacteriumxylanophilumgroup.id.14375.summary<-rename(genus..Eubacteriumxylanophilumgroup.id.14375.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Eubacteriumxylanophilumgroup.id.14375.summary,file = "genus..Eubacteriumxylanophilumgroup.id.14375.summary.csv")
genus..Eubacteriumxylanophilumgroup.id.14375.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Eubacteriumxylanophilumgroup.id.14375.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Eubacteriumxylanophilumgroup.id.14375.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Eubacteriumxylanophilumgroup.id.14375.summary
)
genus..Eubacteriumxylanophilumgroup.id.14375.summary_res_1 <- mr(genus..Eubacteriumxylanophilumgroup.id.14375.summary_data_1)
fwrite(genus..Eubacteriumxylanophilumgroup.id.14375.summary_res_1,file = "genus..Eubacteriumxylanophilumgroup.id.14375.summary_res_1.csv")
fwrite(genus..Eubacteriumxylanophilumgroup.id.14375.summary_data_1,file = "genus..Eubacteriumxylanophilumgroup.id.14375.summary_data_1.csv")
#genes..Ruminococcusgauvreauiigroup.id.11342.summary.txt
genus..Ruminococcusgauvreauiigroup.id.11342.summary<-fread("genus..Ruminococcusgauvreauiigroup.id.11342.summary.txt")
genus..Ruminococcusgauvreauiigroup.id.11342.summary<-rename(genus..Ruminococcusgauvreauiigroup.id.11342.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Ruminococcusgauvreauiigroup.id.11342.summary,file = "genus..Ruminococcusgauvreauiigroup.id.11342.summary.csv")
genus..Ruminococcusgauvreauiigroup.id.11342.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Ruminococcusgauvreauiigroup.id.11342.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Ruminococcusgauvreauiigroup.id.11342.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Ruminococcusgauvreauiigroup.id.11342.summary
)
genus..Ruminococcusgauvreauiigroup.id.11342.summary_res_1 <- mr(genus..Ruminococcusgauvreauiigroup.id.11342.summary_data_1)
fwrite(genus..Ruminococcusgauvreauiigroup.id.11342.summary_res_1,file = "genus..Ruminococcusgauvreauiigroup.id.11342.summary_res_1.csv")
fwrite(genus..Ruminococcusgauvreauiigroup.id.11342.summary_data_1,file = "genus..Ruminococcusgauvreauiigroup.id.11342.summary_data_1.csv")
#genes..Ruminococcusgnavusgroup.id.14376.summary.txt
genus..Ruminococcusgnavusgroup.id.14376.summary<-fread("genus..Ruminococcusgnavusgroup.id.14376.summary.txt")
genus..Ruminococcusgnavusgroup.id.14376.summary<-rename(genus..Ruminococcusgnavusgroup.id.14376.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Ruminococcusgnavusgroup.id.14376.summary,file = "genus..Ruminococcusgnavusgroup.id.14376.summary.csv")
genus..Ruminococcusgnavusgroup.id.14376.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Ruminococcusgnavusgroup.id.14376.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Ruminococcusgnavusgroup.id.14376.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Ruminococcusgnavusgroup.id.14376.summary
)
genus..Ruminococcusgnavusgroup.id.14376.summary_res_1 <- mr(genus..Ruminococcusgnavusgroup.id.14376.summary_data_1)
fwrite(genus..Ruminococcusgnavusgroup.id.14376.summary_res_1,file = "genus..Ruminococcusgnavusgroup.id.14376.summary_res_1.csv")
fwrite(genus..Ruminococcusgnavusgroup.id.14376.summary_data_1,file = "genus..Ruminococcusgnavusgroup.id.14376.summary_data_1.csv")
#genes..Ruminococcustorquesgroup.id.14377.summary.txt
genus..Ruminococcustorquesgroup.id.14377.summary<-fread("genus..Ruminococcustorquesgroup.id.14377.summary.txt")
genus..Ruminococcustorquesgroup.id.14377.summary<-rename(genus..Ruminococcustorquesgroup.id.14377.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus..Ruminococcustorquesgroup.id.14377.summary,file = "genus..Ruminococcustorquesgroup.id.14377.summary.csv")
genus..Ruminococcustorquesgroup.id.14377.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus..Ruminococcustorquesgroup.id.14377.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus..Ruminococcustorquesgroup.id.14377.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus..Ruminococcustorquesgroup.id.14377.summary
)
genus..Ruminococcustorquesgroup.id.14377.summary_res_1 <- mr(genus..Ruminococcustorquesgroup.id.14377.summary_data_1)
fwrite(genus..Ruminococcustorquesgroup.id.14377.summary_res_1,file = "genus..Ruminococcustorquesgroup.id.14377.summary_res_1.csv")
fwrite(genus..Ruminococcustorquesgroup.id.14377.summary_data_1,file = "genus..Ruminococcustorquesgroup.id.14377.summary_data_1.csv")
#genus.Actinomyces.id.423.summary.txt
genus.Actinomyces.id.423.summary<-fread("genus.Actinomyces.id.423.summary.txt")
genus.Actinomyces.id.423.summary<-rename(genus.Actinomyces.id.423.summary, 
                                              id.outcome = bac,
                                              beta.outcome = beta,
                                              se.outcome = SE,
                                              effect_allele.outcome = eff.allele,
                                              other_allele.outcome = ref.allele,
                                              SNP = rsID)
fwrite(genus.Actinomyces.id.423.summary,file = "genus.Actinomyces.id.423.summary.csv")
genus.Actinomyces.id.423.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Actinomyces.id.423.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Actinomyces.id.423.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Actinomyces.id.423.summary
)
genus.Actinomyces.id.423.summary_res_1 <- mr(genus.Actinomyces.id.423.summary_data_1)
fwrite(genus.Actinomyces.id.423.summary_res_1,file = "genus.Actinomyces.id.423.summary_res_1.csv")
fwrite(genus.Actinomyces.id.423.summary_data_1,file = "genus.Actinomyces.id.423.summary_data_1.csv")
#genes.Adlercreutzia.id.812.summary.txt
genus.Adlercreutzia.id.812.summary<-fread("genus.Adlercreutzia.id.812.summary.txt")
genus.Adlercreutzia.id.812.summary<-rename(genus.Adlercreutzia.id.812.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Adlercreutzia.id.812.summary,file = "genus.Adlercreutzia.id.812.summary.csv")
genus.Adlercreutzia.id.812.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Adlercreutzia.id.812.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Adlercreutzia.id.812.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Adlercreutzia.id.812.summary
)
genus.Adlercreutzia.id.812.summary_res_1 <- mr(genus.Adlercreutzia.id.812.summary_data_1)
fwrite(genus.Adlercreutzia.id.812.summary_res_1,file = "genus.Adlercreutzia.id.812.summary_res_1.csv")
fwrite(genus.Adlercreutzia.id.812.summary_data_1,file = "genus.Adlercreutzia.id.812.summary_data_1.csv")
#genes.Akkermansia.id.4037.summary.txt
genus.Akkermansia.id.4037.summary<-fread("genus.Akkermansia.id.4037.summary.txt")
genus.Akkermansia.id.4037.summary<-rename(genus.Akkermansia.id.4037.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Akkermansia.id.4037.summary,file = "genus.Akkermansia.id.4037.summary.csv")
genus.Akkermansia.id.4037.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Akkermansia.id.4037.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Akkermansia.id.4037.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Akkermansia.id.4037.summary
)
genus.Akkermansia.id.4037.summary_res_1 <- mr(genus.Akkermansia.id.4037.summary_data_1)
fwrite(genus.Akkermansia.id.4037.summary_res_1,file = "genus.Akkermansia.id.4037.summary_res_1.csv")
fwrite(genus.Akkermansia.id.4037.summary_data_1,file = "genus.Akkermansia.id.4037.summary_data_1.csv")
#genes.Alistipes.id.968.summary.txt
genus.Alistipes.id.968.summary<-fread("genus.Alistipes.id.968.summary.txt")
genus.Alistipes.id.968.summary<-rename(genus.Alistipes.id.968.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Alistipes.id.968.summary,file = "genus.Alistipes.id.968.summary.csv")
genus.Alistipes.id.968.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Alistipes.id.968.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Alistipes.id.968.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Alistipes.id.968.summary
)
genus.Alistipes.id.968.summary_res_1 <- mr(genus.Alistipes.id.968.summary_data_1)
fwrite(genus.Alistipes.id.968.summary_res_1,file = "genus.Alistipes.id.968.summary_res_1.csv")
fwrite(genus.Alistipes.id.968.summary_data_1,file = "genus.Alistipes.id.968.summary_data_1.csv")
#genes.Allisonella.id.2174.summary.txt
genus.Allisonella.id.2174.summary<-fread("genus.Allisonella.id.2174.summary.txt")
genus.Allisonella.id.2174.summary<-rename(genus.Allisonella.id.2174.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Allisonella.id.2174.summary,file = "genus.Allisonella.id.2174.summary.csv")
genus.Allisonella.id.2174.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Allisonella.id.2174.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Allisonella.id.2174.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Allisonella.id.2174.summary
)
genus.Allisonella.id.2174.summary_res_1 <- mr(genus.Allisonella.id.2174.summary_data_1)
fwrite(genus.Allisonella.id.2174.summary_res_1,file = "genus.Allisonella.id.2174.summary_res_1.csv")
fwrite(genus.Allisonella.id.2174.summary_data_1,file = "genus.Allisonella.id.2174.summary_data_1.csv")
#genes.Alloprevotella.id.961.summary.txt
genus.Alloprevotella.id.961.summary<-fread("genus.Alloprevotella.id.961.summary.txt")
genus.Alloprevotella.id.961.summary<-rename(genus.Alloprevotella.id.961.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Alloprevotella.id.961.summary,file = "genus.Alloprevotella.id.961.summary.csv")
genus.Alloprevotella.id.961.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Alloprevotella.id.961.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Alloprevotella.id.961.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Alloprevotella.id.961.summary
)
genus.Alloprevotella.id.961.summary_res_1 <- mr(genus.Alloprevotella.id.961.summary_data_1)
fwrite(genus.Alloprevotella.id.961.summary_res_1,file = "genus.Alloprevotella.id.961.summary_res_1.csv")
fwrite(genus.Alloprevotella.id.961.summary_data_1,file = "genus.Alloprevotella.id.961.summary_data_1.csv")
#genes.Anaerofilum.id.2053.summary.txt
genus.Anaerofilum.id.2053.summary<-fread("genus.Anaerofilum.id.2053.summary.txt")
genus.Anaerofilum.id.2053.summary<-rename(genus.Anaerofilum.id.2053.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Anaerofilum.id.2053.summary,file = "genus.Anaerofilum.id.2053.summary.csv")
genus.Anaerofilum.id.2053.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Anaerofilum.id.2053.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Anaerofilum.id.2053.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Anaerofilum.id.2053.summary
)
genus.Anaerofilum.id.2053.summary_res_1 <- mr(genus.Anaerofilum.id.2053.summary_data_1)
or_genus.Anaerofilum.id.2053.summary_res_1 <- generate_odds_ratios(genus.Anaerofilum.id.2053.summary_res_1)
fwrite(genus.Anaerofilum.id.2053.summary_res_1,file = "genus.Anaerofilum.id.2053.summary_res_1.csv")
fwrite(genus.Anaerofilum.id.2053.summary_data_1,file = "genus.Anaerofilum.id.2053.summary_data_1.csv")
ps_genus.Anaerofilum.id.2053.summary_data_1 <- mr_scatter_plot(genus.Anaerofilum.id.2053.summary_res_1, genus.Anaerofilum.id.2053.summary_data_1)
ps_genus.Anaerofilum.id.2053.summary_data_1[[1]]
heterogeneity_genus.Anaerofilum.id.2053.summary_data_1<-mr_heterogeneity(genus.Anaerofilum.id.2053.summary_data_1)
fwrite(heterogeneity_genus.Anaerofilum.id.2053.summary_data_1,file = "heterogeneity_genus.Anaerofilum.id.2053.summary_data_1.csv")
genus.Anaerofilum.id.2053.summary_data_1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                                   SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                                   data = genus.Anaerofilum.id.2053.summary_data_1, 
                                                                   NbDistribution = 1000,  SignifThreshold = 0.05)
leaveoneout_genus.Anaerofilum.id.2053.summary_data_1 <- mr_leaveoneout(genus.Anaerofilum.id.2053.summary_data_1)
p_genus.Anaerofilum.id.2053.summary_data_1 <- mr_leaveoneout_plot(leaveoneout_genus.Anaerofilum.id.2053.summary_data_1)
p_genus.Anaerofilum.id.2053.summary_data_1[[1]]
genus.Anaerofilum.id.2053.summary_data_1_presso
#genes.Anaerostipes.id.1991.summary.txt
genus.Anaerostipes.id.1991.summary<-fread("genus.Anaerostipes.id.1991.summary.txt")
genus.Anaerostipes.id.1991.summary<-rename(genus.Anaerostipes.id.1991.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Anaerostipes.id.1991.summary,file = "genus.Anaerostipes.id.1991.summary.csv")
genus.Anaerostipes.id.1991.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Anaerostipes.id.1991.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Anaerostipes.id.1991.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Anaerostipes.id.1991.summary
)
genus.Anaerostipes.id.1991.summary_res_1 <- mr(genus.Anaerostipes.id.1991.summary_data_1)
fwrite(genus.Anaerostipes.id.1991.summary_res_1,file = "genus.Anaerostipes.id.1991.summary_res_1.csv")
fwrite(genus.Anaerostipes.id.1991.summary_data_1,file = "genus.Anaerostipes.id.1991.summary_data_1.csv")
#genes.Anaerotruncus.id.2054.summary.txt
genus.Anaerotruncus.id.2054.summary<-fread("genus.Anaerotruncus.id.2054.summary.txt")
genus.Anaerotruncus.id.2054.summary<-rename(genus.Anaerotruncus.id.2054.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Anaerotruncus.id.2054.summary,file = "genus.Anaerotruncus.id.2054.summary.csv")
genus.Anaerotruncus.id.2054.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Anaerotruncus.id.2054.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Anaerotruncus.id.2054.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Anaerotruncus.id.2054.summary
)
genus.Anaerotruncus.id.2054.summary_res_1 <- mr(genus.Anaerotruncus.id.2054.summary_data_1)
fwrite(genus.Anaerotruncus.id.2054.summary_res_1,file = "genus.Anaerotruncus.id.2054.summary_res_1.csv")
fwrite(genus.Anaerotruncus.id.2054.summary_data_1,file = "genus.Anaerotruncus.id.2054.summary_data_1.csv")
#genes.Bacteroides.id.918.summary.txt
genus.Bacteroides.id.918.summary<-fread("genus.Bacteroides.id.918.summary.txt")
genus.Bacteroides.id.918.summary<-rename(genus.Bacteroides.id.918.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Bacteroides.id.918.summary,file = "genus.Bacteroides.id.918.summary.csv")
genus.Bacteroides.id.918.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Bacteroides.id.918.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Bacteroides.id.918.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Bacteroides.id.918.summary
)
genus.Bacteroides.id.918.summary_res_1 <- mr(genus.Bacteroides.id.918.summary_data_1)
fwrite(genus.Bacteroides.id.918.summary_res_1,file = "genus.Bacteroides.id.918.summary_res_1.csv")
fwrite(genus.Bacteroides.id.918.summary_data_1,file = "genus.Bacteroides.id.918.summary_data_1.csv")
#genes.Barnesiella.id.944.summary.txt
genus.Barnesiella.id.944.summary<-fread("genus.Barnesiella.id.944.summary.txt")
genus.Barnesiella.id.944.summary<-rename(genus.Barnesiella.id.944.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Barnesiella.id.944.summary,file = "genus.Barnesiella.id.944.summary.csv")
genus.Barnesiella.id.944.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Barnesiella.id.944.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Barnesiella.id.944.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Barnesiella.id.944.summary
)
genus.Barnesiella.id.944.summary_res_1 <- mr(genus.Barnesiella.id.944.summary_data_1)
fwrite(genus.Barnesiella.id.944.summary_res_1,file = "genus.Barnesiella.id.944.summary_res_1.csv")
fwrite(genus.Barnesiella.id.944.summary_data_1,file = "genus.Barnesiella.id.944.summary_data_1.csv")
#genes.Bifidobacterium.id.436.summary.txt
genus.Bifidobacterium.id.436.summary<-fread("genus.Bifidobacterium.id.436.summary.txt")
genus.Bifidobacterium.id.436.summary<-rename(genus.Bifidobacterium.id.436.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Bifidobacterium.id.436.summary,file = "genus.Bifidobacterium.id.436.summary.csv")
genus.Bifidobacterium.id.436.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Bifidobacterium.id.436.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Bifidobacterium.id.436.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Bifidobacterium.id.436.summary
)
genus.Bifidobacterium.id.436.summary_res_1 <- mr(genus.Bifidobacterium.id.436.summary_data_1)
fwrite(genus.Bifidobacterium.id.436.summary_res_1,file = "genus.Bifidobacterium.id.436.summary_res_1.csv")
fwrite(genus.Bifidobacterium.id.436.summary_data_1,file = "genus.Bifidobacterium.id.436.summary_data_1.csv")
#genes.Bilophila.id.3170.summary.txt
genus.Bilophila.id.3170.summary<-fread("genus.Bilophila.id.3170.summary.txt")
genus.Bilophila.id.3170.summary<-rename(genus.Bilophila.id.3170.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Bilophila.id.3170.summary,file = "genus.Bilophila.id.3170.summary.csv")
genus.Bilophila.id.3170.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Bilophila.id.3170.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Bilophila.id.3170.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Bilophila.id.3170.summary
)
genus.Bilophila.id.3170.summary_res_1 <- mr(genus.Bilophila.id.3170.summary_data_1)
fwrite(genus.Bilophila.id.3170.summary_res_1,file = "genus.Bilophila.id.3170.summary_res_1.csv")
fwrite(genus.Bilophila.id.3170.summary_data_1,file = "genus.Bilophila.id.3170.summary_data_1.csv")
#genes.Blautia.id.1992.summary.txt
genus.Blautia.id.1992.summary<-fread("genus.Blautia.id.1992.summary.txt")
genus.Blautia.id.1992.summary<-rename(genus.Blautia.id.1992.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Blautia.id.1992.summary,file = "genus.Blautia.id.1992.summary.csv")
genus.Blautia.id.1992.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Blautia.id.1992.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Blautia.id.1992.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Blautia.id.1992.summary
)
genus.Blautia.id.1992.summary_res_1 <- mr(genus.Blautia.id.1992.summary_data_1)
fwrite(genus.Blautia.id.1992.summary_res_1,file = "genus.Blautia.id.1992.summary_res_1.csv")
fwrite(genus.Blautia.id.1992.summary_data_1,file = "genus.Blautia.id.1992.summary_data_1.csv")
#genes.Butyricicoccus.id.2055.summary.txt
genus.Butyricicoccus.id.2055.summary<-fread("genus.Butyricicoccus.id.2055.summary.txt")
genus.Butyricicoccus.id.2055.summary<-rename(genus.Butyricicoccus.id.2055.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Butyricicoccus.id.2055.summary,file = "genus.Butyricicoccus.id.2055.summary.csv")
genus.Butyricicoccus.id.2055.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Butyricicoccus.id.2055.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Butyricicoccus.id.2055.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Butyricicoccus.id.2055.summary
)
genus.Butyricicoccus.id.2055.summary_res_1 <- mr(genus.Butyricicoccus.id.2055.summary_data_1)
fwrite(genus.Butyricicoccus.id.2055.summary_res_1,file = "genus.Butyricicoccus.id.2055.summary_res_1.csv")
fwrite(genus.Butyricicoccus.id.2055.summary_data_1,file = "genus.Butyricicoccus.id.2055.summary_data_1.csv")
#genes.Butyricimonas.id.945.summary.txt
genus.Butyricimonas.id.945.summary<-fread("genus.Butyricimonas.id.945.summary.txt")
genus.Butyricimonas.id.945.summary<-rename(genus.Butyricimonas.id.945.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Butyricimonas.id.945.summary,file = "genus.Butyricimonas.id.945.summary.csv")
genus.Butyricimonas.id.945.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Butyricimonas.id.945.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Butyricimonas.id.945.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Butyricimonas.id.945.summary
)
genus.Butyricimonas.id.945.summary_res_1 <- mr(genus.Butyricimonas.id.945.summary_data_1)
fwrite(genus.Butyricimonas.id.945.summary_res_1,file = "genus.Butyricimonas.id.945.summary_res_1.csv")
fwrite(genus.Butyricimonas.id.945.summary_data_1,file = "genus.Butyricimonas.id.945.summary_data_1.csv")
#genes.Butyrivibrio.id.1993.summary.txt
genus.Butyrivibrio.id.1993.summary<-fread("genus.Butyrivibrio.id.1993.summary.txt")
genus.Butyrivibrio.id.1993.summary<-rename(genus.Butyrivibrio.id.1993.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Butyrivibrio.id.1993.summary,file = "genus.Butyrivibrio.id.1993.summary.csv")
genus.Butyrivibrio.id.1993.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Butyrivibrio.id.1993.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Butyrivibrio.id.1993.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Butyrivibrio.id.1993.summary
)
genus.Butyrivibrio.id.1993.summary_res_1 <- mr(genus.Butyrivibrio.id.1993.summary_data_1)
fwrite(genus.Butyrivibrio.id.1993.summary_res_1,file = "genus.Butyrivibrio.id.1993.summary_res_1.csv")
fwrite(genus.Butyrivibrio.id.1993.summary_data_1,file = "genus.Butyrivibrio.id.1993.summary_data_1.csv")
#genes.CandidatusSoleaferrea.id.11350.summary.txt
genus.CandidatusSoleaferrea.id.11350.summary<-fread("genus.CandidatusSoleaferrea.id.11350.summary.txt")
genus.CandidatusSoleaferrea.id.11350.summary<-rename(genus.CandidatusSoleaferrea.id.11350.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.CandidatusSoleaferrea.id.11350.summary,file = "genus.CandidatusSoleaferrea.id.11350.summary.csv")
genus.CandidatusSoleaferrea.id.11350.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.CandidatusSoleaferrea.id.11350.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.CandidatusSoleaferrea.id.11350.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.CandidatusSoleaferrea.id.11350.summary
)
genus.CandidatusSoleaferrea.id.11350.summary_res_1 <- mr(genus.CandidatusSoleaferrea.id.11350.summary_data_1)
fwrite(genus.CandidatusSoleaferrea.id.11350.summary_res_1,file = "genus.CandidatusSoleaferrea.id.11350.summary_res_1.csv")
fwrite(genus.CandidatusSoleaferrea.id.11350.summary_data_1,file = "genus.CandidatusSoleaferrea.id.11350.summary_data_1.csv")
#genes.Catenibacterium.id.2153.summary.txt
genus.Catenibacterium.id.2153.summary<-fread("genus.Catenibacterium.id.2153.summary.txt")
genus.Catenibacterium.id.2153.summary<-rename(genus.Catenibacterium.id.2153.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Catenibacterium.id.2153.summary,file = "genus.Catenibacterium.id.2153.summary.csv")
genus.Catenibacterium.id.2153.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Catenibacterium.id.2153.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Catenibacterium.id.2153.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Catenibacterium.id.2153.summary
)
genus.Catenibacterium.id.2153.summary_res_1 <- mr(genus.Catenibacterium.id.2153.summary_data_1)
fwrite(genus.Catenibacterium.id.2153.summary_res_1,file = "genus.Catenibacterium.id.2153.summary_res_1.csv")
fwrite(genus.Catenibacterium.id.2153.summary_data_1,file = "genus.Catenibacterium.id.2153.summary_data_1.csv")
#genes.ChristensenellaceaeR.7group.id.11283.summary.txt
genus.ChristensenellaceaeR.7group.id.11283.summary<-fread("genus.ChristensenellaceaeR.7group.id.11283.summary.txt")
genus.ChristensenellaceaeR.7group.id.11283.summary<-rename(genus.ChristensenellaceaeR.7group.id.11283.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.ChristensenellaceaeR.7group.id.11283.summary,file = "genus.ChristensenellaceaeR.7group.id.11283.summary.csv")
genus.ChristensenellaceaeR.7group.id.11283.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.ChristensenellaceaeR.7group.id.11283.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.ChristensenellaceaeR.7group.id.11283.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.ChristensenellaceaeR.7group.id.11283.summary
)
genus.ChristensenellaceaeR.7group.id.11283.summary_res_1 <- mr(genus.ChristensenellaceaeR.7group.id.11283.summary_data_1)
fwrite(genus.ChristensenellaceaeR.7group.id.11283.summary_res_1,file = "genus.ChristensenellaceaeR.7group.id.11283.summary_res_1.csv")
fwrite(genus.ChristensenellaceaeR.7group.id.11283.summary_data_1,file = "genus.ChristensenellaceaeR.7group.id.11283.summary_data_1.csv")
#genes.Clostridiumsensustricto1.id.1873.summary.txt
genus.Clostridiumsensustricto1.id.1873.summary<-fread("genus.Clostridiumsensustricto1.id.1873.summary.txt")
genus.Clostridiumsensustricto1.id.1873.summary<-rename(genus.Clostridiumsensustricto1.id.1873.summary, 
                                         id.outcome = bac,
                                         beta.outcome = beta,
                                         se.outcome = SE,
                                         effect_allele.outcome = eff.allele,
                                         other_allele.outcome = ref.allele,
                                         SNP = rsID)
fwrite(genus.Clostridiumsensustricto1.id.1873.summary,file = "genus.Clostridiumsensustricto1.id.1873.summary.csv")
genus.Clostridiumsensustricto1.id.1873.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Clostridiumsensustricto1.id.1873.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Clostridiumsensustricto1.id.1873.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Clostridiumsensustricto1.id.1873.summary
)
genus.Clostridiumsensustricto1.id.1873.summary_res_1 <- mr(genus.Clostridiumsensustricto1.id.1873.summary_data_1)
fwrite(genus.Clostridiumsensustricto1.id.1873.summary_res_1,file = "genus.Clostridiumsensustricto1.id.1873.summary_res_1.csv")
fwrite(genus.Clostridiumsensustricto1.id.1873.summary_data_1,file = "genus.Clostridiumsensustricto1.id.1873.summary_data_1.csv")
#genus.Collinsella.id.815.summary.txt
genus.Collinsella.id.815.summary<-fread("genus.Collinsella.id.815.summary.txt")
genus.Collinsella.id.815.summary<-rename(genus.Collinsella.id.815.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Collinsella.id.815.summary,file = "genus.Collinsella.id.815.summary.csv")
genus.Collinsella.id.815.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Collinsella.id.815.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Collinsella.id.815.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Collinsella.id.815.summary
)
genus.Collinsella.id.815.summary_res_1 <- mr(genus.Collinsella.id.815.summary_data_1)
fwrite(genus.Collinsella.id.815.summary_res_1,file = "genus.Collinsella.id.815.summary_res_1.csv")
fwrite(genus.Collinsella.id.815.summary_data_1,file = "genus.Collinsella.id.815.summary_data_1.csv")
#genus.Coprobacter.id.949.summary.txt
genus.Coprobacter.id.949.summary<-fread("genus.Coprobacter.id.949.summary.txt")
genus.Coprobacter.id.949.summary<-rename(genus.Coprobacter.id.949.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Coprobacter.id.949.summary,file = "genus.Coprobacter.id.949.summary.csv")
genus.Coprobacter.id.949.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Coprobacter.id.949.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Coprobacter.id.949.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Coprobacter.id.949.summary
)
genus.Coprobacter.id.949.summary_res_1 <- mr(genus.Coprobacter.id.949.summary_data_1)
fwrite(genus.Coprobacter.id.949.summary_res_1,file = "genus.Coprobacter.id.949.summary_res_1.csv")
fwrite(genus.Coprobacter.id.949.summary_data_1,file = "genus.Coprobacter.id.949.summary_data_1.csv")
#genus.Coprococcus1.id.11301.summary.txt
genus.Coprococcus1.id.11301.summary<-fread("genus.Coprococcus1.id.11301.summary.txt")
genus.Coprococcus1.id.11301.summary<-rename(genus.Coprococcus1.id.11301.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Coprococcus1.id.11301.summary,file = "genus.Coprococcus1.id.11301.summary.csv")
genus.Coprococcus1.id.11301.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Coprococcus1.id.11301.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Coprococcus1.id.11301.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Coprococcus1.id.11301.summary
)
genus.Coprococcus1.id.11301.summary_res_1 <- mr(genus.Coprococcus1.id.11301.summary_data_1)
fwrite(genus.Coprococcus1.id.11301.summary_res_1,file = "genus.Coprococcus1.id.11301.summary_res_1.csv")
fwrite(genus.Coprococcus1.id.11301.summary_data_1,file = "genus.Coprococcus1.id.11301.summary_data_1.csv")
#genus.Coprococcus2.id.11302.summary.txt
genus.Coprococcus2.id.11302.summary<-fread("genus.Coprococcus2.id.11302.summary.txt")
genus.Coprococcus2.id.11302.summary<-rename(genus.Coprococcus2.id.11302.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Coprococcus2.id.11302.summary,file = "genus.Coprococcus2.id.11302.summary.csv")
genus.Coprococcus2.id.11302.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Coprococcus2.id.11302.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Coprococcus2.id.11302.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Coprococcus2.id.11302.summary
)
genus.Coprococcus2.id.11302.summary_res_1 <- mr(genus.Coprococcus2.id.11302.summary_data_1)
fwrite(genus.Coprococcus2.id.11302.summary_res_1,file = "genus.Coprococcus2.id.11302.summary_res_1.csv")
fwrite(genus.Coprococcus2.id.11302.summary_data_1,file = "genus.Coprococcus2.id.11302.summary_data_1.csv")
#genus.Coprococcus3.id.11303.summary.txt
genus.Coprococcus3.id.11303.summary<-fread("genus.Coprococcus3.id.11303.summary.txt")
genus.Coprococcus3.id.11303.summary<-rename(genus.Coprococcus3.id.11303.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Coprococcus3.id.11303.summary,file = "genus.Coprococcus3.id.11303.summary.csv")
genus.Coprococcus3.id.11303.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Coprococcus3.id.11303.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Coprococcus3.id.11303.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Coprococcus3.id.11303.summary
)
genus.Coprococcus3.id.11303.summary_res_1 <- mr(genus.Coprococcus3.id.11303.summary_data_1)
fwrite(genus.Coprococcus3.id.11303.summary_res_1,file = "genus.Coprococcus3.id.11303.summary_res_1.csv")
fwrite(genus.Coprococcus3.id.11303.summary_data_1,file = "genus.Coprococcus3.id.11303.summary_data_1.csv")
#genus.DefluviitaleaceaeUCG011.id.11287.summary.txt
genus.DefluviitaleaceaeUCG011.id.11287.summary<-fread("genus.DefluviitaleaceaeUCG011.id.11287.summary.txt")
genus.DefluviitaleaceaeUCG011.id.11287.summary<-rename(genus.DefluviitaleaceaeUCG011.id.11287.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.DefluviitaleaceaeUCG011.id.11287.summary,file = "genus.DefluviitaleaceaeUCG011.id.11287.summary.csv")
genus.DefluviitaleaceaeUCG011.id.11287.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.DefluviitaleaceaeUCG011.id.11287.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.DefluviitaleaceaeUCG011.id.11287.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.DefluviitaleaceaeUCG011.id.11287.summary
)
genus.DefluviitaleaceaeUCG011.id.11287.summary_res_1 <- mr(genus.DefluviitaleaceaeUCG011.id.11287.summary_data_1)
fwrite(genus.DefluviitaleaceaeUCG011.id.11287.summary_res_1,file = "genus.DefluviitaleaceaeUCG011.id.11287.summary_res_1.csv")
fwrite(genus.DefluviitaleaceaeUCG011.id.11287.summary_data_1,file = "genus.DefluviitaleaceaeUCG011.id.11287.summary_data_1.csv")
#genus.Desulfovibrio.id.3173.summary.txt
genus.Desulfovibrio.id.3173.summary<-fread("genus.Desulfovibrio.id.3173.summary.txt")
genus.Desulfovibrio.id.3173.summary<-rename(genus.Desulfovibrio.id.3173.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Desulfovibrio.id.3173.summary,file = "genus.Desulfovibrio.id.3173.summary.csv")
genus.Desulfovibrio.id.3173.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Desulfovibrio.id.3173.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Desulfovibrio.id.3173.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Desulfovibrio.id.3173.summary
)
genus.Desulfovibrio.id.3173.summary_res_1 <- mr(genus.Desulfovibrio.id.3173.summary_data_1)
fwrite(genus.Desulfovibrio.id.3173.summary_res_1,file = "genus.Desulfovibrio.id.3173.summary_res_1.csv")
fwrite(genus.Desulfovibrio.id.3173.summary_data_1,file = "genus.Desulfovibrio.id.3173.summary_data_1.csv")
#genus.Dialister.id.2183.summary.txt
genus.Dialister.id.2183.summary<-fread("genus.Dialister.id.2183.summary.txt")
genus.Dialister.id.2183.summary<-rename(genus.Dialister.id.2183.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Dialister.id.2183.summary,file = "genus.Dialister.id.2183.summary.csv")
genus.Dialister.id.2183.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Dialister.id.2183.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Dialister.id.2183.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Dialister.id.2183.summary
)
genus.Dialister.id.2183.summary_res_1 <- mr(genus.Dialister.id.2183.summary_data_1)
fwrite(genus.Dialister.id.2183.summary_res_1,file = "genus.Dialister.id.2183.summary_res_1.csv")
fwrite(genus.Dialister.id.2183.summary_data_1,file = "genus.Dialister.id.2183.summary_data_1.csv")
#genus.Dorea.id.1997.summary.txt
genus.Dorea.id.1997.summary<-fread("genus.Dorea.id.1997.summary.txt")
genus.Dorea.id.1997.summary<-rename(genus.Dorea.id.1997.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Dorea.id.1997.summary,file = "genus.Dorea.id.1997.summary.csv")
genus.Dorea.id.1997.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Dorea.id.1997.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Dorea.id.1997.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Dorea.id.1997.summary
)
genus.Dorea.id.1997.summary_res_1 <- mr(genus.Dorea.id.1997.summary_data_1)
fwrite(genus.Dorea.id.1997.summary_res_1,file = "genus.Dorea.id.1997.summary_res_1.csv")
fwrite(genus.Dorea.id.1997.summary_data_1,file = "genus.Dorea.id.1997.summary_data_1.csv")
#genus.Eggerthella.id.819.summary.txt
genus.Eggerthella.id.819.summary<-fread("genus.Eggerthella.id.819.summary.txt")
genus.Eggerthella.id.819.summary<-rename(genus.Eggerthella.id.819.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Eggerthella.id.819.summary,file = "genus.Eggerthella.id.819.summary.csv")
genus.Eggerthella.id.819.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Eggerthella.id.819.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Eggerthella.id.819.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Eggerthella.id.819.summary
)
genus.Eggerthella.id.819.summary_res_1 <- mr(genus.Eggerthella.id.819.summary_data_1)
fwrite(genus.Eggerthella.id.819.summary_res_1,file = "genus.Eggerthella.id.819.summary_res_1.csv")
fwrite(genus.Eggerthella.id.819.summary_data_1,file = "genus.Eggerthella.id.819.summary_data_1.csv")
#genus.Eisenbergiella.id.11304.summary.txt
genus.Eisenbergiella.id.11304.summary<-fread("genus.Eisenbergiella.id.11304.summary.txt")
genus.Eisenbergiella.id.11304.summary<-rename(genus.Eisenbergiella.id.11304.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Eisenbergiella.id.11304.summary,file = "genus.Eisenbergiella.id.11304.summary.csv")
genus.Eisenbergiella.id.11304.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Eisenbergiella.id.11304.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Eisenbergiella.id.11304.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Eisenbergiella.id.11304.summary
)
genus.Eisenbergiella.id.11304.summary_res_1 <- mr(genus.Eisenbergiella.id.11304.summary_data_1)
fwrite(genus.Eisenbergiella.id.11304.summary_res_1,file = "genus.Eisenbergiella.id.11304.summary_res_1.csv")
fwrite(genus.Eisenbergiella.id.11304.summary_data_1,file = "genus.Eisenbergiella.id.11304.summary_data_1.csv")
#genus.Enterorhabdus.id.820.summary.txt
genus.Enterorhabdus.id.820.summary<-fread("genus.Enterorhabdus.id.820.summary.txt")
genus.Enterorhabdus.id.820.summary<-rename(genus.Enterorhabdus.id.820.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Enterorhabdus.id.820.summary,file = "genus.Enterorhabdus.id.820.summary.csv")
genus.Enterorhabdus.id.820.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Enterorhabdus.id.820.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Enterorhabdus.id.820.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Enterorhabdus.id.820.summary
)
genus.Enterorhabdus.id.820.summary_res_1 <- mr(genus.Enterorhabdus.id.820.summary_data_1)
fwrite(genus.Enterorhabdus.id.820.summary_res_1,file = "genus.Enterorhabdus.id.820.summary_res_1.csv")
fwrite(genus.Enterorhabdus.id.820.summary_data_1,file = "genus.Enterorhabdus.id.820.summary_data_1.csv")
#genus.Erysipelatoclostridium.id.11381.summary.txt
genus.Erysipelatoclostridium.id.11381.summary<-fread("genus.Erysipelatoclostridium.id.11381.summary.txt")
genus.Erysipelatoclostridium.id.11381.summary<-rename(genus.Erysipelatoclostridium.id.11381.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Erysipelatoclostridium.id.11381.summary,file = "genus.Erysipelatoclostridium.id.11381.summary.csv")
genus.Erysipelatoclostridium.id.11381.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Erysipelatoclostridium.id.11381.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Erysipelatoclostridium.id.11381.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Erysipelatoclostridium.id.11381.summary
)
genus.Erysipelatoclostridium.id.11381.summary_res_1 <- mr(genus.Erysipelatoclostridium.id.11381.summary_data_1)
fwrite(genus.Erysipelatoclostridium.id.11381.summary_res_1,file = "genus.Erysipelatoclostridium.id.11381.summary_res_1.csv")
fwrite(genus.Erysipelatoclostridium.id.11381.summary_data_1,file = "genus.Erysipelatoclostridium.id.11381.summary_data_1.csv")
#genus.ErysipelotrichaceaeUCG003.id.11384.summary.txt
genus.ErysipelotrichaceaeUCG003.id.11384.summary<-fread("genus.ErysipelotrichaceaeUCG003.id.11384.summary.txt")
genus.ErysipelotrichaceaeUCG003.id.11384.summary<-rename(genus.ErysipelotrichaceaeUCG003.id.11384.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.ErysipelotrichaceaeUCG003.id.11384.summary,file = "genus.ErysipelotrichaceaeUCG003.id.11384.summary.csv")
genus.ErysipelotrichaceaeUCG003.id.11384.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.ErysipelotrichaceaeUCG003.id.11384.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.ErysipelotrichaceaeUCG003.id.11384.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.ErysipelotrichaceaeUCG003.id.11384.summary
)
genus.ErysipelotrichaceaeUCG003.id.11384.summary_res_1 <- mr(genus.ErysipelotrichaceaeUCG003.id.11384.summary_data_1)
fwrite(genus.ErysipelotrichaceaeUCG003.id.11384.summary_res_1,file = "genus.ErysipelotrichaceaeUCG003.id.11384.summary_res_1.csv")
fwrite(genus.ErysipelotrichaceaeUCG003.id.11384.summary_data_1,file = "genus.ErysipelotrichaceaeUCG003.id.11384.summary_data_1.csv")
#genus.Escherichia.Shigella.id.3504.summary.txt
genus.Escherichia.Shigella.id.3504.summary<-fread("genus.Escherichia.Shigella.id.3504.summary.txt")
genus.Escherichia.Shigella.id.3504.summary<-rename(genus.Escherichia.Shigella.id.3504.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Escherichia.Shigella.id.3504.summary,file = "genus.Escherichia.Shigella.id.3504.summary.csv")
genus.Escherichia.Shigella.id.3504.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Escherichia.Shigella.id.3504.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Escherichia.Shigella.id.3504.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Escherichia.Shigella.id.3504.summary
)
genus.Escherichia.Shigella.id.3504.summary_res_1 <- mr(genus.Escherichia.Shigella.id.3504.summary_data_1)
fwrite(genus.Escherichia.Shigella.id.3504.summary_res_1,file = "genus.Escherichia.Shigella.id.3504.summary_res_1.csv")
fwrite(genus.Escherichia.Shigella.id.3504.summary_data_1,file = "genus.Escherichia.Shigella.id.3504.summary_data_1.csv")
#genus.Faecalibacterium.id.2057.summary.txt
genus.Faecalibacterium.id.2057.summary<-fread("genus.Faecalibacterium.id.2057.summary.txt")
genus.Faecalibacterium.id.2057.summary<-rename(genus.Faecalibacterium.id.2057.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Faecalibacterium.id.2057.summary,file = "genus.Faecalibacterium.id.2057.summary.csv")
genus.Faecalibacterium.id.2057.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Faecalibacterium.id.2057.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Faecalibacterium.id.2057.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Faecalibacterium.id.2057.summary
)
genus.Faecalibacterium.id.2057.summary_res_1 <- mr(genus.Faecalibacterium.id.2057.summary_data_1)
fwrite(genus.Faecalibacterium.id.2057.summary_res_1,file = "genus.Faecalibacterium.id.2057.summary_res_1.csv")
fwrite(genus.Faecalibacterium.id.2057.summary_data_1,file = "genus.Faecalibacterium.id.2057.summary_data_1.csv")
#genus.FamilyXIIIAD3011group.id.11293.summary.txt
genus.FamilyXIIIAD3011group.id.11293.summary<-fread("genus.FamilyXIIIAD3011group.id.11293.summary.txt")
genus.FamilyXIIIAD3011group.id.11293.summary<-rename(genus.FamilyXIIIAD3011group.id.11293.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.FamilyXIIIAD3011group.id.11293.summary,file = "genus.FamilyXIIIAD3011group.id.11293.summary.csv")
genus.FamilyXIIIAD3011group.id.11293.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.FamilyXIIIAD3011group.id.11293.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.FamilyXIIIAD3011group.id.11293.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.FamilyXIIIAD3011group.id.11293.summary
)
genus.FamilyXIIIAD3011group.id.11293.summary_res_1 <- mr(genus.FamilyXIIIAD3011group.id.11293.summary_data_1)
fwrite(genus.FamilyXIIIAD3011group.id.11293.summary_res_1,file = "genus.FamilyXIIIAD3011group.id.11293.summary_res_1.csv")
fwrite(genus.FamilyXIIIAD3011group.id.11293.summary_data_1,file = "genus.FamilyXIIIAD3011group.id.11293.summary_data_1.csv")
#genus.FamilyXIIIUCG001.id.11294.summary.txt
genus.FamilyXIIIUCG001.id.11294.summary<-fread("genus.FamilyXIIIUCG001.id.11294.summary.txt")
genus.FamilyXIIIUCG001.id.11294.summary<-rename(genus.FamilyXIIIUCG001.id.11294.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.FamilyXIIIUCG001.id.11294.summary,file = "genus.FamilyXIIIUCG001.id.11294.summary.csv")
genus.FamilyXIIIUCG001.id.11294.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.FamilyXIIIUCG001.id.11294.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.FamilyXIIIUCG001.id.11294.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.FamilyXIIIUCG001.id.11294.summary
)
genus.FamilyXIIIUCG001.id.11294.summary_res_1 <- mr(genus.FamilyXIIIUCG001.id.11294.summary_data_1)
fwrite(genus.FamilyXIIIUCG001.id.11294.summary_res_1,file = "genus.FamilyXIIIUCG001.id.11294.summary_res_1.csv")
fwrite(genus.FamilyXIIIUCG001.id.11294.summary_data_1,file = "genus.FamilyXIIIUCG001.id.11294.summary_data_1.csv")
#genus.Flavonifractor.id.2059.summary.txt
genus.Flavonifractor.id.2059.summary<-fread("genus.Flavonifractor.id.2059.summary.txt")
genus.Flavonifractor.id.2059.summary<-rename(genus.Flavonifractor.id.2059.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Flavonifractor.id.2059.summary,file = "genus.Flavonifractor.id.2059.summary.csv")
genus.Flavonifractor.id.2059.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Flavonifractor.id.2059.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Flavonifractor.id.2059.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Flavonifractor.id.2059.summary
)
genus.Flavonifractor.id.2059.summary_res_1 <- mr(genus.Flavonifractor.id.2059.summary_data_1)
fwrite(genus.Flavonifractor.id.2059.summary_res_1,file = "genus.Flavonifractor.id.2059.summary_res_1.csv")
fwrite(genus.Flavonifractor.id.2059.summary_data_1,file = "genus.Flavonifractor.id.2059.summary_data_1.csv")
#genus.Fusicatenibacter.id.11305.summary.txt
genus.Fusicatenibacter.id.11305.summary<-fread("genus.Fusicatenibacter.id.11305.summary.txt")
genus.Fusicatenibacter.id.11305.summary<-rename(genus.Fusicatenibacter.id.11305.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Fusicatenibacter.id.11305.summary,file = "genus.Fusicatenibacter.id.11305.summary.csv")
genus.Fusicatenibacter.id.11305.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Fusicatenibacter.id.11305.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Fusicatenibacter.id.11305.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Fusicatenibacter.id.11305.summary
)
genus.Fusicatenibacter.id.11305.summary_res_1 <- mr(genus.Fusicatenibacter.id.11305.summary_data_1)
fwrite(genus.Fusicatenibacter.id.11305.summary_res_1,file = "genus.Fusicatenibacter.id.11305.summary_res_1.csv")
fwrite(genus.Fusicatenibacter.id.11305.summary_data_1,file = "genus.Fusicatenibacter.id.11305.summary_data_1.csv")
#genus.Gordonibacter.id.821.summary.txt
genus.Gordonibacter.id.821.summary<-fread("genus.Gordonibacter.id.821.summary.txt")
genus.Gordonibacter.id.821.summary<-rename(genus.Gordonibacter.id.821.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Gordonibacter.id.821.summary,file = "genus.Gordonibacter.id.821.summary.csv")
genus.Gordonibacter.id.821.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Gordonibacter.id.821.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Gordonibacter.id.821.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Gordonibacter.id.821.summary
)
genus.Gordonibacter.id.821.summary_res_1 <- mr(genus.Gordonibacter.id.821.summary_data_1)
fwrite(genus.Gordonibacter.id.821.summary_res_1,file = "genus.Gordonibacter.id.821.summary_res_1.csv")
fwrite(genus.Gordonibacter.id.821.summary_data_1,file = "genus.Gordonibacter.id.821.summary_data_1.csv")
#genus.Haemophilus.id.3698.summary.txt
genus.Haemophilus.id.3698.summary<-fread("genus.Haemophilus.id.3698.summary.txt")
genus.Haemophilus.id.3698.summary<-rename(genus.Haemophilus.id.3698.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Haemophilus.id.3698.summary,file = "genus.Haemophilus.id.3698.summary.csv")
genus.Haemophilus.id.3698.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Haemophilus.id.3698.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Haemophilus.id.3698.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Haemophilus.id.3698.summary
)
genus.Haemophilus.id.3698.summary_res_1 <- mr(genus.Haemophilus.id.3698.summary_data_1)
fwrite(genus.Haemophilus.id.3698.summary_res_1,file = "genus.Haemophilus.id.3698.summary_res_1.csv")
fwrite(genus.Haemophilus.id.3698.summary_data_1,file = "genus.Haemophilus.id.3698.summary_data_1.csv")
#genus.Holdemanella.id.11393.summary.txt
genus.Holdemanella.id.11393.summary<-fread("genus.Holdemanella.id.11393.summary.txt")
genus.Holdemanella.id.11393.summary<-rename(genus.Holdemanella.id.11393.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Holdemanella.id.11393.summary,file = "genus.Holdemanella.id.11393.summary.csv")
genus.Holdemanella.id.11393.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Holdemanella.id.11393.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Holdemanella.id.11393.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Holdemanella.id.11393.summary
)
genus.Holdemanella.id.11393.summary_res_1 <- mr(genus.Holdemanella.id.11393.summary_data_1)
fwrite(genus.Holdemanella.id.11393.summary_res_1,file = "genus.Holdemanella.id.11393.summary_res_1.csv")
fwrite(genus.Holdemanella.id.11393.summary_data_1,file = "genus.Holdemanella.id.11393.summary_data_1.csv")
#genus.Holdemania.id.2157.summary.txt
genus.Holdemania.id.2157.summary<-fread("genus.Holdemania.id.2157.summary.txt")
genus.Holdemania.id.2157.summary<-rename(genus.Holdemania.id.2157.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Holdemania.id.2157.summary,file = "genus.Holdemania.id.2157.summary.csv")
genus.Holdemania.id.2157.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Holdemania.id.2157.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Holdemania.id.2157.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Holdemania.id.2157.summary
)
genus.Holdemania.id.2157.summary_res_1 <- mr(genus.Holdemania.id.2157.summary_data_1)
fwrite(genus.Holdemania.id.2157.summary_res_1,file = "genus.Holdemania.id.2157.summary_res_1.csv")
fwrite(genus.Holdemania.id.2157.summary_data_1,file = "genus.Holdemania.id.2157.summary_data_1.csv")
#genus.Howardella.id.2000.summary.txt
genus.Howardella.id.2000.summary<-fread("genus.Howardella.id.2000.summary.txt")
genus.Howardella.id.2000.summary<-rename(genus.Howardella.id.2000.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Howardella.id.2000.summary,file = "genus.Howardella.id.2000.summary.csv")
genus.Howardella.id.2000.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Howardella.id.2000.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Howardella.id.2000.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Howardella.id.2000.summary
)
genus.Howardella.id.2000.summary_res_1 <- mr(genus.Howardella.id.2000.summary_data_1)
fwrite(genus.Howardella.id.2000.summary_res_1,file = "genus.Howardella.id.2000.summary_res_1.csv")
fwrite(genus.Howardella.id.2000.summary_data_1,file = "genus.Howardella.id.2000.summary_data_1.csv")
#genus.Hungatella.id.11306.summary.txt
genus.Hungatella.id.11306.summary<-fread("genus.Hungatella.id.11306.summary.txt")
genus.Hungatella.id.11306.summary<-rename(genus.Hungatella.id.11306.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Hungatella.id.11306.summary,file = "genus.Hungatella.id.11306.summary.csv")
genus.Hungatella.id.11306.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Hungatella.id.11306.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Hungatella.id.11306.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Hungatella.id.11306.summary
)
genus.Hungatella.id.11306.summary_res_1 <- mr(genus.Hungatella.id.11306.summary_data_1)
fwrite(genus.Hungatella.id.11306.summary_res_1,file = "genus.Hungatella.id.11306.summary_res_1.csv")
fwrite(genus.Hungatella.id.11306.summary_data_1,file = "genus.Hungatella.id.11306.summary_data_1.csv")
#genus.Intestinibacter.id.11345.summary.txt
genus.Intestinibacter.id.11345.summary<-fread("genus.Intestinibacter.id.11345.summary.txt")
genus.Intestinibacter.id.11345.summary<-rename(genus.Intestinibacter.id.11345.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Intestinibacter.id.11345.summary,file = "genus.Intestinibacter.id.11345.summary.csv")
genus.Intestinibacter.id.11345.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Intestinibacter.id.11345.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Intestinibacter.id.11345.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Intestinibacter.id.11345.summary
)
genus.Intestinibacter.id.11345.summary_res_1 <- mr(genus.Intestinibacter.id.11345.summary_data_1)
fwrite(genus.Intestinibacter.id.11345.summary_res_1,file = "genus.Intestinibacter.id.11345.summary_res_1.csv")
fwrite(genus.Intestinibacter.id.11345.summary_data_1,file = "genus.Intestinibacter.id.11345.summary_data_1.csv")
#genus.Intestinimonas.id.2062.summary.txt
genus.Intestinimonas.id.2062.summary<-fread("genus.Intestinimonas.id.2062.summary.txt")
genus.Intestinimonas.id.2062.summary<-rename(genus.Intestinimonas.id.2062.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Intestinimonas.id.2062.summary,file = "genus.Intestinimonas.id.2062.summary.csv")
genus.Intestinimonas.id.2062.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Intestinimonas.id.2062.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Intestinimonas.id.2062.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Intestinimonas.id.2062.summary
)
genus.Intestinimonas.id.2062.summary_res_1 <- mr(genus.Intestinimonas.id.2062.summary_data_1)
or_genus.Intestinimonas.id.2062.summary_res_1 <- generate_odds_ratios(genus.Intestinimonas.id.2062.summary_res_1)
fwrite(genus.Intestinimonas.id.2062.summary_res_1,file = "genus.Intestinimonas.id.2062.summary_res_1.csv")
fwrite(genus.Intestinimonas.id.2062.summary_data_1,file = "genus.Intestinimonas.id.2062.summary_data_1.csv")
ps_genus.Intestinimonas.id.2062.summary_data_1 <- mr_scatter_plot(genus.Intestinimonas.id.2062.summary_res_1, genus.Intestinimonas.id.2062.summary_data_1)
ps_genus.Intestinimonas.id.2062.summary_data_1[[1]]
heterogeneity_genus.Intestinimonas.id.2062.summary_data_1<-mr_heterogeneity(genus.Intestinimonas.id.2062.summary_data_1)
fwrite(heterogeneity_genus.Intestinimonas.id.2062.summary_data_1,file = "heterogeneity_genus.Intestinimonas.id.2062.summary_data_1.csv")
genus.Intestinimonas.id.2062.summary_data_1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                                   SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                                   data = genus.Intestinimonas.id.2062.summary_data_1, 
                                                                   NbDistribution = 1000,  SignifThreshold = 0.05)
leaveoneout_genus.Intestinimonas.id.2062.summary_data_1 <- mr_leaveoneout(genus.Intestinimonas.id.2062.summary_data_1)
p_genus.Intestinimonas.id.2062.summary_data_1 <- mr_leaveoneout_plot(leaveoneout_genus.Intestinimonas.id.2062.summary_data_1)
p_genus.Intestinimonas.id.2062.summary_data_1[[1]]
genus.Intestinimonas.id.2062.summary_data_1_presso
#genus.Lachnoclostridium.id.11308.summary.txt
genus.Lachnoclostridium.id.11308.summary<-fread("genus.Lachnoclostridium.id.11308.summary.txt")
genus.Lachnoclostridium.id.11308.summary<-rename(genus.Lachnoclostridium.id.11308.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Lachnoclostridium.id.11308.summary,file = "genus.Lachnoclostridium.id.11308.summary.csv")
genus.Lachnoclostridium.id.11308.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Lachnoclostridium.id.11308.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Lachnoclostridium.id.11308.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Lachnoclostridium.id.11308.summary
)
genus.Lachnoclostridium.id.11308.summary_res_1 <- mr(genus.Lachnoclostridium.id.11308.summary_data_1)
fwrite(genus.Lachnoclostridium.id.11308.summary_res_1,file = "genus.Lachnoclostridium.id.11308.summary_res_1.csv")
fwrite(genus.Lachnoclostridium.id.11308.summary_data_1,file = "genus.Lachnoclostridium.id.11308.summary_data_1.csv")
#genus.Lachnospira.id.2004.summary.txt
genus.Lachnospira.id.2004.summary<-fread("genus.Lachnospira.id.2004.summary.txt")
genus.Lachnospira.id.2004.summary<-rename(genus.Lachnospira.id.2004.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.Lachnospira.id.2004.summary,file = "genus.Lachnospira.id.2004.summary.csv")
genus.Lachnospira.id.2004.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Lachnospira.id.2004.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Lachnospira.id.2004.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Lachnospira.id.2004.summary
)
genus.Lachnospira.id.2004.summary_res_1 <- mr(genus.Lachnospira.id.2004.summary_data_1)
fwrite(genus.Lachnospira.id.2004.summary_res_1,file = "genus.Lachnospira.id.2004.summary_res_1.csv")
fwrite(genus.Lachnospira.id.2004.summary_data_1,file = "genus.Lachnospira.id.2004.summary_data_1.csv")
#genus.LachnospiraceaeFCS020group.id.11314.summary.txt
genus.LachnospiraceaeFCS020group.id.11314.summary<-fread("genus.LachnospiraceaeFCS020group.id.11314.summary.txt")
genus.LachnospiraceaeFCS020group.id.11314.summary<-rename(genus.LachnospiraceaeFCS020group.id.11314.summary, 
                                                       id.outcome = bac,
                                                       beta.outcome = beta,
                                                       se.outcome = SE,
                                                       effect_allele.outcome = eff.allele,
                                                       other_allele.outcome = ref.allele,
                                                       SNP = rsID)
fwrite(genus.LachnospiraceaeFCS020group.id.11314.summary,file = "genus.LachnospiraceaeFCS020group.id.11314.summary.csv")
genus.LachnospiraceaeFCS020group.id.11314.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.LachnospiraceaeFCS020group.id.11314.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.LachnospiraceaeFCS020group.id.11314.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.LachnospiraceaeFCS020group.id.11314.summary
)
genus.LachnospiraceaeFCS020group.id.11314.summary_res_1 <- mr(genus.LachnospiraceaeFCS020group.id.11314.summary_data_1)
fwrite(genus.LachnospiraceaeFCS020group.id.11314.summary_res_1,file = "genus.LachnospiraceaeFCS020group.id.11314.summary_res_1.csv")
fwrite(genus.LachnospiraceaeFCS020group.id.11314.summary_data_1,file = "genus.LachnospiraceaeFCS020group.id.11314.summary_data_1.csv")
#genus.LachnospiraceaeNC2004group.id.11316.summary.txt
genus.LachnospiraceaeNC2004group.id.11316.summary<-fread("genus.LachnospiraceaeNC2004group.id.11316.summary.txt")
genus.LachnospiraceaeNC2004group.id.11316.summary<-rename(genus.LachnospiraceaeNC2004group.id.11316.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.LachnospiraceaeNC2004group.id.11316.summary,file = "genus.LachnospiraceaeNC2004group.id.11316.summary.csv")
genus.LachnospiraceaeNC2004group.id.11316.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.LachnospiraceaeNC2004group.id.11316.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.LachnospiraceaeNC2004group.id.11316.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.LachnospiraceaeNC2004group.id.11316.summary
)
genus.LachnospiraceaeNC2004group.id.11316.summary_res_1 <- mr(genus.LachnospiraceaeNC2004group.id.11316.summary_data_1)
fwrite(genus.LachnospiraceaeNC2004group.id.11316.summary_res_1,file = "genus.LachnospiraceaeNC2004group.id.11316.summary_res_1.csv")
fwrite(genus.LachnospiraceaeNC2004group.id.11316.summary_data_1,file = "genus.LachnospiraceaeNC2004group.id.11316.summary_data_1.csv")
#genus.LachnospiraceaeND3007group.id.11317.summary.txt
genus.LachnospiraceaeND3007group.id.11317.summary<-fread("genus.LachnospiraceaeND3007group.id.11317.summary.txt")
genus.LachnospiraceaeND3007group.id.11317.summary<-rename(genus.LachnospiraceaeND3007group.id.11317.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.LachnospiraceaeND3007group.id.11317.summary,file = "genus.LachnospiraceaeND3007group.id.11317.summary.csv")
genus.LachnospiraceaeND3007group.id.11317.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.LachnospiraceaeND3007group.id.11317.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.LachnospiraceaeND3007group.id.11317.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.LachnospiraceaeND3007group.id.11317.summary
)
genus.LachnospiraceaeND3007group.id.11317.summary_res_1 <- mr(genus.LachnospiraceaeND3007group.id.11317.summary_data_1)
fwrite(genus.LachnospiraceaeND3007group.id.11317.summary_res_1,file = "genus.LachnospiraceaeND3007group.id.11317.summary_res_1.csv")
fwrite(genus.LachnospiraceaeND3007group.id.11317.summary_data_1,file = "genus.LachnospiraceaeND3007group.id.11317.summary_data_1.csv")
#genus.LachnospiraceaeNK4A136group.id.11319.summary.txt
genus.LachnospiraceaeNK4A136group.id.11319.summary<-fread("genus.LachnospiraceaeNK4A136group.id.11319.summary.txt")
genus.LachnospiraceaeNK4A136group.id.11319.summary<-rename(genus.LachnospiraceaeNK4A136group.id.11319.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.LachnospiraceaeNK4A136group.id.11319.summary,file = "genus.LachnospiraceaeNK4A136group.id.11319.summary.csv")
genus.LachnospiraceaeNK4A136group.id.11319.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.LachnospiraceaeNK4A136group.id.11319.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.LachnospiraceaeNK4A136group.id.11319.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.LachnospiraceaeNK4A136group.id.11319.summary
)
genus.LachnospiraceaeNK4A136group.id.11319.summary_res_1 <- mr(genus.LachnospiraceaeNK4A136group.id.11319.summary_data_1)
fwrite(genus.LachnospiraceaeNK4A136group.id.11319.summary_res_1,file = "genus.LachnospiraceaeNK4A136group.id.11319.summary_res_1.csv")
fwrite(genus.LachnospiraceaeNK4A136group.id.11319.summary_data_1,file = "genus.LachnospiraceaeNK4A136group.id.11319.summary_data_1.csv")
#genus.LachnospiraceaeUCG001.id.11321.summary.txt
genus.LachnospiraceaeUCG001.id.11321.summary<-fread("genus.LachnospiraceaeUCG001.id.11321.summary.txt")
genus.LachnospiraceaeUCG001.id.11321.summary<-rename(genus.LachnospiraceaeUCG001.id.11321.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.LachnospiraceaeUCG001.id.11321.summary,file = "genus.LachnospiraceaeUCG001.id.11321.summary.csv")
genus.LachnospiraceaeUCG001.id.11321.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.LachnospiraceaeUCG001.id.11321.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.LachnospiraceaeUCG001.id.11321.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.LachnospiraceaeUCG001.id.11321.summary
)
genus.LachnospiraceaeUCG001.id.11321.summary_res_1 <- mr(genus.LachnospiraceaeUCG001.id.11321.summary_data_1)
fwrite(genus.LachnospiraceaeUCG001.id.11321.summary_res_1,file = "genus.LachnospiraceaeUCG001.id.11321.summary_res_1.csv")
fwrite(genus.LachnospiraceaeUCG001.id.11321.summary_data_1,file = "genus.LachnospiraceaeUCG001.id.11321.summary_data_1.csv")
#genus.LachnospiraceaeUCG004.id.11324.summary.txt
genus.LachnospiraceaeUCG004.id.11324.summary<-fread("genus.LachnospiraceaeUCG004.id.11324.summary.txt")
genus.LachnospiraceaeUCG004.id.11324.summary<-rename(genus.LachnospiraceaeUCG004.id.11324.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.LachnospiraceaeUCG004.id.11324.summary,file = "genus.LachnospiraceaeUCG004.id.11324.summary.csv")
genus.LachnospiraceaeUCG004.id.11324.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.LachnospiraceaeUCG004.id.11324.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.LachnospiraceaeUCG004.id.11324.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.LachnospiraceaeUCG004.id.11324.summary
)
genus.LachnospiraceaeUCG004.id.11324.summary_res_1 <- mr(genus.LachnospiraceaeUCG004.id.11324.summary_data_1)
fwrite(genus.LachnospiraceaeUCG004.id.11324.summary_res_1,file = "genus.LachnospiraceaeUCG004.id.11324.summary_res_1.csv")
fwrite(genus.LachnospiraceaeUCG004.id.11324.summary_data_1,file = "genus.LachnospiraceaeUCG004.id.11324.summary_data_1.csv")
#genus.LachnospiraceaeUCG008.id.11328.summary.txt
genus.LachnospiraceaeUCG008.id.11328.summary<-fread("genus.LachnospiraceaeUCG008.id.11328.summary.txt")
genus.LachnospiraceaeUCG008.id.11328.summary<-rename(genus.LachnospiraceaeUCG008.id.11328.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.LachnospiraceaeUCG008.id.11328.summary,file = "genus.LachnospiraceaeUCG008.id.11328.summary.csv")
genus.LachnospiraceaeUCG008.id.11328.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.LachnospiraceaeUCG008.id.11328.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.LachnospiraceaeUCG008.id.11328.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.LachnospiraceaeUCG008.id.11328.summary
)
genus.LachnospiraceaeUCG008.id.11328.summary_res_1 <- mr(genus.LachnospiraceaeUCG008.id.11328.summary_data_1)
fwrite(genus.LachnospiraceaeUCG008.id.11328.summary_res_1,file = "genus.LachnospiraceaeUCG008.id.11328.summary_res_1.csv")
fwrite(genus.LachnospiraceaeUCG008.id.11328.summary_data_1,file = "genus.LachnospiraceaeUCG008.id.11328.summary_data_1.csv")
#genus.LachnospiraceaeUCG010.id.11330.summary.txt
genus.LachnospiraceaeUCG010.id.11330.summary<-fread("genus.LachnospiraceaeUCG010.id.11330.summary.txt")
genus.LachnospiraceaeUCG010.id.11330.summary<-rename(genus.LachnospiraceaeUCG010.id.11330.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.LachnospiraceaeUCG010.id.11330.summary,file = "genus.LachnospiraceaeUCG010.id.11330.summary.csv")
genus.LachnospiraceaeUCG010.id.11330.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.LachnospiraceaeUCG010.id.11330.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.LachnospiraceaeUCG010.id.11330.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.LachnospiraceaeUCG010.id.11330.summary
)
genus.LachnospiraceaeUCG010.id.11330.summary_res_1 <- mr(genus.LachnospiraceaeUCG010.id.11330.summary_data_1)
fwrite(genus.LachnospiraceaeUCG010.id.11330.summary_res_1,file = "genus.LachnospiraceaeUCG010.id.11330.summary_res_1.csv")
fwrite(genus.LachnospiraceaeUCG010.id.11330.summary_data_1,file = "genus.LachnospiraceaeUCG010.id.11330.summary_data_1.csv")
#genus.Lactobacillus.id.1837.summary.txt
genus.Lactobacillus.id.1837.summary<-fread("genus.Lactobacillus.id.1837.summary.txt")
genus.Lactobacillus.id.1837.summary<-rename(genus.Lactobacillus.id.1837.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Lactobacillus.id.1837.summary,file = "genus.Lactobacillus.id.1837.summary.csv")
genus.Lactobacillus.id.1837.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Lactobacillus.id.1837.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Lactobacillus.id.1837.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Lactobacillus.id.1837.summary
)
genus.Lactobacillus.id.1837.summary_res_1 <- mr(genus.Lactobacillus.id.1837.summary_data_1)
fwrite(genus.Lactobacillus.id.1837.summary_res_1,file = "genus.Lactobacillus.id.1837.summary_res_1.csv")
fwrite(genus.Lactobacillus.id.1837.summary_data_1,file = "genus.Lactobacillus.id.1837.summary_data_1.csv")
#genus.Lactococcus.id.1851.summary.txt
genus.Lactococcus.id.1851.summary<-fread("genus.Lactococcus.id.1851.summary.txt")
genus.Lactococcus.id.1851.summary<-rename(genus.Lactococcus.id.1851.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Lactococcus.id.1851.summary,file = "genus.Lactococcus.id.1851.summary.csv")
genus.Lactococcus.id.1851.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Lactococcus.id.1851.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Lactococcus.id.1851.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Lactococcus.id.1851.summary
)
genus.Lactococcus.id.1851.summary_res_1 <- mr(genus.Lactococcus.id.1851.summary_data_1)
fwrite(genus.Lactococcus.id.1851.summary_res_1,file = "genus.Lactococcus.id.1851.summary_res_1.csv")
fwrite(genus.Lactococcus.id.1851.summary_data_1,file = "genus.Lactococcus.id.1851.summary_data_1.csv")
#genes.Marvinbryantia.id.2005.summary.txt
genus.Marvinbryantia.id.2005.summary<-fread("genus.Marvinbryantia.id.2005.summary.txt")
genus.Marvinbryantia.id.2005.summary<-rename(genus.Marvinbryantia.id.2005.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Marvinbryantia.id.2005.summary,file = "genus.Marvinbryantia.id.2005.summary.csv")
genus.Marvinbryantia.id.2005.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Marvinbryantia.id.2005.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Marvinbryantia.id.2005.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Marvinbryantia.id.2005.summary
)
genus.Marvinbryantia.id.2005.summary_res_1 <- mr(genus.Marvinbryantia.id.2005.summary_data_1)
fwrite(genus.Marvinbryantia.id.2005.summary_res_1,file = "genus.Marvinbryantia.id.2005.summary_res_1.csv")
fwrite(genus.Marvinbryantia.id.2005.summary_data_1,file = "genus.Marvinbryantia.id.2005.summary_data_1.csv")
#genes.Methanobrevibacter.id.123.summary.txt
genus.Methanobrevibacter.id.123.summary<-fread("genus.Methanobrevibacter.id.123.summary.txt")
genus.Methanobrevibacter.id.123.summary<-rename(genus.Methanobrevibacter.id.123.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Methanobrevibacter.id.123.summary,file = "genus.Methanobrevibacter.id.123.summary.csv")
genus.Methanobrevibacter.id.123.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Methanobrevibacter.id.123.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Methanobrevibacter.id.123.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Methanobrevibacter.id.123.summary
)
genus.Methanobrevibacter.id.123.summary_res_1 <- mr(genus.Methanobrevibacter.id.123.summary_data_1)
fwrite(genus.Methanobrevibacter.id.123.summary_res_1,file = "genus.Methanobrevibacter.id.123.summary_res_1.csv")
fwrite(genus.Methanobrevibacter.id.123.summary_data_1,file = "genus.Methanobrevibacter.id.123.summary_data_1.csv")
#genes.Odoribacter.id.952.summary.txt
genus.Odoribacter.id.952.summary<-fread("genus.Odoribacter.id.952.summary.txt")
genus.Odoribacter.id.952.summary<-rename(genus.Odoribacter.id.952.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Odoribacter.id.952.summary,file = "genus.Odoribacter.id.952.summary.csv")
genus.Odoribacter.id.952.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Odoribacter.id.952.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Odoribacter.id.952.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Odoribacter.id.952.summary
)
genus.Odoribacter.id.952.summary_res_1 <- mr(genus.Odoribacter.id.952.summary_data_1)
fwrite(genus.Odoribacter.id.952.summary_res_1,file = "genus.Odoribacter.id.952.summary_res_1.csv")
fwrite(genus.Odoribacter.id.952.summary_data_1,file = "genus.Odoribacter.id.952.summary_data_1.csv")
#genes.Olsenella.id.822.summary.txt
genus.Olsenella.id.822.summary<-fread("genus.Olsenella.id.822.summary.txt")
genus.Olsenella.id.822.summary<-rename(genus.Olsenella.id.822.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Olsenella.id.822.summary,file = "genus.Olsenella.id.822.summary.csv")
genus.Olsenella.id.822.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Olsenella.id.822.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Olsenella.id.822.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Olsenella.id.822.summary
)
genus.Olsenella.id.822.summary_res_1 <- mr(genus.Olsenella.id.822.summary_data_1)
fwrite(genus.Olsenella.id.822.summary_res_1,file = "genus.Olsenella.id.822.summary_res_1.csv")
fwrite(genus.Olsenella.id.822.summary_data_1,file = "genus.Olsenella.id.822.summary_data_1.csv")
#genes.Oscillibacter.id.2063.summary.txt
genus.Oscillibacter.id.2063.summary<-fread("genus.Oscillibacter.id.2063.summary.txt")
genus.Oscillibacter.id.2063.summary<-rename(genus.Oscillibacter.id.2063.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Oscillibacter.id.2063.summary,file = "genus.Oscillibacter.id.2063.summary.csv")
genus.Oscillibacter.id.2063.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Oscillibacter.id.2063.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Oscillibacter.id.2063.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Oscillibacter.id.2063.summary
)
genus.Oscillibacter.id.2063.summary_res_1 <- mr(genus.Oscillibacter.id.2063.summary_data_1)
fwrite(genus.Oscillibacter.id.2063.summary_res_1,file = "genus.Oscillibacter.id.2063.summary_res_1.csv")
fwrite(genus.Oscillibacter.id.2063.summary_data_1,file = "genus.Oscillibacter.id.2063.summary_data_1.csv")
#genes.Oscillospira.id.2064.summary.txt
genus.Oscillospira.id.2064.summary<-fread("genus.Oscillospira.id.2064.summary.txt")
genus.Oscillospira.id.2064.summary<-rename(genus.Oscillospira.id.2064.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Oscillospira.id.2064.summary,file = "genus.Oscillospira.id.2064.summary.csv")
genus.Oscillospira.id.2064.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Oscillospira.id.2064.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Oscillospira.id.2064.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Oscillospira.id.2064.summary
)
genus.Oscillospira.id.2064.summary_res_1 <- mr(genus.Oscillospira.id.2064.summary_data_1)
fwrite(genus.Oscillospira.id.2064.summary_res_1,file = "genus.Oscillospira.id.2064.summary_res_1.csv")
fwrite(genus.Oscillospira.id.2064.summary_data_1,file = "genus.Oscillospira.id.2064.summary_data_1.csv")
#genes.Oxalobacter.id.2978.summary.txt
genus.Oxalobacter.id.2978.summary<-fread("genus.Oxalobacter.id.2978.summary.txt")
genus.Oxalobacter.id.2978.summary<-rename(genus.Oxalobacter.id.2978.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Oxalobacter.id.2978.summary,file = "genus.Oxalobacter.id.2978.summary.csv")
genus.Oxalobacter.id.2978.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Oxalobacter.id.2978.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Oxalobacter.id.2978.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Oxalobacter.id.2978.summary
)
genus.Oxalobacter.id.2978.summary_res_1 <- mr(genus.Oxalobacter.id.2978.summary_data_1)
or_genus.Oxalobacter.id.2978.summary_res_1 <- generate_odds_ratios(genus.Oxalobacter.id.2978.summary_res_1)
fwrite(genus.Oxalobacter.id.2978.summary_res_1,file = "genus.Oxalobacter.id.2978.summary_res_1.csv")
fwrite(genus.Oxalobacter.id.2978.summary_data_1,file = "genus.Oxalobacter.id.2978.summary_data_1.csv")
ps_genus.Oxalobacter.id.2978.summary_data_1 <- mr_scatter_plot(genus.Oxalobacter.id.2978.summary_res_1, genus.Oxalobacter.id.2978.summary_data_1)
ps_genus.Oxalobacter.id.2978.summary_data_1[[1]]
heterogeneity_genus.Oxalobacter.id.2978.summary_data_1<-mr_heterogeneity(genus.Oxalobacter.id.2978.summary_data_1)
fwrite(heterogeneity_genus.Oxalobacter.id.2978.summary_data_1,file = "heterogeneity_genus.Oxalobacter.id.2978.summary_data_1.csv")
genus.Oxalobacter.id.2978.summary_data_1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                                   SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                                   data = genus.Oxalobacter.id.2978.summary_data_1, 
                                                                   NbDistribution = 1000,  SignifThreshold = 0.05)
leaveoneout_genus.Oxalobacter.id.2978.summary_data_1 <- mr_leaveoneout(genus.Oxalobacter.id.2978.summary_data_1)
p_genus.Oxalobacter.id.2978.summary_data_1 <- mr_leaveoneout_plot(leaveoneout_genus.Oxalobacter.id.2978.summary_data_1)
p_genus.Oxalobacter.id.2978.summary_data_1[[1]]
genus.Oxalobacter.id.2978.summary_data_1_presso
#genes.Parabacteroides.id.954.summary.txt
genus.Parabacteroides.id.954.summary<-fread("genus.Parabacteroides.id.954.summary.txt")
genus.Parabacteroides.id.954.summary<-rename(genus.Parabacteroides.id.954.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Parabacteroides.id.954.summary,file = "genus.Parabacteroides.id.954.summary.csv")
genus.Parabacteroides.id.954.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Parabacteroides.id.954.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Parabacteroides.id.954.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Parabacteroides.id.954.summary
)
genus.Parabacteroides.id.954.summary_res_1 <- mr(genus.Parabacteroides.id.954.summary_data_1)
fwrite(genus.Parabacteroides.id.954.summary_res_1,file = "genus.Parabacteroides.id.954.summary_res_1.csv")
fwrite(genus.Parabacteroides.id.954.summary_data_1,file = "genus.Parabacteroides.id.954.summary_data_1.csv")
#genes.Paraprevotella.id.962.summary.txt
genus.Paraprevotella.id.962.summary<-fread("genus.Paraprevotella.id.962.summary.txt")
genus.Paraprevotella.id.962.summary<-rename(genus.Paraprevotella.id.962.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Paraprevotella.id.962.summary,file = "genus.Paraprevotella.id.962.summary.csv")
genus.Paraprevotella.id.962.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Paraprevotella.id.962.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Paraprevotella.id.962.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Paraprevotella.id.962.summary
)
genus.Paraprevotella.id.962.summary_res_1 <- mr(genus.Paraprevotella.id.962.summary_data_1)
fwrite(genus.Paraprevotella.id.962.summary_res_1,file = "genus.Paraprevotella.id.962.summary_res_1.csv")
fwrite(genus.Paraprevotella.id.962.summary_data_1,file = "genus.Paraprevotella.id.962.summary_data_1.csv")
#genes.Parasutterella.id.2892.summary.txt
genus.Parasutterella.id.2892.summary<-fread("genus.Parasutterella.id.2892.summary.txt")
genus.Parasutterella.id.2892.summary<-rename(genus.Parasutterella.id.2892.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Parasutterella.id.2892.summary,file = "genus.Parasutterella.id.2892.summary.csv")
genus.Parasutterella.id.2892.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Parasutterella.id.2892.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Parasutterella.id.2892.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Parasutterella.id.2892.summary
)
genus.Parasutterella.id.2892.summary_res_1 <- mr(genus.Parasutterella.id.2892.summary_data_1)
fwrite(genus.Parasutterella.id.2892.summary_res_1,file = "genus.Parasutterella.id.2892.summary_res_1.csv")
fwrite(genus.Parasutterella.id.2892.summary_data_1,file = "genus.Parasutterella.id.2892.summary_data_1.csv")
#genes.Peptococcus.id.2037.summary.txt
genus.Peptococcus.id.2037.summary<-fread("genus.Peptococcus.id.2037.summary.txt")
genus.Peptococcus.id.2037.summary<-rename(genus.Peptococcus.id.2037.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Peptococcus.id.2037.summary,file = "genus.Peptococcus.id.2037.summary.csv")
genus.Peptococcus.id.2037.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Peptococcus.id.2037.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Peptococcus.id.2037.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Peptococcus.id.2037.summary
)
genus.Peptococcus.id.2037.summary_res_1 <- mr(genus.Peptococcus.id.2037.summary_data_1)
or_genus.Peptococcus.id.2037.summary_res_1 <- generate_odds_ratios(genus.Peptococcus.id.2037.summary_res_1)
fwrite(genus.Peptococcus.id.2037.summary_res_1,file = "genus.Peptococcus.id.2037.summary_res_1.csv")
fwrite(genus.Peptococcus.id.2037.summary_data_1,file = "genus.Peptococcus.id.2037.summary_data_1.csv")
ps_genus.Peptococcus.id.2037.summary_data_1 <- mr_scatter_plot(genus.Peptococcus.id.2037.summary_res_1, genus.Peptococcus.id.2037.summary_data_1)
ps_genus.Peptococcus.id.2037.summary_data_1[[1]]
heterogeneity_genus.Peptococcus.id.2037.summary_data_1<-mr_heterogeneity(genus.Peptococcus.id.2037.summary_data_1)
fwrite(heterogeneity_genus.Peptococcus.id.2037.summary_data_1,file = "heterogeneity_genus.Peptococcus.id.2037.summary_data_1.csv")
genus.Peptococcus.id.2037.summary_data_1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                                   SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                                   data = genus.Peptococcus.id.2037.summary_data_1, 
                                                                   NbDistribution = 1000,  SignifThreshold = 0.05)
leaveoneout_genus.Peptococcus.id.2037.summary_data_1 <- mr_leaveoneout(genus.Peptococcus.id.2037.summary_data_1)
p_genus.Peptococcus.id.2037.summary_data_1 <- mr_leaveoneout_plot(leaveoneout_genus.Peptococcus.id.2037.summary_data_1)
p_genus.Peptococcus.id.2037.summary_data_1[[1]]
genus.Peptococcus.id.2037.summary_data_1_presso
#genes.Phascolarctobacterium.id.2168.summary.txt
genus.Phascolarctobacterium.id.2168.summary<-fread("genus.Phascolarctobacterium.id.2168.summary.txt")
genus.Phascolarctobacterium.id.2168.summary<-rename(genus.Phascolarctobacterium.id.2168.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Phascolarctobacterium.id.2168.summary,file = "genus.Phascolarctobacterium.id.2168.summary.csv")
genus.Phascolarctobacterium.id.2168.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Phascolarctobacterium.id.2168.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Phascolarctobacterium.id.2168.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Phascolarctobacterium.id.2168.summary
)
genus.Phascolarctobacterium.id.2168.summary_res_1 <- mr(genus.Phascolarctobacterium.id.2168.summary_data_1)
fwrite(genus.Phascolarctobacterium.id.2168.summary_res_1,file = "genus.Phascolarctobacterium.id.2168.summary_res_1.csv")
fwrite(genus.Phascolarctobacterium.id.2168.summary_data_1,file = "genus.Phascolarctobacterium.id.2168.summary_data_1.csv")
#genes.Prevotella7.id.11182.summary.txt
genus.Prevotella7.id.11182.summary<-fread("genus.Prevotella7.id.11182.summary.txt")
genus.Prevotella7.id.11182.summary<-rename(genus.Prevotella7.id.11182.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Prevotella7.id.11182.summary,file = "genus.Prevotella7.id.11182.summary.csv")
genus.Prevotella7.id.11182.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Prevotella7.id.11182.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Prevotella7.id.11182.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Prevotella7.id.11182.summary
)
genus.Prevotella7.id.11182.summary_res_1 <- mr(genus.Prevotella7.id.11182.summary_data_1)
fwrite(genus.Prevotella7.id.11182.summary_res_1,file = "genus.Prevotella7.id.11182.summary_res_1.csv")
fwrite(genus.Prevotella7.id.11182.summary_data_1,file = "genus.Prevotella7.id.11182.summary_data_1.csv")
#genes.Prevotella9.id.11183.summary.txt
genus.Prevotella9.id.11183.summary<-fread("genus.Prevotella9.id.11183.summary.txt")
genus.Prevotella9.id.11183.summary<-rename(genus.Prevotella9.id.11183.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Prevotella9.id.11183.summary,file = "genus.Prevotella9.id.11183.summary.csv")
genus.Prevotella9.id.11183.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Prevotella9.id.11183.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Prevotella9.id.11183.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Prevotella9.id.11183.summary
)
genus.Prevotella9.id.11183.summary_res_1 <- mr(genus.Prevotella9.id.11183.summary_data_1)
fwrite(genus.Prevotella9.id.11183.summary_res_1,file = "genus.Prevotella9.id.11183.summary_res_1.csv")
fwrite(genus.Prevotella9.id.11183.summary_data_1,file = "genus.Prevotella9.id.11183.summary_data_1.csv")
#genes.RikenellaceaeRC9gutgroup.id.11191.summary.txt
genus.RikenellaceaeRC9gutgroup.id.11191.summary<-fread("genus.RikenellaceaeRC9gutgroup.id.11191.summary.txt")
genus.RikenellaceaeRC9gutgroup.id.11191.summary<-rename(genus.RikenellaceaeRC9gutgroup.id.11191.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.RikenellaceaeRC9gutgroup.id.11191.summary,file = "genus.RikenellaceaeRC9gutgroup.id.11191.summary.csv")
genus.RikenellaceaeRC9gutgroup.id.11191.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.RikenellaceaeRC9gutgroup.id.11191.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.RikenellaceaeRC9gutgroup.id.11191.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.RikenellaceaeRC9gutgroup.id.11191.summary
)
genus.RikenellaceaeRC9gutgroup.id.11191.summary_res_1 <- mr(genus.RikenellaceaeRC9gutgroup.id.11191.summary_data_1)
fwrite(genus.RikenellaceaeRC9gutgroup.id.11191.summary_res_1,file = "genus.RikenellaceaeRC9gutgroup.id.11191.summary_res_1.csv")
fwrite(genus.RikenellaceaeRC9gutgroup.id.11191.summary_data_1,file = "genus.RikenellaceaeRC9gutgroup.id.11191.summary_data_1.csv")
#genes.Romboutsia.id.11347.summary.txt
genus.Romboutsia.id.11347.summary<-fread("genus.Romboutsia.id.11347.summary.txt")
genus.Romboutsia.id.11347.summary<-rename(genus.Romboutsia.id.11347.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Romboutsia.id.11347.summary,file = "genus.Romboutsia.id.11347.summary.csv")
genus.Romboutsia.id.11347.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Romboutsia.id.11347.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Romboutsia.id.11347.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Romboutsia.id.11347.summary
)
genus.Romboutsia.id.11347.summary_res_1 <- mr(genus.Romboutsia.id.11347.summary_data_1)
fwrite(genus.Romboutsia.id.11347.summary_res_1,file = "genus.Romboutsia.id.11347.summary_res_1.csv")
fwrite(genus.Romboutsia.id.11347.summary_data_1,file = "genus.Romboutsia.id.11347.summary_data_1.csv")
#genes.Roseburia.id.2012.summary.txt
genus.Roseburia.id.2012.summary<-fread("genus.Roseburia.id.2012.summary.txt")
genus.Roseburia.id.2012.summary<-rename(genus.Roseburia.id.2012.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Roseburia.id.2012.summary,file = "genus.Roseburia.id.2012.summary.csv")
genus.Roseburia.id.2012.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Roseburia.id.2012.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Roseburia.id.2012.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Roseburia.id.2012.summary
)
genus.Roseburia.id.2012.summary_res_1 <- mr(genus.Roseburia.id.2012.summary_data_1)
fwrite(genus.Roseburia.id.2012.summary_res_1,file = "genus.Roseburia.id.2012.summary_res_1.csv")
fwrite(genus.Roseburia.id.2012.summary_data_1,file = "genus.Roseburia.id.2012.summary_data_1.csv")
#genes.Ruminiclostridium5.id.11355.summary.txt
genus.Ruminiclostridium5.id.11355.summary<-fread("genus.Ruminiclostridium5.id.11355.summary.txt")
genus.Ruminiclostridium5.id.11355.summary<-rename(genus.Ruminiclostridium5.id.11355.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Ruminiclostridium5.id.11355.summary,file = "genus.Ruminiclostridium5.id.11355.summary.csv")
genus.Ruminiclostridium5.id.11355.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Ruminiclostridium5.id.11355.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Ruminiclostridium5.id.11355.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Ruminiclostridium5.id.11355.summary
)
genus.Ruminiclostridium5.id.11355.summary_res_1 <- mr(genus.Ruminiclostridium5.id.11355.summary_data_1)
fwrite(genus.Ruminiclostridium5.id.11355.summary_res_1,file = "genus.Ruminiclostridium5.id.11355.summary_res_1.csv")
fwrite(genus.Ruminiclostridium5.id.11355.summary_data_1,file = "genus.Ruminiclostridium5.id.11355.summary_data_1.csv")
#genes.Ruminiclostridium6.id.11356.summary.txt
genus.Ruminiclostridium6.id.11356.summary<-fread("genus.Ruminiclostridium6.id.11356.summary.txt")
genus.Ruminiclostridium6.id.11356.summary<-rename(genus.Ruminiclostridium6.id.11356.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Ruminiclostridium6.id.11356.summary,file = "genus.Ruminiclostridium6.id.11356.summary.csv")
genus.Ruminiclostridium6.id.11356.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Ruminiclostridium6.id.11356.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Ruminiclostridium6.id.11356.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Ruminiclostridium6.id.11356.summary
)
genus.Ruminiclostridium6.id.11356.summary_res_1 <- mr(genus.Ruminiclostridium6.id.11356.summary_data_1)
fwrite(genus.Ruminiclostridium6.id.11356.summary_res_1,file = "genus.Ruminiclostridium6.id.11356.summary_res_1.csv")
fwrite(genus.Ruminiclostridium6.id.11356.summary_data_1,file = "genus.Ruminiclostridium6.id.11356.summary_data_1.csv")
#genes.Ruminiclostridium9.id.11357.summary.txt
genus.Ruminiclostridium9.id.11357.summary<-fread("genus.Ruminiclostridium9.id.11357.summary.txt")
genus.Ruminiclostridium9.id.11357.summary<-rename(genus.Ruminiclostridium9.id.11357.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Ruminiclostridium9.id.11357.summary,file = "genus.Ruminiclostridium9.id.11357.summary.csv")
genus.Ruminiclostridium9.id.11357.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Ruminiclostridium9.id.11357.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Ruminiclostridium9.id.11357.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Ruminiclostridium9.id.11357.summary
)
genus.Ruminiclostridium9.id.11357.summary_res_1 <- mr(genus.Ruminiclostridium9.id.11357.summary_data_1)
fwrite(genus.Ruminiclostridium9.id.11357.summary_res_1,file = "genus.Ruminiclostridium9.id.11357.summary_res_1.csv")
fwrite(genus.Ruminiclostridium9.id.11357.summary_data_1,file = "genus.Ruminiclostridium9.id.11357.summary_data_1.csv")
#genes.RuminococcaceaeNK4A214group.id.11358.summary.txt
genus.RuminococcaceaeNK4A214group.id.11358.summary<-fread("genus.RuminococcaceaeNK4A214group.id.11358.summary.txt")
genus.RuminococcaceaeNK4A214group.id.11358.summary<-rename(genus.RuminococcaceaeNK4A214group.id.11358.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.RuminococcaceaeNK4A214group.id.11358.summary,file = "genus.RuminococcaceaeNK4A214group.id.11358.summary.csv")
genus.RuminococcaceaeNK4A214group.id.11358.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.RuminococcaceaeNK4A214group.id.11358.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.RuminococcaceaeNK4A214group.id.11358.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.RuminococcaceaeNK4A214group.id.11358.summary
)
genus.RuminococcaceaeNK4A214group.id.11358.summary_res_1 <- mr(genus.RuminococcaceaeNK4A214group.id.11358.summary_data_1)
fwrite(genus.RuminococcaceaeNK4A214group.id.11358.summary_res_1,file = "genus.RuminococcaceaeNK4A214group.id.11358.summary_res_1.csv")
fwrite(genus.RuminococcaceaeNK4A214group.id.11358.summary_data_1,file = "genus.RuminococcaceaeNK4A214group.id.11358.summary_data_1.csv")
#genes.RuminococcaceaeUCG002.id.11360.summary.txt
genus.RuminococcaceaeUCG002.id.11360.summary<-fread("genus.RuminococcaceaeUCG002.id.11360.summary.txt")
genus.RuminococcaceaeUCG002.id.11360.summary<-rename(genus.RuminococcaceaeUCG002.id.11360.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.RuminococcaceaeUCG002.id.11360.summary,file = "genus.RuminococcaceaeUCG002.id.11360.summary.csv")
genus.RuminococcaceaeUCG002.id.11360.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.RuminococcaceaeUCG002.id.11360.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.RuminococcaceaeUCG002.id.11360.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.RuminococcaceaeUCG002.id.11360.summary
)
genus.RuminococcaceaeUCG002.id.11360.summary_res_1 <- mr(genus.RuminococcaceaeUCG002.id.11360.summary_data_1)
fwrite(genus.RuminococcaceaeUCG002.id.11360.summary_res_1,file = "genus.RuminococcaceaeUCG002.id.11360.summary_res_1.csv")
fwrite(genus.RuminococcaceaeUCG002.id.11360.summary_data_1,file = "genus.RuminococcaceaeUCG002.id.11360.summary_data_1.csv")
#genes.RuminococcaceaeUCG003.id.11361.summary.txt
genus.RuminococcaceaeUCG003.id.11361.summary<-fread("genus.RuminococcaceaeUCG003.id.11361.summary.txt")
genus.RuminococcaceaeUCG003.id.11361.summary<-rename(genus.RuminococcaceaeUCG003.id.11361.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.RuminococcaceaeUCG003.id.11361.summary,file = "genus.RuminococcaceaeUCG003.id.11361.summary.csv")
genus.RuminococcaceaeUCG003.id.11361.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.RuminococcaceaeUCG003.id.11361.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.RuminococcaceaeUCG003.id.11361.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.RuminococcaceaeUCG003.id.11361.summary
)
genus.RuminococcaceaeUCG003.id.11361.summary_res_1 <- mr(genus.RuminococcaceaeUCG003.id.11361.summary_data_1)
fwrite(genus.RuminococcaceaeUCG003.id.11361.summary_res_1,file = "genus.RuminococcaceaeUCG003.id.11361.summary_res_1.csv")
fwrite(genus.RuminococcaceaeUCG003.id.11361.summary_data_1,file = "genus.RuminococcaceaeUCG003.id.11361.summary_data_1.csv")
#genes.RuminococcaceaeUCG004.id.11362.summary.txt
genus.RuminococcaceaeUCG004.id.11362.summary<-fread("genus.RuminococcaceaeUCG004.id.11362.summary.txt")
genus.RuminococcaceaeUCG004.id.11362.summary<-rename(genus.RuminococcaceaeUCG004.id.11362.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.RuminococcaceaeUCG004.id.11362.summary,file = "genus.RuminococcaceaeUCG004.id.11362.summary.csv")
genus.RuminococcaceaeUCG004.id.11362.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.RuminococcaceaeUCG004.id.11362.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.RuminococcaceaeUCG004.id.11362.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.RuminococcaceaeUCG004.id.11362.summary
)
genus.RuminococcaceaeUCG004.id.11362.summary_res_1 <- mr(genus.RuminococcaceaeUCG004.id.11362.summary_data_1)
fwrite(genus.RuminococcaceaeUCG004.id.11362.summary_res_1,file = "genus.RuminococcaceaeUCG004.id.11362.summary_res_1.csv")
fwrite(genus.RuminococcaceaeUCG004.id.11362.summary_data_1,file = "genus.RuminococcaceaeUCG004.id.11362.summary_data_1.csv")
#genes.RuminococcaceaeUCG005.id.11363.summary.txt
genus.RuminococcaceaeUCG005.id.11363.summary<-fread("genus.RuminococcaceaeUCG005.id.11363.summary.txt")
genus.RuminococcaceaeUCG005.id.11363.summary<-rename(genus.RuminococcaceaeUCG005.id.11363.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.RuminococcaceaeUCG005.id.11363.summary,file = "genus.RuminococcaceaeUCG005.id.11363.summary.csv")
genus.RuminococcaceaeUCG005.id.11363.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.RuminococcaceaeUCG005.id.11363.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.RuminococcaceaeUCG005.id.11363.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.RuminococcaceaeUCG005.id.11363.summary
)
genus.RuminococcaceaeUCG005.id.11363.summary_res_1 <- mr(genus.RuminococcaceaeUCG005.id.11363.summary_data_1)
or_genus.RuminococcaceaeUCG005.id.11363.summary_res_1 <- generate_odds_ratios(genus.RuminococcaceaeUCG005.id.11363.summary_res_1)
fwrite(genus.RuminococcaceaeUCG005.id.11363.summary_res_1,file = "genus.RuminococcaceaeUCG005.id.11363.summary_res_1.csv")
fwrite(genus.RuminococcaceaeUCG005.id.11363.summary_data_1,file = "genus.RuminococcaceaeUCG005.id.11363.summary_data_1.csv")
ps_genus.RuminococcaceaeUCG005.id.11363.summary_data_1 <- mr_scatter_plot(genus.RuminococcaceaeUCG005.id.11363.summary_res_1, genus.RuminococcaceaeUCG005.id.11363.summary_data_1)
ps_genus.RuminococcaceaeUCG005.id.11363.summary_data_1[[1]]
heterogeneity_genus.RuminococcaceaeUCG005.id.11363.summary_data_1<-mr_heterogeneity(genus.RuminococcaceaeUCG005.id.11363.summary_data_1)
fwrite(heterogeneity_genus.RuminococcaceaeUCG005.id.11363.summary_data_1,file = "heterogeneity_genus.RuminococcaceaeUCG005.id.11363.summary_data_1.csv")
genus.RuminococcaceaeUCG005.id.11363.summary_data_1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                                   SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                                   data = genus.RuminococcaceaeUCG005.id.11363.summary_data_1, 
                                                                   NbDistribution = 1000,  SignifThreshold = 0.05)
leaveoneout_genus.RuminococcaceaeUCG005.id.11363.summary_data_1 <- mr_leaveoneout(genus.RuminococcaceaeUCG005.id.11363.summary_data_1)
p_genus.RuminococcaceaeUCG005.id.11363.summary_data_1 <- mr_leaveoneout_plot(leaveoneout_genus.RuminococcaceaeUCG005.id.11363.summary_data_1)
p_genus.RuminococcaceaeUCG005.id.11363.summary_data_1[[1]]
genus.RuminococcaceaeUCG005.id.11363.summary_data_1_presso
#genes.RuminococcaceaeUCG009.id.11366.summary.txt
genus.RuminococcaceaeUCG009.id.11366.summary<-fread("genus.RuminococcaceaeUCG009.id.11366.summary.txt")
genus.RuminococcaceaeUCG009.id.11366.summary<-rename(genus.RuminococcaceaeUCG009.id.11366.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.RuminococcaceaeUCG009.id.11366.summary,file = "genus.RuminococcaceaeUCG009.id.11366.summary.csv")
genus.RuminococcaceaeUCG009.id.11366.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.RuminococcaceaeUCG009.id.11366.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.RuminococcaceaeUCG009.id.11366.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.RuminococcaceaeUCG009.id.11366.summary
)
genus.RuminococcaceaeUCG009.id.11366.summary_res_1 <- mr(genus.RuminococcaceaeUCG009.id.11366.summary_data_1)
fwrite(genus.RuminococcaceaeUCG009.id.11366.summary_res_1,file = "genus.RuminococcaceaeUCG009.id.11366.summary_res_1.csv")
fwrite(genus.RuminococcaceaeUCG009.id.11366.summary_data_1,file = "genus.RuminococcaceaeUCG009.id.11366.summary_data_1.csv")
#genes.RuminococcaceaeUCG010.id.11367.summary.txt
genus.RuminococcaceaeUCG010.id.11367.summary<-fread("genus.RuminococcaceaeUCG010.id.11367.summary.txt")
genus.RuminococcaceaeUCG010.id.11367.summary<-rename(genus.RuminococcaceaeUCG010.id.11367.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.RuminococcaceaeUCG010.id.11367.summary,file = "genus.RuminococcaceaeUCG010.id.11367.summary.csv")
genus.RuminococcaceaeUCG010.id.11367.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.RuminococcaceaeUCG010.id.11367.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.RuminococcaceaeUCG010.id.11367.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.RuminococcaceaeUCG010.id.11367.summary
)
genus.RuminococcaceaeUCG010.id.11367.summary_res_1 <- mr(genus.RuminococcaceaeUCG010.id.11367.summary_data_1)
fwrite(genus.RuminococcaceaeUCG010.id.11367.summary_res_1,file = "genus.RuminococcaceaeUCG010.id.11367.summary_res_1.csv")
fwrite(genus.RuminococcaceaeUCG010.id.11367.summary_data_1,file = "genus.RuminococcaceaeUCG010.id.11367.summary_data_1.csv")
#genes.RuminococcaceaeUCG011.id.11368.summary.txt
genus.RuminococcaceaeUCG011.id.11368.summary<-fread("genus.RuminococcaceaeUCG011.id.11368.summary.txt")
genus.RuminococcaceaeUCG011.id.11368.summary<-rename(genus.RuminococcaceaeUCG011.id.11368.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.RuminococcaceaeUCG011.id.11368.summary,file = "genus.RuminococcaceaeUCG011.id.11368.summary.csv")
genus.RuminococcaceaeUCG011.id.11368.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.RuminococcaceaeUCG011.id.11368.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.RuminococcaceaeUCG011.id.11368.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.RuminococcaceaeUCG011.id.11368.summary
)
genus.RuminococcaceaeUCG011.id.11368.summary_res_1 <- mr(genus.RuminococcaceaeUCG011.id.11368.summary_data_1)
fwrite(genus.RuminococcaceaeUCG011.id.11368.summary_res_1,file = "genus.RuminococcaceaeUCG011.id.11368.summary_res_1.csv")
fwrite(genus.RuminococcaceaeUCG011.id.11368.summary_data_1,file = "genus.RuminococcaceaeUCG011.id.11368.summary_data_1.csv")
#genes.RuminococcaceaeUCG013.id.11370.summary.txt
genus.RuminococcaceaeUCG013.id.11370.summary<-fread("genus.RuminococcaceaeUCG013.id.11370.summary.txt")
genus.RuminococcaceaeUCG013.id.11370.summary<-rename(genus.RuminococcaceaeUCG013.id.11370.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.RuminococcaceaeUCG013.id.11370.summary,file = "genus.RuminococcaceaeUCG013.id.11370.summary.csv")
genus.RuminococcaceaeUCG013.id.11370.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.RuminococcaceaeUCG013.id.11370.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.RuminococcaceaeUCG013.id.11370.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.RuminococcaceaeUCG013.id.11370.summary
)
genus.RuminococcaceaeUCG013.id.11370.summary_res_1 <- mr(genus.RuminococcaceaeUCG013.id.11370.summary_data_1)
fwrite(genus.RuminococcaceaeUCG013.id.11370.summary_res_1,file = "genus.RuminococcaceaeUCG013.id.11370.summary_res_1.csv")
fwrite(genus.RuminococcaceaeUCG013.id.11370.summary_data_1,file = "genus.RuminococcaceaeUCG013.id.11370.summary_data_1.csv")
#genes.RuminococcaceaeUCG014.id.11371.summary.txt
genus.RuminococcaceaeUCG014.id.11371.summary<-fread("genus.RuminococcaceaeUCG014.id.11371.summary.txt")
genus.RuminococcaceaeUCG014.id.11371.summary<-rename(genus.RuminococcaceaeUCG014.id.11371.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.RuminococcaceaeUCG014.id.11371.summary,file = "genus.RuminococcaceaeUCG014.id.11371.summary.csv")
genus.RuminococcaceaeUCG014.id.11371.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.RuminococcaceaeUCG014.id.11371.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.RuminococcaceaeUCG014.id.11371.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.RuminococcaceaeUCG014.id.11371.summary
)
genus.RuminococcaceaeUCG014.id.11371.summary_res_1 <- mr(genus.RuminococcaceaeUCG014.id.11371.summary_data_1)
fwrite(genus.RuminococcaceaeUCG014.id.11371.summary_res_1,file = "genus.RuminococcaceaeUCG014.id.11371.summary_res_1.csv")
fwrite(genus.RuminococcaceaeUCG014.id.11371.summary_data_1,file = "genus.RuminococcaceaeUCG014.id.11371.summary_data_1.csv")
#genes.Ruminococcus1.id.11373.summary.txt
genus.Ruminococcus1.id.11373.summary<-fread("genus.Ruminococcus1.id.11373.summary.txt")
genus.Ruminococcus1.id.11373.summary<-rename(genus.Ruminococcus1.id.11373.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Ruminococcus1.id.11373.summary,file = "genus.Ruminococcus1.id.11373.summary.csv")
genus.Ruminococcus1.id.11373.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Ruminococcus1.id.11373.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Ruminococcus1.id.11373.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Ruminococcus1.id.11373.summary
)
genus.Ruminococcus1.id.11373.summary_res_1 <- mr(genus.Ruminococcus1.id.11373.summary_data_1)
fwrite(genus.Ruminococcus1.id.11373.summary_res_1,file = "genus.Ruminococcus1.id.11373.summary_res_1.csv")
fwrite(genus.Ruminococcus1.id.11373.summary_data_1,file = "genus.Ruminococcus1.id.11373.summary_data_1.csv")
#genes.Ruminococcus2.id.11374.summary.txt
genus.Ruminococcus2.id.11374.summary<-fread("genus.Ruminococcus2.id.11374.summary.txt")
genus.Ruminococcus2.id.11374.summary<-rename(genus.Ruminococcus2.id.11374.summary, 
                                                          id.outcome = bac,
                                                          beta.outcome = beta,
                                                          se.outcome = SE,
                                                          effect_allele.outcome = eff.allele,
                                                          other_allele.outcome = ref.allele,
                                                          SNP = rsID)
fwrite(genus.Ruminococcus2.id.11374.summary,file = "genus.Ruminococcus2.id.11374.summary.csv")
genus.Ruminococcus2.id.11374.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Ruminococcus2.id.11374.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Ruminococcus2.id.11374.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Ruminococcus2.id.11374.summary
)
genus.Ruminococcus2.id.11374.summary_res_1 <- mr(genus.Ruminococcus2.id.11374.summary_data_1)
fwrite(genus.Ruminococcus2.id.11374.summary_res_1,file = "genus.Ruminococcus2.id.11374.summary_res_1.csv")
fwrite(genus.Ruminococcus2.id.11374.summary_data_1,file = "genus.Ruminococcus2.id.11374.summary_data_1.csv")
#genes.Sellimonas.id.14369.summary.txt
genus.Sellimonas.id.14369.summary<-fread("genus.Sellimonas.id.14369.summary.txt")
genus.Sellimonas.id.14369.summary<-rename(genus.Sellimonas.id.14369.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.Sellimonas.id.14369.summary,file = "genus.Sellimonas.id.14369.summary.csv")
genus.Sellimonas.id.14369.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Sellimonas.id.14369.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Sellimonas.id.14369.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Sellimonas.id.14369.summary
)
genus.Sellimonas.id.14369.summary_res_1 <- mr(genus.Sellimonas.id.14369.summary_data_1)
fwrite(genus.Sellimonas.id.14369.summary_res_1,file = "genus.Sellimonas.id.14369.summary_res_1.csv")
fwrite(genus.Sellimonas.id.14369.summary_data_1,file = "genus.Sellimonas.id.14369.summary_data_1.csv")
#genes.Senegalimassilia.id.11160.summary.txt
genus.Senegalimassilia.id.11160.summary<-fread("genus.Senegalimassilia.id.11160.summary.txt")
genus.Senegalimassilia.id.11160.summary<-rename(genus.Senegalimassilia.id.11160.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.Senegalimassilia.id.11160.summary,file = "genus.Senegalimassilia.id.11160.summary.csv")
genus.Senegalimassilia.id.11160.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Senegalimassilia.id.11160.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Senegalimassilia.id.11160.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Senegalimassilia.id.11160.summary
)
genus.Senegalimassilia.id.11160.summary_res_1 <- mr(genus.Senegalimassilia.id.11160.summary_data_1)
fwrite(genus.Senegalimassilia.id.11160.summary_res_1,file = "genus.Senegalimassilia.id.11160.summary_res_1.csv")
fwrite(genus.Senegalimassilia.id.11160.summary_data_1,file = "genus.Senegalimassilia.id.11160.summary_data_1.csv")
#genes.Slackia.id.825.summary.txt
genus.Slackia.id.825.summary<-fread("genus.Slackia.id.825.summary.txt")
genus.Slackia.id.825.summary<-rename(genus.Slackia.id.825.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.Slackia.id.825.summary,file = "genus.Slackia.id.825.summary.csv")
genus.Slackia.id.825.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Slackia.id.825.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Slackia.id.825.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Slackia.id.825.summary
)
genus.Slackia.id.825.summary_res_1 <- mr(genus.Slackia.id.825.summary_data_1)
fwrite(genus.Slackia.id.825.summary_res_1,file = "genus.Slackia.id.825.summary_res_1.csv")
fwrite(genus.Slackia.id.825.summary_data_1,file = "genus.Slackia.id.825.summary_data_1.csv")
#genes.Streptococcus.id.1853.summary.txt
genus.Streptococcus.id.1853.summary<-fread("genus.Streptococcus.id.1853.summary.txt")
genus.Streptococcus.id.1853.summary<-rename(genus.Streptococcus.id.1853.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.Streptococcus.id.1853.summary,file = "genus.Streptococcus.id.1853.summary.csv")
genus.Streptococcus.id.1853.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Streptococcus.id.1853.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Streptococcus.id.1853.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Streptococcus.id.1853.summary
)
genus.Streptococcus.id.1853.summary_res_1 <- mr(genus.Streptococcus.id.1853.summary_data_1)
fwrite(genus.Streptococcus.id.1853.summary_res_1,file = "genus.Streptococcus.id.1853.summary_res_1.csv")
fwrite(genus.Streptococcus.id.1853.summary_data_1,file = "genus.Streptococcus.id.1853.summary_data_1.csv")
#genes.Subdoligranulum.id.2070.summary.txt
genus.Subdoligranulum.id.2070.summary<-fread("genus.Subdoligranulum.id.2070.summary.txt")
genus.Subdoligranulum.id.2070.summary<-rename(genus.Subdoligranulum.id.2070.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.Subdoligranulum.id.2070.summary,file = "genus.Subdoligranulum.id.2070.summary.csv")
genus.Subdoligranulum.id.2070.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Subdoligranulum.id.2070.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Subdoligranulum.id.2070.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Subdoligranulum.id.2070.summary
)
genus.Subdoligranulum.id.2070.summary_res_1 <- mr(genus.Subdoligranulum.id.2070.summary_data_1)
fwrite(genus.Subdoligranulum.id.2070.summary_res_1,file = "genus.Subdoligranulum.id.2070.summary_res_1.csv")
fwrite(genus.Subdoligranulum.id.2070.summary_data_1,file = "genus.Subdoligranulum.id.2070.summary_data_1.csv")
#genes.Sutterella.id.2896.summary.txt
genus.Sutterella.id.2896.summary<-fread("genus.Sutterella.id.2896.summary.txt")
genus.Sutterella.id.2896.summary<-rename(genus.Sutterella.id.2896.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.Sutterella.id.2896.summary,file = "genus.Sutterella.id.2896.summary.csv")
genus.Sutterella.id.2896.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Sutterella.id.2896.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Sutterella.id.2896.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Sutterella.id.2896.summary
)
genus.Sutterella.id.2896.summary_res_1 <- mr(genus.Sutterella.id.2896.summary_data_1)
or_genus.Sutterella.id.2896.summary_res_1 <- generate_odds_ratios(genus.Sutterella.id.2896.summary_res_1)
fwrite(genus.Sutterella.id.2896.summary_res_1,file = "genus.Sutterella.id.2896.summary_res_1.csv")
fwrite(genus.Sutterella.id.2896.summary_data_1,file = "genus.Sutterella.id.2896.summary_data_1.csv")
ps_genus.Sutterella.id.2896.summary_data_1 <- mr_scatter_plot(genus.Sutterella.id.2896.summary_res_1, genus.Sutterella.id.2896.summary_data_1)
ps_genus.Sutterella.id.2896.summary_data_1[[1]]
heterogeneity_genus.Sutterella.id.2896.summary_data_1<-mr_heterogeneity(genus.Sutterella.id.2896.summary_data_1)
fwrite(heterogeneity_genus.Sutterella.id.2896.summary_data_1,file = "heterogeneity_genus.Sutterella.id.2896.summary_data_1.csv")
genus.Sutterella.id.2896.summary_data_1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                                   SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                                   data = genus.Sutterella.id.2896.summary_data_1, 
                                                                   NbDistribution = 1000,  SignifThreshold = 0.05)
leaveoneout_genus.Sutterella.id.2896.summary_data_1 <- mr_leaveoneout(genus.Sutterella.id.2896.summary_data_1)
p_genus.Sutterella.id.2896.summary_data_1 <- mr_leaveoneout_plot(leaveoneout_genus.Sutterella.id.2896.summary_data_1)
p_genus.Sutterella.id.2896.summary_data_1[[1]]
genus.Sutterella.id.2896.summary_data_1_presso
#genes.Terrisporobacter.id.11348.summary.txt
genus.Terrisporobacter.id.11348.summary<-fread("genus.Terrisporobacter.id.11348.summary.txt")
genus.Terrisporobacter.id.11348.summary<-rename(genus.Terrisporobacter.id.11348.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.Terrisporobacter.id.11348.summary,file = "genus.Terrisporobacter.id.11348.summary.csv")
genus.Terrisporobacter.id.11348.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Terrisporobacter.id.11348.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Terrisporobacter.id.11348.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Terrisporobacter.id.11348.summary
)
genus.Terrisporobacter.id.11348.summary_res_1 <- mr(genus.Terrisporobacter.id.11348.summary_data_1)
or_genus.Terrisporobacter.id.11348.summary_res_1 <- generate_odds_ratios(genus.Terrisporobacter.id.11348.summary_res_1)
fwrite(genus.Terrisporobacter.id.11348.summary_res_1,file = "genus.Terrisporobacter.id.11348.summary_res_1.csv")
fwrite(genus.Terrisporobacter.id.11348.summary_data_1,file = "genus.Terrisporobacter.id.11348.summary_data_1.csv")
ps_genus.Terrisporobacter.id.11348.summary_data_1 <- mr_scatter_plot(genus.Terrisporobacter.id.11348.summary_res_1, genus.Terrisporobacter.id.11348.summary_data_1)
ps_genus.Terrisporobacter.id.11348.summary_data_1[[1]]
heterogeneity_genus.Terrisporobacter.id.11348.summary_data_1<-mr_heterogeneity(genus.Terrisporobacter.id.11348.summary_data_1)
fwrite(heterogeneity_genus.Terrisporobacter.id.11348.summary_data_1,file = "heterogeneity_genus.Terrisporobacter.id.11348.summary_data_1.csv")
genus.Terrisporobacter.id.11348.summary_data_1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                                                   SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                                                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                                                   data = genus.Terrisporobacter.id.11348.summary_data_1, 
                                                                   NbDistribution = 1000,  SignifThreshold = 0.05)
leaveoneout_genus.Terrisporobacter.id.11348.summary_data_1 <- mr_leaveoneout(genus.Terrisporobacter.id.11348.summary_data_1)
p_genus.Terrisporobacter.id.11348.summary_data_1 <- mr_leaveoneout_plot(leaveoneout_genus.Terrisporobacter.id.11348.summary_data_1)
p_genus.Terrisporobacter.id.11348.summary_data_1[[1]]
genus.Terrisporobacter.id.11348.summary_data_1_presso
#genes.Turicibacter.id.2162.summary.txt
genus.Turicibacter.id.2162.summary<-fread("genus.Turicibacter.id.2162.summary.txt")
genus.Turicibacter.id.2162.summary<-rename(genus.Turicibacter.id.2162.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.Turicibacter.id.2162.summary,file = "genus.Turicibacter.id.2162.summary.csv")
genus.Turicibacter.id.2162.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Turicibacter.id.2162.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Turicibacter.id.2162.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Turicibacter.id.2162.summary
)
genus.Turicibacter.id.2162.summary_res_1 <- mr(genus.Turicibacter.id.2162.summary_data_1)
fwrite(genus.Turicibacter.id.2162.summary_res_1,file = "genus.Turicibacter.id.2162.summary_res_1.csv")
fwrite(genus.Turicibacter.id.2162.summary_data_1,file = "genus.Turicibacter.id.2162.summary_data_1.csv")
#genes.Tyzzerella3.id.11335.summary.txt
genus.Tyzzerella3.id.11335.summary<-fread("genus.Tyzzerella3.id.11335.summary.txt")
genus.Tyzzerella3.id.11335.summary<-rename(genus.Tyzzerella3.id.11335.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.Tyzzerella3.id.11335.summary,file = "genus.Tyzzerella3.id.11335.summary.csv")
genus.Tyzzerella3.id.11335.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Tyzzerella3.id.11335.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Tyzzerella3.id.11335.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Tyzzerella3.id.11335.summary
)
genus.Tyzzerella3.id.11335.summary_res_1 <- mr(genus.Tyzzerella3.id.11335.summary_data_1)
fwrite(genus.Tyzzerella3.id.11335.summary_res_1,file = "genus.Tyzzerella3.id.11335.summary_res_1.csv")
fwrite(genus.Tyzzerella3.id.11335.summary_data_1,file = "genus.Tyzzerella3.id.11335.summary_data_1.csv")
#genes.unknowngenus.id.1000000073.summary.txt
genus.unknowngenus.id.1000000073.summary<-fread("genus.unknowngenus.id.1000000073.summary.txt")
genus.unknowngenus.id.1000000073.summary<-rename(genus.unknowngenus.id.1000000073.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.1000000073.summary,file = "genus.unknowngenus.id.1000000073.summary.csv")
genus.unknowngenus.id.1000000073.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.1000000073.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.1000000073.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.1000000073.summary
)
genus.unknowngenus.id.1000000073.summary_res_1 <- mr(genus.unknowngenus.id.1000000073.summary_data_1)
fwrite(genus.unknowngenus.id.1000000073.summary_res_1,file = "genus.unknowngenus.id.1000000073.summary_res_1.csv")
fwrite(genus.unknowngenus.id.1000000073.summary_data_1,file = "genus.unknowngenus.id.1000000073.summary_data_1.csv")
#genes.unknowngenus.id.1000001215.summary.txt
genus.unknowngenus.id.1000001215.summary<-fread("genus.unknowngenus.id.1000001215.summary.txt")
genus.unknowngenus.id.1000001215.summary<-rename(genus.unknowngenus.id.1000001215.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.1000001215.summary,file = "genus.unknowngenus.id.1000001215.summary.csv")
genus.unknowngenus.id.1000001215.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.1000001215.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.1000001215.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.1000001215.summary
)
genus.unknowngenus.id.1000001215.summary_res_1 <- mr(genus.unknowngenus.id.1000001215.summary_data_1)
fwrite(genus.unknowngenus.id.1000001215.summary_res_1,file = "genus.unknowngenus.id.1000001215.summary_res_1.csv")
fwrite(genus.unknowngenus.id.1000001215.summary_data_1,file = "genus.unknowngenus.id.1000001215.summary_data_1.csv")
#genes.unknowngenus.id.1000005472.summary.txt
genus.unknowngenus.id.1000005472.summary<-fread("genus.unknowngenus.id.1000005472.summary.txt")
genus.unknowngenus.id.1000005472.summary<-rename(genus.unknowngenus.id.1000005472.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.1000005472.summary,file = "genus.unknowngenus.id.1000005472.summary.csv")
genus.unknowngenus.id.1000005472.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.1000005472.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.1000005472.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.1000005472.summary
)
genus.unknowngenus.id.1000005472.summary_res_1 <- mr(genus.unknowngenus.id.1000005472.summary_data_1)
fwrite(genus.unknowngenus.id.1000005472.summary_res_1,file = "genus.unknowngenus.id.1000005472.summary_res_1.csv")
fwrite(genus.unknowngenus.id.1000005472.summary_data_1,file = "genus.unknowngenus.id.1000005472.summary_data_1.csv")
#genes.unknowngenus.id.1000005479.summary.txt
genus.unknowngenus.id.1000005479.summary<-fread("genus.unknowngenus.id.1000005479.summary.txt")
genus.unknowngenus.id.1000005479.summary<-rename(genus.unknowngenus.id.1000005479.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.1000005479.summary,file = "genus.unknowngenus.id.1000005479.summary.csv")
genus.unknowngenus.id.1000005479.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.1000005479.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.1000005479.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.1000005479.summary
)
genus.unknowngenus.id.1000005479.summary_res_1 <- mr(genus.unknowngenus.id.1000005479.summary_data_1)
fwrite(genus.unknowngenus.id.1000005479.summary_res_1,file = "genus.unknowngenus.id.1000005479.summary_res_1.csv")
fwrite(genus.unknowngenus.id.1000005479.summary_data_1,file = "genus.unknowngenus.id.1000005479.summary_data_1.csv")
#genes.unknowngenus.id.1000006162.summary.txt
genus.unknowngenus.id.1000006162.summary<-fread("genus.unknowngenus.id.1000006162.summary.txt")
genus.unknowngenus.id.1000006162.summary<-rename(genus.unknowngenus.id.1000006162.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.1000006162.summary,file = "genus.unknowngenus.id.1000006162.summary.csv")
genus.unknowngenus.id.1000006162.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.1000006162.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.1000006162.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.1000006162.summary
)
genus.unknowngenus.id.1000006162.summary_res_1 <- mr(genus.unknowngenus.id.1000006162.summary_data_1)
fwrite(genus.unknowngenus.id.1000006162.summary_res_1,file = "genus.unknowngenus.id.1000006162.summary_res_1.csv")
fwrite(genus.unknowngenus.id.1000006162.summary_data_1,file = "genus.unknowngenus.id.1000006162.summary_data_1.csv")
#genes.unknowngenus.id.1868.summary.txt
genus.unknowngenus.id.1868.summary<-fread("genus.unknowngenus.id.1868.summary.txt")
genus.unknowngenus.id.1868.summary<-rename(genus.unknowngenus.id.1868.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.1868.summary,file = "genus.unknowngenus.id.1868.summary.csv")
genus.unknowngenus.id.1868.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.1868.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.1868.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.1868.summary
)
genus.unknowngenus.id.1868.summary_res_1 <- mr(genus.unknowngenus.id.1868.summary_data_1)
fwrite(genus.unknowngenus.id.1868.summary_res_1,file = "genus.unknowngenus.id.1868.summary_res_1.csv")
fwrite(genus.unknowngenus.id.1868.summary_data_1,file = "genus.unknowngenus.id.1868.summary_data_1.csv")
#genes.unknowngenus.id.2001.summary.txt
genus.unknowngenus.id.2001.summary<-fread("genus.unknowngenus.id.2001.summary.txt")
genus.unknowngenus.id.2001.summary<-rename(genus.unknowngenus.id.2001.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.2001.summary,file = "genus.unknowngenus.id.2001.summary.csv")
genus.unknowngenus.id.2001.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.2001.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.2001.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.2001.summary
)
genus.unknowngenus.id.2001.summary_res_1 <- mr(genus.unknowngenus.id.2001.summary_data_1)
fwrite(genus.unknowngenus.id.2001.summary_res_1,file = "genus.unknowngenus.id.2001.summary_res_1.csv")
fwrite(genus.unknowngenus.id.2001.summary_data_1,file = "genus.unknowngenus.id.2001.summary_data_1.csv")
#genes.unknowngenus.id.2041.summary.txt
genus.unknowngenus.id.2041.summary<-fread("genus.unknowngenus.id.2041.summary.txt")
genus.unknowngenus.id.2041.summary<-rename(genus.unknowngenus.id.2041.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.2041.summary,file = "genus.unknowngenus.id.2041.summary.csv")
genus.unknowngenus.id.2041.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.2041.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.2041.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.2041.summary
)
genus.unknowngenus.id.2041.summary_res_1 <- mr(genus.unknowngenus.id.2041.summary_data_1)
fwrite(genus.unknowngenus.id.2041.summary_res_1,file = "genus.unknowngenus.id.2041.summary_res_1.csv")
fwrite(genus.unknowngenus.id.2041.summary_data_1,file = "genus.unknowngenus.id.2041.summary_data_1.csv")
#genes.unknowngenus.id.2071.summary.txt
genus.unknowngenus.id.2071.summary<-fread("genus.unknowngenus.id.2071.summary.txt")
genus.unknowngenus.id.2071.summary<-rename(genus.unknowngenus.id.2071.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.2071.summary,file = "genus.unknowngenus.id.2071.summary.csv")
genus.unknowngenus.id.2071.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.2071.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.2071.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.2071.summary
)
genus.unknowngenus.id.2071.summary_res_1 <- mr(genus.unknowngenus.id.2071.summary_data_1)
fwrite(genus.unknowngenus.id.2071.summary_res_1,file = "genus.unknowngenus.id.2071.summary_res_1.csv")
fwrite(genus.unknowngenus.id.2071.summary_data_1,file = "genus.unknowngenus.id.2071.summary_data_1.csv")
#genes.unknowngenus.id.2755.summary.txt
genus.unknowngenus.id.2755.summary<-fread("genus.unknowngenus.id.2755.summary.txt")
genus.unknowngenus.id.2755.summary<-rename(genus.unknowngenus.id.2755.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.2755.summary,file = "genus.unknowngenus.id.2755.summary.csv")
genus.unknowngenus.id.2755.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.2755.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.2755.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.2755.summary
)
genus.unknowngenus.id.2755.summary_res_1 <- mr(genus.unknowngenus.id.2755.summary_data_1)
fwrite(genus.unknowngenus.id.2755.summary_res_1,file = "genus.unknowngenus.id.2755.summary_res_1.csv")
fwrite(genus.unknowngenus.id.2755.summary_data_1,file = "genus.unknowngenus.id.2755.summary_data_1.csv")
#genes.unknowngenus.id.826.summary.txt
genus.unknowngenus.id.826.summary<-fread("genus.unknowngenus.id.826.summary.txt")
genus.unknowngenus.id.826.summary<-rename(genus.unknowngenus.id.826.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.826.summary,file = "genus.unknowngenus.id.826.summary.csv")
genus.unknowngenus.id.826.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.826.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.826.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.826.summary
)
genus.unknowngenus.id.826.summary_res_1 <- mr(genus.unknowngenus.id.826.summary_data_1)
fwrite(genus.unknowngenus.id.826.summary_res_1,file = "genus.unknowngenus.id.826.summary_res_1.csv")
fwrite(genus.unknowngenus.id.826.summary_data_1,file = "genus.unknowngenus.id.826.summary_data_1.csv")
#genes.unknowngenus.id.959.summary.txt
genus.unknowngenus.id.959.summary<-fread("genus.unknowngenus.id.959.summary.txt")
genus.unknowngenus.id.959.summary<-rename(genus.unknowngenus.id.959.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.unknowngenus.id.959.summary,file = "genus.unknowngenus.id.959.summary.csv")
genus.unknowngenus.id.959.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.unknowngenus.id.959.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.unknowngenus.id.959.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.unknowngenus.id.959.summary
)
genus.unknowngenus.id.959.summary_res_1 <- mr(genus.unknowngenus.id.959.summary_data_1)
fwrite(genus.unknowngenus.id.959.summary_res_1,file = "genus.unknowngenus.id.959.summary_res_1.csv")
fwrite(genus.unknowngenus.id.959.summary_data_1,file = "genus.unknowngenus.id.959.summary_data_1.csv")
#genes.Veillonella.id.2198.summary.txt
genus.Veillonella.id.2198.summary<-fread("genus.Veillonella.id.2198.summary.txt")
genus.Veillonella.id.2198.summary<-rename(genus.Veillonella.id.2198.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.Veillonella.id.2198.summary,file = "genus.Veillonella.id.2198.summary.csv")
genus.Veillonella.id.2198.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Veillonella.id.2198.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Veillonella.id.2198.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Veillonella.id.2198.summary
)
genus.Veillonella.id.2198.summary_res_1 <- mr(genus.Veillonella.id.2198.summary_data_1)
fwrite(genus.Veillonella.id.2198.summary_res_1,file = "genus.Veillonella.id.2198.summary_res_1.csv")
fwrite(genus.Veillonella.id.2198.summary_data_1,file = "genus.Veillonella.id.2198.summary_data_1.csv")
#genes.Victivallis.id.2256.summary.txt
genus.Victivallis.id.2256.summary<-fread("genus.Victivallis.id.2256.summary.txt")
genus.Victivallis.id.2256.summary<-rename(genus.Victivallis.id.2256.summary, 
                                             id.outcome = bac,
                                             beta.outcome = beta,
                                             se.outcome = SE,
                                             effect_allele.outcome = eff.allele,
                                             other_allele.outcome = ref.allele,
                                             SNP = rsID)
fwrite(genus.Victivallis.id.2256.summary,file = "genus.Victivallis.id.2256.summary.csv")
genus.Victivallis.id.2256.summary <- read_outcome_data(
  snps = iv_graves_1$SNP,
  filename = "genus.Victivallis.id.2256.summary.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "P.weighetdSumZ",
  samplesize_col = "N"
)
genus.Victivallis.id.2256.summary_data_1 <- harmonise_data(
  exposure_dat = iv_graves_1, 
  outcome_dat = genus.Victivallis.id.2256.summary
)
genus.Victivallis.id.2256.summary_res_1 <- mr(genus.Victivallis.id.2256.summary_data_1)
fwrite(genus.Victivallis.id.2256.summary_res_1,file = "genus.Victivallis.id.2256.summary_res_1.csv")
fwrite(genus.Victivallis.id.2256.summary_data_1,file = "genus.Victivallis.id.2256.summary_data_1.csv")
fwrite(or_family.Oxalobacteraceae.id.2966.summary,file = "or_family.Oxalobacteraceae.id.2966.summary.csv")
fwrite(or_genus..Clostridiuminnocuumgroup.id.14397.summary_res_1,file = "or_genus..Clostridiuminnocuumgroup.id.14397.summary_res_1.csv")
fwrite(or_genus.Anaerofilum.id.2053.summary_res_1,file = "or_genus.Anaerofilum.id.2053.summary_res_1.csv")
fwrite(or_genus.Intestinimonas.id.2062.summary_res_1,file = "or_genus.Intestinimonas.id.2062.summary_res_1.csv")
fwrite(or_genus.Oxalobacter.id.2978.summary_res_1,file = "or_genus.Oxalobacter.id.2978.summary_res_1.csv")
fwrite(or_genus.Peptococcus.id.2037.summary_res_1,file = "or_genus.Peptococcus.id.2037.summary_res_1.csv")
fwrite(or_genus.RuminococcaceaeUCG005.id.11363.summary_res_1,file = "or_genus.RuminococcaceaeUCG005.id.11363.summary_res_1.csv")
fwrite(or_genus.Sutterella.id.2896.summary_res_1,file = "or_genus.Sutterella.id.2896.summary_res_1.csv")
fwrite(or_genus.Terrisporobacter.id.11348.summary_res_1,file = "or_genus.Terrisporobacter.id.11348.summary_res_1.csv")