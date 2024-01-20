Software or program: R 4.2.3
1.Source Code for MR Analysis of Anxiety and Stress-Related Disorders and Endometrial Cancer
library(MRPracticals)
library(TwoSampleMR)
library(MRInstruments)
library(MendelianRandomization)
data(gwas_catalog)
exposure_gwas <- subset(gwas_catalog, grepl("31116379", PubmedID)) #最大GWAS
anxiety<-anxiety[anxiety$`P-VALUE`<1*10^-6,]
exposure_data<-format_data(exposure_gwas)
anxiety<-clump_data(anxiety, clump_r2 = 0.001)
ao<-available_outcomes()

#outcome来源
outcome_gwas <- subset(ao, grepl("O'Mara TA", author))
outcome1<- extract_outcome_data(snps = SNP, 
                      outcomes = "ebi-a-GCST006465",
                      proxies = TRUE,
                      rsq = 0.6,
                      align_alleles = 1,
                      palindromes = 1,
                      maf_threshold = 0.3,
                      access_token = ieugwasr::check_access_token(),
                      splitsize = 10000,
                      proxy_splitsize = 500) #子宫内膜样组织学
outcome2<- extract_outcome_data(snps = anxiety$SNP, 
                      outcomes = "ebi-a-GCST006466",
                      proxies = TRUE,
                      rsq = 0.8,
                      align_alleles = 1,
                      palindromes = 1,
                      maf_threshold = 0.3,
                      access_token = ieugwasr::check_access_token(),
                      splitsize = 10000,
                      proxy_splitsize = 500) #非子宫内膜样组织学
outcome3<- extract_outcome_data(snps = anxiety$SNP, 
                      outcomes = "ebi-a-GCST006464",
                      proxies = TRUE,
                      rsq = 0.8,
                      align_alleles = 1,
                      palindromes = 1,
                      maf_threshold = 0.3,
                      access_token = ieugwasr::check_access_token(),
                      splitsize = 10000,
                      proxy_splitsize = 500) #子宫内膜癌
H_data1 <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome1
)
H_data2 <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome2
)
H_data3 <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome3
)
mrresu1<- mr(H_data1)
mrresu2<- mr(H_data2)
mrresu3<- mr(H_data3)
bx<- as.vector(as.matrix(anxiety$`OR or BETA`)) #anxiety
bxse<- get_se(bx,anxiety$`P-VALUE`) #anxiety
pval<- as.vector(as.matrix(anxietySNP$`P-VALUE`))
anxietyse<- get_se(bx, pval)
SNP<- as.vector(as.matrix(exposure_data$`SNP`))
by1<- as.vector(as.matrix(outcome1$beta.outcome))
byse1<- as.vector(as.vector(as.matrix(outcome1$se.outcome)))#子宫内膜样组织学
by2<- as.vector(as.matrix(outcome2$beta.outcome))
byse2<- as.vector(as.vector(as.matrix(outcome2$se.outcome)))#非子宫内膜样组织学
by3<- as.vector(as.matrix(outcome3$beta.outcome))
byse3<- as.vector(as.vector(as.matrix(outcome3$se.outcome)))#子宫内膜癌
object1<- - mr_input(bx, bxse, by1, byse1, exposure = "anxiety", outcome="EDC", snps= SNP)
mr_ivw(object1, robust=TRUE)
mr_egger(object1, robust=TRUE)
mr_median(object1,weighting = "weighted",distribution = "normal",
alpha = 0.05,iterations = 10000,seed = 314159265)
object2<- mr_input(bx, bxse, by2, byse2, exposure = "anxiety", outcome="EDC", snps= SNP)
mr_ivw(object2, robust=TRUE)
mr_egger(object2, robust=TRUE)
mr_median(object2,weighting = "weighted",distribution = "normal",
alpha = 0.05,iterations = 10000,seed = 314159265)
object<- mr_input(bx, bxse, by3, byse3, exposure = "anxiety", outcome="EDC", snps= SNP)
mr_ivw(object3, robust=TRUE)
mr_egger(object3, robust=TRUE)
mr_median(object3,weighting = "weighted",distribution = "normal",
alpha = 0.05,iterations = 10000,seed = 314159265)

2.Source code for MR analysis of depression on endometrial cancer

library(MendelianRandomization)
library(readxl)
depressionSNP <- read_excel("D:/我的文件/情绪和妇产科疾病
                            （约稿/抑郁GWAS欧洲/depressionSNP.xlsx")
outcome1<- extract_outcome_data(snps = exposure_data$SNP, 
                                outcomes = "ebi-a-GCST006465") #子宫内膜样组织学
outcome2<- extract_outcome_data(snps = exposure_data$SNP, 
                                outcomes = "ebi-a-GCST006466") #非子宫内膜样组织学
outcome3<- extract_outcome_data(snps = exposure_data$SNP, 
                                outcomes = "ebi-a-GCST006464") #子宫内膜癌
bx<- as.vector(as.matrix(exposure_data$beta.exposure))
pval<- as.vector(as.matrix(exposure_data$pval.exposure))
depressionse<- get_se(bx, pval)
by1<- as.vector(as.matrix(outcome1$beta.outcome))
byse1<- as.vector(as.vector(as.matrix(outcome1$se.outcome)))
by2<- as.vector(as.matrix(outcome2$beta.outcome))
byse2<- as.vector(as.vector(as.matrix(outcome2$se.outcome)))
by3<- as.vector(as.matrix(outcome3$beta.outcome))
byse3<- as.vector(as.vector(as.matrix(outcome3$se.outcome)))
object1<- mr_input(bx, depressionse, by1, byse1, exposure = "MDD", outcome="EDC") #snps= SNP)
mr_ivw(object1, robust=TRUE)
mr_egger(object1, robust=TRUE)
mr_median(object1,weighting = "weighted",distribution = "normal",
alpha = 0.05,iterations = 10000,seed = 314159265)
object2<- mr_input(bx, depressionse, by2, byse2, exposure = "MDD", outcome="EDC") #snps= SNP)
mr_ivw(object2, robust=TRUE)
mr_egger(object2, robust=TRUE)
mr_median(object2,weighting = "weighted",distribution = "normal",
alpha = 0.05,iterations = 10000,seed = 314159265)
object3<- mr_input(bx, depressionse, by3, byse3, exposure = "MDD", outcome="EDC") #snps= SNP)
mr_ivw(object3, robust=TRUE)
mr_egger(object3, robust=TRUE)
mr_median(object3,weighting = "weighted",distribution = "normal",
alpha = 0.05,iterations = 10000,seed = 314159265)

3. Source code for multivariate MR analysis
Library(TwoSampleMR)
mv<- c("ieu-b-102", "ukb-a-248","ukb-b-17422","ieu-a-1090", "finn-b-E4_POCS")
MVMRexposure<- mv_extract_exposures(mv)
mvoutcome_data1<- extract_outcome_data(snps = MVMRexposure$SNP, 
                                      outcomes = "ebi-a-GCST006464")
mvoutcome_data2<- extract_outcome_data(snps = MVMRexposure$SNP, 
                                      outcomes = "ebi-a-GCST006465")
mvoutcome_data3<- extract_outcome_data(snps = MVMRexposure$SNP, 
                                      outcomes = "ebi-a-GCST006466")
mvH.data1<-mv_harmonise_data(MVMRexposure, mvoutcome_data1, harmonise_strictness = 2)
mvH.data2<-mv_harmonise_data(MVMRexposure, mvoutcome_data2, harmonise_strictness = 2)
mvH.data3<-mv_harmonise_data(MVMRexposure, mvoutcome_data3, harmonise_strictness = 2)

mv_multiple(mvH.data1)
mv_multiple(mvH.data2)
mv_multiple(mvH.data3)
