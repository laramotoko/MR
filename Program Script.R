##Project infomation####
#HU yuan, 2021-05-04, Capstone Project: Urinary Sodium with Type 2 Diabetes and Its Risk Factors: A Two-Sample Mendelian Randomization
#Exposure: Urinary sodium
#Outcome: Type 2 diabetes

getwd()
setwd("C:/Users/lara5/Documents/MPH/PBHT6900F-Master Capstone Project/Data")

##Labraries####
library("dplyr")
library("dtplyr")
library("tidyr")
library("MendelianRandomization")
library("TwoSampleMR")
library("MRInstruments")
library("RadialMR")
library("MRPRESSO")
library("MRPracticals")
library("MRMix")
library("mr.raps")
library("ggplot2")


#@UNa, T2DM unadjusted for BMI####
#Input data####
#Exposure
library(readxl)
E1_UNa_SNPs_datainput <- read_excel("C:/Users/lara5/Documents/MPH/PBHT6900F-Master Capstone Project/Data/E1-UNa SNPs-datainput.xlsx", 
                                    col_types = c("text", "numeric", "text", 
                                                  "text", "text", "numeric", "text", 
                                                  "numeric", "numeric", "numeric"))

#Outcome
library(readxl)
O2_T2DM_2018 <- read.delim("~/MPH/PBHT6900F-Master Capstone Project/Data/O2-Mahajan.NatGenet2018b.T2D.European.txt")
colnames(O2_T2DM_2018)


#Merge####
#New variable for exposure SNP, NSNP=CHR:BP
E1_UNa_SNPs_datainput$NSNP <- paste0(E1_UNa_SNPs_datainput$CHR, ":", E1_UNa_SNPs_datainput$BP)
#Change the name of outcome variable
O2_T2DM_2018$NSNP <- O2_T2DM_2018$SNP


#New dataset SNP_50_UNa_T2DM
SNP_50_UNa_T2DM <- merge(x=E1_UNa_SNPs_datainput, y=O2_T2DM_2018, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Harmonization####
#Variable names£¬adjust beta, effective allele£¬change the valuence of beta according to effective allele
#Identify variants that do not share the same allele pair between datasets, and either correct this if possible or eliminate such variants
SNP_50_UNa_T2DM <- SNP_50_UNa_T2DM %>%
  mutate(effective_allele = ALLELE1,
         other_allele = ALLELE0,
         beta_UNa = BETA,
         se_UNa = SE.x,
         beta_T2DM = ifelse(EA == effective_allele, Beta, 0-Beta),
         se_T2DM = SE.y)


#Analysis####
#Input
MR_SNP_50 <- mr_input(bx = SNP_50_UNa_T2DM$beta_UNa, bxse = SNP_50_UNa_T2DM$se_UNa,
                      by = SNP_50_UNa_T2DM$beta_T2DM, byse = SNP_50_UNa_T2DM$se_T2DM,
                      exposure = "UNa", outcome = "T2DM", snps = SNP_50_UNa_T2DM$SNP.x,
                      effect_allele = SNP_50_UNa_T2DM$effective_allele,
                      other_allele = SNP_50_UNa_T2DM$other_allele,
                      eaf = SNP_50_UNa_T2DM$A1FREQ)

#Median, IVW, Egger
mr_allmethods(MR_SNP_50)
#P value: Simple median=0.000 IVW=0.019, MR-Egger=0.528, intercept=0.208


#@T2DM adjusted for BMI####
#Input####
library(readxl)
O1_T2DM_2018_BMI <- read.delim("~/MPH/PBHT6900F-Master Capstone Project/Data/O1-Mahajan.NatGenet2018b.T2Dbmiadj.European.txt")
colnames(O1_T2DM_2018_BMI)


#Merge####
#Change the name of outcome variable
O1_T2DM_2018_BMI$NSNP <- O1_T2DM_2018_BMI$SNP

#New dataset SNP_50_UNa_T2DM_BMI
SNP_50_UNa_T2DM_BMI <- merge(x=E1_UNa_SNPs_datainput, y=O1_T2DM_2018_BMI, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Harmonization####
SNP_50_UNa_T2DM_BMI <- SNP_50_UNa_T2DM_BMI %>%
  mutate(effective_allele = ALLELE1,
         other_allele = ALLELE0,
         beta_UNa = BETA,
         se_UNa = SE.x,
         beta_T2DM_BMI = ifelse(EA == effective_allele, Beta, 0-Beta),
         se_T2DM_BMI = SE.y)


#Analysis####
MR_SNP_50_BMI <- mr_input(bx = SNP_50_UNa_T2DM_BMI$beta_UNa, bxse = SNP_50_UNa_T2DM_BMI$se_UNa,
                          by = SNP_50_UNa_T2DM_BMI$beta_T2DM_BMI, byse = SNP_50_UNa_T2DM_BMI$se_T2DM_BMI,
                          exposure = "UNa", outcome = "T2DM_BMI", snps = SNP_50_UNa_T2DM_BMI$SNP.x,
                          effect_allele = SNP_50_UNa_T2DM_BMI$effective_allele,
                          other_allele = SNP_50_UNa_T2DM_BMI$other_allele,
                          eaf = SNP_50_UNa_T2DM_BMI$A1FREQ)


#Median, IVW, Egger
mr_allmethods(MR_SNP_50_BMI)
#P value: Simple median=0.001 IVW=0.369, MR-Egger=0.079, intercept=0.040



#@T2DM without UKB####
#Input data####
#Exposure
library(readxl)
E1_UNa_SNPs_datainput <- read_excel("C:/Users/lara5/Documents/MPH/PBHT6900F-Master Capstone Project/Data/E1-UNa SNPs-datainput.xlsx", 
                                    col_types = c("text", "numeric", "text", 
                                                  "text", "text", "numeric", "text", 
                                                  "numeric", "numeric", "numeric"))

#Outcome
library(readxl)
O3_T2DM_NUKB <- read.delim("~/MPH/PBHT6900F-Master Capstone Project/Data/METAANALYSIS_DIAGRAM_SE1.txt")
colnames(O3_T2DM_NUKB)


#Merge####
#New variable for exposure SNP, NSNP=CHR:BP
E1_UNa_SNPs_datainput$NSNP<-paste0(E1_UNa_SNPs_datainput$CHR, ":", E1_UNa_SNPs_datainput$BP)
#Change the name of outcome variable
O3_T2DM_NUKB$NSNP<-O3_T2DM_NUKB$Chr.Position

#New data set SNP_50_UNa_T2DM_NUKB
SNP_50_UNa_T2DM_NUKB<-merge(x=E1_UNa_SNPs_datainput, y=O3_T2DM_NUKB, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Harmonization####
SNP_50_UNa_T2DM_NUKB<-SNP_50_UNa_T2DM_NUKB %>%
  mutate(effective_allele = ALLELE1,
         other_allele = ALLELE0,
         beta_UNa = BETA,
         se_UNa = SE,
         beta_T2DM = ifelse(Allele1 == effective_allele, Effect, 0-Effect),
         se_T2DM = StdErr)


#Analysis####
#Input
MR_SNP_50_NUKB<-mr_input(bx = SNP_50_UNa_T2DM_NUKB$beta_UNa, bxse = SNP_50_UNa_T2DM_NUKB$se_UNa,
                         by = SNP_50_UNa_T2DM_NUKB$beta_T2DM, byse = SNP_50_UNa_T2DM_NUKB$se_T2DM,
                         exposure = "UNa", outcome = "T2DM", snps = SNP_50_UNa_T2DM_NUKB$SNP,
                         effect_allele = SNP_50_UNa_T2DM_NUKB$effective_allele,
                         other_allele = SNP_50_UNa_T2DM_NUKB$other_allele,
                         eaf = SNP_50_UNa_T2DM_NUKB$A1FREQ)

#Median, IVW, Egger
mr_allmethods(MR_SNP_50_NUKB)



#@T2DM without UKB adjusted with BMI####
#Input data####
#Exposure
library(readxl)
E1_UNa_SNPs_datainput <- read_excel("C:/Users/lara5/Documents/MPH/PBHT6900F-Master Capstone Project/Data/E1-UNa SNPs-datainput.xlsx", 
                                    col_types = c("text", "numeric", "text", 
                                                  "text", "text", "numeric", "text", 
                                                  "numeric", "numeric", "numeric"))

#Outcome
library(readxl)
O3_T2DM_NUKB_BMI <- read.delim("~/MPH/PBHT6900F-Master Capstone Project/Data/METAANALYSIS_DIAGRAM_SE1.BMI.txt")
colnames(O3_T2DM_NUKB_BMI)


#Merge####
#New variable for exposure SNP, NSNP=CHR:BP
E1_UNa_SNPs_datainput$NSNP<-paste0(E1_UNa_SNPs_datainput$CHR, ":", E1_UNa_SNPs_datainput$BP)
#Change the name of outcome variable
O3_T2DM_NUKB_BMI$NSNP<-O3_T2DM_NUKB_BMI$Chr.Position

#New data set SNP_50_UNa_T2DM_NUKB
SNP_50_UNa_T2DM_NUKB_BMI<-merge(x=E1_UNa_SNPs_datainput, y=O3_T2DM_NUKB_BMI, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Harmonization####
SNP_50_UNa_T2DM_NUKB_BMI<-SNP_50_UNa_T2DM_NUKB_BMI %>%
  mutate(effective_allele = ALLELE1,
         other_allele = ALLELE0,
         beta_UNa = BETA,
         se_UNa = SE,
         beta_T2DM = ifelse(Allele1 == effective_allele, Effect, 0-Effect),
         se_T2DM = StdErr)


#Analysis####
#Input
MR_SNP_50_NUKB_BMI<-mr_input(bx = SNP_50_UNa_T2DM_NUKB_BMI$beta_UNa, bxse = SNP_50_UNa_T2DM_NUKB_BMI$se_UNa,
                             by = SNP_50_UNa_T2DM_NUKB_BMI$beta_T2DM, byse = SNP_50_UNa_T2DM_NUKB_BMI$se_T2DM,
                             exposure = "UNa", outcome = "T2DM", snps = SNP_50_UNa_T2DM_NUKB_BMI$SNP,
                             effect_allele = SNP_50_UNa_T2DM_NUKB_BMI$effective_allele,
                             other_allele = SNP_50_UNa_T2DM_NUKB_BMI$other_allele,
                             eaf = SNP_50_UNa_T2DM_NUKB_BMI$A1FREQ)

#Median, IVW, Egger
mr_allmethods(MR_SNP_50_NUKB_BMI)

#@Multivariable-BMI#### 
#BMI Input####
library(readxl)
M1_BMI <- read_excel("C:/Users/lara5/Documents/MPH/PBHT6900F-Master Capstone Project/Data/M1_BMI_SNP_2015_Euro.xlsx", 
                     col_types = c("text", "numeric", "text", "numeric", 
                                   "text", "text", "numeric", "numeric", 
                                   "numeric", "numeric", "text", "numeric"))


#UNa GWAS Input####  
E1_UNa_GWAS <- read_excel("C:/Users/lara5/Documents/MPH/PBHT6900F-Master Capstone Project/Data/PazokiR_prePMID_Sodium.GWAS.xlsx")


#Merge####
#Merge UNa and BMI####
#UNa_BMI=SNPs UNa+SNPsBMI
Mer_UNa <- E1_UNa_SNPs_datainput %>%
  select(SNP, CHR, BP, "Gene nearest to Locus")
colnames(Mer_UNa) < -c("SNP", "CHR", "BP", "GEN")

Mer_BMI < -M1_BMI %>%
  select(SNP, CHR, BP_37, Nearest)
colnames(Mer_BMI) <- c("SNP", "CHR", "BP", "GEN")

#change variables type
Mer_BMI$BP <- as.character(Mer_BMI$BP)

#Bind
SNP_UNa_BMI <- rbind(Mer_UNa, Mer_BMI)


#Merge SNP_UNa_BMI with UNa GWAS####
Mer_SNP_UNa_BMI_GUNa <- merge(x=SNP_UNa_BMI, y=E1_UNa_GWAS, by.x = "SNP", by.y = "SNP", all.x = TRUE)


#Merge with BMI GWAS####
#Import
M1_BMI_GWAS <- read.delim("~/MPH/PBHT6900F-Master Capstone Project/Data/BMI_2015_Euro.txt")

Mer_SNP_UNa_BMI_GUNa_GBMI <- merge(x=Mer_SNP_UNa_BMI_GUNa, y=M1_BMI_GWAS, by.x = "SNP", by.y = "SNP", all.x = TRUE)


#1.Merge with T2DM GWAS####
#CHR:BP
Mer_SNP_UNa_BMI_GUNa_GBMI$NSNP <- paste0(as.character(Mer_SNP_UNa_BMI_GUNa_GBMI$CHR.x), ":", Mer_SNP_UNa_BMI_GUNa_GBMI$BP.x)

Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM <- merge(x=Mer_SNP_UNa_BMI_GUNa_GBMI, y=O2_T2DM_2018, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Looking fo missing data####
table(is.na(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM$ALLELE1))


#Harmonization####
Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM <- Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM %>%
  mutate(SNP=SNP.x, 
         effect_allele=ALLELE1,
         other_allele=ALLELE0,
         beta_UNa=BETA,
         beta_BMI=ifelse(A1 == effect_allele, b, 0-b),
         beta_T2DM=ifelse(EA == effect_allele, Beta, 0-Beta),
         se_UNa=SE.x,
         se_BMI=se,
         se_T2DM=SE.y)


#Analysis####
#Input####
Input_MR <- mr_mvinput(
  bx = as.matrix(cbind(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM$beta_UNa, Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM$beta_BMI)),
  bxse = as.matrix(cbind(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM$se_UNa, Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM$se_BMI)),
  by = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM$beta_T2DM,
  byse = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM$se_T2DM,
  exposure = c("UNa", "BMI"),
  outcome = "T2DM",
  snps = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM$SNP,
  effect_allele = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM$effect_allele,
  other_allele = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM$other_allele,
  eaf = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM$A1FREQ)


#MR####
mr_mvmedian(Input_MR)
mr_mvegger(Input_MR)
mr_mvivw(Input_MR)


#2. with T2DM_BMI GWAS####
#CHR:BP
Mer_SNP_UNa_BMI_GUNa_GBMI$NSNP <- paste0(as.character(Mer_SNP_UNa_BMI_GUNa_GBMI$CHR.x), ":", Mer_SNP_UNa_BMI_GUNa_GBMI$BP.x)

Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI <- merge(x=Mer_SNP_UNa_BMI_GUNa_GBMI, y=O1_T2DM_2018_BMI, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Looking fo missing data####
table(is.na(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI$ALLELE1))


#Harmonization####
Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI <- Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI %>%
  mutate(SNP=SNP.x, 
         effect_allele=ALLELE1,
         other_allele=ALLELE0,
         beta_UNa=BETA,
         beta_BMI=ifelse(A1 == effect_allele, b, 0-b),
         beta_T2DM=ifelse(EA == effect_allele, Beta, 0-Beta),
         se_UNa=SE.x,
         se_BMI=se,
         se_T2DM=SE.y)


#Analysis####
#Input####
Input_MR <- mr_mvinput(
  bx = as.matrix(cbind(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI$beta_UNa, Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI$beta_BMI)),
  bxse = as.matrix(cbind(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI$se_UNa, Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI$se_BMI)),
  by = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI$beta_T2DM,
  byse = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI$se_T2DM,
  exposure = c("UNa", "BMI"),
  outcome = "T2DM",
  snps = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI$SNP,
  effect_allele = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI$effect_allele,
  other_allele = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI$other_allele,
  eaf = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_BMI$A1FREQ)


#MR####
mr_mvmedian(Input_MR)
mr_mvegger(Input_MR)
mr_mvivw(Input_MR)


#3.Merge with T2DM_NUKB GWAS####
#CHR:BP
Mer_SNP_UNa_BMI_GUNa_GBMI$NSNP <- paste0(as.character(Mer_SNP_UNa_BMI_GUNa_GBMI$CHR.x), ":", Mer_SNP_UNa_BMI_GUNa_GBMI$BP.x)

Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB <- merge(x=Mer_SNP_UNa_BMI_GUNa_GBMI, y=O3_T2DM_NUKB, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Looking fo missing data####
table(is.na(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB$ALLELE1))


#Harmonization####
Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB <- Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB %>%
  mutate(SNP=SNP, 
         effect_allele=ALLELE1,
         other_allele=ALLELE0,
         beta_UNa=BETA,
         beta_BMI=ifelse(A1 == effect_allele, b, 0-b),
         beta_T2DM=ifelse(Allele1 == effect_allele, Effect, 0-Effect),
         se_UNa=SE,
         se_BMI=se,
         se_T2DM=StdErr)


#Analysis####
#Input####
Input_MR <- mr_mvinput(
  bx = as.matrix(cbind(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB$beta_UNa, Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB$beta_BMI)),
  bxse = as.matrix(cbind(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB$se_UNa, Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB$se_BMI)),
  by = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB$beta_T2DM,
  byse = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB$se_T2DM,
  exposure = c("UNa", "BMI"),
  outcome = "T2DM",
  snps = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB$SNP,
  effect_allele = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB$effect_allele,
  other_allele = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB$other_allele,
  eaf = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB$A1FREQ)


#MR####
mr_mvmedian(Input_MR)
mr_mvegger(Input_MR)
mr_mvivw(Input_MR)


#4.Merge with T2DM_NUKB_BMI GWAS####
#CHR:BP
Mer_SNP_UNa_BMI_GUNa_GBMI$NSNP <- paste0(as.character(Mer_SNP_UNa_BMI_GUNa_GBMI$CHR.x), ":", Mer_SNP_UNa_BMI_GUNa_GBMI$BP.x)

Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI <- merge(x=Mer_SNP_UNa_BMI_GUNa_GBMI, y=O3_T2DM_NUKB_BMI, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Looking fo missing data####
table(is.na(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI$ALLELE1))


#Harmonization####
Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI <- Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI %>%
  mutate(SNP=SNP, 
         effect_allele=ALLELE1,
         other_allele=ALLELE0,
         beta_UNa=BETA,
         beta_BMI=ifelse(A1 == effect_allele, b, 0-b),
         beta_T2DM=ifelse(Allele1 == effect_allele, Effect, 0-Effect),
         se_UNa=SE,
         se_BMI=se,
         se_T2DM=StdErr)


#Analysis####
#Input####
Input_MR <- mr_mvinput(
  bx = as.matrix(cbind(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI$beta_UNa, Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI$beta_BMI)),
  bxse = as.matrix(cbind(Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI$se_UNa, Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI$se_BMI)),
  by = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI$beta_T2DM,
  byse = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI$se_T2DM,
  exposure = c("UNa", "BMI"),
  outcome = "T2DM",
  snps = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI$SNP,
  effect_allele = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI$effect_allele,
  other_allele = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI$other_allele,
  eaf = Mer_SNP_UNa_BMI_GUNa_GBMI_GT2DM_NUKB_BMI$A1FREQ)


#MR####
mr_mvmedian(Input_MR)
mr_mvegger(Input_MR)
mr_mvivw(Input_MR)


#@Multivariable-eGFR####
#eGFR Input####
M2_eGFR_SNPs <- read_excel("~/MPH/PBHT6900F-Master Capstone Project/Data/M2_eGFR_SNPs.xlsx", 
                           col_types = c("text", "text", "text", 
                                         "text", "text", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric"))

M2_eGFR_GWAS <- read.csv("~/MPH/PBHT6900F-Master Capstone Project/Data/M2_eGFR_GWAS.txt", sep="")


#Merge####
#Merge UNa and eGFR####
Mer_eGFR <- M2_eGFR_SNPs

Mer_eGFR <- Mer_eGFR %>% separate(NSNP, into = c("CHR", "BP"), sep = ":")

Mer_eGFR <- Mer_eGFR %>%
  select(SNP, CHR, BP, Nearest)
colnames(Mer_eGFR) <- c("SNP", "CHR", "BP", "GEN")
Mer_eGFR$BP <- as.character(Mer_eGFR$BP)

#bind
SNP_UNa_eGFR <- rbind(Mer_UNa, Mer_eGFR)


#Merge SNP_UNa_eGFR with UNa GWAS####
Mer_SNP_UNa_eGFR_GUNa <- merge(x=SNP_UNa_eGFR, y=E1_UNa_GWAS, by.x = "SNP", by.y = "SNP", all.x = TRUE)


#Merge with eGFR GWAS####
colnames(M2_eGFR_GWAS[3])
colnames(M2_eGFR_GWAS)[3] <- "SNP"

Mer_SNP_UNa_eGFR_GUNa_GeGFR <- merge(x=Mer_SNP_UNa_eGFR_GUNa, y=M2_eGFR_GWAS, by.x= "SNP", by.y= "SNP", all.x=TRUE)


#1.Merge with T2DM GWAS####
#CHR:BP
Mer_SNP_UNa_eGFR_GUNa_GeGFR$NSNP <- paste0(as.character(Mer_SNP_UNa_eGFR_GUNa_GeGFR$CHR.x), ":", Mer_SNP_UNa_eGFR_GUNa_GeGFR$BP.x)
Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM <- merge(x=Mer_SNP_UNa_eGFR_GUNa_GeGFR, y=O2_T2DM_2018, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Looking for missing data####
table(is.na(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM$ALLELE1))


#Harmonization####
Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM <- Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM %>%
  mutate(SNP=SNP.x, 
         effect_allele=ALLELE1,
         other_allele=ALLELE0,
         beta_UNa=BETA,
         beta_eGFR=ifelse(Allele1 == effect_allele, Effect, 0-Effect),
         beta_T2DM=ifelse(EA == effect_allele, Beta, 0-Beta),
         se_UNa=SE.x,
         se_eGFR=StdErr,
         se_T2DM=SE.y)


#Analysis####
#Input####
Input_MR <- mr_mvinput(
  bx = as.matrix(cbind(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM$beta_UNa, Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM$beta_eGFR)),
  bxse = as.matrix(cbind(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM$se_UNa, Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM$se_eGFR)),
  by = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM$beta_T2DM,
  byse = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM$se_T2DM,
  exposure = c("UNa", "eGFR"),
  outcome = "T2DM",
  snps = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM$SNP,
  effect_allele = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM$effect_allele,
  other_allele = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM$other_allele,
  eaf = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM$A1FREQ)


#MR####
mr_mvmedian(Input_MR)
mr_mvegger(Input_MR)
mr_mvivw(Input_MR)


#2.Merge with T2DM_BMI GWAS####
#CHR:BP
Mer_SNP_UNa_eGFR_GUNa_GeGFR$NSNP <- paste0(as.character(Mer_SNP_UNa_eGFR_GUNa_GeGFR$CHR.x), ":", Mer_SNP_UNa_eGFR_GUNa_GeGFR$BP.x)
Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI <- merge(x=Mer_SNP_UNa_eGFR_GUNa_GeGFR, y=O1_T2DM_2018_BMI, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Looking for missing data####
table(is.na(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI$ALLELE1))


#Harmonization####
Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI <- Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI %>%
  mutate(SNP=SNP.x, 
         effect_allele=ALLELE1,
         other_allele=ALLELE0,
         beta_UNa=BETA,
         beta_eGFR=ifelse(Allele1 == effect_allele, Effect, 0-Effect),
         beta_T2DM=ifelse(EA == effect_allele, Beta, 0-Beta),
         se_UNa=SE.x,
         se_eGFR=StdErr,
         se_T2DM=SE.y)


#Analysis####
#Input####
Input_MR <- mr_mvinput(
  bx = as.matrix(cbind(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI$beta_UNa, Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI$beta_eGFR)),
  bxse = as.matrix(cbind(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI$se_UNa, Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI$se_eGFR)),
  by = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI$beta_T2DM,
  byse = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI$se_T2DM,
  exposure = c("UNa", "eGFR"),
  outcome = "T2DM",
  snps = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI$SNP,
  effect_allele = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI$effect_allele,
  other_allele = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI$other_allele,
  eaf = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_BMI$A1FREQ)


#MR####
mr_mvmedian(Input_MR)
mr_mvegger(Input_MR)
mr_mvivw(Input_MR)


#3.Merge with T2DM_NUKB GWAS####
#CHR:BP
Mer_SNP_UNa_eGFR_GUNa_GeGFR$NSNP <- paste0(as.character(Mer_SNP_UNa_eGFR_GUNa_GeGFR$CHR.x), ":", Mer_SNP_UNa_eGFR_GUNa_GeGFR$BP.x)
Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB <- merge(x=Mer_SNP_UNa_eGFR_GUNa_GeGFR, y=O3_T2DM_NUKB, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Looking for missing data####
table(is.na(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB$ALLELE1))


#Harmonization####
Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB <- Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB %>%
  mutate(SNP=SNP, 
         effect_allele=ALLELE1,
         other_allele=ALLELE0,
         beta_UNa=BETA,
         beta_eGFR=ifelse(Allele1.x == effect_allele, Effect.x, 0-Effect.x),
         beta_T2DM=ifelse(Allele1.y == effect_allele, Effect.y, 0-Effect.y),
         se_UNa=SE,
         se_eGFR=StdErr.x,
         se_T2DM=StdErr.y)


#Analysis####
#Input####
Input_MR <- mr_mvinput(
  bx = as.matrix(cbind(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB$beta_UNa, Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB$beta_eGFR)),
  bxse = as.matrix(cbind(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB$se_UNa, Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB$se_eGFR)),
  by = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB$beta_T2DM,
  byse = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB$se_T2DM,
  exposure = c("UNa", "eGFR"),
  outcome = "T2DM",
  snps = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB$SNP,
  effect_allele = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB$effect_allele,
  other_allele = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB$other_allele,
  eaf = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB$A1FREQ)


#MR####
mr_mvmedian(Input_MR)
mr_mvegger(Input_MR)
mr_mvivw(Input_MR)


#4.Merge with T2DM_NUKB_BMI GWAS####
#CHR:BP
Mer_SNP_UNa_eGFR_GUNa_GeGFR$NSNP <- paste0(as.character(Mer_SNP_UNa_eGFR_GUNa_GeGFR$CHR.x), ":", Mer_SNP_UNa_eGFR_GUNa_GeGFR$BP.x)
Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI <- merge(x=Mer_SNP_UNa_eGFR_GUNa_GeGFR, y=O3_T2DM_NUKB_BMI, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)


#Looking for missing data####
table(is.na(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI$ALLELE1))


#Harmonization####
Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI <- Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI %>%
  mutate(SNP=SNP, 
         effect_allele=ALLELE1,
         other_allele=ALLELE0,
         beta_UNa=BETA,
         beta_eGFR=ifelse(Allele1.x == effect_allele, Effect.x, 0-Effect.x),
         beta_T2DM=ifelse(Allele1.y == effect_allele, Effect.y, 0-Effect.y),
         se_UNa=SE,
         se_eGFR=StdErr.x,
         se_T2DM=StdErr.y)


#Analysis####
#Input####
Input_MR <- mr_mvinput(
  bx = as.matrix(cbind(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI$beta_UNa, Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI$beta_eGFR)),
  bxse = as.matrix(cbind(Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI$se_UNa, Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI$se_eGFR)),
  by = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI$beta_T2DM,
  byse = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI$se_T2DM,
  exposure = c("UNa", "eGFR"),
  outcome = "T2DM",
  snps = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI$SNP,
  effect_allele = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI$effect_allele,
  other_allele = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI$other_allele,
  eaf = Mer_SNP_UNa_eGFR_GUNa_GeGFR_GT2DM_NUKB_BMI$A1FREQ)


#MR####
mr_mvmedian(Input_MR)
mr_mvegger(Input_MR)
mr_mvivw(Input_MR)


#@Multivariable-eGFR+BMI####
#Merge####
#Merge UNa and BMI####
Mer_UNa <- E1_UNa_SNPs_datainput %>%
  select(SNP, CHR, BP, "Gene nearest to Locus")
colnames(Mer_UNa) < -c("SNP", "CHR", "BP", "GEN")

Mer_BMI < -M1_BMI %>%
  select(SNP, CHR, BP_37, Nearest)
colnames(Mer_BMI) <- c("SNP", "CHR", "BP", "GEN")

#bind####
SNP_UNa_BMI_eGFR <- rbind(SNP_UNa_BMI, Mer_eGFR)


#Merge SNAP_UNa_BMI_eGFR with UNa GWAS####
Mer_SNP_UNa_BMI_eGFR_GUNa <- merge(x=SNP_UNa_BMI_eGFR, y=E1_UNa_GWAS, by.x = "SNP", by.y = "SNP", all.x = TRUE)

#Merge with BMI GWAS####
Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI <- merge(x=Mer_SNP_UNa_BMI_eGFR_GUNa, y=M1_BMI_GWAS, by.x = "SNP", by.y = "SNP", all.x = TRUE)

#Merge with eGFR GWAS####
Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR <- merge(x=Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI, y=M2_eGFR_GWAS, by.x = "SNP", by.y = "SNP", all.x = TRUE)

#1.Merge with T2DM GWAS####
Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR$NSNP <- paste0(as.character(Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR$CHR.x), ":", Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR$BP.x)
Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM <- merge(x=Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR, y=O2_T2DM_2018, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)

#Looking for missing data####
table(is.na(Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$A1))


#Harmonization####
Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM <- Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM %>%
  mutate(SNP=NSNP, 
         effect_allele=ALLELE1,
         other_allele=ALLELE0,
         beta_UNa=BETA,
         beta_BMI=ifelse(A1 == effect_allele, b, 0-b),
         beta_eGFR=ifelse(Allele1 == effect_allele, Effect, 0-Effect),
         beta_T2DM=ifelse(EA == effect_allele, Beta, 0-Beta),
         se_UNa=SE.x,
         se_BMI=se,
         se_eGFR=StdErr,
         se_T2DM=SE.y)


#Analysis####
#Input####
Input_MR <- mr_mvinput(
  bX = as.matrix(cbind(Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$beta_UNa, Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$beta_BMI, Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$beta_eGFR)),
  bXse = as.matrix(cbind(Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$se_UNa, Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$se_BMI, Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$se_eGFR)),
  bY = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$beta_T2DM,
  bYse = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$se_T2DM,
  exposure = c("UNa", "BMI","eGFR"),
  outcome = "T2DM",
  snps = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$SNP,
  effect_allele = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$effect_allele,
  other_allele = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$other_allele,
  eaf = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$A1FREQ)

#MR####
mr_mvmedian(Input_MR)
mr_mvegger(Input_MR)
mr_mvivw(Input_MR)

r_input <- format_mvmr(
  BXGs = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM[, c("beta_UNa", "beta_BMI", "beta_eGFR")],
  BYG = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$beta_T2DM,
  seBXGs = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM[, c("se_UNa", "se_BMI", "se_eGFR")],
  seBYG = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$se_T2DM,
  RSID = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM$SNP
)
ivw_mvmr(r_input)
pleiotropy_mvmr(r_input, 0)




#2.Merge with T2DM_NUKB GWAS#### 
Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR$NSNP <- paste0(as.character(Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR$CHR.x), ":", Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR$BP.x)
Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB <- merge(x=Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR, y=O3_T2DM_NUKB, by.x = "NSNP", by.y = "NSNP", all.x = TRUE)

#Looking for missing data####
table(is.na(Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$ALLELE1))


#Harmonization####
Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB <- Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB %>%
  mutate(SNP=NSNP, 
         effect_allele=ALLELE1,
         other_allele=ALLELE0,
         beta_UNa=BETA,
         beta_BMI=ifelse(A1 == effect_allele, b, 0-b),
         beta_eGFR=ifelse(Allele1.x == effect_allele, Effect.x, 0-Effect.x),
         beta_T2DM=ifelse(Allele1.y == effect_allele, Effect.y, 0-Effect.y),
         se_UNa=SE,
         se_BMI=se,
         se_eGFR=StdErr.x,
         se_T2DM=StdErr.y)


#Analysis####
#Input####
Input_MR <- mr_mvinput(
  bX = as.matrix(cbind(Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$beta_UNa, Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$beta_BMI, Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$beta_eGFR)),
  bXse = as.matrix(cbind(Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$se_UNa, Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$se_BMI, Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$se_eGFR)),
  bY = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$beta_T2DM,
  bYse = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$se_T2DM,
  exposure = c("UNa", "BMI","eGFR"),
  outcome = "T2DM_NUKB",
  snps = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$SNP,
  effect_allele = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$effect_allele,
  other_allele = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$other_allele,
  eaf = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$A1FREQ)


#MR####
mr_mvmedian(Input_MR)
mr_mvegger(Input_MR)
mr_mvivw(Input_MR)

r_input <- format_mvmr(
  BXGs = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB[, c("beta_UNa", "beta_BMI", "beta_eGFR")],
  BYG = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$beta_T2DM,
  seBXGs = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB[, c("se_UNa", "se_BMI", "se_eGFR")],
  seBYG = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$se_T2DM,
  RSID = Mer_SNP_UNa_BMI_eGFR_GUNa_GBMI_GeGFR_GT2DM_NUKB$SNP
)
ivw_mvmr(r_input)
pleiotropy_mvmr(r_input, 0)
snpcov_mvmr(r_input)

#@Bidirectional####
#T2DM SNPs
#Merge####

#Change the name of outcome variable
T2DM_SNP$SNP <- T2DM_SNP$RSID
T2DM_SNP$beta <- log(T2DM_SNP$OR)
T2DM_SNP <- T2DM_SNP %>% separate(`CI 95%`, into = c("upper", "lower"), sep = "-")
T2DM_SNP$upper <- as.numeric(T2DM_SNP$upper)
T2DM_SNP$lower <- as.numeric(T2DM_SNP$lower)
T2DM_SNP$upper <- log(T2DM_SNP$upper)
T2DM_SNP$lower <- log(T2DM_SNP$lower)
T2DM_SNP$SE <- (T2DM_SNP$lower- T2DM_SNP$upper)/3.92

#New dataset SNP_50_UNa_T2DM
SNP_T2DM_GUNa <- merge(x=T2DM_SNP, y=E1_UNa_GWAS, by.x = "SNP", by.y = "SNP", all.x = TRUE)


#Harmonization####
#Variable names£¬adjust beta, effective allele£¬change the valuence of beta according to effective allele
#Identify variants that do not share the same allele pair between datasets, and either correct this if possible or eliminate such variants
SNP_T2DM_GUNa <- SNP_T2DM_GUNa %>%
  mutate(effective_allele =allele1,
         other_allele = allele2,
         beta_T2DM = beta,
         se_T2DM = SE.x,
         beta_UNa = ifelse(ALLELE1 == effective_allele, BETA, 0-BETA),
         se_UNa = SE.y)


#Analysis####
#Input
MR_SNP_T2DM_GUNa <- mr_input(bx = SNP_T2DM_GUNa$beta_T2DM, bxse = SNP_T2DM_GUNa$se_T2DM,
                             by = SNP_T2DM_GUNa$beta_UNa, byse = SNP_T2DM_GUNa$se_UNa,
                             exposure = "T2DM", outcome = "UNa", snps = SNP_T2DM_GUNa$NSNP,
                             effect_allele = SNP_T2DM_GUNa$effective_allele,
                             other_allele = SNP_T2DM_GUNa$other_allele,
                             eaf = SNP_T2DM_GUNa$freq1)

#Median, IVW, Egger
mr_allmethods(MR_SNP_T2DM_GUNa)
mr_ivw(MR_SNP_T2DM_GUNa)
mr_median(MR_SNP_T2DM_GUNa)
mr_egger(MR_SNP_T2DM_GUNa)
