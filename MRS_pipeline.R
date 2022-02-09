#####################################
##############STEP ONE###############
#####################################
### Generate co-methylated regions ###
# install.packages("/projects/huels_lab/MRS/02_code/03-CTMethod/01-CoMeBack/comeback_0.1.0.tar.gz", 
# repos = NULL, type="source")
library(comeback)
source(MRS_func.R)

# Load phenotypes file
pheno = read.csv("/projects/huels_lab/MRS/03_results/01-Gestational_age/GA_Combined_updated.csv")
#Need two columns, ID and ethnicity/race
pheno = pheno[,c("SampleID","ethnicity")]
colnames(pheno) = c("ID", "race")

# Load DNA methylation (beta) file
betas = get(load("/projects/huels_lab/DCHS/PACE/01_data/01_DNAm/04_QCed_data_combined/DRAKNC.RData"))
betas = slot(betas,"assayData")[["betas"]]
betas = data.frame(t(betas))
betas$ID = rownames(betas)  

# Correct for ancestry first (regress out ancestry) 
# Could skip this step and just use betas to calculate co-methylation regions
betas_final = merge(betas, pheno, by = "ID")
res <- matrix(NA, ncol = ncol(betas) - 1, nrow = nrow(betas_final))
rownames(res) <- betas_final$ID
colnames(res) <- colnames(betas_final)[2:ncol(betas)]

for(i in 1:(ncol(betas)-1)){
  cpg = colnames(betas_final)[(i+1)]
  dat2 = betas_final[,c(cpg, "race")]
  colnames(dat2) = c("cpg", "race")
  lm1 <- lm(cpg ~ race, data = dat2)
  res[,i] = resid(lm1)
}
save(res, file="/projects/huels_lab/MRS/03_results/03-CTMethod/00-CMR/02-Drakenstein_Newborn/Newborn_betas_controlled_for_ancestry.RData")

#Run cmr function to generate Co-methylation regions
set.seed(1234)
cmrs <- cmr(Mdata = res, corlo = 0.3)
#cmrs is the output from cmr(), and saved all the co-methylation regions in a list
save(cmrs, file="/projects/huels_lab/MRS/03_results/03-CTMethod/00-CMR/02-Drakenstein_Newborn/Newborn_CMRs_controlled_for_ancestry.RData")

#read ref cmr
cmrs_ref <- readRDS("/projects/huels_lab/MRS/03_results/03-CTMethod/00-CMR/03-Reference/cmr_B_noCTC_Scor30_crldst400_mxd1K_mxld2_medNoCTC_disFree_age25_80.rds")

# To generate Co-methylation regions dataframe for later analysis, 
# you can either input cmrs calculated from your own dataset,
CoMeRegion = GenCoMeRegion(cmrs = cmrs, beta = res, Overlap = F)
# or a cmrs calculated from a reference dataset,
CoMeRegion = GenCoMeRegion(beta = res, reference = cmrs_ref, Overlap = F)
# or both and then overlap the cmrs
CoMeRegion = GenCoMeRegion(cmrs = cmrs, beta = res, reference = cmrs_ref, Overlap = T)

#CoMeRegion is a matrix that assigned a unique number to each co-methylation region
save(CoMeRegion, file = "/projects/huels_lab/MRS/03_results/03-CTMethod/00-CMR/01-Adults/CoMeRegion_Newborn_All.rda")

#####################################
##############STEP TWO###############
#####################################
### Calculate MRS
library(dplyr)    
library(robustHD)

# Load real phenotype and methylation data
Betas_all = get(load("/projects/huels_lab/DCHS/PACE/01_data/01_DNAm/04_QCed_data_combined/DRAKNC.RData"))
Betas_all = slot(Betas_all,"assayData")[["betas"]]
Betas_all = data.frame(Betas_all)
beta = data.frame(t(Betas_all))
beta[1:6,1:6] #check beta file
#        cg00000029 cg00000109 cg00000165 cg00000236 cg00000289 cg00000292
#UUU_11  0.11111111 0.11111111 0.11111111 0.11111111 0.11111111 0.11111111 
#UUU_21  0.11111111 0.11111111 0.11111111 0.11111111 0.11111111 0.11111111
#UUU_31  0.11111111 0.11111111 0.11111111 0.11111111 0.11111111 0.11111111
#UUU_41  0.11111111 0.11111111 0.11111111 0.11111111 0.11111111 0.11111111
#UUU_51  0.11111111 0.11111111 0.11111111 0.11111111 0.11111111 0.11111111
#UUU_61  0.11111111 0.11111111 0.11111111 0.11111111 0.11111111 0.11111111

#Load phenotype
pheno = read.csv("/projects/huels_lab/MRS/03_results/01-Gestational_age/GA_Combined_updated.csv")
#pheno needs to have a column called "ID" that matches beta file's rownames
colnames(pheno)[1] = "ID"

###Load smoking summary statistics (newborn)
SS_newborn <- read.csv("/projects/huels_lab/MRS/01_data/03_smokingSS/SmokingSS_newborn_pace.csv")
#Your final summary statistics dataframe should have at least three columns: CpGs names, coefficient and p-value
SS = SS_newborn[,c(1,2,4)] #subset the first, second and fourth columns of SS_newborn
#change the columns according to your dataset, 1st: CpGs, 2nd: coefficients, 3rd: p-value
#Rename the column names as "Marker", "BETA" and "Pvalue"
colnames(SS) = c("Marker", "BETA", "Pvalue")
head(SS)
#      Marker    BETA    Pvalue
#1 cg05575921 -0.1 1.00e-2
#2 cg12803068  0.1 1.00e-2
#3 cg04180046  0.1 1.00e-2
#4 cg09935388 -0.1 1.00e-2
#5 cg25949550 -0.1 1.00e-2
#6 cg12876356 -0.1 1.00e-2

#Get the smallest p-value
minpvalue = min(SS$Pvalue)
minpvalue = sapply(strsplit(as.character(minpvalue), "-"), "[", 2)
###Load Co-methylation regions for newborns
load("/projects/huels_lab/MRS/03_results/03-CTMethod/00-CMR/02-Drakenstein_Newborn/CoMeRegion_Newborn_All.rda")
#Specify how p-value threshold, for example, if you want 5 * 10 ^ (-2), specify pthread to be 2
Pthred = 2:minpvalue
MRS = GenMRS(beta, SS_newborn, Pthred, pheno, CoMeRegion, CoMeBack = T)
MRS = MRS$MRS
write.csv(MRS, "/projects/huels_lab/MRS/03_results/03-CTMethod/02-Smoking_MRS/MRS_smoking_newborn_PT.csv")
