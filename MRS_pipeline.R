#####################################
##############STEP ONE###############
#####################################

#####################################
# Load phenotypes file
pheno = read.csv("GA_Combined_updated.csv")
#Need two columns, ID and ethnicity/race
pheno = pheno[,c("SampleID","ethnicity")]
colnames(pheno) = c("ID", "race")

# Load DNA methylation (beta) file
betas = get(load("DRAKNC.RData"))
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
save(res, file="betas_test.RData")

##########################################################
###Run cmr function to generate Co-methylation regions####
##########################################################
### Generate co-methylated regions ###
# install.packages("comeback_0.1.0.tar.gz", repos = NULL, type="source")
#comeback_0.1.0.tar.gz can be found at https://bitbucket.org/flopflip/comeback.
library(comeback)
source("MRS_func.R")

load("betas_test.RData")
res = t(res) #column are CpG sites, rows are for samples
cmrs <- cmr(Mdata = res, corlo = 0.3)
#cmrs is the output from cmr(), and saved all the co-methylation regions in a list
save(cmrs, file="CMR.test.RData")

#read ref cmr
#cmrs_ref <- readRDS("ref.rds")

# To generate Co-methylation regions dataframe for later analysis, 
# you can either input cmrs calculated from your own dataset,
CoMeRegion = GenCoMeRegion(cmrs = cmrs, beta = res, Overlap = F)
# or a cmrs calculated from a reference dataset,
CoMeRegion = GenCoMeRegion(beta = res, reference = cmrs_ref, Overlap = F)
# or both and then overlap the cmrs
CoMeRegion = GenCoMeRegion(cmrs = cmrs, beta = res, reference = cmrs_ref, Overlap = T)

#CoMeRegion is a matrix that assigned a unique number to each co-methylation region
save(CoMeRegion, file = "CoMeRegion.rda")

#####################################
##############STEP TWO###############
#####################################
### Calculate MRS
library(dplyr)    
library(robustHD)

# Load real phenotype and methylation data
DNAm = get(load("betas_test.RData"))
DNAm = t(DNAm)
DNAm[1:6,1:6] #check DNAm file

###Load smoking summary statistics (newborn)
SS_newborn <- read.csv("SmokingSS_Adults_CHARGE.csv")
#The summary statistics file should have at least four columns: CpGs names, beta coefficient, 
#standand errors and p-values
SS = SS_newborn[,c(1,2,3,4)]
#change the columns according to your dataset, 1st: CpGs, 2nd: coefficients, 3rd: SE and 4th: p-value
#Rename the column names as "Marker", "BETA" and "Pvalue"
colnames(SS) = c("Marker", "BETA", "SE", "Pvalue")
head(SS)
#      Marker    BETA    SE Pvalue
#1 cg05575921 -0.1 1  1.00e-2
#2 cg12803068  0.1 1  1.00e-2
#3 cg04180046  0.1 1  1.00e-2
#4 cg09935388 -0.1 1  1.00e-2
#5 cg25949550 -0.1 1  1.00e-2
#6 cg12876356 -0.1 1  1.00e-2

#Get the smallest p-value
minpvalue = min(SS$Pvalue[SS$Pvalue != 0])
minpvalue = sapply(strsplit(as.character(minpvalue), "-"), "[", 2)
###Load Co-methylation regions for newborns -> CoMeRegion
load("CoMeRegion.rda")
#Specify how p-value threshold, for example, if you want 5 * 10 ^ (-2), specify pthread to be 2
Pthred = 2:minpvalue
MRS = GenMRS(DNAm, SS, Pthred, CoMeRegion, CoMeBack = T, weightSE = F)
#if weightSE = T, weights = BETA/SE, where BETA is the effect size
#Basic information of MRS
write.csv(MRS$pvalueinfo, "MRS_pvalueinfo.csv", row.names = F)
write.csv(MRS$MRS, "MRS.csv", row.names = F)

MRS$ID
#####################################
##############STEP THREE#############
#####################################
#Read MRS file
MRS = read.csv("MRS.csv")

#Compare prediction performance of MRS and select the optimized MRS
#Load phenotype
pheno = read.csv("pheno.csv")
pheno = pheno[,c("ID","pheno")] 
head(pheno)
#pheno needs to have at least one column called "ID" that matches beta file's rownames, 
#and another column called "pheno" for phenotype (e.g. prenatal smoking)
colnames(pheno) = c("ID", "pheno")
#Merge MRS data and Phenotype data
MRS_Pheno = merge(MRS, pheno, by = "ID")

#Prediction performance of all MRS would be saved in CorResults
CorResults = matrix(NA, ncol(MRS)-1, 1)
rownames(CorResults) = colnames(MRS)[-1]
for (j in 2:ncol(MRS)){
  CorResults[(j-1),1] = cor(MRS_Pheno[,j], MRS_Pheno$pheno, use = "complete")^2
}
CorResults = data.frame(CorResults)
colnames(CorResults) = "R2"
#Maximum prediction
max(CorResults)
#P-value threshold that leads to maximum prediction
MaxPNumber = which(CorResults == max(CorResults))
MaxPNumber
rownames(CorResults)[MaxPNumber] #"P5e-23" and "P5e-22"

#Plot figures
library(ggplot2)
CorResults[,"Pvalue"] = as.numeric(sapply(strsplit(rownames(CorResults), "P5e."), "[", 2))
CorResults[is.na(CorResults$Pvalue),"Pvalue"] = -log10(as.numeric(sapply(strsplit(rownames(CorResults)[is.na(CorResults$Pvalue)], "P"), "[", 2))/5)
ggplot(aes(x = Pvalue, y = R2), data = CorResults) +
  geom_line(size = 1.3)+
  geom_point() +
  scale_y_continuous(name="Prediction accuracy") +
  scale_x_continuous(name = "-log10(P-value/5)") +
  scale_color_brewer(palette="Paired") + 
  ggtitle("") +
  theme( plot.title = element_text(family = "Times", face = "bold", size = 28, hjust = 0.5), #Adjusts the text properties of the plot title
         axis.title = element_text(family = "Times", face = "bold", size = 24), #Adjusts the text properties of the axis titles
         text = element_text(size = 20), #Adjusts the size of the x-axis items
         axis.line = element_line(colour = "black"), #Adds axis lines in black color
         panel.background = element_blank(), #Removes gray background
         axis.ticks = element_blank(), #Removes axis tick marks
         legend.title = element_blank(),
         legend.position = "bottom"
  )
ggsave("CorResults.png", plot = last_plot(), height = 6, width = 6) 
