rm(list = ls())
library(openxlsx)
library(BacArena)
source(file.path('.','R','BN-BacArena.R'))

# Define arena
arena <- Arena(m = 20, n = 20, seed = 10)

# Add organism
model <- Bac(loadMAT(file.path('.','test','models','Alistipes_putredinis_DSM_17216.mat')), limit_growth = T)
arena <- addOrg(arena, model, amount = 3, biomass = 0.5)

model <- Bac(loadMAT(file.path('.','test','models','Bacteroides_uniformis_ATCC_8492.mat')), limit_growth = T)
arena <- addOrg(arena, model, amount = 10, biomass = 0.5)

model <- Bac(loadMAT(file.path('.','test','models','Barnesiella_intestinihominis_YIT_11860.mat')), limit_growth = T)
arena <- addOrg(arena, model, amount = 10, biomass = 0.5)

# Add media
media <- read.xlsx(file.path('.','test','media','M1.xlsx'))
mets <- media$REACTION
cant <- media$mM
arena <- addSubs(arena, smax = cant, mediac = mets, unit = 'mM', difunc = 'r', addAnyway = T)

# Load coefficient matrix
bacMat <- read.xlsx(file.path('.','test','coeffMat','cellCoeff.xlsx'), colNames = F)
nutMat <- read.xlsx(file.path('.','test','coeffMat','nutCoeff.xlsx'), colNames = F)
colnames(nutMat) <- mets

# Simulate
data <- simEnvBN(arena, 3, bacCoeff = bacMat, nutCoeff = nutMat)
