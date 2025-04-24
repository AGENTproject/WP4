library(ASRgenomics)
library(data.table)
library(rlist) # used version: rlist_0.4.6.1
library(BGLR)  # used version: BGLR_1.0.8
library(caret) # used version: caret_6.0-86
library(parallel)
library(tryCatchLog)
library(adegenet)
library(R.utils)

# remove everything in the working environment
rm(list = ls())
gc(reset = TRUE)

# set working directory
setwd('~/AGENT')

Sys.setenv(OPENBLAS_NUM_THREADS="1")
# Sys.setenv(R_MAX_VSIZE = "100G")

# load the EG-BLUP GS modeling/prediction generic function
source('00_EG-BLUP_Function.R')

### Script Input Configuration #################################################

crop  <- 'Wheat' # Wheat or Barley (case sensitive!)
trait <- 'DTH'
env   <- 'CREA'

acc_id <- 'Biosample'  # column represents the AGENT ids
env_id <- 'Institute'  # column represents the institute FAO code

geno.file  <- paste0('./AGENT_', crop, '_VCF_Assay/AGENT_', crop, '.csv.gz')
map.file   <- paste0('./AGENT_', crop, '_VCF_Assay/AGENT_', crop, '_map.csv')
blues.file <- paste0('./WP3_BLUEs_Inputs/Final_', crop, '_with_Biosample_ID.csv')
sig.marker <- paste0('./WP4_Outputs/', crop, '/sig_markers/', env, '_', trait, '_minimum_group_snps.csv')

#' list of significant markers (e.g., GWAS or/and QTLome output)
# sel_markers <- c('Chr6B:90649005:T:C', 'Chr3B:785438778:C:G', 'Chr4A:594007637:C:G')
sel_markers <- read.csv(sig.marker)[[1]]

#' define the number of the population structure clusters 
cluster_no <- 4

#' BGLR testing models parameters
test.num.crossvalid <- 2
test.num.folds      <- 5
test.nIter          <- 10000
test.burnIn         <- 3000
test.MAF            <- 0.01

#' BGLR best model parameters
best.nIter          <- 10000
best.burnIn         <- 3000
best.MAF            <- 0.01

### Reading the input files ####################################################

models_accuracy <- data.frame(ncol(3))

#' 2 cols data.frame (1st column is your genotype id, and 2nd column is trait)
#' NOTE: only one line per entry (i.e., genotype id)!
BLUE <- read.csv(blues.file)
BLUE <- BLUE[BLUE[, env_id] == env, c(acc_id, trait)]

metadata <- read.csv(map.file)
metadata$alleles <- NULL
colnames(metadata) <- c('Marker', 'CHR', 'LOC', 'REF', 'ALT')

geno.data <- as.data.frame(data.table::fread(file = geno.file))
rownames(geno.data) <- geno.data[,1]
geno.data <- geno.data[,-1]

# impute genotypic data (mean imputation)
qc <- qc.filtering(M = as.matrix(geno.data), impute = TRUE, message = TRUE)

geno.data <- cbind(metadata, t(qc$M.clean))
rownames(geno.data) <- NULL


#' subset phenotypic dataset to exclude entries with no genotypic data
pheno.data <- merge(colnames(geno.data)[-c(1:5)], BLUE, by.x = 1, by.y = 1, all.x = TRUE, sort = FALSE)

colnames(pheno.data) <- c('Genotype_Matrix', 'Trait')

#' 2D matrix, genotype id in row names, and climate variable in col names
acc.climate <- read.csv(paste0('./WP3_BLUEs_Inputs/', crop, '_climate_matrix.csv'), row.names = 1)


### Climate matrix #############################################################

#' scale climate indices, calculate distance, then similarity matrix between acc.
acc.climate <- scale(acc.climate, center = TRUE, scale = TRUE)
penv <- as.matrix(dist(acc.climate, method = "euclidean", diag = TRUE, upper = TRUE))
penv <- 1 - penv/max(penv, na.rm = TRUE)

#' invert a symmetric, positive definite square matrix
penv1 <- chol2inv(penv)

#' computes the nearest positive definite (if it is not)
if(!lqmm::is.positive.definite(penv1)) penv1 <- lqmm::make.positive.definite(penv1)

penv1 <- as(solve(penv1, tol = 1e-19), "dgCMatrix")

rownames(penv1) <- rownames(penv)
colnames(penv1) <- colnames(penv)

penv1 <- as.matrix(penv1)

#' garbage collecting
rm(penv)

### Kinship split ##############################################################

acc.groups <- read.csv(file = paste0('./AGENT_', crop, '_VCF_Assay/AGENT_', crop, '_kinship_groups.csv'), row.names = 1)

save.image(file = paste0('gs_', crop, '_inputs.RData'))

#' restrict the pheno, geno, and groups to the accessions have coordinates/climatic data
geno.data  <- geno.data[, c(1:5, which(colnames(geno.data) %in% colnames(penv1)))]
pheno.data <- pheno.data[pheno.data[,1] %in% colnames(penv1),]
acc.groups <- acc.groups[acc.groups[,1] %in% colnames(penv1),]

acc.climate <- penv1


### GS (Testing Models) ########################################################
start.time <- Sys.time()

### Model 0 (EGBLUP) ###########################################################
model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                               groups = NULL, sig_markers = NULL, climate = NULL,
                               num.folds = test.num.folds, nIter = test.nIter, 
                               burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)

models_accuracy <- rbind(models_accuracy, model_accuracy)


### Model 1 (Kinship split) ####################################################
model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                               groups = acc.groups, sig_markers = NULL, climate = NULL,
                               num.folds = test.num.folds, nIter = test.nIter, 
                               burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)

models_accuracy <- rbind(models_accuracy, model_accuracy)


### Model 2 (Fixing major markers) #############################################
model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                               groups = NULL, sig_markers = sel_markers, climate = NULL,
                               num.folds = test.num.folds, nIter = test.nIter, 
                               burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)

models_accuracy <- rbind(models_accuracy, model_accuracy)


### Model 3 (Climate matrix) ###################################################
model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                               groups = NULL, sig_markers = NULL, climate = acc.climate,
                               num.folds = test.num.folds, nIter = test.nIter, 
                               burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)

models_accuracy <- rbind(models_accuracy, model_accuracy)


### Model 4 (Kinship split + Fixing major markers) #############################
model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                               groups = acc.groups, sig_markers = sel_markers, climate = NULL,
                               num.folds = test.num.folds, nIter = test.nIter, 
                               burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)

models_accuracy <- rbind(models_accuracy, model_accuracy)


### Model 5 (Kinship split + Climate matrix) ###################################
model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                               groups = acc.groups, sig_markers = NULL, climate = acc.climate,
                               num.folds = test.num.folds, nIter = test.nIter, 
                               burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)

models_accuracy <- rbind(models_accuracy, model_accuracy)


### Model 6 (Fixing major markers + Climate matrix) ############################
model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                               groups = NULL, sig_markers = sel_markers, climate = acc.climate,
                               num.folds = test.num.folds, nIter = test.nIter, 
                               burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)

models_accuracy <- rbind(models_accuracy, model_accuracy)


### Model 7 (Kinship split + Fixing major markers + Climate matrix) ############
model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                               groups = acc.groups, sig_markers = sel_markers, climate = acc.climate,
                               num.folds = test.num.folds, nIter = test.nIter, 
                               burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)

models_accuracy <- rbind(models_accuracy, model_accuracy)

end.time <- Sys.time()
round(end.time - start.time, 2)

### Select Best Model ##########################################################

sink(paste0('./WP4_Outputs/', crop, '/', env, '_', trait, '_select_best_model.txt'))
models_accuracy

fit <- lm(cor ~ model, data = models_accuracy)
anova(fit)

test <- agricolae::LSD.test(fit, 'model')
test$groups

best_model <- rownames(test$groups[1,])
writeLines(paste("\nBest Model is:", best_model))
sink()


### Get Best Model Predictions #################################################

all.acc <- colnames(geno.data)[-c(1:5)]
all.pheno.data <- data.frame(Genotype_Matrix = all.acc)
all.pheno.data <- merge(all.pheno.data, pheno.data, by = 1, all.x = TRUE, sort = FALSE)


if (startsWith(best_model, "Kinship")){
  #' reshape the data.frame of allele data to match the required format for the df2genind function
  geno.GD <- as.data.frame(t(geno.data[, -c(1:5)]))
  colnames(geno.GD) <- geno.data[,1]
  
  #' convert the data.frame of allele data to a genind object (exclude non genetic data)
  genind_obj <- adegenet::df2genind(geno.GD, ploidy = 1, ind.names = rownames(geno.GD), 
                                    loc.names = colnames(geno.GD))
  
  #' cluster identification using successive K-means
  #' you may need to increase your memory limit to avoid allocate vector error because of size
  #' e.g., memory.limit(size = 32000)
  cluster_obj <- adegenet::find.clusters(genind_obj, max.n.clust = cluster_no, 
                                         n.pca = cluster_no + 1, n.clust = cluster_no)
  
  #' get the factor that giving group membership for each individual
  acc.groups <- as.data.frame(cluster_obj$grp)
  acc.groups <- cbind(rownames(acc.groups), acc.groups)
  colnames(acc.groups) <- c("AGENT_ID", "group")
  rownames(acc.groups) <- NULL
}

### Model 0 (EGBLUP) ###########################################################
if (best_model == 'EGBLUP') {
  predictions <- agend_egblup(all.pheno.data, geno.data, MAF = best.MAF, pred_mode = TRUE, 
                              groups = NULL, sig_markers = NULL, climate = NULL,
                              nIter = best.nIter, burnIn = best.burnIn)
}

### Model 1 (Kinship split) ####################################################
if (best_model == 'Kinship') {
  predictions <- agend_egblup(all.pheno.data, geno.data, MAF = best.MAF, pred_mode = TRUE, 
                              groups = acc.groups, sig_markers = NULL, climate = NULL,
                              nIter = best.nIter, burnIn = best.burnIn)
}

### Model 2 (Fixing major markers) #############################################
if (best_model == 'Fixing') {
  predictions <- agend_egblup(all.pheno.data, geno.data, MAF = best.MAF, pred_mode = TRUE, 
                              groups = NULL, sig_markers = sel_markers, climate = NULL,
                              nIter = best.nIter, burnIn = best.burnIn)
}

### Model 3 (Climate matrix) ###################################################
if (best_model == 'Climatic') {
  predictions <- agend_egblup(all.pheno.data, geno.data, MAF = best.MAF, pred_mode = TRUE, 
                              groups = NULL, sig_markers = NULL, climate = acc.climate,
                              nIter = best.nIter, burnIn = best.burnIn)
}

### Model 4 (Kinship split + Fixing major markers) #############################
if (best_model == 'Kinship+Fixing') {
  predictions <- agend_egblup(all.pheno.data, geno.data, MAF = best.MAF, pred_mode = TRUE, 
                              groups = acc.groups, sig_markers = sel_markers, climate = NULL,
                              nIter = best.nIter, burnIn = best.burnIn)
}

### Model 5 (Kinship split + Climate matrix) ###################################
if (best_model == 'Kinship+Climatic') {
  predictions <- agend_egblup(all.pheno.data, geno.data, MAF = best.MAF, pred_mode = TRUE, 
                              groups = acc.groups, sig_markers = NULL, climate = acc.climate,
                              nIter = best.nIter, burnIn = best.burnIn)
}

### Model 6 (Fixing major markers + Climate matrix) ############################
if (best_model == 'Fixing+Climatic') {
  predictions <- agend_egblup(all.pheno.data, geno.data, MAF = best.MAF, pred_mode = TRUE, 
                              groups = NULL, sig_markers = sel_markers, climate = acc.climate,
                              nIter = best.nIter, burnIn = best.burnIn)
}

### Model 7 ((Kinship split + Fixing major markers + Climate matrix) ###########
if (best_model == 'Kinship+Fixing+Climatic') {
  predictions <- agend_egblup(all.pheno.data, geno.data, MAF = best.MAF, pred_mode = TRUE, 
                              groups = acc.groups, sig_markers = sel_markers, climate = acc.climate,
                              nIter = best.nIter, burnIn = best.burnIn)
}

write.csv(predictions, file = paste0('./WP4_Outputs/', crop, '/', env, '_', trait, '_gs_predictions.csv'))
