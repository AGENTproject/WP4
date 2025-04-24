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

gs.output.folder <- paste0('WP4_Outputs/', crop, '/GS_Output')
if (!file.exists(gs.output.folder)) dir.create(gs.output.folder)

# imputed genotypic markers in hapmap file format
marker_matrix <- paste0('./AGENT_', crop, '_VCF_Assay/AGENT_', crop, '.csv.gz')
marker_map    <- paste0('./AGENT_', crop, '_VCF_Assay/AGENT_', crop, '_map.csv')

pheno_file <- paste0('./WP3_BLUEs_Inputs/Final_', crop, '_with_Biosample_ID.csv')

acc_id <- 'Biosample'  # column represents the AGENT ids
env_id <- 'Institute'  # column represents the institute FAO code

traits <- c("DTH", "PLH", "TKW")
# traits <- c("HT", "FT", "PH", "TGW", "leaf_rust", "yellow_rust", "Yield")

# 2D matrix, genotype id in row names, and climate variable in col names
climate_file <- paste0('./WP3_BLUEs_Inputs/', crop, '_climate_matrix.csv')

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

SNPs.clean  <- as.data.frame(data.table::fread(marker_matrix))
SNPs.map    <- read.csv(marker_map)
BLUE        <- read.csv(pheno_file)
acc.climate <- read.csv(climate_file)

acc.climate <- acc.climate[!duplicated(acc.climate$Biosample),]
rownames(acc.climate) <- acc.climate$Biosample
acc.climate <- acc.climate[, -c(1:4)]
acc.climate <- acc.climate[complete.cases(acc.climate),]

### Reshape and impute markers data ############################################

metadata <- SNPs.map[, 1:4]
metadata$ref <- substr(metadata$alleles, 1, 1)
metadata$minor <- substr(metadata$alleles, 3, 3)
metadata$alleles <- NULL
colnames(metadata) <- c('Marker', 'CHR', 'LOC', 'REF', 'ALT')

rownames(SNPs.clean) <- SNPs.clean[, 1]
SNPs.clean <- SNPs.clean[, -1]

qc <- ASRgenomics::qc.filtering(M = as.matrix(SNPs.clean), impute = TRUE, message = TRUE)

SNPs.clean <- as.data.frame(t(qc$M.clean))

SNPs.clean$Marker <- rownames(SNPs.clean)
geno.data <- merge(metadata, SNPs.clean, by = "Marker", all.y = TRUE)

rm(SNPs.clean, qc)
gc()


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

#' save intermediate results
save.image(file = paste0('gs_', crop, '_inputs.RData'))

### loop over env and traits ###################################################

BLUE <- read.csv(paste0('./WP3_BLUEs_Inputs/Final_', crop, '_with_Biosample_ID.csv'))

for (env in unique(BLUE[, env_id])) {
  for (trait in traits) {
    
    # env   <- 'ICARDA' # unique(BLUE[, env_id])
    # trait <- 'TKW'    # DTH, PLH, TKW
    
    #' load intermediate results
    load(paste0('gs_', crop, '_inputs.RData'))
    
    all.geno.data <- geno.data
    
    gwas.file <- paste0('./WP4_Outputs/', crop, '/', env, '_', trait, '_GWAS/GAPIT.Association.Filter_GWAS_results.csv')
    
    if(!file.exists(gwas.file)) {
      next
    } else {
      gwas_results <- read.csv(gwas.file, row.names = 1)
    }
    
    #' list of significant markers (e.g., GWAS or/and QTLome output)
    sel_markers <- gwas_results$SNP
    
    #' All entries with a given genotypic data:
    pheno.data <- merge(colnames(geno.data)[-c(1:5)], BLUE[BLUE$Institute == env, c('Biosample', trait)],
                        by.x = 1, by.y = 1, all.x = TRUE, sort = FALSE)
    
    colnames(pheno.data) <- c('Genotype_Matrix', 'Trait')
    
    #' restrict the pheno, geno, and groups to the accessions have coordinates/climatic data
    geno.data  <- geno.data[, c(1:5, which(colnames(geno.data) %in% colnames(penv1)))]
    pheno.data <- pheno.data[pheno.data[,1] %in% colnames(penv1),]
    acc.groups <- acc.groups[acc.groups[,1] %in% colnames(penv1),]
    
    acc.climate <- penv1
    
    if (sum(!is.na(pheno.data$Trait)) == 0) {
      next
    }
    
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
    if (length(sel_markers) > 0) {
      model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                                     groups = NULL, sig_markers = sel_markers, climate = NULL,
                                     num.folds = test.num.folds, nIter = test.nIter, 
                                     burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)
      
      models_accuracy <- rbind(models_accuracy, model_accuracy)
    }
    
    
    ### Model 3 (Climate matrix) ###################################################
    model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                                   groups = NULL, sig_markers = NULL, climate = acc.climate,
                                   num.folds = test.num.folds, nIter = test.nIter, 
                                   burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)
    
    models_accuracy <- rbind(models_accuracy, model_accuracy)
    
    
    ### Model 4 (Kinship split + Fixing major markers) #############################
    if (length(sel_markers) > 0) {
      model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                                     groups = acc.groups, sig_markers = sel_markers, climate = NULL,
                                     num.folds = test.num.folds, nIter = test.nIter, 
                                     burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)
      
      models_accuracy <- rbind(models_accuracy, model_accuracy)
    }
    
    
    ### Model 5 (Kinship split + Climate matrix) ###################################
    model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                                   groups = acc.groups, sig_markers = NULL, climate = acc.climate,
                                   num.folds = test.num.folds, nIter = test.nIter, 
                                   burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)
    
    models_accuracy <- rbind(models_accuracy, model_accuracy)
    
    
    ### Model 6 (Fixing major markers + Climate matrix) ############################
    if (length(sel_markers) > 0) {
      model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                                     groups = NULL, sig_markers = sel_markers, climate = acc.climate,
                                     num.folds = test.num.folds, nIter = test.nIter, 
                                     burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)
      
      models_accuracy <- rbind(models_accuracy, model_accuracy)
    }
    
    
    ### Model 7 (Kinship split + Fixing major markers + Climate matrix) ############
    if (length(sel_markers) > 0) {
      model_accuracy <- agend_egblup(pheno.data, geno.data, MAF = test.MAF, pred_mode = FALSE, 
                                     groups = acc.groups, sig_markers = sel_markers, climate = acc.climate,
                                     num.folds = test.num.folds, nIter = test.nIter, 
                                     burnIn = test.burnIn, num.crossvalid = test.num.crossvalid)
      
      models_accuracy <- rbind(models_accuracy, model_accuracy)
    }

    ### Select Best Model ##########################################################
    
    output_file <- paste0(gs.output.folder, '/Best_Model_', crop, '_', env, '_', trait, '.txt')
    write('Select Best Model:', file = output_file)
    
    capture.output(models_accuracy, file = output_file, append = TRUE)
    
    fit <- lm(cor ~ model, data = models_accuracy)
    capture.output(anova(fit), file = output_file, append = TRUE)
    
    test <- agricolae::LSD.test(fit, 'model')
    capture.output(test$groups, file = output_file, append = TRUE)
    
    best_model <- rownames(test$groups[1,])
    write(paste("\nBest Model is:", best_model), file = output_file, append = TRUE)
    
    other_model <- rownames(test$groups)[which(!grepl("Climatic", rownames(test$groups)))[1]]
    write(paste("\nThe other Model is:", other_model), file = output_file, append = TRUE)
    
    end.time <- Sys.time()
    capture.output(round(end.time - start.time, 2), file = output_file, append = TRUE)
    
    write(paste("\nTotal number of BLUE values:", sum(!is.na(pheno.data$Trait))), file = output_file, append = TRUE)
    
    
    ### Get Best Model Predictions #################################################
    
    start.time <- Sys.time()
    
    geno.data <- all.geno.data
    all.acc <- colnames(geno.data)[-c(1:5)]
    all.pheno.data <- data.frame(Genotype_Matrix = all.acc)
    all.pheno.data <- merge(all.pheno.data, pheno.data, by = 1, all.x = TRUE, sort = FALSE)
    
    
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
    
    end.time <- Sys.time()
    round(end.time - start.time, 2)
    
    write.csv(predictions, file = paste0(gs.output.folder, '/Predictions_', crop, '_', trait, '_', env, '.csv'))
  }
}