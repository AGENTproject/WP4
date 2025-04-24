# Event: Datathon WP4 AGENT Project
# License: GPLv3

# load required libraries
library(data.table)
library(bigstep)
library(tidyr)
library(dplyr)

# remotes::install_github("jiabowang/GAPIT")
library(GAPIT)
# source('http://zzlab.net/GAPIT/gapit_functions.txt')

# https://vsni.co.uk/free-software/asrgenomics
library(ASRgenomics)

# remove everything in the working environment
rm(list = ls())
gc(reset = TRUE)

# set working directory
setwd('~/AGENT')

# load input data imported via the Shiny interface (deprecated)
# load('/home/rstudio/ShinyApps/agent/AGENT_T4_JOB_20230523224926.RData')

# e.g., 4 for fixed value, or NA to use Gao et al. 2008 estimation
gwas_lod <- 5.25

run_stepwise  <- TRUE
QTLome_submit <- TRUE

crop <- 'Wheat' # Wheat or Barley (case sensitive!)

traits <- c("DTH", "PLH", "TKW")
#traits <- c("HT","FT","PH","TGW","leaf_rust","yellow_rust","Yield")

outputs.folder <- paste0('WP4_Outputs/', crop)
if (!file.exists(outputs.folder)) dir.create(outputs.folder)

sig.markers.folder <- paste0('./', outputs.folder, '/sig_markers')
if (!file.exists(sig.markers.folder)) dir.create(sig.markers.folder)

geno.file  <- paste0('./AGENT_', crop, '_VCF_Assay/AGENT_', crop, '.csv.gz')
map.file   <- paste0('./AGENT_', crop, '_VCF_Assay/AGENT_', crop, '_map.csv')
blues.file <- paste0('./WP3_BLUEs_Inputs/Final_', crop, '_with_Biosample_ID.csv')

qtlome.file   <- paste0('./WP4_QTLome/QTLome_', crop, '_1000bp_range.csv')
qtlome.margin <- 1000

acc_id <- 'Biosample'  # column represents the AGENT ids
env_id <- 'Institute'  # column represents the institute FAO code

# read phenotypic input data
BLUEs  <- read.csv(blues.file)

# read genotypic input data
geno.map  <- read.csv(file = map.file)
geno.data <- as.data.frame(data.table::fread(file = geno.file))

# reshape required hapmap format
rownames(geno.data) <- geno.data[,1]
geno.data <- geno.data[,-1]
geno.data <- t(geno.data)

# rebuild hapmap metadata columns
SNPs <- geno.map[,1:4]

SNPs$strand <- '+'
SNPs[,6:11] <- NA

colnames(SNPs)  <- c('rs#', 'alleles', 'chrom', 'pos', 'strand', 'assembly#', 
                     'center', 'protLSID', 'assayLSID', 'panelLSID', 'QCcode')

SNPs <- cbind(SNPs, geno.data)

# get the SNPs metadata (rs#, alleles, chrom, pos)
SNPs_metadata <- SNPs[, 1:4]

# free memory in the working space
rm(geno.map, geno.data)

# impute genotypic data (mean imputation)
qc <- qc.filtering(M = t(SNPs[, -c(1:11)]), impute = TRUE, message = TRUE)

SNPs.clean <- cbind(SNPs[, 1:11], t(qc$M.clean))

# free memory in the working space
rm(SNPs, qc)

################################################################################
# define the significantly associated markers and the explained variance
# significance of marker-trait association is defined based on the simpleM
# method as described by Gao et al. 2008
#
# Reference:
# Gao, X.; Starmer, J.; Martin, E.R. A multiple testing correction method for
# genetic association studies using correlated single nucleotide polymorphisms.
# Genet. Epidemiol. 2008, 32, 361–369.
# https://onlinelibrary.wiley.com/doi/10.1002/gepi.20310

if (is.na(gwas_lod)) {
  # reshape SNPs data to fit in the code snippet shared by Marcel from IPK
  dta.G.GWAS <- SNPs.clean[, -c(2, 5:11)]
  
  # the matrix has columns Marker,CHR, Position and one for every SNP (0,1,2)
  colnames(dta.G.GWAS) <- c('Marker', 'CHR', 'Position', colnames(dta.G.GWAS)[-c(1:3)])
  
  # significance level for GWAS
  q <- 0.05
  
  rownames(dta.G.GWAS) <- dta.G.GWAS$Marker
  
  chromosomes <- unique(dta.G.GWAS$CHR)
  M.eff.all   <- c()
  
  pb <- txtProgressBar(min = 0, max = length(chromosomes), initial = 0, style = 3)
  
  for (i in  1:length(chromosomes)) {
    marker.interest  <- dta.G.GWAS[which(dta.G.GWAS$CHR == chromosomes[i]), 1]
    
    if (length(marker.interest) == 1) {
      M.eff <- 1
    } else {
      dta.G.chromosome <- t(dta.G.GWAS[which(dta.G.GWAS$Marker %in% marker.interest), -c(1:3)])
      
      CLD <- matrix(nrow = (ncol(dta.G.chromosome)-1), ncol = (ncol(dta.G.chromosome)-1))
      CLD <- cor(dta.G.chromosome)
      
      eig.CLD     <- eigen(CLD)
      eig.CLD     <- eig.CLD$values
      sum.eig.CLD <- sum(eig.CLD)
      
      for (l in 1:length(eig.CLD)) {
        C <- (sum(eig.CLD[1:l]) / sum.eig.CLD)
        if (C < 0.995) { M.eff <- l }
      }
      
      M.eff <- M.eff + 1
    }
    
    M.eff.all <- rbind(M.eff.all, M.eff)
    
    rm(M.eff)
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  rm(dta.G.GWAS)
  
  p <- q / sum(M.eff.all)
  
  gwas_lod <- -log10(p) 
}

start.time <- Sys.time()

for (env in unique(BLUEs[, env_id])) {
  for (trait in traits) {

    # env   <- 'IPK'
    # trait <- 'DTH'
    
    # filter to select the required records and columns only
    BLUE <- BLUEs[BLUEs[, env_id] == env, c(acc_id, trait)]
    
    # omit missing phenotypic records
    BLUE <- na.omit(BLUE)
    
    # rename accession id as Taxa to fit with the GAPIT format
    colnames(BLUE) <- c('Taxa', trait)

    ### GWAS (GAPIT) ###############################################################
    #
    # https://www.zzlab.net/GAPIT/
    # 
    # In GAPIT3, the Bonferroni multiple test threshold is used by default to 
    # determine significance, it only affect the visualization of results (Manhattan)
    # but not change the association results, users can adjust the output p-values 
    # by using the base function p.adjust() in R
    # 
    # By default, GAPIT3 conducts GWAS using the BLINK method, which has the 
    # highest statistical power & computation efficiency among all methods 
    # implemented (e.g., MLM)
    # 
    # All methods (GWAS and GS) have linear computing time to number of markers.
    # However, they have mixed computing complexity to the number of individuals.
    # Most of these methods have computing time complexity that are cubic to number
    # of individuals.
    
    # make sure that all accessions produced to gwas are exists 
    # in both geno and pheno datasets (exists in pheno even if the trait value is NA)
    geno.gwas  <- SNPs.clean[, c(1:11, which(colnames(SNPs.clean) %in% BLUE[,1]))]
    pheno.gwas <- BLUE[BLUE[,1] %in% colnames(geno.gwas),]
    
    # skip if no enough data for GWAS analysis
    if (nrow(pheno.gwas) < 10) next
    
    geno.GD <- as.data.frame(t(geno.gwas[, -c(1:11)]))
    geno.GD <- cbind(colnames(geno.gwas)[-c(1:11)], geno.GD)
    colnames(geno.GD) <- c('taxa', geno.gwas[,1])
    rownames(geno.GD) <- NULL
    
    geno.GM <- geno.gwas[, c(1,3,4)]
    colnames(geno.GM) <- c('SNP', 'Chromosome', 'Position')
    
    # https://rdrr.io/github/jiabowang/GAPIT3/man/GAPIT.html
    myGAPIT <- GAPIT(
      Y  = pheno.gwas,
      GD = geno.GD,
      GM = geno.GM,
      model = 'MLM',
      PCA.total = 3,
      SNP.MAF = 0,
      output.numerical = FALSE,
      file.output = TRUE,
      SNP.P3D = TRUE
    )
    
    # move GWAS output files into associated folder
    folder.name <- paste0(outputs.folder, '/', env, '_', trait, '_GWAS')
    if (!file.exists(folder.name)) dir.create(folder.name)
    
    gwas_outputs <- list.files(pattern = 'GAPIT\\..*')
    file.rename(gwas_outputs, paste0(folder.name, '/', gwas_outputs))
    
    # get GWAS results/output
    GWAS <- myGAPIT$GWAS[,c(1:4,9)]
    colnames(GWAS) <- c('Marker', 'Chr', 'Pos', 'p', 'effect') 
    
    # calculate the LOD score (logarithm of the odds) = -log10(pValue)
    GWAS$LOD <- -1 * log10(GWAS$p)
    
    # select SNPs markers exceed the given LOD threshold
    sel_markers <- GWAS[GWAS$LOD > gwas_lod & !is.na(GWAS$p), 'Marker']
    
    sel_SNPs <- SNPs.clean[which(SNPs.clean[,1] %in% sel_markers), ]
    sel_SNPs <- as.data.frame(t(sel_SNPs[,-c(1:11)]))
    
    # merge genetic matrix & phenotypic data using accession id as a matching key
    data <- merge(pheno.gwas, sel_SNPs, by.x = 'Taxa', by.y = 0)
    
    # remove accessions with missing trait value
    data <- data[complete.cases(data),]
    
    ### QTLome Reporting ###########################################################
    # check if the candidate markers are already listed in QTLome or not yet
    
    if (QTLome_submit && length(sel_markers) > 0) {
      # load QTLome file
      QTLome <- read.csv(qtlome.file)
      
      # data.frame forQTLome reporting
      QTLome_report <- data.frame(sel_marker_id = numeric(),
                                  sel_marker    = character(),
                                  sel_chrom     = character(),
                                  sel_position  = numeric(),
                                  sel_trait     = character(),
                                  sel_env       = character(),
                                  qtlome_status = character(),
                                  qtlome_marker = character(),
                                  position      = numeric(),
                                  qtlome_ref    = character(),
                                  qtlome_trait  = character())
      
      # loop over selected markers 
      for (i in 1:length(sel_markers)) {
        sel_marker <- strsplit(sel_markers[i], '\\:')[[1]]
        sel_chrom  <- substring(sel_marker[1], 4)
        sel_pos    <- as.numeric(sel_marker[2])
        
        sel_qtlome <- QTLome[QTLome$Position >= sel_pos - qtlome.margin & 
                               QTLome$Position <= sel_pos + qtlome.margin &
                               QTLome$Chromosome == sel_chrom,]
        
        if (nrow(sel_qtlome) > 0) {
          for (j in 1:nrow(sel_qtlome)) {
            QTLome_report[nrow(QTLome_report) + 1,] <- c(i,
                                                         sel_markers[i],
                                                         sel_chrom,
                                                         sel_pos,
                                                         trait,
                                                         env,
                                                         '-',
                                                         sel_qtlome[j, 'Name'],
                                                         sel_qtlome[j, 'Position'],
                                                         sel_qtlome[j, 'Marker_on_QTLome'],
                                                         sel_qtlome[j, 'trait_on_QTLome'])
          }
        } else {
          QTLome_report[nrow(QTLome_report) + 1,] <- c(i, sel_markers[i], sel_chrom, sel_pos, trait, env, 'NEW', '', '', '', '')
        }
      }
      
      write.csv(QTLome_report, paste0(sig.markers.folder, '/', env, '_', trait, '_qtlome_report.csv'))
    }
    
    
    ### Stepwise Regression ########################################################
    
    if (run_stepwise & length(sel_markers) > 1) {
      # Stepwise Selection for Large Data Sets
      # Selecting linear and generalized linear models for large data sets using 
      # modified stepwise procedure and modern selection criteria (like modifications 
      # of Bayesian Information Criterion). Selection can be performed on data which 
      # exceed RAM capacity.
      
      # the stepwise procedure for big data
      # first, we have to convert our data to a proper format, an object of class big
      reg_data <- bigstep::prepare_data(y = data[,trait], X = data[,-c(1:2)])
      
      # then, we can use the stepwise procedure with, for example, aic or bic
      results <- bigstep::stepwise(reg_data, crit = bic)
      
      # get the list of the markers in the final model
      sel_markers <- results$model
    }
    
    # sub-genotypic matrix includes selected markers in the final stepwise model
    data <- data[!duplicated(data$Taxa),]
    rownames(data) <- data[,1]
    data <- subset(data, select = sel_markers)
    
    # transpose sub-genotypic matrix back (SNPs in rows and accessions in columns)
    data <- as.data.frame(t(data))
    
    # merge SNPs metadata for selected markers
    data <- merge(SNPs_metadata, data, by.x = 1, by.y = 0)
    
    # save the output of sub-genotypic matrix of the selected markers
    write.csv(data, paste0(sig.markers.folder, '/', env, '_', trait, '_sub-genotyping_matrix.csv'))
    
    
    ### Minimum Subset #############################################################
    
    if (length(sel_markers) > 0) {
      # minor (alternative) allele notation
      alt_code <- 0
      
      # calculate the minor allele frequency (MAF) and observation (MAO) per Acc.
      MAF_acc <- apply(data[,-c(1:4)], MARGIN = 2, function(x) sum(x == alt_code)/length(x))
      MAO_acc <- apply(data[,-c(1:4)], MARGIN = 2, function(x) sum(x == alt_code))
      
      # reshape calculated stats in the data.frame format
      snps_stats <- as.data.frame(cbind(MAF_acc, MAO_acc))
      
      # save the stats (i.e., MAF, MAO) for the next pipeline step
      write.csv(snps_stats, paste0(sig.markers.folder, '/', env, '_', trait, '_snps_stats_for_each_landrace.csv'))
      
      # find the minimum subset of accessions that covers all candidate markers
      SNP_effect <- GWAS[GWAS$Marker %in% data$`rs#`, ]
      SNP_effect <- merge(data[, 1:4], SNP_effect[, -c(2:3)], by = 1)
      
      minimum_group <- data.frame(Accession = character(), 
                                  Cumulative_MAO = integer(), 
                                  MAF = numeric())
      
      # calculate the minor allele frequency (MAF) for each SNP
      MAF_snp <- apply(data[,-c(1:4)], MARGIN = 1, function(x) sum(x == alt_code)/length(x))
      
      # calculate weighted SNPs scores to increase the chance of including accessions 
      # represents less frequent SNPs, account for SNP effect and magnitude 
      rest_wSNPs <- cbind(data[,1:4], data[,-c(1:4)] * SNP_effect[, 'effect'] * sqrt((1 / MAF_snp)))
      
      # if SNP has negative effect, then the target SNP value is 0 (reference allele)
      rest_wSNPs  <- rest_wSNPs[SNP_effect[, 'effect'] > 0, ]
      accumulated <- sum(SNP_effect[, 'effect'] < 0)
      
      #' loop while we still have minor allele not represented 
      #' in the minimum group of accessions
      while (nrow(rest_wSNPs) > 0) {
        #' calculate the scores of remaining markers not represented yet 
        #' in the minimum group of accessions 
        rest_score <- apply(rest_wSNPs[,-c(1:4)], MARGIN = 2, sum)
        
        #' get the accession name that has the max score, if there are more than one, 
        #' on tie, select the one with the highest MAF all over the whole set of 
        #' selected markers
        a <- names(which.max(MAF_acc[names(rest_score[rest_score == max(rest_score)])]))
        
        #' calculate the accumulated unique minor alleles represents in the current 
        #' minimum group of accessions
        accumulated <- accumulated + sum(data[data[,1] %in% rest_wSNPs[,1], a] == 2)
        
        #' add the selected accession name and the accumulated number of minor alleles 
        #' represented to the minimum group of accessions
        minimum_group <- rbind(minimum_group, 
                               data.frame(Accession = a, 
                                          Accumulated_MAO = accumulated, 
                                          MAF_acc = MAF_acc[a]))
        
        #' keep all SNPs markers not represented in the last accession added to the 
        #' minimum group
        rest_wSNPs <- rest_wSNPs[rest_wSNPs[,a] == 0, ]
      }
      
      write.csv(minimum_group, paste0(sig.markers.folder, '/', env, '_', trait, '_minimum_group_accessions.csv'), row.names = FALSE)
      
      #' subset the genotyping matrix to have only the selected markers and the optimum 
      #' group of accessions
      minimum_group_snps <- data[, c(colnames(data)[1:4], minimum_group$Accession)]
      
      effect <- SNP_effect[, 'effect']
      occurs <- apply(minimum_group_snps, MARGIN = 1, function(x) sum(x == 0))
      occurs <- ifelse(effect > 0, length(minimum_group$Accession) - occurs, occurs)
      
      minimum_group_snps <- cbind(occurs, effect, minimum_group_snps)
      minimum_group_snps <- minimum_group_snps[, c(3, 1:2, 4:ncol(minimum_group_snps))]
      
      write.csv(minimum_group_snps, paste0(sig.markers.folder, '/', env, '_', trait, '_minimum_group_snps.csv'), row.names = FALSE)
    }
  }
}

end.time <- Sys.time()
(time.taken <- round(end.time - start.time, 2))

### Summary Reporting ##########################################################

outputs.folder <- paste0('./WP4_Outputs/', crop, '/')

# get summary with all markers and their respective LOD for all genebanks
# to be used to create a general Manhattan plot
for (trait in traits) {
  report <- data.frame()
  for (env in unique(BLUEs[, env_id])) {
    gwas.file <- paste0(outputs.folder, env, '_', trait, '_GWAS/GAPIT.Association.GWAS_Results.MLM.', trait, '.csv')
    if (!file.exists(gwas.file)) next
    temp <- read.csv(gwas.file)
    if (nrow(report) == 0) { report <- temp[,1:3] }
    temp$LOD <- -log10(temp$P.value)
    temp <- temp[, c('LOD', 'Effect', 'MAF')]
    colnames(temp) <- paste(env, trait, colnames(temp), sep = '.')
    report <- cbind(report, temp)
  }

  marker_map <- report[, c(1:3)]
  report <- round(report[, -c(1:3)], 3)
  lod_cols <- seq(1, ncol(report), by = 3)
  report <- cbind(marker_map, sig_marker = apply(report[, lod_cols], 1, function(row) ifelse(any(row > gwas_lod, na.rm = TRUE),1,0)), report)
  write.csv(report, paste0('./WP4_Outputs/', crop, '/',crop, '_',trait,'_GWAS_Markers_Summary.csv'))
}

################################################################################

markers <- as.data.frame(data.table::fread(file = geno.file))
rownames(markers) <- markers[, 'BIOSAMPLE_ID']

markers <- t(markers[rownames(markers) %in% BLUEs$Biosample, -1])

df <- data.frame()

for (env in unique(BLUEs[, env_id])) {
  for (trait in traits) {
    # env <- 'IPGR'
    # trait <- 'DTH'

    file_path <- paste0('./WP4_Outputs/', crop, '/sig_markers/', env, '_', trait, '_minimum_group_snps.csv')
    if (!file.exists(file_path)) next

    qtl <- read.csv(file_path)
    qtl <- qtl[, c(1,4:6)]

    qtl$env   <- env
    qtl$trait <- trait

    file_path <- paste0('./WP4_Outputs/', crop, '/', crop, '_', trait, '_GWAS_Markers_Summary.csv')
    if (!file.exists(file_path)) next

    gwas <- read.csv(file_path, row.names = 1)

    params <- paste0(env, '.', trait, '.', c('LOD', 'Effect', 'MAF'))
    if (!all(params %in% colnames(gwas))) next

    gwas <- gwas[, c('SNP', params)]
    colnames(gwas) <- c('rs.', 'LOD', 'Effect', 'MAF')

    qtl <- merge(qtl, gwas, by = 'rs.')

    if (nrow(qtl) == 1) {
      temp <- cbind(qtl, t(markers[rownames(markers) == qtl$rs.,]))
    } else {
      temp <- merge(qtl, markers[rownames(markers) %in% qtl$rs.,], by.x = 'rs.', by.y = 0)
    }

    df <- rbind(df, temp)
  }
}

write.csv(df, paste0('./WP4_Outputs/', crop, '/', crop, '_sig_markers_subset.csv'))

accession_cols <- colnames(df)[-c(1:9)]

summary_df <- df %>%
  select(trait, Effect, all_of(accession_cols)) %>%
  pivot_longer(cols = all_of(accession_cols), names_to = "Accession", values_to = "Value") %>%
  filter(Value == 2) %>%
  group_by(Accession, trait) %>%
  summarise(Sum_Effect = sum(Effect, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = trait, values_from = Sum_Effect, values_fill = 0)

write.csv(summary_df, paste0('./WP4_Outputs/', crop, '/', crop, '_QTL_effects_per_accession.csv'))

