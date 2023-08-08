#' Event: Datathon WP4 AGENT Project
#' License: GPLv3

#' load required libraries
library(data.table)
library(bigstep)
library(inline)

library(rlist) # used version: rlist_0.4.6.1
library(BGLR)  # used version: BGLR_1.0.8
library(caret) # used version: caret_6.0-86
library(parallel)
library(tryCatchLog)

#' https://github.com/jiabowang/GAPIT
library(GAPIT3)
# source('http://zzlab.net/GAPIT/gapit_functions.txt')

#' https://vsni.co.uk/free-software/asrgenomics
library(ASRgenomics)

#' set working directory
setwd('.')

rm(list = ls())

#' Shiny app path for the frontend interface
shiny_app <- './'

#' Poland Datathon datasets: 
job_id <- 'AGENT_T4_JOB_20230508060723'

#' load input data sent from the Shiny job requester application
load(paste0(shiny_app, job_id, '.RData'))

config$egblup_all  <- FALSE
config$calc_gao    <- TRUE
config$QTLome_file <- 'QTLome_report_template.csv'

#' get the SNPs metadata (rs#, alleles, chrom, pos)
SNPs_metadata <- SNPs[, 1:4]

#' get the pheno input data
#' NOTE: need to get weights or SE for these estimates
BLUE <- BLUE[, c(config$acc_id, config$trait)]
BLUE <- na.omit(BLUE)

#' rename accession id as Taxa to fit with the GAPIT format
colnames(BLUE) <- c('Taxa', config$trait)

################################################################################
# Pre-processing/filtering SNPs & Acc.
hapmap_snp_attr <- c('rs#', 'alleles', 'chrom', 'pos', 'strand', 'assembly#', 
                     'center', 'protLSID', 'assayLSID', 'panelLSID', 'QCcode')

geno.M <- SNPs[, -c(1:11)]
rownames(geno.M) <- SNPs[, 1]
geno.M <- t(geno.M)

qc <- qc.filtering(
  M = geno.M, 
  marker.callrate = ifelse(config$snp_filtering, config$snp_rate/100, 1), 
  ind.callrate = ifelse(config$acc_filtering, config$acc_rate/100, 1), 
  maf = ifelse(config$maf_filtering, config$maf_rate/100, 0), 
  impute = (config$impute == 'TRUE'),
  message = TRUE
)

SNPs.clean <- SNPs[SNPs[,1] %in% colnames(qc$M.clean), 
                   colnames(SNPs) %in% c(hapmap_snp_attr, rownames(qc$M.clean))]

SNPs.clean <- cbind(SNPs.clean[, 1:11], t(qc$M.clean))

################################################################################
# SNPrune
if (config$snprune) {
  #' Reference:
  #' Calus, M.P.L., Vandenplas, J. SNPrune: an efficient algorithm to prune large 
  #' SNP array and sequence datasets based on high linkage disequilibrium. Genet 
  #' Sel Evol 50, 34 (2018). https://doi.org/10.1186/s12711-018-0404-z
  #' 
  #' GitHub Repository: 
  #' https://github.com/mariocalus/SNPrune
  #' 
  #' Installation instructions: (only first time, Linux machine required)
  #' chmod  +x SNPrune
  #' sudo cp SNPrune /sbin/SNPrune
  #' sudo cp seed_used.prv /sbin/seed_used.prv
  #' 
  #' SNPrune call code was developed by Alsamman, Alsamman M. (ICARDA)
  
  dataFolder <- 'SNPrune/Data'
  SNPInfo    <- SNPs.clean[, 1:11]
  hapMapNum  <- t(SNPs.clean[, -c(1:11)])
  
  #' delete data folder if it exists before staring new process
  if (file.exists(dataFolder)) unlink(dataFolder, recursive = TRUE)
  
  #' create new data folder
  dir.create(dataFolder)
  
  #' add row names in the "n:" formate
  rownames(hapMapNum) <- paste(1:nrow(hapMapNum), ':', sep = '')
  
  #' save the hapMapNumeric to a file with no header and no spaces
  write.table(hapMapNum,
              file.path(dataFolder, 'hapMapNumeric.txt'),
              quote = FALSE,
              sep = '',
              col.names = FALSE)
  
  #' use sed in system to replace ":" with " " in the file
  system(paste("sed -i 's/:/ /g' ",
               file.path(dataFolder, 'hapMapNumeric.txt'),
               sep = ''))
  
  #' create the SNPrune parameters file in temp folder in DataFolder
  sink(file.path(dataFolder, 'SNPrune.inp'))
  
  #' SNPn:
  #' total number of SNPs
  cat(ncol(hapMapNum), sep = "\n")
  
  #' genoType:
  #' types of data genotypes means 0,1,2, while alleles means A,C,G,T
  cat('genotypes', sep = "\n")
  
  #' inFile:
  cat('hapMapNumeric.txt', sep = "\n")
  
  #' outFile:
  cat('SNPruneOut.txt', sep = "\n")
  
  #' genoFormat:
  #' format of the genotype file: either with ("sparse") 
  #' or without spaces ("dense") between columns with genotypes
  cat('dense', sep = "\n")
  
  #' pruneLD:
  #' prune for SNP pairs with LD above a pre-define R2-threshold ("high_LD"; default), 
  #' or for SNP pairs in complete LD ("identical")
  cat('high_LD', sep = "\n")
  
  #' threads:
  cat('1', sep = "\n")
  
  #' R2threshold:
  cat('0.99', sep = "\n")
  
  #' lowestAllele:
  #' when using "alleles": lowest allele code used. 
  #' when using "genotypes": put "0"
  cat('0', sep = "\n")
  
  #' highestAllele:
  #' when using "alleles": highest allele code used. 
  #' when using "genotypes": put "1"
  cat('1', sep = "\n")
  
  #' end of the SNPrune parameters file
  sink()
  
  #' save the current working directory
  currentDir <- getwd()
  
  #' change the working directory to the temp folder in DataFolder
  setwd(dataFolder)
  
  #' run the SNPrune program
  system('SNPrune')
  
  #' change the working directory to the current directory
  #' fixing some issues with the output file
  #' remove the extra spaces in the output file of removed_SNPs.txt and 
  #' retained_SNPs.txt using perl -pi -e 's/\s+$//'
  system("perl -pi -e 's/[ ]+/ /g' removed_SNPs.txt")
  system("perl -pi -e 's/[ ]+/ /g' retained_SNPs.txt")
  
  #' remove first space in the output file of removed_SNPs.txt and
  #' retained_SNPs.txt using perl -pi -e 's/^\s+//'
  system("perl -pi -e 's/^[ ]+//g' removed_SNPs.txt")
  system("perl -pi -e 's/^[ ]+//g' retained_SNPs.txt")
  
  #' return the original working directory!
  setwd(currentDir)
  
  #' check summary stats
  cat(readLines(file.path(dataFolder, 'summary_stats.txt')), sep = "\n")
  
  #' read the SNPrune output file and ignore extra spaces
  retainedSNPs <- read.table(file.path(dataFolder, 'retained_SNPs.txt'),
                             sep = ' ',
                             stringsAsFactors = FALSE,
                             strip.white = TRUE,
                             row.names = 1)
  
  #' remove SNPs from the hapMap according to the removed_SNPs.txt file
  SNPs.clean <- SNPs.clean[retainedSNPs[, 1], ]
}

################################################################################
#' define the significantly associated markers and the explained variance
#' significance of marker-trait association is defined based on the simpleM 
#' method as described by Gao et al. 2008
#' 
#' Reference:
#' Gao, X.; Starmer, J.; Martin, E.R. A multiple testing correction method for 
#' genetic association studies using correlated single nucleotide polymorphisms. 
#' Genet. Epidemiol. 2008, 32, 361–369.
#' https://onlinelibrary.wiley.com/doi/10.1002/gepi.20310
if (config$calc_gao) {
  #' reshape SNPs data to fit in the code snippet shared by Marcel from IPK
  dta.G.GWAS <- SNPs.clean[, -c(2, 5:11)]
  
  #' the matrix has columns Marker,CHR, Position and one for every SNP (0,1,2)
  colnames(dta.G.GWAS) <- c('Marker', 'CHR', 'Position', colnames(dta.G.GWAS)[-c(1:3)])
  
  #' significance level for GWAS
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
  
  config$gwas_lod <- -log10(p)
}

################################################################################
#' 1. GAPIT (do GWAS)
if (config$perform_gwas) {
  #' https://www.zzlab.net/GAPIT/
  #' 
  #' In GAPIT3, the Bonferroni multiple test threshold is used by default to 
  #' determine significance, it only affect the visualization of results (Manhattan)
  #' but not change the association results, users can adjust the output p-values 
  #' by using the base function p.adjust() in R
  #' 
  #' By default, GAPIT3 conducts GWAS using the BLINK method, which has the 
  #' highest statistical power & computation efficiency among all methods 
  #' implemented (e.g., MLM)
  #' 
  #' All methods (GWAS and GS) have linear computing time to number of markers.
  #' However, they have mixed computing complexity to the number of individuals.
  #' Most of these methods have computing time complexity that are cubic to number
  #' of individuals.
  
  #' make sure that all accessions produced to gwas are exists 
  #' in both geno and pheno datasets (exists in pheno even if the trait value is NA)
  geno.gwas  <- SNPs.clean[, c(1:11, which(colnames(SNPs.clean) %in% BLUE[,1]))]
  pheno.gwas <- BLUE[BLUE[,1] %in% colnames(geno.gwas),]
  
  geno.GD <- t(geno.gwas[, -c(1:11)])
  geno.GD <- as.data.frame(geno.GD)
  geno.GD <- cbind(colnames(geno.gwas)[-c(1:11)], geno.GD)
  colnames(geno.GD) <- c('taxa', geno.gwas[,1])
  rownames(geno.GD) <- NULL
  
  geno.GM <- geno.gwas[, c(1,3,4)]
  colnames(geno.GM) <- c('SNP', 'Chromosome', 'Position')
  
  #' https://rdrr.io/github/jiabowang/GAPIT3/man/GAPIT.html
  tryLog(
    myGAPIT <- GAPIT(
      Y = pheno.gwas,
      # G = rbind(colnames(geno.gwas), geno.gwas),
      GD = geno.GD,
      GM = geno.GM,
      model = "MLM",
      PCA.total = 3,
      SNP.MAF = ifelse(config$maf_filtering, config$maf_rate/100, 0),
      output.numerical = TRUE,
      file.output = TRUE
    ), 
    write.error.dump.file = TRUE
  )
  
  #' move GWAS output files into the GWAS folder
  if (!file.exists('GWAS')) dir.create('GWAS')
  gwas_outputs <- list.files(pattern = 'GAPIT\\..*')
  file.rename(gwas_outputs, paste0('GWAS/', gwas_outputs))
  output_files <- list.files('GWAS', full.names = TRUE)
  
  #' get GWAS results/output
  GWAS <- myGAPIT$GWAS[,c(1:4,9)]
  colnames(GWAS) <- c('Marker', 'Chr', 'Pos', 'p', 'effect') 
  
  #' calculate the LOD score (logarithm of the odds) = -log10(pValue)
  GWAS$LOD <- -1 * log10(GWAS$p)
  
  #' select SNPs markers exceed the given LOD threshold
  sel_markers <- GWAS[GWAS$LOD > config$gwas_lod & !is.na(GWAS$p), 'Marker']
  # sel_markers <- read.delim("clipboard")[[1]]
  
  sel_SNPs <- SNPs.clean[which(SNPs.clean[,1] %in% sel_markers), c(1, which(colnames(SNPs.clean) %in% myGAPIT$GD[,1]))]
  sel_SNPs <- SNPs.clean[which(SNPs.clean[,1] %in% sel_markers), ]
  rownames(sel_SNPs) <- sel_SNPs[,1]
  sel_SNPs <- sel_SNPs[,-1]
  sel_SNPs <- as.data.frame(t(sel_SNPs[,-c(1:10)]))
  
  #' merge genetic matrix & phenotypic data using accession id as a matching key
  #' NOTE: not all BLUE has SNPs!
  data <- merge(pheno.gwas, sel_SNPs, by.x = 'Taxa', by.y = 0)
  
  #' remove accessions with missing trait value
  data <- data[complete.cases(data),]
  
  ################################################################################
  #' 2. stepwise regression
  if (config$stepwise && !is.null(sel_markers)) {
    #' Stepwise Selection for Large Data Sets
    #' Selecting linear and generalized linear models for large data sets using 
    #' modified stepwise procedure and modern selection criteria (like modifications 
    #' of Bayesian Information Criterion). Selection can be performed on data which 
    #' exceed RAM capacity.
    
    #' the stepwise procedure for big data
    #' first, we have to convert our data to a proper format, an object of class big
    reg_data <- bigstep::prepare_data(y = data[,config$trait], X = data[,-c(1:2)])
    
    #' then, we can use the stepwise procedure with, for example, aic or bic
    results <- bigstep::stepwise(reg_data, crit = bic)
    
    #' get the list of the markers in the final model
    sel_markers <- results$model
  } else {
    sel_markers <- colnames(data)[-c(1:2)]
  }
  
  #' sub-genotyping matrix includes only markers selected in the final stepwise 
  #' regression model
  rownames(data) <- data[,1]
  data <- subset(data, select = sel_markers)
  
  #' transpose sub-genotyping matrix back (SNPs in rows and accessions in columns)
  data <- as.data.frame(t(data))
  
  #' merge SNPs metadata for selected markers
  data <- merge(SNPs_metadata, data, by.x = 1, by.y = 0)
  
  #' save the output of sub-genotyping matrix of the selected markers for the next 
  #' pipeline step
  write.csv(data, config$sub_matrix)
  output_files <- c(output_files, config$sub_matrix)
  
  ################################################################################
  #' QTLome
  #' check if the candidate markers are already listed in QTLome or not yet
  if (config$QTLome_submit && !is.null(sel_markers)) {
    #' TODO: we have to agree on the QTLome file format/template
    #' 
    #' Position reference: 
    #' POPSEQ_2012, POPSEQ_2017, IBSC_2012, MorexV1_2017, MorexV2_2019, 
    #' MorexV3_2021, MorexGenome_floresta, and Igartua_2019
    #' 
    #' Trait vocabulary (as agreed in WP5): 
    #' - Thousand grain weight: Plant grain weight (13)? Kernel weight (7)
    #' - Plant height: Height (32)
    #' - Flowering time: Days to Heading (10)?
    
    QTLome <- read.csv('Tidy 220525 Barley QTLome_LK.csv')
    
    #' data.frame for QTLome reporting
    QTLome_report <- cbind(QTLome_add = 0, 
                           QTLome_name = NA, 
                           trait_affected = config$trait, 
                           data[,1:4])
    QTLome_report[5,]
    #' starting id number used for naming new markers
    id <- 1
    
    for(i in 1:nrow(QTLome_report)) {
      #' check if current marker match any existing QTLome entry 
      #' (i.e., same chrom and pos in the start-end range)
      QTLome_match <- which(QTLome$Chr   == QTLome_report[i, 'chrom'] & 
                              QTLome$START <= QTLome_report[i, 'pos'] & 
                              QTLome$STOP  >= QTLome_report[i, 'pos'])
      
      #' if find any, report it in the QTLome_name column 
      #' (use comma if we have more than one match)
      #' else, give it a specific naming style
      QTLome_name <- ifelse(length(QTLome_match) > 0, 
                            paste(unique(QTLome[QTLome_match, 'Number']), collapse = ', '),
                            paste0('Q.AGENT.', sprintf('%03d', id)))
      
      QTLome_report[i, 'QTLome_name'] <- QTLome_name
      
      if (length(QTLome_match) == 0) { 
        id <- id + 1
        QTLome_report[i, 'QTLome_add'] <- 1
      } 
    }
    
    write.csv(QTLome_report, config$QTLome_file)
    output_files <- c(output_files, config$QTLome_file)
  }
  
  if (config$min_subset && !is.null(sel_markers)) {
    ################################################################################
    #' 3. calculate MAF for each trait-landrace.
    
    #' calculate the minor allele frequency (MAF) and observation (MAO) 
    #' for each accession
    MAF_acc <- apply(data[,-c(1:4)], MARGIN = 2, function(x) sum(x == 2)/length(x))
    MAO_acc <- apply(data[,-c(1:4)], MARGIN = 2, function(x) sum(x == 2))
    
    #' formulate calculated stats in the data.frame format
    snps_stats <- as.data.frame(cbind(MAF_acc, MAO_acc))
    
    #' save the stats (e.g., MAF, MAO) for the next pipeline step
    write.csv(snps_stats, config$stats_file)
    output_files <- c(output_files, config$stats_file)
    
    ################################################################################
    #' 4. find the minimum subset of accessions that covers all candidate markers
    SNP_effect <- GWAS[GWAS$Marker %in% data$`rs#`, ]
    SNP_effect <- merge(data[,1:4], SNP_effect[,-c(2:3)], by = 1)
    
    minimum_group <- data.frame(Accession = character(), 
                                Cumulative_MAO = integer(), 
                                MAF = numeric())
    
    #' calculate the minor allele frequency (MAF) for each SNP
    MAF_snp <- apply(data[,-c(1:4)], MARGIN = 1, function(x) sum(x == 2)/length(x))
    
    #' calculate weighted SNPs scores to increase the chance of including accessions 
    #' represents less frequent SNPs, account for SNP effect and magnitude 
    # rest_wSNPs  <- cbind(data[,1:4], data[,-c(1:4)] * sqrt((1 / MAF_snp)))
    rest_wSNPs <- cbind(data[,1:4], data[,-c(1:4)] * SNP_effect[, 'effect'] * sqrt((1 / MAF_snp)))
    
    #' if SNP has negative effect, then the target SNP value is 0 (reference allele)
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
    
    write.csv(minimum_group, config$min_grp_acc_file, row.names = FALSE)
    output_files <- c(output_files, config$min_grp_acc_file)
    
    #' subset the genotyping matrix to have only the selected markers and the optimum 
    #' group of accessions
    minimum_group_snps <- data[, c(colnames(data)[1:4], minimum_group$Accession)]
    
    effect <- SNP_effect[, 'effect']
    occurs <- apply(minimum_group_snps, MARGIN = 1, function(x) sum(x == 0))
    occurs <- ifelse(effect > 0, length(minimum_group$Accession) - occurs, occurs)
    
    minimum_group_snps <- cbind(occurs, effect, minimum_group_snps)
    minimum_group_snps <- minimum_group_snps[, c(3, 1:2, 4:ncol(minimum_group_snps))]
    
    write.csv(minimum_group_snps, config$min_grp_snp_file, row.names = FALSE)
    output_files <- c(output_files, config$min_grp_snp_file)
  }
}

if (config$perform_gs) {
  if (config$eg_blup) {
    ################################################################################
    #' reshapping the genotype data for the EG-BLUP section
    metadata <- SNPs.clean[, 1:4]
    metadata$ref <- substr(metadata$alleles, 1, 1)
    metadata$minor <- substr(metadata$alleles, 3, 3)
    metadata$alleles <- NULL
    colnames(metadata) <- c('Marker', 'CHR', 'LOC', 'REF', 'ALT')
    
    geno.egblup <- cbind(metadata, SNPs.clean[, -c(1:11)])
    
    #' reshapping the phenotype data for the EG-BLUP section
    
    #' Only for entries with a given phenotypic data:
    # pheno.egblup <- BLUE
    # rownames(pheno.egblup) <- NULL
    
    #' All entries with a given genotypic data:
    pheno.egblup <- merge(colnames(geno.egblup)[-c(1:5)], BLUE, 
                          by.x = 1, by.y = 1, all.x = TRUE, sort = FALSE)
    
    colnames(pheno.egblup) <- c('Genotype_Matrix', 'Trait')
    
    ################################################################################
    #
    #         Sample R code for EG-BLUP
    #
    #         Code is a derivative of the publication: 
    #         Choosing the right tool: Leveraging of plant genetic resources in wheat (Triticum aestivum L.) 
    #         benefits from selection of a suitable genomic prediction model
    #----------------------------------------------------------------------------------------------------------
    #         Marcel O. Berkner1, Albert W. Schulthess1, Yusheng Zhao1, Yong Jiang1, Markus Oppermann1 and Jochen C. Reif1
    #----------------------------------------------------------------------------------------------------------
    #         1 Leibniz Institute of Plant Genetics and Crop Plant Research (IPK), 06466 Stadt Seeland, Germany
    #         
    #         Corresponding author: Jochen C. Reif, Email: reif@ipk-gatersleben.de
    #----------------------------------------------------------------------------------------------------------
    #This code was developed with R version 4.0.2
    #
    
    #' Load genotypic data and pure based on minor-allele-frequency
    Geno <- geno.egblup
    
    #' Select Marker on minor allele frequency threshold
    #' Set minor allele frequency threshold (example of 1%)
    # MAF <- 0.01
    MAF <- ifelse(config$maf_filtering, config$maf_rate/100, 0)
    
    #' Omit additional information on SNP markers
    Freq <- data.frame(c(1:nrow(Geno)), Geno[,'Marker'], (Geno[,c(6:ncol(Geno))]))
    
    #' Transform SNP matrix into numeric
    Freq[, c(3:ncol(Freq))] <- sapply(Freq[, c(3:ncol(Freq))], as.numeric)
    
    #' Calculate frequencies of both alleles
    Freq <- data.frame(Freq[, c(1,2)],
                       (rowMeans(Freq[, c(3:ncol(Freq))]))/2,
                       1 - ((rowMeans(Freq[, c(3:ncol(Freq))]))/2))
    
    colnames(Freq) <- c('RowNo', 'Marker', 'freq.allele1', 'freq.allele2')
    f <- data.frame(Freq[,c(1,2)],c(1:nrow(Freq)))
    
    for(i in 1:nrow(Freq)) {
      #' Select the minor allele
      if (Freq[i,3] <= 0.5) {
        f[i,3] <- Freq[i,3]
      } else {
        f[i,3] <- Freq[i,4]
      }
    }
    
    colnames(f) <- c('RowNo', 'Marker', 'MAF')
    f <- f[f[,3] >= MAF,]
    
    #' Selected SNP according to minor allele frequency threshold
    selmarker <- as.character(f[,2])
    
    #' Create G matrix, marker annotation of -1,0,1
    G <- data.frame(Geno[,c('Marker', 'CHR', 'LOC')], Geno[, c(6:ncol(Geno))]-1)
    
    #' Reduction of matrix to selected markers
    G <- G[c(which(G$Marker %in% selmarker)),]
    
    #' Load phenotypic data
    #' Load phenotype  matrix; phenotype information is given as BLUE
    Pheno <- pheno.egblup
    accession.order <- colnames(G)
    accession.order <- accession.order[4:length(colnames(G))]
    
    #' Bring phenotypic and genotypic information in same order
    Pheno <- Pheno[(match(accession.order, Pheno$Genotype_Matrix)),] 
    
    #' Select trait if multiple traits exist
    dta.P <- Pheno
    
    #' NOTE: stop this line to get prediction for all entries with genotypic data
    if (!config$egblup_all) { 
      dta.P <- dta.P[!(is.na(dta.P$Trait)),]
    }
    
    p <- as.matrix(dta.P[,'Trait'])
    
    #' Index-Matrix for runs of five-fold-cross-validation 
    #' Randomly assigned index-matrix
    
    #' Select number of five-fold cross-validations
    num.crossvalid <- 2
    
    #' random split of each round of cross-validation into five parts
    idx.M <- data.frame(as.numeric(row.names.data.frame(dta.P)),
                        c(createFolds(as.numeric(row.names.data.frame(dta.P)), k = 5, list = FALSE)))
    
    for (i in 1:(num.crossvalid - 1)) {
      idx.M <- data.frame(idx.M, c(createFolds(as.numeric(row.names.data.frame(dta.P)), k = 5, list = FALSE)))
    }
    
    #' column name according to the round of cross-validation
    colnames(idx.M) <- c('RowNo', 1:num.crossvalid)
    
    #' Counterpart to index-matrix for running the models
    x    <- c(1:(num.crossvalid * 5))
    run  <- c(rep(1:num.crossvalid, each = 5, times = 1))
    fold <- c(rep(1:5, each = 1, times = length(x)/5))
    x.M  <- data.frame(x, run, fold)
    
    rm(run,fold)
    
    #' Create common vector for the evaluation of runs
    run <- c(1:num.crossvalid)
    run <- as.numeric(run)
    cor <- c(1:num.crossvalid)
    cor <- as.numeric(cor)
    
    #' EG-BLUP
    #' Addititive Relationship Matrix
    #' First method based on VanRaden, 2008
    
    #' Create marker matrix M with dimension n x m for n genotypes and m markers
    M <- Geno[c(which(Geno$Marker %in% selmarker)), 6:ncol(Geno)]
    rownames(M) <- Geno[c(which(Geno$Marker %in% selmarker)), 'Marker']
    
    M <- M[, c(dta.P$Genotype_Matrix)]
    M <- M[!apply(M, 1, sd) == 0,]
    M <- M - 1
    M <- t(M)
    
    freq.p <- (colSums(M == 1) + 1/2 * colSums(M == 0))/nrow(M)
    P <- 2*(t(matrix(rep(freq.p,nrow(M)), ncol = nrow(M))) - 0.5)
    Z <- M - P
    ZZ <- Z %*% t(Z)
    vp <- 2*sum(freq.p*(1 - freq.p))
    A <- ZZ/vp
    
    rm(Z, ZZ, P, vp, freq.p)
    
    #' Covariance Matrix H
    #' Hadamard Product of A as described by Henderson 1985
    H <- A * A
    rownames(H) <- dta.P$Genotype_Matrix
    colnames(H) <- dta.P$Genotype_Matrix
    
    #' Function of the model
    para.EGBLUP <- function(i) {
      #' Link the current run and fold with information which is assigned in the index-matrix
      idx <- idx.M[, c(1, (x.M[i, 'run'] + 1))]
      
      p.train <- dta.P$Trait
      
      #' Replace the phenotypic values of the test set with NA
      p.train[idx[idx[, 2] == x.M[i, 'fold'], 1]] <- NA
      
      # set.seed(1002)
      EGBLUP <- BGLR::BGLR(y = p.train,
                           ETA = ETA2,
                           nIter = 10000,  # set >= 10000 for large data
                           burnIn = 3000,  # set >= 3000 for large data
                           verbose = FALSE)  
      
      result[[i]] <- list(list(data.frame(dta.P$Genotype_Matrix,EGBLUP[1],EGBLUP[14],EGBLUP[15])),
                          list(data.frame(EGBLUP[6],EGBLUP[7],EGBLUP[16],EGBLUP[17],EGBLUP[18],EGBLUP[19])))
    }
    
    result <- list()
    ETA2   <- list(list(K1 = A, model = 'RKHS'), list(K2 = H, model = 'RKHS'))
    
    #' Select number of cores based on the available resources
    no_cores <- detectCores() - 2
    
    cl <- makeCluster(no_cores)
    clusterExport(cl = cl, varlist = list('dta.P', 'ETA2', 'idx.M', 'x.M', 'result'))
    
    tryLog(result <- parLapply(cl, x, para.EGBLUP), write.error.dump.file = TRUE)
    stopCluster(cl)
    
    #' Evaluation of model
    Pred.abil.EGBLUP <- data.frame(run, cor)
    PredEB <- data.frame(ncol(5))
    
    #' Combine the results for all folds of the five-fold cross-validations with 
    #' the information of the respective fold of the cross-validation
    for(j in 1:length(x)){
      result.loop <- as.data.frame(result[[j]][1])
      
      #' NOTE: stop this line to get prediction for all entries with genotypic data
      if (!config$egblup_all) { 
        test.set <- c(which(is.na(result.loop$y)))
      }
      
      p.test <- data.frame(dta.P,
                           result.loop$yHat,
                           c(rep(j, each = 1, times = nrow(dta.P))),
                           c(rep(x.M[j, 'run'], each = 1, times = nrow(dta.P))),
                           c(rep(x.M[j, 'fold'], each = 1, times = nrow(dta.P))))
      
      #' NOTE: stop this line to get prediction for all entries with genotypic data
      if (!config$egblup_all) { 
        p.test <- p.test[test.set,]
      }
      
      PredEB <- rbind(PredEB, p.test)
    }
    
    colnames(PredEB) <- c('Genotype_Matrix', 'Trait', 'pred.Trait', 'x', 'run', 'fold')
    
    #' Calculate correlation between all folds of each five-fold cross-validations
    for(i in 1:length(unique(x.M$run))) {
      cor.set <- PredEB[c(which(PredEB$run == i)),]
      Pred.abil.EB <- cor(cor.set[, c('Trait', 'pred.Trait')])
      Pred.abil.EGBLUP[Pred.abil.EGBLUP[, 'run'] == i, 2] <- Pred.abil.EB[1, 2]
    }
    
    #' NOTE: update this line to get prediction for all entries with genotypic data
    if (!config$egblup_all) { 
      crossvalid_result <- aggregate(cbind(Trait, pred.Trait) ~ Genotype_Matrix, data = PredEB, FUN = mean)
    } else {
      crossvalid_result <- aggregate(cbind(pred.Trait) ~ Genotype_Matrix, data = PredEB, FUN = mean)
    }
    
    write.csv(crossvalid_result, config$egblup_results_file, row.names = FALSE)
    output_files <- c(output_files, config$egblup_results_file)
    
    unlink(c('mu.dat', 'varE.dat', 'ETA_1_varU.dat', 'ETA_2_varU.dat'), force = TRUE)
  }
}

zip(paste0(shiny_app, 'www/output/', job_id, '.zip'), output_files)
unlink(output_files, force = TRUE, recursive = TRUE)

#' send an email from R script
#' https://pkgs.rstudio.com/blastula/
