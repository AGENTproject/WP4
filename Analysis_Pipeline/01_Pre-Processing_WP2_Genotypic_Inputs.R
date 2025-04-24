library(ASRgenomics)
library(data.table)
library(R.utils)
library(inline)
library(vcfR)

rm(list=ls())
gc(reset = TRUE)
setwd('~/AGENT')

# C code to convert markers data into numeric notation effectively
source('00_SNPs_Numeric_Recoding.R')

# ?ASRgenomics::qc.filtering
my.marker.callrate <- 0.2
my.ind.callrate    <- 0.2
my.maf             <- 0.01
my.heterozygosity  <- 1
my.Fis             <- 1

crop <- 'Wheat' # Wheat or Barley (case sensitive!)

# original VCF on FAIRDOM filtered by WP2 for presence rate >= 80% and MAF > 1%
# this step was critical to reduce the size and enable this script to handle it

if (crop == 'Wheat') {
  # https://urgi.versailles.inrae.fr/fairdom/assays/72
  src.url <- 'https://files.ipk-gatersleben.de/file/lrvo2KjNXy9t7OU4/feh32gcUtTiAu5mX/gbs_diversity_AGENT_wheat_ChineseSpringV2_1_biosampleIDs_MAF0o01_PP0o8_ENA.vcf.gz'

  # source: https://urgi.versailles.inrae.fr/fairdom/data_files/423
  # includes AGENT IDs and IPK ACCENUMBs alongside the BioSample IDs
  ids.file <- 'IDs_231012.csv'
} else {
  # https://urgi.versailles.inrae.fr/fairdom/assays/63
  src.url <- 'https://files.ipk-gatersleben.de/file/LnsLjnr9jEZeWT8b/IXXwJGSlpZVDiWH5/gbs_diversity_AGENT_barley_MorexV3_biosampleIDs_MAF0o01_PP0o8_ENA.vcf.gz'

  # source: https://urgi.versailles.inrae.fr/fairdom/documents/1885
  # includes AGENT IDs and IPK ACCENUMBs alongside the BioSample IDs
  ids.file <- 'IDs_230717.csv'
}

dst.file <- paste0('./AGENT_', crop, '_VCF_Assay/', basename(src.url))
ids.file <- paste0('./AGENT_', crop, '_VCF_Assay/', ids.file)
snp.file <- paste0('./AGENT_', crop, '_VCF_Assay/AGENT_', crop, '.csv.gz')
map.file <- paste0('./AGENT_', crop, '_VCF_Assay/AGENT_', crop, '_map.csv')
kin.file <- paste0('./AGENT_', crop, '_VCF_Assay/AGENT_', crop, '_kinship.csv.gz')

# download the filtered VCF file if it is not exists
if (!file.exists(dst.file)) {
  options(timeout = max(300, getOption('timeout')))
  utils::download.file(src.url, dst.file)
}

# read the vcf file
vcf.data  <- vcfR::read.vcfR(dst.file)

# convert vcf file format into the required hapmap format
hapmap <- vcfR::vcfR2hapmap(vcf.data)
hapmap <- hapmap[-1,]

# free memory in the working space
rm(vcf.data)

# remove extra hapmap columns and get the SNP matrix
geno.data <- hapmap[, -c(1:11)]

# assign the marker names to the matrix row names
rownames(geno.data) <- hapmap$rs

# transpose the SNP matrix to have acc. in rows and markers in columns
geno.data <- t(geno.data)

################################################################################
#' use NA to code missing values instead of 'NN' string, then
#' recode markers using numeric notation 0, 1, 2 
#' (1 for heterozygous and 2 for reference allele)
#'
# geno.data[geno.data == 'NN'] <- NA
# geno.data <- ASRgenomics::snp.recode(M = geno.data)$Mrecode
################################################################################

# encode using single letter representation (heterozygous as generic X)
geno.data[geno.data == 'NN'] <- 'N'
geno.data[geno.data == 'AA'] <- 'A'
geno.data[geno.data == 'TT'] <- 'T'
geno.data[geno.data == 'CC'] <- 'C'
geno.data[geno.data == 'GG'] <- 'G'

geno.data[!geno.data %in% c('N', 'A', 'T', 'C', 'G')] <- 'X'

# split the SNPs rs to get the marker metadata (e.g., Chr1A:1235919:A:G)
metadata <- data.table::tstrsplit(hapmap$rs, '\\:')

# convert it into geno map data.frame
geno.map <- as.data.frame(do.call(cbind, metadata))

# add the original snp name and allele columns
geno.map <- cbind(hapmap$rs, paste0(geno.map[, 3], '/', geno.map[, 4]), geno.map)

# get the standard column names
colnames(geno.map) <- c('rs', 'alleles', 'chrom', 'pos', 'ref', 'alt')

# call the C function to convert markers data into numeric notation effectively 
# geno.data <- snps_numeric_recoding(geno.data, geno.map[,5])

numeric.geno.data <- NULL

# loop through the data frame in batches
for (start_row in seq(1, nrow(geno.data), by = 15000)) {
  
  # define the end row for the current batch
  end_row <- min(start_row + 15000 - 1, nrow(geno.data))
  
  numeric.geno.data <- rbind(numeric.geno.data, snps_numeric_recoding(geno.data[start_row:end_row,], geno.map[,5]))
}

geno.data <- numeric.geno.data
rm(numeric.geno.data)

rownames(geno.data) <- colnames(hapmap)[-c(1:11)]

# filtering the genomic dataset
M_filter <- ASRgenomics::qc.filtering(M      = as.matrix(geno.data), 
                                      maf    = my.maf,
                                      Fis    = my.Fis,
                                      impute = FALSE, 
                                      marker.callrate = my.marker.callrate,
                                      ind.callrate    = my.ind.callrate,
                                      heterozygosity  = my.heterozygosity)

################################################################################
# # read the id mapping metadata
# ids <- read.csv(ids.file)
# 
# # use AGENT id when exists and IPK acc. number elsewhere (NA)
# ids$ID <- ifelse(!is.na(ids$AGENT_ID), ids$AGENT_ID, ids$IPK_ACCENUMB)
# 
# # get the AGENT id mapped to the existing biosample id
# AGENT_ID <- merge(rownames(M_filter$M.clean), ids, by.x = 1, by.y = 1, all.x = TRUE, sort = FALSE)$ID
# 
# # bind the agent ids as the first column in the filtered data.frame
# geno.data <- cbind(AGENT_ID, M_filter$M.clean)
# 
# # remove duplicated samples for the same acc. id (first match)
# geno.data <- geno.data[!duplicated(geno.data[,1]),]
################################################################################

# bind the BIOSAMPLE ids as the first column in the data.frame
BIOSAMPLE_ID <- rownames(M_filter$M.clean)
geno.data <- cbind(BIOSAMPLE_ID, M_filter$M.clean)

# save the filtered genotypic matrix in a gz csv file format
data.table::fwrite(geno.data, snp.file)

# save the map in a normal csv file
geno.map <- geno.map[geno.map$rs %in% colnames(geno.data),]
write.csv(geno.map, map.file, row.names = FALSE)

### Kinship matrix #############################################################

geno.data <- as.data.frame(data.table::fread(file = snp.file))
rownames(geno.data) <- geno.data[,1]
geno.data <- geno.data[,-1]
# View(geno.data[1:50,1:50])

# calculate the kinship matrix
kinship <- ASRgenomics::G.matrix(M = as.matrix(geno.data), method = 'VanRaden')$G
write.csv(kinship, gzfile(kin.file))
