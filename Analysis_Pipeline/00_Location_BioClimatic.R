library(QBMS)

# remove everything in the working environment
rm(list = ls())

# set working directory
setwd('~/AGENT')

# Wheat or Barley (case sensitive!)
crop  <- 'Wheat'

# read phenotypic input data
pheno_acc <- read.csv(paste0('WP3_BLUEs_Inputs/Final_', crop, '_with_Biosample_ID.csv'))

# keep records have coordinates values
pheno_acc <- pheno_acc[!is.na(pheno_acc$Latitude) & !is.na(pheno_acc$Longitude), ]

# keep the unique list of coordinates
unique_loc <- pheno_acc[!duplicated(pheno_acc[c('Latitude', 'Longitude')]), c('Latitude', 'Longitude')]

# download TerraClimate netCDF data files
QBMS::ini_terraclimate('1980-01-01', '2023-12-31', c('ppt', 'tmin', 'tmax'), data_path = './Climate_Data/')

# extract climate data from TerraClimat and calculate biovars using QBMS package (calculations takes ~12 hours for ~5000 locations)
data <- QBMS::get_terraclimate(unique_loc$Latitude, unique_loc$Longitude, '1980-01-01', '2023-12-31', c('ppt', 'tmin', 'tmax'), offline = TRUE, data_path = './Climate_Data/')

# create empty data.frame 
biovars <- data.frame()

for (i in 1:nrow(unique_loc)) {
  # skip to the next record if there is no data
  if (all(is.na(data$biovars[[i]][[1]]))) next
  
  # get the biovars per year for current coordinats/location
  # https://www.worldclim.org/data/bioclim.html
  temp <- data$biovars[[i]]
  
  # add current location coordinates info
  temp$Latitude  <- unique_loc$Latitude[i]
  temp$Longitude <- unique_loc$Longitude[i]
  
  # bind current location biovars to the global data.frame
  # one record for each year/latitude/longitude
  biovars <- rbind(biovars, temp)
}

avg_biovars <- aggregate(. ~ Latitude + Longitude, subset(biovars, select = -c(year)), mean)
acc_coords  <- pheno_acc[, c('Biosample', 'Latitude', 'Longitude')]
acc_climate <- merge(acc_coords, avg_biovars, by = c('Latitude', 'Longitude'), all.x = TRUE, sort = FALSE)

# rownames(acc_climate) <- acc_climate$Biosample
# acc_climate <- acc_climate[, -c(1:3)]

write.csv(acc_climate, paste0('./WP3_BLUEs_Inputs/', crop, '_climate_matrix.csv'))
