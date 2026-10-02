#Dust exposure matching for SOMI data
#Using the Grid-based approximation and vectorized extraction (GAVE) method

rm(list=ls())

library(dplyr)
library(tictoc)
library(data.table)

#Read in SOMI dataset
#Note: this dataset is not provided under a data use agreement with CDPH
somi_dust<-readRDS("data/somi_dust.RDS")

#Variables: 
#VS_unique: ID
#DOB: birth date
#X_year: birth year
#lat: latitude for address
#lon: longitude for address
#x_approx: approximated longitude for address (rounded to 0.01*0.01 grid resolution)
#y_approx: approximated latitude for address (rounded to 0.01*0.01 grid resolution)
#period: exposure periods (13-23 biweekly periods per ID, depending on #gestational weeks), in the same format as the dust data variables

#Read in cleaned vD dust data
vD_dust<-readRDS("data/vD_dust_CA_wide.RDS")
#Variables:
#lon: longitude for exposure grid point
#lat: latitude for exposure grid point
#dyyyyddd: average dust PM2.5 during the biweekly period from the dddth day of year yyyy to the (ddd+13)th day of year yyyy


#Read in pre-prepared dataframe that has the points needed for matching with different buffer
point_list<-readRDS("data/buffer_point_list.RDS")
#Note: The code for creating this list is included in the "01.data_preparation.R" document


#########################################################
#Code for GAVE exposure matching

#prepare datasets in data.table format for fast extraction
setDT(somi)
setDT(somi_dust)
setDT(point_list)
setDT(vD_dust)
setkey(vD_dust, lon, lat)
setkey(somi_dust, VS_unique) 

n <- nrow(somi)

#set up empty columns to store the results
somi_dust[, `:=`(
  dust_buff_15 = as.numeric(NA),
  dust_buff_10 = as.numeric(NA),
  dust_buff_5  = as.numeric(NA),
  dust_buff_2  = as.numeric(NA)
)]


# -------- FUNCTION --------

gave_one_id <- function(i) {
  # Extract ID and approximated grid coordinates for this subject
  id <- somi[i, VS_unique]
  x_approx <- somi[i, x_approx]
  y_approx <- somi[i, y_approx]
  
  # Generate target grid points by applying pre-computed relative offsets
  target_points <- copy(point_list)
  target_points[, `:=`(
    lon = round(x_approx + lon_rel, 5),
    lat = round(y_approx + lat_rel, 5)
  )]
  
  # Get exposure time periods of interest for this subject
  time <- somi_dust[VS_unique == id, period]
  
  # Join exposure data at target grid points
  tp <- vD_dust[target_points, on = .(lon, lat)]
  
  # Keep only relevant columns: coordinates, buffer flags, and selected time periods
  tp <- tp[, c("lon", "lat", "buff15", "buff10", "buff5", "buff2", time), with = FALSE]
  
  # Identify exposure value columns (exclude coordinates and buffer indicators)
  value_cols <- setdiff(names(tp), c("lon", "lat", "buff15", "buff10", "buff5", "buff2"))
  
  # Compute mean exposure within each buffer
  somi_dust[VS_unique == id, `:=`(
    dust_buff_15 = unname(colMeans(tp[buff15 == 1, ..value_cols], na.rm = TRUE)),
    dust_buff_10 = unname(colMeans(tp[buff10 == 1, ..value_cols], na.rm = TRUE)),
    dust_buff_5  = unname(colMeans(tp[buff5  == 1, ..value_cols], na.rm = TRUE)),
    dust_buff_2  = unname(colMeans(tp[buff2  == 1, ..value_cols], na.rm = TRUE))
  )]
}

#test for one iteration
tic()
gave_one_id(1)
toc()


#Begin GAVE process for all subjects

# -------- PROGRESS BAR --------
pb <- txtProgressBar(min = 0, max = n, style = 3) # Initialize the progress bar

# -------- LOOP --------

for (i in 1:n) {
  gave_one_id(i) # GAVE computation
  setTxtProgressBar(pb, i) # Update the progress bar
}

close(pb) # Close the progress bar after the loop

# -------- SAVE --------
saveRDS(somi_dust,"data/somi_dust.RDS")
#somi_dust<-readRDS("data/somi_dust.RDS")

