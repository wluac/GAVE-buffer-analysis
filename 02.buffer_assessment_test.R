#Test run: buffer-based spatial exposure assessment
#Compare computational speed of three methods:
# (1) Grid-based approximation and vectorized extraction (GAVE)
# (2) conventional x km circular buffer analysis
# (3) exactextract R

#All intermediate data sets are provided with the Rcode 

rm(list=ls())
set.seed(10092025)

library(dplyr)
library(sf)
library(tictoc)
library(ggplot2)
library(corrplot)
library(data.table)
library(terra)
library(tidyterra)
library(exactextractr)

##################################################################
#(0) Read in intermediate data sets

#read in california boundary map
CAmap<-read_sf("CA_State_Boundary.gpkg")

#read in data for 1000 random addresses in California
sample_address<-readRDS("sample_address.RDS")

#make point object and set projection
sample_address_pt <- sample_address%>%
  mutate(x=lon)%>%mutate(y=lat)%>%
  st_as_sf(coords = c("lon", "lat"),crs = 4326)

#transform the projection for sample address
sample_address_pt <- st_transform(sample_address_pt, st_crs(CAmap))

#check location
ggplot()+geom_sf(data=CAmap)+geom_sf(data = sample_address_pt,color = "red")


#read in van Donkelaar dust data and assign exposure period
vD_dust<-readRDS("vD_dust_CA_wide.RDS")

#make point object and set projection
vD_dust_pts <- vD_dust%>%
  mutate(x=lon)%>%mutate(y=lat)%>%
  st_as_sf(coords = c("lon", "lat"),crs = 4326)

#transform the projection for vD_dust data
vD_dust_pts <- st_transform(vD_dust_pts, st_crs(CAmap))

###################################################################################
#(1) Buffer analysis with Grid-based approximation and vectorized extraction (GAVE)

#First, select a lat and lon address in the center of California as a reference location
center<-data.frame(y_center = 37.16, x_center = -119.53)

#make point object and set projection
center <- st_as_sf(center, coords = c("x_center", "y_center"),crs = 4326)
#transform the projection 
center <- st_transform(center, crs = st_crs(CAmap))

#create a buffer for that approximated center point
buffer_km<-10
center_buffer<-st_buffer(center, dist = buffer_km*1000)

#check location, make sure it's in California and around 37 degree lat
ggplot()+geom_sf(data=CAmap)+geom_sf(data = center_buffer,color = "red")

#intersect the buffer area with the vD data, get the points within the buffer
idx <- st_intersects(center_buffer, vD_dust_pts, sparse = FALSE)[1,]
point_list<-data.frame(x=vD_dust_pts[idx,]$x, y=vD_dust_pts[idx,]$y)

#calculate the relative lon and lat of these points compared to the approximated center point
point_list<-point_list%>%
  mutate(lon_rel = x +119.53)%>%
  mutate(lat_rel = y -37.16)%>%
  mutate(buff10 = 1)

#summarize the relative lat and lons (should sum to 0)
summary(point_list$lat_rel)
summary(point_list$lon_rel)

#keep three variables
point_list<-point_list[,c("lon_rel","lat_rel","buff10")]

#Note: for each sample address, we can find points at these relative locations, 
#and calculate average values at these points

#create an empty vector to store the results
approx_results<-rep(as.numeric(NA), nrow(sample_address))

#prepare objects for faster extraction
setDT(sample_address)
setDT(vD_dust)
setkey(vD_dust, lon, lat)
setDT(point_list)

tic("Approximated buffer exposure assessment")
#for each sample address
for (i in 1:1000) {
  #find the approximated point at the center of four points in the exposure grids
  x_approx = round(sample_address[i,lon], digits = 2)
  y_approx = round(sample_address[i,lat], digits = 2)
  #get the list of coordinates for points that should be included in the calculation
  target_points <- point_list%>%
    filter(buff10==1)%>%
    mutate(lon = round (x_approx + lon_rel, digits = 5))%>%
    mutate(lat = round (y_approx + lat_rel, digits = 5))
  #get the time period of interest
  time = sample_address[i, period]
  #get the values for the time peiod of interest at each of these points on the list
  target_points <- vD_dust[target_points, on = .(lon, lat)]
  target_points <- target_points[, c("lon", "lat", time), with = FALSE]
  #calculate the average of these dust values
  approx_results[i]<-unname(colMeans(target_points, na.rm = T)[3])
}
toc()

#Approximated buffer exposure assessment: 4.69 sec elapsed

###########################################################################
#(1.5) Use grid-based approximation to calculate different buffers together

#read in point_list dataset 
#this contains the relative x- and y- offsets of the grid points
#covered in 2-, 5-, 10-, and 15-km buffers
point_list<-readRDS("buffer_point_list.RDS")


#Use the same method in (1) to calculate average dust levels within different buffers
#create an empty vector to store the results
approx_results_multi<-data.frame(
  id = sample_address$id,
  x = sample_address$lon,
  y = sample_address$lat,
  dust_buff15 = as.numeric(NA),
  dust_buff10 = as.numeric(NA),
  dust_buff5 = as.numeric(NA),
  dust_buff2 = as.numeric(NA)
)

#prepare objects for faster extraction
setDT(approx_results_multi)
setDT(point_list)
setDT(vD_dust)
setDT(sample_address)
setkey(vD_dust, lon, lat)

tic("Approximated buffer exposure assessment")
#for each sample address
for (i in 1:1000) {
  #find the approximated point at the center of four points in the exposure grids
  x_approx = round(sample_address[i, lon], digits = 2)
  y_approx = round(sample_address[i, lat], digits = 2)
  #get the list of coordinates for points that should be included in the calculation
  target_points <- copy(point_list)
  target_points[, `:=`(
    lon = round(x_approx + lon_rel, 5),
    lat = round(y_approx + lat_rel, 5)
  )]
  #get the time period of interest
  time = sample_address[i, period]
  #get the values for the time peiod of interest at each of these points on the list
  target_points <- vD_dust[target_points, on = .(lon, lat)]
  target_points <- target_points[, c("lon","lat","buff15","buff10","buff5","buff2", time), with = FALSE]
  #calculate the average of these dust values
  approx_results_multi[i, `:=`(
    dust_buff15 = unname(colMeans(target_points[buff15 == 1, ][, -(1:6)], na.rm = TRUE)),
    dust_buff10 = unname(colMeans(target_points[buff10 == 1, ][, -(1:6)], na.rm = TRUE)),
    dust_buff5  = unname(colMeans(target_points[buff5  == 1, ][, -(1:6)], na.rm = TRUE)),
    dust_buff2  = unname(colMeans(target_points[buff2  == 1, ][, -(1:6)], na.rm = TRUE))
  )]
}
toc()
#Approximated buffer exposure assessment: 7.58 sec elapsed

#remove intermediate objects to clear environment before the next step
rm(i, time, x_approx, y_approx, target_points)

#####################################################################
# (2) Use conventional buffer analysis to calculate different buffers

#read in un-transformed dust data again
vD_dust<-readRDS("vD_dust_CA_wide.RDS")
#make point object and set projection
vD_dust_pts <- vD_dust%>%
  mutate(x=lon)%>%mutate(y=lat)%>%
  st_as_sf(coords = c("lon", "lat"),crs = 4326)
#transform the projection for vD_dust data
vD_dust_pts <- st_transform(vD_dust_pts, st_crs(CAmap))

#read in un-transformed sample address again
sample_address<-readRDS("sample_address.RDS")
#make point object and set projection
sample_address_pt <- sample_address%>%
  mutate(x=lon)%>%mutate(y=lat)%>%
  st_as_sf(coords = c("lon", "lat"),crs = 4326)
#transform the projection for sample address
sample_address_pt <- st_transform(sample_address_pt, st_crs(CAmap))

#set empty dataframe to store conventional buffer results
buffer_results_multi<-sample_address%>%
  mutate(dust_buff15=as.numeric(NA))%>%
  mutate(dust_buff10=as.numeric(NA))%>%
  mutate(dust_buff5=as.numeric(NA))%>%
  mutate(dust_buff2=as.numeric(NA))

#15 km buffer:
#Set buffer size and create buffer around addresses
buffer_km<-15
address_buffers<-st_buffer(sample_address_pt, dist = buffer_km*1000)

#check location
ggplot()+geom_sf(data=CAmap)+geom_sf(data = address_buffers,color = "red")

#compute mean exposure within each buffer for the corresponding period
tic("Conventional buffer exposure assessment")
buffer_results_multi$dust_buff15 <- mapply(function(buffer, time) {

  #Keep only points within buffer
  idx <- st_intersects(buffer, vD_dust_pts, sparse = FALSE)[1,]
  #If no points in buffer, return NA
  if (!any(idx)) return(NA_real_)
  #Compute mean dust exposure
  mean(vD_dust_pts[idx,time, drop=T], na.rm=T)
  
}, st_geometry(address_buffers), address_buffers$period)
toc()

#Conventional buffer exposure assessment: 1718.89 sec elapsed

#10km buffer
#Set buffer size and create buffer around addresses
buffer_km<-10
address_buffers<-st_buffer(sample_address_pt, dist = buffer_km*1000)

#check location
ggplot()+geom_sf(data=CAmap)+geom_sf(data = address_buffers,color = "red")

#compute mean exposure within each buffer for the corresponding period
tic("Conventional buffer exposure assessment")
buffer_results_multi$dust_buff10 <- mapply(function(buffer, time) {
  
  #Keep only points within buffer
  idx <- st_intersects(buffer, vD_dust_pts, sparse = FALSE)[1,]
  #If no points in buffer, return NA
  if (!any(idx)) return(NA_real_)
  #Compute mean dust exposure
  mean(vD_dust_pts[idx,time, drop=T], na.rm=T)
  
}, st_geometry(address_buffers), address_buffers$period)
toc()
#Conventional buffer exposure assessment: 1742.63 sec elapsed


#5 km buffer:
#Set buffer size and create buffer around addresses
buffer_km<-5
address_buffers<-st_buffer(sample_address_pt, dist = buffer_km*1000)

#check location
ggplot()+geom_sf(data=CAmap)+geom_sf(data = address_buffers,color = "red")

#compute mean exposure within each buffer for the corresponding period
tic("Traditional buffer exposure assessment")
buffer_results_multi$dust_buff5 <- mapply(function(buffer, time) {
  
  #Keep only points within buffer
  idx <- st_intersects(buffer, vD_dust_pts, sparse = FALSE)[1,]
  #If no points in buffer, return NA
  if (!any(idx)) return(NA_real_)
  #Compute mean dust exposure
  mean(vD_dust_pts[idx,time, drop=T], na.rm=T)
  
}, st_geometry(address_buffers), address_buffers$period)
toc()
#Conventional buffer exposure assessment: 1727.23 sec elapsed

#2 km buffer:
#Set buffer size and create buffer around addresses
buffer_km<-2
address_buffers<-st_buffer(sample_address_pt, dist = buffer_km*1000)

#check location
ggplot()+geom_sf(data=CAmap)+geom_sf(data = address_buffers,color = "red")

#compute mean exposure within each buffer for the corresponding period
tic("Traditional buffer exposure assessment")
buffer_results_multi$dust_buff2 <- mapply(function(buffer, time) {
  
  #Keep only points within buffer
  idx <- st_intersects(buffer, vD_dust_pts, sparse = FALSE)[1,]
  #If no points in buffer, return NA
  if (!any(idx)) return(NA_real_)
  #Compute mean dust exposure
  mean(vD_dust_pts[idx,time, drop=T], na.rm=T)
  
}, st_geometry(address_buffers), address_buffers$period)
toc()
#Conventional buffer exposure assessment: 1683.19 sec elapsed


#check correlation across different buffer sizes
cor(buffer_results_multi[,2:5], use = "complete.obs")

#check correlation between conventional buffer and grid-based approximation by buffer size
cor(buffer_results_multi$dust_buff15, approx_results_multi$dust_buff15, use = "complete.obs")
#0.9999627
cor(buffer_results_multi$dust_buff10, approx_results_multi$dust_buff10, use = "complete.obs")
#0.9999373
cor(buffer_results_multi$dust_buff5, approx_results_multi$dust_buff5, use = "complete.obs")
#0.999801
cor(buffer_results_multi$dust_buff2, approx_results_multi$dust_buff2, use = "complete.obs")
#0.9992373


######################################################################
#(3) Conduct buffer analysis with exactextract r package
#exactextractr need raster exposure data

#read in un-transformed dust data again
vD_dust<-readRDS("vD_dust_CA_wide.RDS")
#make point object and set projection
vD_dust_pts <- vD_dust%>%
  mutate(x=lon)%>%mutate(y=lat)%>%
  st_as_sf(coords = c("lon", "lat"),crs = 4326)

#rasterize one layer dust data
template <- rast(
  xmin = min(vD_dust$lon),
  xmax = max(vD_dust$lon),
  ymin = min(vD_dust$lat),
  ymax = max(vD_dust$lat),
  resolution = c(0.01, 0.01),
  crs = st_crs(vD_dust_pts)$wkt
)
template <- rast(
  ncols = 1042, nrows=1000,
  xmin = -124.415, xmax = -114.005, ymin = 32.005, ymax = 41.995,
  crs = st_crs(vD_dust_pts)$wkt
)

tic("rasterize data")
vD_dust_raster <- rasterize(vD_dust_pts, template, field = "dust2000001", fun = mean)
toc()
#rasterize data: 13.06 sec elapsed

#create buffer polygon
buffer_km<-10
address_buffers<-st_buffer(sample_address_pt, dist = buffer_km*1000)
#transform projection
address_buffers <- st_transform(address_buffers, crs = 4326)

#plot raster and buffer polygon
plot(vD_dust_raster)
plot(st_geometry(address_buffers), add = TRUE, border = "red")

#calculate exposure in buffer with exactextractr package for one layer only 
tic("Exactextractr buffer exposure assessment")
exact_results <- exact_extract(vD_dust_raster, address_buffers, 'mean')
toc()

#Exactextractr buffer exposure assessment: 0.21 sec elapsed


#conduct exactextract method by rasterizing one layer at a time

#set empty dataframe to store conventional buffer results
exact_results_multi<-sample_address%>%
  mutate(dust_buff15=as.numeric(NA))%>%
  mutate(dust_buff10=as.numeric(NA))%>%
  mutate(dust_buff5=as.numeric(NA))%>%
  mutate(dust_buff2=as.numeric(NA))

#get a list of the exposure timepoints
t_list = unique(sample_address$period)


#conduct exactextract (10km buffer)
#for each timepoint
tic("exactextract exposure assessment")
for(t in 1:length(t_list)){
  #rasterize the layer of exposure
  vD_dust_raster <- rasterize(vD_dust_pts, template, field = t_list[t])
  #find addresses that are assessed at that time point
  x = sample_address$period==t_list[t]
  #conduct exact extract and save the results
  exact_results_multi[x,"dust_buff10"]<-exact_extract(vD_dust_raster, address_buffers[x, ], 'mean')
}
toc()
#exactextract exposure assessment: 6375.39 sec elapsed

gc()

#conduct exactextract (2, 5, 10, 15km buffer)
#create buffer polygons
address_buffers15<-st_buffer(sample_address_pt, dist = 15*1000)
address_buffers10<-st_buffer(sample_address_pt, dist = 10*1000)
address_buffers5<-st_buffer(sample_address_pt, dist = 5*1000)
address_buffers2<-st_buffer(sample_address_pt, dist = 2*1000)
#transform projection
address_buffers15 <- st_transform(address_buffers15, crs = 4326)
address_buffers10 <- st_transform(address_buffers10, crs = 4326)
address_buffers5 <- st_transform(address_buffers5, crs = 4326)
address_buffers2 <- st_transform(address_buffers2, crs = 4326)


#for each timepoint
tic("exactextract exposure assessment")
for(t in 1:length(t_list)){
  #rasterize the layer of exposure
  vD_dust_raster <- rasterize(vD_dust_pts, template, field = t_list[t])
  #find addresses that are assessed at that time point
  x = sample_address$period==t_list[t]
  #conduct exact extract and save the results
  exact_results_multi[x,"dust_buff15"]<-exact_extract(vD_dust_raster, address_buffers15[x, ], 'mean')
  exact_results_multi[x,"dust_buff10"]<-exact_extract(vD_dust_raster, address_buffers10[x, ], 'mean')
  exact_results_multi[x,"dust_buff5"]<-exact_extract(vD_dust_raster, address_buffers5[x, ], 'mean')
  exact_results_multi[x,"dust_buff2"]<-exact_extract(vD_dust_raster, address_buffers2[x, ], 'mean')
}
toc()

#exactextract exposure assessment: 6543.89 sec elapsed


#check correlation between conventional buffer and exactextractr by buffer size
cor(buffer_results_multi$dust_buff15, exact_results_multi$dust_buff15, use = "complete.obs")
#0.999963
cor(buffer_results_multi$dust_buff10, exact_results_multi$dust_buff10, use = "complete.obs")
#0.9999536
cor(buffer_results_multi$dust_buff5, exact_results_multi$dust_buff5, use = "complete.obs")
#0.9998623
cor(buffer_results_multi$dust_buff2, exact_results_multi$dust_buff2, use = "complete.obs")
#0.9997873


