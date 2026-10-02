#Generation of intermediate data-sets used for the buffer-assessment benchmarking test

#This code is used to generate the two data-sets:
#buffer_point_list.RDS, and
#sample_address.RDS
#used for the benchmarking test for buffer analysis methods
#Both generated data-sets are provided with the code. 

rm(list=ls())
set.seed(10092025)

library(sf)
library(dplyr)
library(ggplot2)

#################################################
#(1) Generate 1000 random addresses in California
#format: lat and lon with 5 decimal places

#read in california boundary map
CAmap<-read_sf("CA_State_Boundary.gpkg")

#generate 2500 random lats and lons, within the range of California
sample_address<-data.frame(
  lat = round(runif(2500, min=32.32, max=42),5),
  lon = round(runif(2500, min=-124.26, max=-114.8),5)
)

#add x and y variable to back up lon and lat before transformation
sample_address<-sample_address%>%
  mutate(x=lon)%>%
  mutate(y=lat)

#make point object and set projection
sample_address_pt <- st_as_sf(sample_address, 
                            coords = c("lon", "lat"),crs = 4326)

#transform the projection for sample address
sample_address_pt <- st_transform(sample_address_pt, st_crs(CAmap))

#For each sample address, see if it falls in CA boundary
withinCA<-st_intersects(sample_address_pt, CAmap, sparse = F)

#keep the first 1000 random addresses within California
sample_address<-sample_address[withinCA,]
sample_address_pt<-sample_address_pt[withinCA,]
sample_address<-sample_address[1:1000,]
sample_address_pt<-sample_address_pt[1:1000,]

#assign IDs to these random addresses
sample_address$id<-1:1000
sample_address_pt$id<-1:1000

#rearrange variables
sample_address<-sample_address[,c("id","lon","lat")]

rm(withinCA, sample_address_pt)


###############################################################
#(2) Read in van Donkelaar dust data and assign exposure period
vD_dust<-readRDS("vD_dust_CA_wide.RDS")

#make point object and set projection
vD_dust_pts <- vD_dust%>%
  mutate(x=lon)%>%mutate(y=lat)%>%
  st_as_sf(coords = c("lon", "lat"),crs = 4326)

#transform the projection for vD_dust data
vD_dust_pts <- st_transform(vD_dust_pts, st_crs(CAmap))

#This dataset is structured to have one grid point per row and one biweekly period per column
#randomly draw one biweekly period and assign it to the sample addresses
#can be interpreted as birth date or other relevant exposure period
sample_address$period<-colnames(vD_dust)[sample(3:600, 1000, replace = T)]

#save sample address (use this set of address for the tests)
#saveRDS(sample_address, "sample_address.RDS")
sample_address<-readRDS("sample_address.RDS")

#make point object and set projection
sample_address_pt <- sample_address%>%
  mutate(x=lon)%>%mutate(y=lat)%>%
  st_as_sf(coords = c("lon", "lat"),crs = 4326)

#transform the projection for sample address
sample_address_pt <- st_transform(sample_address_pt, st_crs(CAmap))

#check location
ggplot()+geom_sf(data=CAmap)+geom_sf(data = sample_address_pt,color = "red")


###############################################################
#(3) Prepare buffer_point_list dataset
#this contains the relative x- and y- offsets of the grid points
#covered in 2-, 5-, 10-, and 15-km buffers

#First, select a lat and lon address in the center of California as a reference location
center<-data.frame(y_center = 37.16, x_center = -119.53)

#make point object and set projection
center <- st_as_sf(center, coords = c("x_center", "y_center"),crs = 4326)

#transform the projection 
center <- st_transform(center, crs = st_crs(CAmap))

#check location, make sure it's in California and around 37 degree lat
ggplot()+geom_sf(data=CAmap)+geom_sf(data = center,color = "red")

#create a buffer for that approximated center point
buffer_km<-15
center_buffer<-st_buffer(center, dist = buffer_km*1000)

#intersect the buffer area with the vD data, get the points within the buffer
idx <- st_intersects(center_buffer, vD_dust_pts, sparse = FALSE)[1,]
point_list<-data.frame(x=vD_dust_pts[idx,]$x, y=vD_dust_pts[idx,]$y)

#calculate the relative lon and lat of these points compared to the approximated center point
point_list<-point_list%>%
  mutate(lon_rel = x +119.53)%>%
  mutate(lat_rel = y -37.16)%>%
  mutate(buff15 = 1)

#summarize the relative lat and lons (should sum to 0)
summary(point_list$lat_rel)
summary(point_list$lon_rel)


#repeat the process for 10km buffer
buffer_km<-10
center_buffer<-st_buffer(center, dist = buffer_km*1000)
ggplot()+geom_sf(data=CAmap)+geom_sf(data = center_buffer,color = "red")
idx <- st_intersects(center_buffer, vD_dust_pts, sparse = FALSE)[1,]
point_list_2<-data_frame(x=vD_dust_pts[idx,]$x, y=vD_dust_pts[idx,]$y)
point_list_2$buff10<-1
point_list<-left_join(point_list, point_list_2, by=c("x", "y"))

#repeat the process for 5km buffer
buffer_km<-5
center_buffer<-st_buffer(center, dist = buffer_km*1000)
ggplot()+geom_sf(data=CAmap)+geom_sf(data = center_buffer,color = "red")
idx <- st_intersects(center_buffer, vD_dust_pts, sparse = FALSE)[1,]
point_list_3<-data_frame(x=vD_dust_pts[idx,]$x, y=vD_dust_pts[idx,]$y)
point_list_3$buff5<-1
point_list<-left_join(point_list, point_list_3, by=c("x", "y"))

#repeat the process for 2km buffer
buffer_km<-2
center_buffer<-st_buffer(center, dist = buffer_km*1000)
ggplot()+geom_sf(data=CAmap)+geom_sf(data = center_buffer,color = "red")
idx <- st_intersects(center_buffer, vD_dust_pts, sparse = FALSE)[1,]
point_list_4<-data_frame(x=vD_dust_pts[idx,]$x, y=vD_dust_pts[idx,]$y)
point_list_4$buff2<-1
point_list<-left_join(point_list, point_list_4, by=c("x", "y"))

#replace all NA with 0
point_list[is.na(point_list)] <- 0

#chevk the relative lat and lons (should sum to 0)
summary(point_list[point_list$buff15==1,]$lat_rel)
summary(point_list[point_list$buff10==1,]$lat_rel)
summary(point_list[point_list$buff5==1,]$lat_rel)
summary(point_list[point_list$buff2==1,]$lat_rel)

#save dataset for future use
point_list<-point_list[,c("lon_rel","lat_rel","buff15","buff10","buff5","buff2")]
saveRDS(point_list,"buffer_point_list.RDS")



