# GAVE-buffer-analysis
Implementation of the grid-based approximation and vectorized extraction (GAVE) method for buffer exposure assessment in environmental epidemiology.

The files in this repository are provided to accompany the publication entitled “Short Communication: A Raster-Free Method for Rapid Buffer-Based Exposure Assessment in Environmental Epidemiology” by Wenxin Lu and Alexandra Heaney


R code files:

01.data_preparation.R

Provides the code used to generate intermediate datasets used for the benchmarking test, including buffer_point_list.RDS and sample_address.RDS. Both generated datasets are also provided. 

02.buffer_assessment_test.R

Provides the code used to conduct the benchmarking test comparing the performances of three buffer analysis methods, as described in Section 2.2 of the publication. 
Note that the computational speed of these tests may differ based on the device used. 

03.SOMI_example.R

Provides the code for a real-life application of the grid-based approximation and vectorized extraction (GAVE) method for buffer analysis in the Study of Outcomes in Mothers and Infants cohort, as described in Section 2.3 of the publication.
The SOMI cohort data are not publicly available under a data use agreement with the California Department of Public Health (CDPH). But this code can be adapted for other applications. 



Intermediate datasets:


CA_State_Boundary.gpkg

The California State Boundary, used in the benchmarking tests for 1000 random addressed in California. 

vD_dust_CA_wide.RDS

The dust-specific PM2.5 data (from SatPM2.5 V5.NA.05/V5.NA.05.02, https://sites.wustl.edu/acag/surface-pm2-5/), defined at 0.01° × 0.01° spatial resolution, used in the benchmarking test. 
Variables “lon” and “lat” are the longitude and latitude for the 0.01° × 0.01° grid cells. 
The other columns are named as “dustyyyyddd”, such as dust2000001: the average dust-specific PM2.5 during the biweekly period from the 1st day of 2000 (Jan 1, 2000) to the (1+13)th day of 2000 (Jan 14, 2000).

sample_address.RDS

The 1000 random addresses in California were generated for the benchmarking test.
Variables include: “id”, “lon”, “lat”, and “period”. The “period” variable is in the same format as the dust PM2.5 variable names in vD_dust_CA_wide.RDS.

buffer_point_list.RDS

The relative longitude and latitude offsets for the 0.01° × 0.01° grid points covered in 15-, 10-, 5-, and 2-km circular buffers in California. This dataset was created in 01.data_preparation.R and used in the benchmarking test and SOMI application.
Variables include: 
“lon_rel” and “lat_rel”: relative longitude and latitude offsets for each grid point.
“buff15”, “buff10”, “buff5”, and “buff2”: binary variables flagging grid points covered in 15-, 10-, 5-, and 2-km circular buffers, respectively. 
