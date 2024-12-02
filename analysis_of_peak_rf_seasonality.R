

rm(list=ls())
library(ncdf4)
library(terra)
library(lmomco)
library(Kendall)
library(ggplot2)
library(devtools)
library(plotly)
library(circular)
#library(CircStats)
library(trend)


#library(devtools) # To install copula from git hub
#install_github("mmaechler/copula")
library(copula) #install.packages("copula",dependencies = TRUE)

## Function to check leap year
# Function to check for a leap year
is_leap_year <- function(year) {
  if ((year %% 4 == 0 && year %% 100 != 0) || year %% 400 == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}




list_of_imd_data_files<- list.files("C:/Users/Indian_Rainfall_0_pnt_25")

# Find the years of records of IMD
years_of_record<- matrix(nrow = length(list_of_imd_data_files), ncol = 1)
number_of_days_each_year<- matrix(nrow = length(list_of_imd_data_files), ncol = 1)
for (i in 1:length(list_of_imd_data_files)){
  
  temp_list_file<- list_of_imd_data_files[i]
  years_of_record[i]<- as.numeric(substr(temp_list_file, 9,12))
  
  
  if (is_leap_year(as.numeric(substr(temp_list_file, 9,12))) ==T){
    number_of_days_each_year[i,1]<- 366
  }else{
    number_of_days_each_year[i,1]<- 365
  }
  
  
  
}









path_IMD_ls<- paste("C:\\Users\\Indian_Rainfall_0_pnt_25\\",list_of_imd_data_files[1],sep = "")
nc_IMD_ls<- nc_open(path_IMD_ls)
name_var<- attributes(nc_IMD_ls$var)
dim_lon_IMD<- ncvar_get(nc_IMD_ls, attributes(nc_IMD_ls$dim)$names[1])
dim_lat_IMD<- ncvar_get(nc_IMD_ls, attributes(nc_IMD_ls$dim)$names[2])

IMD_data<- ncvar_get(nc_IMD_ls, "RAINFALL")

nc_close(nc_IMD_ls)


if (is_leap_year(as.numeric(substr(path_IMD_ls,43,46))) == T){
  n_days<- 366
}else{
  n_days<- 365
}



## Identify all the 
all_long_serially_IMD<- matrix(nrow = 1, ncol=1)
all_lat_serially_IMD<- matrix(nrow = 1, ncol=1)
for (i in 1:length(dim_lon_IMD)){
  
  for (j in 1:length(dim_lat_IMD)){
    
    temp_IMD_long<- dim_lon_IMD[i]
    temp_IMD_lat<- dim_lat_IMD[j]
    
    IMD_data_one_grid<- IMD_data[i,j,]
    
    length_na_IMD<- length(which(is.na(IMD_data_one_grid)==T))
    
    if (length_na_IMD != n_days){
      all_long_serially_IMD<- rbind(all_long_serially_IMD,temp_IMD_long)
      all_lat_serially_IMD<- rbind(all_lat_serially_IMD,temp_IMD_lat)
      
    }
    
    
  }
  
}


all_long_serially_IMD<- all_long_serially_IMD[-1,]
all_lat_serially_IMD<- all_lat_serially_IMD[-1,]





## Extract the IMD data for the grids without NAs___________________________________

IMD_data_removing_na<- matrix(nrow = 1, ncol = length(all_long_serially_IMD))
for (f in 1:length(list_of_imd_data_files)){
  
  path_IMD<- paste("C:\\Users\\Indian_Rainfall_0_pnt_25\\",list_of_imd_data_files[f],sep = "")
  nc_IMD<- nc_open(path_IMD)
  name_var<- attributes(nc_IMD$var)
  dim_lon_IMD<- ncvar_get(nc_IMD, attributes(nc_IMD$dim)$names[1])
  dim_lat_IMD<- ncvar_get(nc_IMD, attributes(nc_IMD$dim)$names[2])
  
  IMD_data<- ncvar_get(nc_IMD, "RAINFALL")
  nc_close(nc_IMD_ls)
  
  numb_days<- number_of_days_each_year[f,1]
  
  data_IMD_all_grid_for_a_year<- matrix(nrow = numb_days, ncol = length(all_long_serially_IMD))
  for (i in 1:length(all_long_serially_IMD)){
    ind_lon<- which(dim_lon_IMD == all_long_serially_IMD[i])
    ind_lat<- which(dim_lat_IMD == all_lat_serially_IMD[i])
    
    IMD_data_for_a_grid<- IMD_data[ind_lon,ind_lat,]
    
    data_IMD_all_grid_for_a_year[,i]<- IMD_data_for_a_grid
    
    #cat(i)
  }
  
  
  IMD_data_removing_na<- rbind(IMD_data_removing_na,data_IMD_all_grid_for_a_year)
  
  print(paste("Data extracted for year ",years_of_record[f,1],"-----------To be extracted till ",years_of_record[nrow(years_of_record),1],sep = ""))
  
}

IMD_data_removing_na<- IMD_data_removing_na[-1,]










### ____________________________________________________ 7 day total rainfall


accumulated_days_RF<- 7 - 1


## To average over for 7 days as Gunter Bloschl paper did it to capture the data occurring from same event
IMD_data_removing_na_7Days_movngAVG<- matrix(nrow = nrow(IMD_data_removing_na), ncol = ncol(IMD_data_removing_na))
for (i in 1:nrow(IMD_data_removing_na)){
  
  end_ind<- i
  start_ind<- end_ind - accumulated_days_RF
  
  indces_to_be_extracted<- c(start_ind:end_ind)
  
  postiv_ind<- which(indces_to_be_extracted > 0)
  
  ind_pos<- indces_to_be_extracted[postiv_ind]
  if (length(ind_pos)<2){
    data_7day<- matrix(IMD_data_removing_na[ind_pos,], nrow = 1)
  }else{
    data_7day<- IMD_data_removing_na[ind_pos,]
    identify_na<- which(is.na(data_7day[,1])==T)
    if (length(identify_na) > 0){
      data_7day<- data_7day[-identify_na,]
    }
    
  }
  
  av_7day<- apply(data_7day,2,FUN = mean)
  
  IMD_data_removing_na_7Days_movngAVG[i,]<- av_7day
  
  print(paste("Average Done for ",accumulated_days_RF," Days","       ",i,".............To be done till ",nrow(IMD_data_removing_na),sep=""))
}






















## Extraction of dates of peak rainfall_______________________________________________


dates_of_peak_rainfall<- matrix(nrow= nrow(years_of_record), ncol = ncol(IMD_data_removing_na))
max_rf_data<- matrix(nrow= nrow(years_of_record), ncol = ncol(IMD_data_removing_na))
for (i in 1:ncol(dates_of_peak_rainfall)){
  
  temp_data<- IMD_data_removing_na[,i]
  #temp_data<- IMD_data_removing_na_7Days_movngAVG[,i] #For 7 day accmulated rainfall
  
  peak_rf_dates<- matrix(nrow = nrow(years_of_record),ncol = 1)
  peak_rf<- matrix(nrow = nrow(years_of_record),ncol = 1)
  for (j in 1:nrow(years_of_record)){
    
    end_ind<- sum(number_of_days_each_year[1:j,1])
    strt_ind<- end_ind - (number_of_days_each_year[j,1]) + 1
    
    one_year_data<- temp_data[strt_ind:end_ind]
    
    one_year_data<- one_year_data[!is.na(one_year_data)]
    
   
    
    
    
    
    ind_max<- which(one_year_data == max(one_year_data))
    
    if (length(ind_max) > 1){
      ind_max<- ind_max[1]
    }
    
    peak_rf[j,1]<- max(one_year_data)
    
    peak_rf_dates[j,1]<- (ind_max/(number_of_days_each_year[j,1]))*2*pi
  }
  
  
  #plot(circular(one_year_data))
  
  
  dates_of_peak_rainfall[,i]<- peak_rf_dates[,1]
  max_rf_data[,i]<- peak_rf[,1]
  
  print(paste("Peakdate for grid", i, " to extracted till ",ncol(IMD_data_removing_na),sep = ""))
}











## Testing for watson two sample test ______________________________________-
total_length_rec<- nrow(years_of_record)
half_rec<- round(total_length_rec/2,0)

two_sample_test_watson<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
diff_mean_two_sample_test<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
diff_resultant_length_two_sample_test<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
for (i in 1:ncol(dates_of_peak_rainfall)){
  
  first_data<- dates_of_peak_rainfall[1:half_rec,i]
  sec_data<- dates_of_peak_rainfall[(half_rec+1):nrow(dates_of_peak_rainfall),i]
  
  first_mean<- mean.circular(first_data)
  sec_mean<- mean.circular(sec_data)
  
  dif_mean<- abs(first_mean - sec_mean)
  
  if (dif_mean > pi){
    diff_mean_two_sample_test[1,i]<- (2*pi) - dif_mean
  }else{
    diff_mean_two_sample_test[1,i]<- dif_mean
  }
  
  
  
  
  
  
  diff_resultant_length_two_sample_test[1,i]<- abs((rho.circular(first_data)) - (rho.circular(sec_data)))
  
  
  watson_test<- watson.two.test(first_data, sec_data, alpha = 0.05)
  
  watson_output<- capture.output(watson_test)[6]
  
  if (watson_output == "Reject Null Hypothesis "){
    two_sample_test_watson[i]<- 2
    
  }else if (watson_output == "Do Not Reject Null Hypothesis "){
    watson_test<- watson.two.test(first_data, sec_data, alpha = 0.1)
    watson_output<- capture.output(watson_test)[6]
    if (watson_output == "Reject Null Hypothesis "){
      two_sample_test_watson[i]<- 1
    }else{
      two_sample_test_watson[i]<- 0
    }
    

  }
  
  cat(i)
}





library(raster)
###Spatial plot_________________________
library(maptools)
india_watershed_di<- readShapeSpatial("D:/Research_IIT_Roorkee/India_shapefile/india_India_Country_Boundary_MAPOG/india_India_Country_Boundary.shp")

library(mapproj)
#as_map<- map(NW_US_di)

#library(gridExtra)

library(plotrix)

#setwd("G:/Kalai/Circular_statistics/Results/figures")
#tiff('seasonality_in_NW_US_211_sites.tif')
#plot(NW_US_di,lwd=2)

min_long<- extent(india_watershed_di)[1]
max_long<- extent(india_watershed_di)[2]
min_lat<- extent(india_watershed_di)[3]
max_lat<- extent(india_watershed_di)[4]

#gaps_long<- 5
long_at<- c(70,75,80,85,90,95)
lat_at<- c(10,15,20,25,30,35)

png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/watson_2_sample_test.png")


map(india_watershed_di,xlim=c(67.5,98), ylim=c(6,38),lwd=2)
title(main="Watson Two Sample Test")



map.axes(cex.axis=1.2,lwd=3,at=c(long_at,lat_at),
         labels=c('70°E','75°E','80°E','85°E','90°E','95°E','10°N','15°N','20°N','25°N',"30°N",'35°N'))


north(xy=cbind(95,35),type=2)

#north.arrow(95,35,len,lab='NORTH',cex.lab=1,tcol='black')

for (i in 1:length(all_lat_serially_IMD)){
  
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  
  if (two_sample_test_watson[1,i] == 2){
    col_pnt<- "tomato4"
  }else if (two_sample_test_watson[1,i] == 1){
    col_pnt<- "red"
  }else{
    col_pnt<- "grey"
  }
  points(temp_long,temp_lat,col = col_pnt, cex = 1,pch=16)
  
  
}



#i<- c(1051,1681,2821) #Sample i

#points ______marking the sig difference in Mean and resultant length
#if (i == 2821){
i<- 2821 #Significant difference in the mean
temp_lat<- all_lat_serially_IMD[i]
temp_long<- all_long_serially_IMD[i]
points(temp_long,temp_lat,col = "green4", cex = 2,pch=3, lwd=2)
points(87,19,col = "green4", cex = 2,pch=3, lwd=2)
text(93, 19.2,paste(temp_long,"°E;",temp_lat,"°N",sep = ""),cex=1)
#}else if (i == 1681){
i<- 1681  #Significant diff in the resultant
temp_lat<- all_lat_serially_IMD[i]
temp_long<- all_long_serially_IMD[i]
points(temp_long,temp_lat,col = "green4", cex = 2,pch=4, lwd=2)
points(87,17.5,col = "green4", cex = 2,pch=4,lwd=2)
text(93, 17.3,paste(temp_long,"°E;",temp_lat,"°N",sep = ""),cex=1)
#}else if (i == 1051){
i<- 1051
temp_lat<- all_lat_serially_IMD[i]
temp_long<- all_long_serially_IMD[i]
points(temp_long,temp_lat,col = "green4", cex = 2,pch=8, lwd=2)
points(87,15.5,col = "green4", cex = 2,pch=8, lwd=2)
text(93, 15.3,paste(temp_long,"°E;",temp_lat,"°N",sep = ""),cex=1)
#}



# For legend
points(82.5,13,col = "tomato4", cex = 1.5,pch=15)
points(82.5,11,col = "red", cex = 1.5,pch=15)
points(82.5,9,col = "grey", cex = 1.5,pch=15)

text(84.6, 13.2,"0.05",cex=1.2)
text(84.6, 11.2,"0.1",cex=1.2)
text(84.6, 9.2,"> 0.1",cex=1.2)


dev.off()








## _____________________________ linear circular association ________________

#Function to test for J-W-M correlation to be significant

R2xtCorrCoeff <- function(lvar, cvar){
  rxc <- cor(lvar, cos(cvar))
  rxs <- cor(lvar, sin(cvar))
  rcs <- cor(cos(cvar), sin(cvar))
  R2xtVal <- ((rxc*rxc)+(rxs*rxs)-(2*rxc*rxs*rcs))/(1-rcs*rcs)
  return(R2xtVal)
}


R2xtIndTestRand <-function(lvar, cvar, NR){
  R2xtObs <- R2xtCorrCoeff(lvar, cvar)
  nxtrm <- 1
  for(r in 1:NR){
    lvarRand <- sample(lvar)
    R2xtRand <- R2xtCorrCoeff(lvarRand,cvar)
    if (R2xtRand >= R2xtObs){nxtrm <- nxtrm+1 }
  }
  pval <- nxtrm/(NR+1)
  return(c(R2xtObs, pval))
}

sample_siz<- 999

p_val_JWM_test<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
cor_coef_JWM_test<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
for (i in 1:ncol(dates_of_peak_rainfall)){
  
  temp_dates<- dates_of_peak_rainfall[,i]
  
  linear_var<- 1:length(temp_dates)
  
  jwm_cor_p_val<- R2xtIndTestRand(linear_var,temp_dates,sample_siz)
  
  p_val_JWM_test[1,i]<- jwm_cor_p_val[2]
  cor_coef_JWM_test[1,i]<- jwm_cor_p_val[1]
  cat(i)
}



## Plot_______________________
#map(india_watershed_di,xlim=c(min_long,max_long), ylim=c(min_lat,max_lat),lwd=2)

png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/j_w_m.png")


map(india_watershed_di,xlim=c(67.5,98), ylim=c(6,38),lwd=2)
title(main="Johnson-Wehrly-Mardia correlation coefficient")


map.axes(cex.axis=1.2,lwd=3,at=c(long_at,lat_at),
         labels=c('70°E','75°E','80°E','85°E','90°E','95°E','10°N','15°N','20°N','25°N',"30°N",'35°N'))


north(xy=cbind(95,35),type=2)


for (i in 1:length(all_lat_serially_IMD)){
  
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  
  if (p_val_JWM_test[1,i] <= 0.05){
    col_pnt<- "darkslateblue"
  }else if (p_val_JWM_test[1,i] <= 0.1){
    col_pnt<- "blue"
  }else{
    col_pnt<- "grey"
  }
  points(temp_long,temp_lat,col = col_pnt, cex = 1,pch=16)
  
  
}

points(82.5,13,col = "darkslateblue", cex = 1.5,pch=15)
points(82.5,11,col = "blue", cex = 1.5,pch=15)
points(82.5,9,col = "grey", cex = 1.5,pch=15)

text(84.6, 13.2,"0.05",cex=1.2)
text(84.6, 11.2,"0.1",cex=1.2)
text(84.6, 9.2,"> 0.1",cex=1.2)

dev.off()




#Function to test for Mardia correlation to be significant _____________________________________
UniformScores <- function(lvar, cvar){
  ranklvar <- rank(lvar, ties.method="random")
  n <- length(cvar)
  cvar2 <- 0
  for (j in 1:n){
    cvar2[ranklvar[j]] <- cvar[j]
  }
  rankcvar <- rank(cvar2, ties.method="random")
  uscores <- rankcvar*2*pi/n
  return(uscores)
}


Ustar <- function(uniscores){
  
  n <- length(uniscores)
  Tc<- 0
  Ts<- 0
  
  for (j in 1:n){
    Tc <- Tc+j*cos(uniscores[j])
    Ts <- Ts+j*sin(uniscores[j])
  }
  
  Ustar <- (Ts*Ts)+(Tc*Tc)
  return(Ustar)
  
}


MardiaRankIndTestRand <- function(uniscores, NR){
  
  UstarObs <- Ustar(uniscores)
  nxtrm <- 1
  for (r in 1:NR){
    uniscoresRand <- sample(uniscores)
    UstarRand <- Ustar(uniscoresRand)
    if (UstarRand >= UstarObs){nxtrm <- nxtrm + 1 }
  }
  pval <- nxtrm/(NR+1)
  return(c(UstarObs, pval))
  
}




p_val_Mardia_test<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
#cor_coef_Mardia_test<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
for (i in 1:ncol(dates_of_peak_rainfall)){
  
  temp_dates<- dates_of_peak_rainfall[,i]
  
  linear_var<- 1:length(temp_dates)
  
  uscores<- UniformScores(linear_var, temp_dates)
  
  Mardia_cor_p_val<- MardiaRankIndTestRand(uscores,sample_siz)
  
  p_val_Mardia_test[1,i]<- Mardia_cor_p_val[2]
  #cor_coef_JWM_test[1,i]<- jwm_cor_p_val[1]
  cat(i)
}



## Plot_______________________
png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/mardia.png")


map(india_watershed_di,xlim=c(67.5,98), ylim=c(6,38),lwd=2)
title(main="Mardia coefficient")


map.axes(cex.axis=1.2,lwd=3,at=c(long_at,lat_at),
         labels=c('70°E','75°E','80°E','85°E','90°E','95°E','10°N','15°N','20°N','25°N',"30°N",'35°N'))


north(xy=cbind(95,35),type=2)



for (i in 1:length(all_lat_serially_IMD)){
  
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  
  if (p_val_Mardia_test[1,i] <= 0.05){
    col_pnt<- "darkgreen"
  }else if (p_val_Mardia_test[1,i] <= 0.1){
    col_pnt<- "green"
  }else{
    col_pnt<- "grey"
  }
  points(temp_long,temp_lat,col = col_pnt, cex = 1,pch=16)
  
  
}

points(82.5,13,col = "darkgreen", cex = 1.5,pch=15)
points(82.5,11,col = "green", cex = 1.5,pch=15)
points(82.5,9,col = "grey", cex = 1.5,pch=15)

text(84.6, 13.2,"0.05",cex=1.2)
text(84.6, 11.2,"0.1",cex=1.2)
text(84.6, 9.2,"> 0.1",cex=1.2)

dev.off()






## _________________________________________ mean date and Resultant length 
mean_date_peak_rf<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
resultant_date_rf<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
p_value_mk_test_rf<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
p_value_petitt_test_rf<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
petitt_test_change_points<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))


#mat_x_bar<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
#mat_y_bar<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))

rao_1_seas_0_unifo<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
raleigh_p_value<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
for (i in 1:ncol(dates_of_peak_rainfall)){
  
  temp_dates<- dates_of_peak_rainfall[,i]
  
  mean_date_peak_rf[1,i]<- mean.circular(temp_dates, na.rm=FALSE)
  resultant_date_rf[1,i]<- rho.circular(temp_dates, na.rm = FALSE)
  
  mean_plus_pi<- mean.circular(temp_dates, na.rm=FALSE) + pi
  mean_minus_pi<- mean.circular(temp_dates, na.rm=FALSE) - pi
  
  
  
  dates_rescaled<- temp_dates
  
  ind_les_zero<- which(dates_rescaled < 0)
  
  if (length(ind_les_zero) != 0){
    dates_rescaled[ind_les_zero]<- (2*pi) + dates_rescaled[ind_les_zero]
  }
  
  
  
  ####________________Recling
  if (mean_plus_pi > pi){

    
    ind_rescal<- which(temp_dates > mean_plus_pi)
    
    dates_rescaled[ind_rescal]<- dates_rescaled[ind_rescal] - (2*pi)
  }else if (mean_plus_pi > ((3/2)*pi)){
    
    ind_rescal<- which(temp_dates > mean_plus_pi)
    
    dates_rescaled[ind_rescal]<- dates_rescaled[ind_rescal] - (2*pi)
    
  }else if (mean_plus_pi > ((2)*pi)){
    
    ind_rescal<- which(temp_dates < mean_minus_pi)
    dates_rescaled[ind_rescal]<- dates_rescaled[ind_rescal] + (2*pi)
  }else if (mean_minus_pi > (pi/2)){
    
    ind_rescal<- which(temp_dates < mean_minus_pi)
    dates_rescaled[ind_rescal]<- dates_rescaled[ind_rescal] + (2*pi)
  }
  
  
  
  ##________________________ Uniformity test
  
  # Rao spacing
  rao_spac<- rao.spacing.test(temp_dates, alpha = 0.05)
  rao_output<- capture.output(rao_spac)[6]
  
  if (rao_output == "Reject null hypothesis of uniformity "){
    rao_1_seas_0_unifo[1,i]<- 2
    
  }else if (rao_output == "Do not reject null hypothesis of uniformity "){
    rao_spac<- rao.spacing.test(temp_dates, alpha = 0.1)
    rao_output<- capture.output(rao_spac)[6]
    if (rao_output == "Reject null hypothesis of uniformity "){
      rao_1_seas_0_unifo[1,i]<- 1
      
    }else{
      rao_1_seas_0_unifo[1,i]<- 0
    }
    

  }
  
  
  # Rayleigh___________
  rayleg_test<- rayleigh.test(temp_dates, mu=mean.circular(temp_dates))
  raleigh_p_value[1,i]<- rayleg_test$p.value
  
  
  
  
  ##_______________________________ MK and Petit test
  mk_dates<- MannKendall(dates_rescaled)
  petit_dates<- pettitt.test(dates_rescaled)
  
  p_value_mk_test_rf[1,i]<- mk_dates$sl
  p_value_petitt_test_rf[1,i]<- petit_dates$p.value
  petitt_test_change_points[1,i]<- petit_dates$estimate[1]
  
  
  # tem_dat_na_rm<- temp_dates
  # 
  # cos_data<- cos(tem_dat_na_rm)
  # sin_data<- sin(tem_dat_na_rm)
  # x_bar<- (sum(cos_data))/length(tem_dat_na_rm)
  # y_bar<- (sum(sin_data))/length(tem_dat_na_rm)
  # mean_data<- mean(as.circular(tem_dat_na_rm))
  # 
  # 
  # mat_x_bar[i]<- x_bar
  # mat_y_bar[i]<- y_bar
  #mat_mean_date[i]<- mean_data
  
  #mat_resultant_len[i]<- rho.circular(tem_dat_na_rm)
  
  
  cat(i)
  
}


##________________________________________________________Raleigh and Rao test for uniformity


a<- which(rao_1_seas_0_unifo == 0)
b<- which(raleigh_p_value > 0.1)
intersec_rao_and_rayleigh_uniformity<- intersect(a,b)

#i<- 2022

#ind_dnsty_resul<- i

for (i in 1:length(intersec_rao_and_rayleigh_uniformity)){ #activate this for the search
  
  set.seed(100)
  simulated_data <- rvonmises(n=nrow(peak_rf_dates), mu=circular(pi), kappa=0)
  
  ind_dnsty_resul<- intersec_rao_and_rayleigh_uniformity[i]
  
  data<- dates_of_peak_rainfall[,ind_dnsty_resul]
  
  
  data_na_rm_simulated<- simulated_data[!is.na(simulated_data)]
  data_na_rm<- data[!is.na(data)]
  
  #mle_estimates<- mle.vonmises(data_na_rm)
  mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
  mle_estimates_simulated<- mle.vonmises(data_na_rm_simulated)
  
  
  
  #density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
  density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
  density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated$kappa)
  
  
  cos_data<- cos(data_na_rm)
  sin_data<- sin(data_na_rm)
  x_bar<- (sum(cos_data))/length(data_na_rm)
  y_bar<- (sum(sin_data))/length(data_na_rm)
  #plot(as.circular(rho.circular(data_na_rm)))
  
  
  cos_data_simulation<- cos(data_na_rm_simulated)
  sin_data_simulation<- sin(data_na_rm_simulated)
  x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
  y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)
  
  
  ####PLots_____________________________________________________________________________________________
  ##Density and resultant length
  
  
  
  #png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/simulated_unimodal_data_density_leg.png",sep = ""))
  png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/rao_raleigh_uniform/uniform_seasonality",ind_dnsty_resul,".png",sep = ""))
  
  
  op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
  #plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
  plot(density_site, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
  par(new=T)
  plot(density_site_silulated, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 22,lty=1, lwd=3,axes = F,col="red",points.bg="red")
  arrows(0,0,x_bar,y_bar, lwd = 3, lty = 1)
  par(new=T)
  arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "red")
  par(new=T)
  axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
                labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
                modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
  par(new=F)
  #legend("topright",legend = c("Dates of peak RF","Simulated dates from UD","Density of peak RF", "Density of simulated dates"),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","black","red"))
  dev.off()
  
  
}







# Pettitt--------------------------------------------------------------------Plot of petit test 
a<- which(p_value_petitt_test_rf <= 0.05)
b<- which(petitt_test_change_points > 45 & petitt_test_change_points < 75)
#c<- intersect(a,b)

inter_betn_pettt_chgePNT_with_5per_sig<- intersect(a,b)
#i<- c(551, 552,623,624,625,703,704,705,706,785,786,787,871,872,966,967,1066,1159,1160)
i<- 705

#for (i in 1:length(inter_betn_pettt_chgePNT_with_5per_sig)){ #activate this for the search

ind_dnsty_resul<- i # deactivate this for the search
#ind_dnsty_resul<- inter_betn_pettt_chgePNT_with_5per_sig[i] #activate this for the search

change_pnt<- petitt_test_change_points[ind_dnsty_resul]

first_data<- dates_of_peak_rainfall[1:change_pnt,ind_dnsty_resul]
sec_data<- dates_of_peak_rainfall[(change_pnt+1):nrow(dates_of_peak_rainfall),ind_dnsty_resul]

data_na_rm<- first_data[!is.na(first_data)]
data_na_rm_simulated<- sec_data[!is.na(sec_data)]


fst_md<- mean.circular(first_data)
if (fst_md < 0){
  fst_md<- (2*pi) + fst_md
}

fst_md_J<- round(((fst_md/(2*pi))*365),0)
fst_rl<- round(rho.circular(first_data),2)

sec_md<- mean.circular(sec_data)
if (sec_md < 0){
  sec_md<- (2*pi) + sec_md
}

sec_md_J<- round(((sec_md/(2*pi))*365),0)
sec_rl<- round(rho.circular(sec_data),2)




#mle_estimates<- mle.vonmises(data_na_rm)
mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
mle_estimates_simulated<- bw.cv.ml.circular(data_na_rm_simulated,kernel = c("vonmises"))



#density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated)


cos_data<- cos(data_na_rm)
sin_data<- sin(data_na_rm)
x_bar<- (sum(cos_data))/length(data_na_rm)
y_bar<- (sum(sin_data))/length(data_na_rm)
#plot(as.circular(rho.circular(data_na_rm)))


cos_data_simulation<- cos(data_na_rm_simulated)
sin_data_simulation<- sin(data_na_rm_simulated)
x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)




#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/pettitt_test/pettit_two_sample_le_",ind_dnsty_resul,".png",sep = ""))
png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/pettitt_test/pettit_two_sample",ind_dnsty_resul,".png",sep = ""))

op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
#plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
plot(density_site, main = paste("WT=5%,","E=",all_long_serially_IMD[ind_dnsty_resul],",N=",all_lat_serially_IMD[ind_dnsty_resul],
                                #",DM=",round(diff_mean_two_sample_test[1,ind_dnsty_resul],3),",DRL=",round(diff_resultant_length_two_sample_test[1,ind_dnsty_resul],3),
                                ",D_1=",fst_md_J, ",D_2=",sec_md_J,",R_1=",fst_rl, ",R_2=",sec_rl,
                                sep = ""),ticks = F,points.plot=T, xlim=c(-1.1,1.1),ylim = c(-1.2,1.3), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F,col="darkslategray3")
par(new=T)
plot(density_site_silulated, main = paste("WT=5%,","E=",all_long_serially_IMD[ind_dnsty_resul],",N=",all_lat_serially_IMD[ind_dnsty_resul],
                                          #",DM=",round(diff_mean_two_sample_test[1,ind_dnsty_resul],3),",DRL=",round(diff_resultant_length_two_sample_test[1,ind_dnsty_resul],3),
                                          ",D_1=",fst_md_J, ",D_2=",sec_md_J,",R_1=",fst_rl, ",R_2=",sec_rl,
                                          sep = ""),ticks = F,points.plot=T, xlim=c(-1.1,1.1),ylim = c(-1.2,1.3), sub="", xlab="", ylab="", points.pch = 22,lty=1, lwd=3,axes = F,col="darkred")

arrows(0,0,x_bar,y_bar, lwd = 3, lty = 1,col="darkslategray3")
par(new=T)
arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "darkred")
par(new=T)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)
#legend("topright",legend = c(paste("Dates of PRF 1901-",years_of_record[change_pnt,1],sep = ""),paste("Dates of PRF ",(years_of_record[change_pnt,1] + 1),"-2023",sep=""),paste("Density of PRF 1901-",years_of_record[change_pnt,1],sep = ""), paste("Density of PRF ",(years_of_record[change_pnt,1] + 1),"-2023",sep="")),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","darkslategray3","darkred"))
dev.off()

#cat(i) #activate this for the search
#} #activate this for the search
























###__________________________________________________________Plots_______Rayleigh, Rao, MK test, and Pettitt test

png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/mk_test.png")
map(india_watershed_di,xlim=c(67.5,98), ylim=c(6,38),lwd=2)
title(main="MK Test")
map.axes(cex.axis=1.2,lwd=3,at=c(long_at,lat_at),
         labels=c('70°E','75°E','80°E','85°E','90°E','95°E','10°N','15°N','20°N','25°N',"30°N",'35°N'))
north(xy=cbind(95,35),type=2)

for (i in 1:length(all_lat_serially_IMD)){
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  
  if (p_value_mk_test_rf[1,i] <= 0.05){
    col_pnt<- "pink4"
  }else if (p_value_mk_test_rf[1,i] <= 0.1){
    col_pnt<- "pink"
  }else{
    col_pnt<- "grey"
  }
  points(temp_long,temp_lat,col = col_pnt, cex = 1,pch=16)
}


i<- 2724
temp_lat<- all_lat_serially_IMD[i]
temp_long<- all_long_serially_IMD[i]
points(temp_long,temp_lat,col = "black", cex = 2,pch=8, lwd=2)
points(87,15.5,col = "black", cex = 2,pch=8, lwd=2)
text(93, 15.3,paste(temp_long,"°E;",temp_lat,"°N",sep = ""),cex=1)


points(82.5,13,col = "pink4", cex = 1.5,pch=15)
points(82.5,11,col = "pink", cex = 1.5,pch=15)
points(82.5,9,col = "grey", cex = 1.5,pch=15)

text(84.6, 13.2,"0.05",cex=1.2)
text(84.6, 11.2,"0.1",cex=1.2)
text(84.6, 9.2,"> 0.1",cex=1.2)

dev.off()




#____ Pettit
png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/pettitt.png")
map(india_watershed_di,xlim=c(67.5,98), ylim=c(6,38),lwd=2)
#map(india_watershed_di,xlim=c(min_long,max_long), ylim=c(min_lat,max_lat),lwd=2)
title(main="Pettitt Test")
map.axes(cex.axis=1.2,lwd=3,at=c(long_at,lat_at),
         labels=c('70°E','75°E','80°E','85°E','90°E','95°E','10°N','15°N','20°N','25°N',"30°N",'35°N'))
north(xy=cbind(95,35),type=2)

for (i in 1:length(all_lat_serially_IMD)){
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  
  if (p_value_petitt_test_rf[1,i] <= 0.05){
    col_pnt<- "darkcyan"
  }else if (p_value_petitt_test_rf[1,i] <= 0.1){
    col_pnt<- "cyan"
  }else{
    col_pnt<- "grey"
  }
  points(temp_long,temp_lat,col = col_pnt, cex = 1,pch=16)
}



i<- 705
temp_lat<- all_lat_serially_IMD[i]
temp_long<- all_long_serially_IMD[i]
points(temp_long,temp_lat,col = "orange", cex = 2,pch=4, lwd=2)
points(87,15.5,col = "orange", cex = 2,pch=4, lwd=2)
text(93, 15.3,paste(temp_long,"°E;",temp_lat,"°N",sep = ""),cex=1)

# For legends
points(82.5,13,col = "darkcyan", cex = 1.5,pch=15)
points(82.5,11,col = "cyan", cex = 1.5,pch=15)
points(82.5,9,col = "grey", cex = 1.5,pch=15)

text(84.6, 13.2,"0.05",cex=1.2)
text(84.6, 11.2,"0.1",cex=1.2)
text(84.6, 9.2,"> 0.1",cex=1.2)

dev.off()




# Rao__________
png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/rao_spacing_test.png")
map(india_watershed_di,xlim=c(67.5,98), ylim=c(6,38),lwd=2)
#map(india_watershed_di,xlim=c(min_long,max_long), ylim=c(min_lat,max_lat),lwd=2)
title(main="Rao Spacing Test")
map.axes(cex.axis=1.2,lwd=3,at=c(long_at,lat_at),
         labels=c('70°E','75°E','80°E','85°E','90°E','95°E','10°N','15°N','20°N','25°N',"30°N",'35°N'))
north(xy=cbind(95,35),type=2)

for (i in 1:length(all_lat_serially_IMD)){
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  
  if (rao_1_seas_0_unifo[1,i] == 2){
    col_pnt<- "sienna4"
  }else if (rao_1_seas_0_unifo[1,i] == 1){
    col_pnt<- "sienna1"
  }else{
    col_pnt<- "grey"
  }
  points(temp_long,temp_lat,col = col_pnt, cex = 1,pch=16)
}

points(82.5,13,col = "sienna4", cex = 1.5,pch=15)
points(82.5,11,col = "sienna1", cex = 1.5,pch=15)
points(82.5,9,col = "grey", cex = 1.5,pch=15)

text(84.6, 13.2,"0.05",cex=1.2)
text(84.6, 11.2,"0.1",cex=1.2)
text(84.6, 9.2,"> 0.1",cex=1.2)

dev.off()




# Raleigh________________-
png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/rayleigh_uniformity_test.png")
map(india_watershed_di,xlim=c(67.5,98), ylim=c(6,38),lwd=2)
#map(india_watershed_di,xlim=c(min_long,max_long), ylim=c(min_lat,max_lat),lwd=2)
title(main="Rayleigh Uniformity Test")
map.axes(cex.axis=1.2,lwd=3,at=c(long_at,lat_at),
         labels=c('70°E','75°E','80°E','85°E','90°E','95°E','10°N','15°N','20°N','25°N',"30°N",'35°N'))
north(xy=cbind(95,35),type=2)

for (i in 1:length(all_lat_serially_IMD)){
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  
  if (raleigh_p_value[1,i] <= 0.05){
    col_pnt<- "skyblue4"
  }else if (raleigh_p_value[1,i] <= 0.1){
    col_pnt<- "dodgerblue"
  }else{
    col_pnt<- "grey"
  }
  points(temp_long,temp_lat,col = col_pnt, cex = 1,pch=16)
}

points(82.5,13,col = "skyblue4", cex = 1.5,pch=15)
points(82.5,11,col = "dodgerblue", cex = 1.5,pch=15)
points(82.5,9,col = "grey", cex = 1.5,pch=15)

#text(paste(""))
text(84.6, 13.2,"0.05",cex=1.2)
text(84.6, 11.2,"0.1",cex=1.2)
text(84.6, 9.2,"> 0.1",cex=1.2)

dev.off()




## Plot for mean date and resultant__________________________________________________
mean_date_in_julian<- (mean_date_peak_rf/(2*pi))*365

mean_date_in_julian<- round(mean_date_in_julian,0)

library(RColorBrewer)
#display.brewer.pal(n = 12, name = 'RdBu')
names_colours<- brewer.pal(n = 11, name = "Spectral")#Maximum it can give 11 colours

#names_colours_ne<- c(names_colours[1:9],"blue","violet","darkslateblue")
names_colours_ne<- c(names_colours[1:9],"blue","blue4","darkslateblue")

days_month<- c(31,28,31,30,31,30,31,31,30,31,30,31)


start_month<- matrix(nrow = 1, ncol = 1)
end_month<- matrix(nrow = 1, ncol = 1)
for (i in 1:length(days_month)){
  
  en_moth<- sum(days_month[1:i])
  
  end_month<- rbind(end_month,en_moth)
  start_month<- rbind(start_month,(en_moth - (days_month[i] - 1)))
  
}

end_month<- end_month[-1,]
start_month<- start_month[-1,]




# Mean date________________-
png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/mean_date_rf.png")
map(india_watershed_di,xlim=c(67.5,98), ylim=c(6,38),lwd=2)
#map(india_watershed_di,xlim=c(min_long,max_long), ylim=c(min_lat,max_lat),lwd=2)
title(main="Mean Date of peak RF")
map.axes(cex.axis=1.2,lwd=3,at=c(long_at,lat_at),
         labels=c('70°E','75°E','80°E','85°E','90°E','95°E','10°N','15°N','20°N','25°N',"30°N",'35°N'))
north(xy=cbind(95,35),type=2)

names_month<- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")


for (i in 1:ncol(mean_date_in_julian)){
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  
  
  temp_mean_juldate<- mean_date_in_julian[1,i]
  
  if (temp_mean_juldate <= 0){
    temp_mean_juldate<- 365 + temp_mean_juldate
  }
  
  st_date<- which(start_month <= temp_mean_juldate)
  en_date<- which(end_month >= temp_mean_juldate)
  
  if (st_date[length(st_date)] == en_date[1] ){
    ind_col<- en_date[1]
    
    colr<- names_colours_ne[ind_col]
    
    points(temp_long,temp_lat,col = colr, cex = 1,pch=16)
  }else{
    break
  }
  
  
  
}

dif_lat_ind<- 1
lat_st_ind<- 19

for (i in 1:length(names_colours_ne)){
  lat_pnt<- lat_st_ind - (dif_lat_ind*i)
  points(86.5,lat_pnt,col = names_colours_ne[i], cex = 1,pch=15)
  text(88, lat_pnt,names_month[i],cex=0.8)
}

#points(82.5,14,col = "grey", cex = 1,pch=15)

dev.off()







# Resultant Length________________-
#display.brewer.pal(n = 11, name = 'Spectral')
names_colours_rl<- brewer.pal(n = 11, name = "Spectral")
samp_rl<- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)


png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/resultant_length_rf.png")
map(india_watershed_di,xlim=c(67.5,98), ylim=c(6,38),lwd=2)
#map(india_watershed_di,xlim=c(min_long,max_long), ylim=c(min_lat,max_lat),lwd=2)
title(main="Resultant length of dates of peak RF")
map.axes(cex.axis=1.2,lwd=3,at=c(long_at,lat_at),
         labels=c('70°E','75°E','80°E','85°E','90°E','95°E','10°N','15°N','20°N','25°N',"30°N",'35°N'))
north(xy=cbind(95,35),type=2)

for (i in 1:ncol(resultant_date_rf)){
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  
  temp_rl<- round(resultant_date_rf[1,i],1)
  
  ind_samp<- which(samp_rl == temp_rl)
  
  col_rl<- names_colours_rl[ind_samp]
  
  points(temp_long,temp_lat,col = col_rl, cex = 1,pch=16)
  
}


dif_lat_ind<- 1
lat_st_ind<- 18

for (i in 1:length(names_colours_rl)){
  lat_pnt<- lat_st_ind - (dif_lat_ind*i)
  points(86,lat_pnt,col = names_colours_rl[i], cex = 1,pch=15)
  text(87.5, lat_pnt, as.character(samp_rl[i]),cex=0.7)
}

#points(82.5,14,col = "grey", cex = 1,pch=15)

dev.off()









## Find the ____________________________ modality

#Unimodal_plot_________________________________________________________________________________
## Select the common plots with high resultant length
a<- which(all_long_serially_IMD > 72.5 & all_long_serially_IMD < 75)
b<- which(all_lat_serially_IMD > 10 & all_lat_serially_IMD < 15)
intersect(a,b)
resultant_date_rf[1,970]
ind_dnsty_resul<- 970




data<- dates_of_peak_rainfall[,ind_dnsty_resul]
data_na_rm<- data[!is.na(data)]

#mle_estimates<- mle.vonmises(data_na_rm)
mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))

#density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
density_site<- density(as.circular(data_na_rm),bw=mle_estimates)

cos_data<- cos(data_na_rm)
sin_data<- sin(data_na_rm)
x_bar<- (sum(cos_data))/length(data_na_rm)
y_bar<- (sum(sin_data))/length(data_na_rm)
#plot(as.circular(rho.circular(data_na_rm)))

####PLots_____________________________________________________________________________________________
##Density and resultant length


png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/density_resultant_length.png")

op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
plot(density_site, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
arrows(0,0,x_bar,y_bar, lwd = 3)
#mtext(site_number, font = 2)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
legend("topright",legend = c("Dates of peak RF","Density of peak RF"),lty = c(0,1),lwd = c(0,3),pch = c(16,NA))
dev.off()




#library(circular)
#library(CircStats)
#####Rose diagram
png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/density_rose_dia.png")


op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
rose.diag(as.circular(data_na_rm), main = "",bins = 12,axes=FALSE,points.plot=T, xlim=c(-1,1),ylim = c(-1,1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3)
par(new=T)
plot(as.circular(data_na_rm), main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1),ylim = c(-1,1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
#mtext(site_number, font = 2)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)
legend("topright",legend = c("Dates of peak RF"),pch = c(16))

dev.off()






#Bimodal_plot_________________________________________________________________________________
## Select the common plots with high resultant length
a<- which(all_long_serially_IMD > 77.5 & all_long_serially_IMD < 81)
b<- which(all_lat_serially_IMD > 8 & all_lat_serially_IMD < 15)
intersect(a,b)

all_intersec<- intersect(a,b)
resultant_date_rf[1,2271]

#for (i in 1:length(all_intersec)){ ###____________________To be activated for search

  
ind_dnsty_resul<- 2271

#ind_dnsty_resul<- all_intersec[i]


data<- dates_of_peak_rainfall[,ind_dnsty_resul]
data_na_rm<- data[!is.na(data)]

#mle_estimates<- mle.vonmises(data_na_rm)
mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))

#density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
density_site<- density(as.circular(data_na_rm),bw=mle_estimates)


cos_data<- cos(data_na_rm)
sin_data<- sin(data_na_rm)
x_bar<- (sum(cos_data))/length(data_na_rm)
y_bar<- (sum(sin_data))/length(data_na_rm)
#plot(as.circular(rho.circular(data_na_rm)))

####PLots_____________________________________________________________________________________________
##Density and resultant length


#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/bimodal_density_resultant_length_",all_intersec[i],".png",sep = ""))
png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/bimodal_density_resultant_length.png",sep = ""))

op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
#plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
plot(density_site, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.5,1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)

arrows(0,0,x_bar,y_bar, lwd = 3)
#mtext(site_number, font = 2)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
legend("topright",legend = c("Dates of peak RF","Density of peak RF"),lty = c(0,1),lwd = c(0,3),pch = c(16,NA))
dev.off()




#library(circular)
#library(CircStats)
#####Rose diagram
#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/bimodal_density_rose_dia_",all_intersec[i],".png",sep = ""))
png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/bimodal_density_rose_dia.png",sep = ""))


op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
rose.diag(as.circular(data_na_rm), main = "",bins = 12,axes=FALSE,points.plot=T, xlim=c(-1,1),ylim = c(-1,1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3)
par(new=T)
plot(as.circular(data_na_rm), main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1),ylim = c(-1,1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
#plot(as.circular(data_na_rm), main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1),ylim = c(-1,1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
#mtext(site_number, font = 2)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)
legend("topright",legend = c("Dates of peak RF"),pch = c(16))

dev.off()

#cat(i) ###____________________To be activated for search
#} ###____________________To be activated for search








### Identify the number of modes of seasonality of peak rf___________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________-
## __________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________-

# Unimodal__________________________________________________________
set.seed(100)
simulated_data <- rvonmises(n=nrow(peak_rf_dates), mu=circular(pi), kappa=0)

ind_dnsty_resul<- 970

data<- dates_of_peak_rainfall[,ind_dnsty_resul]


data_na_rm_simulated<- simulated_data[!is.na(simulated_data)]
data_na_rm<- data[!is.na(data)]

#mle_estimates<- mle.vonmises(data_na_rm)
mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
mle_estimates_simulated<- mle.vonmises(data_na_rm_simulated)



#density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated$kappa)


cos_data<- cos(data_na_rm)
sin_data<- sin(data_na_rm)
x_bar<- (sum(cos_data))/length(data_na_rm)
y_bar<- (sum(sin_data))/length(data_na_rm)
#plot(as.circular(rho.circular(data_na_rm)))


cos_data_simulation<- cos(data_na_rm_simulated)
sin_data_simulation<- sin(data_na_rm_simulated)
x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)


####PLots_____________________________________________________________________________________________
##Density and resultant length



#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/simulated_unimodal_data_density_leg.png",sep = ""))
png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/simulated_unimodal_data_density",ind_dnsty_resul,".png",sep = ""))


op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
#plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
plot(density_site, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
par(new=T)
plot(density_site_silulated, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 22,lty=1, lwd=3,axes = F,col="red",points.bg="red")
arrows(0,0,x_bar,y_bar, lwd = 3, lty = 1)
par(new=T)
arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "red")
par(new=T)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)
#legend("topright",legend = c("Dates of peak RF","Simulated dates from UD","Density of peak RF", "Density of simulated dates"),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","black","red"))
dev.off()

















#Bimodal_____________________________________________________________
#Generate data from a von Mises distribution
set.seed(100)
simulated_data <- rvonmises(n=nrow(peak_rf_dates), mu=circular(pi), kappa=0)


all_intersec<- intersect(a,b)
resultant_date_rf[1,2271]

#for (i in 1:length(all_intersec)){ ###____________________To be activated for search
ind_dnsty_resul<- 2271

data<- dates_of_peak_rainfall[,ind_dnsty_resul]


data_na_rm_simulated<- simulated_data[!is.na(simulated_data)]
data_na_rm<- data[!is.na(data)]

#mle_estimates<- mle.vonmises(data_na_rm)
mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
mle_estimates_simulated<- mle.vonmises(data_na_rm_simulated)



#density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated$kappa)


cos_data<- cos(data_na_rm)
sin_data<- sin(data_na_rm)
x_bar<- (sum(cos_data))/length(data_na_rm)
y_bar<- (sum(sin_data))/length(data_na_rm)
#plot(as.circular(rho.circular(data_na_rm)))


cos_data_simulation<- cos(data_na_rm_simulated)
sin_data_simulation<- sin(data_na_rm_simulated)
x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)


####PLots_____________________________________________________________________________________________
##Density and resultant length


#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/bimodal_density_resultant_length_",all_intersec[i],".png",sep = ""))
#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/simulated_bimodal_data_density_leg.png",sep = ""))
png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/simulated_bimodal_data_density",ind_dnsty_resul,".png",sep = ""))


op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
#plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
plot(density_site, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
par(new=T)
plot(density_site_silulated, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 22,lty=1, lwd=3,axes = F,col="red",points.bg="red")
arrows(0,0,x_bar,y_bar, lwd = 3, lty = 1)
par(new=T)
arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "red")
par(new=T)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)
#legend("topright",legend = c("Dates of peak RF","Simulated dates from UD","Density of peak RF", "Density of simulated dates"),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","black","red"))
dev.off()





### ____________________________________________________ Another plot of bimodality


#Bimodal_____________________________________________________________
#Generate data from a von Mises distribution
set.seed(100)
simulated_data <- rvonmises(n=nrow(peak_rf_dates), mu=circular(pi), kappa=0)


all_intersec<- intersect(a,b)
resultant_date_rf[1,1051]

#for (i in 1:length(all_intersec)){ ###____________________To be activated for search

ind_dnsty_resul<- 1051
data<- dates_of_peak_rainfall[,ind_dnsty_resul]


data_na_rm_simulated<- simulated_data[!is.na(simulated_data)]
data_na_rm<- data[!is.na(data)]

#mle_estimates<- mle.vonmises(data_na_rm)
mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
mle_estimates_simulated<- mle.vonmises(data_na_rm_simulated)



#density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated$kappa)


cos_data<- cos(data_na_rm)
sin_data<- sin(data_na_rm)
x_bar<- (sum(cos_data))/length(data_na_rm)
y_bar<- (sum(sin_data))/length(data_na_rm)
#plot(as.circular(rho.circular(data_na_rm)))


cos_data_simulation<- cos(data_na_rm_simulated)
sin_data_simulation<- sin(data_na_rm_simulated)
x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)


####PLots_____________________________________________________________________________________________
##Density and resultant length


#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/bimodal_density_resultant_length_",all_intersec[i],".png",sep = ""))
#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/simulated_bimodal_data_density_leg.png",sep = ""))
png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/simulated_bimodal_data_density",ind_dnsty_resul,".png",sep = ""))


op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
#plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
plot(density_site, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
par(new=T)
plot(density_site_silulated, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 22,lty=1, lwd=3,axes = F,col="red",points.bg="red")
arrows(0,0,x_bar,y_bar, lwd = 3, lty = 1)
par(new=T)
arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "red")
par(new=T)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)
#legend("topright",legend = c("Dates of peak RF","Simulated dates from UD","Density of peak RF", "Density of simulated dates"),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","black","red"))
dev.off()







## Bimodalaty reason search__________________________________
a<- which(all_lat_serially_IMD < 35 & all_lat_serially_IMD > 32.5)
b<- which(all_long_serially_IMD < 80 & all_long_serially_IMD > 77.5)

c<- intersect(a,b)


for (i in 1:length(c)){
  set.seed(100)
  simulated_data <- rvonmises(n=nrow(peak_rf_dates), mu=circular(pi), kappa=0)
  
  ind_dnsty_resul<- c[i]
  data<- dates_of_peak_rainfall[,ind_dnsty_resul]
  
  
  data_na_rm_simulated<- simulated_data[!is.na(simulated_data)]
  data_na_rm<- data[!is.na(data)]
  
  #mle_estimates<- mle.vonmises(data_na_rm)
  mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
  mle_estimates_simulated<- mle.vonmises(data_na_rm_simulated)
  
  
  
  #density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
  density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
  density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated$kappa)
  
  
  cos_data<- cos(data_na_rm)
  sin_data<- sin(data_na_rm)
  x_bar<- (sum(cos_data))/length(data_na_rm)
  y_bar<- (sum(sin_data))/length(data_na_rm)
  #plot(as.circular(rho.circular(data_na_rm)))
  
  
  cos_data_simulation<- cos(data_na_rm_simulated)
  sin_data_simulation<- sin(data_na_rm_simulated)
  x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
  y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)
  
  
  ####PLots_____________________________________________________________________________________________
  ##Density and resultant length
  
  
  #png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/bimodal_density_resultant_length_",all_intersec[i],".png",sep = ""))
  #png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/simulated_bimodal_data_density_leg.png",sep = ""))
  png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/insight_inves_bimodal_in_ladak/",ind_dnsty_resul,".png",sep = ""))
  
  
  op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
  #plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
  plot(density_site, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
  par(new=T)
  plot(density_site_silulated, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 22,lty=1, lwd=3,axes = F,col="red",points.bg="red")
  arrows(0,0,x_bar,y_bar, lwd = 3, lty = 1)
  par(new=T)
  arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "red")
  par(new=T)
  axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
                labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
                modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
  par(new=F)
  #legend("topright",legend = c("Dates of peak RF","Simulated dates from UD","Density of peak RF", "Density of simulated dates"),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","black","red"))
  dev.off()
  
  
  
}

# Greater than two modes search_________________________
## 
a<- which(all_lat_serially_IMD < 31 & all_lat_serially_IMD > 30)
b<- which(all_long_serially_IMD < 80 & all_long_serially_IMD > 79)

c<- intersect(a,b)


for (i in 1:length(c)){
  set.seed(100)
  simulated_data <- rvonmises(n=nrow(peak_rf_dates), mu=circular(pi), kappa=0)
  
  ind_dnsty_resul<- c[i]
  data<- dates_of_peak_rainfall[,ind_dnsty_resul]
  
  
  data_na_rm_simulated<- simulated_data[!is.na(simulated_data)]
  data_na_rm<- data[!is.na(data)]
  
  #mle_estimates<- mle.vonmises(data_na_rm)
  mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
  mle_estimates_simulated<- mle.vonmises(data_na_rm_simulated)
  
  
  
  #density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
  density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
  density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated$kappa)
  
  
  cos_data<- cos(data_na_rm)
  sin_data<- sin(data_na_rm)
  x_bar<- (sum(cos_data))/length(data_na_rm)
  y_bar<- (sum(sin_data))/length(data_na_rm)
  #plot(as.circular(rho.circular(data_na_rm)))
  
  
  cos_data_simulation<- cos(data_na_rm_simulated)
  sin_data_simulation<- sin(data_na_rm_simulated)
  x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
  y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)
  
  
  ####PLots_____________________________________________________________________________________________
  ##Density and resultant length
  
  
  #png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/bimodal_density_resultant_length_",all_intersec[i],".png",sep = ""))
  #png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/simulated_bimodal_data_density_leg.png",sep = ""))
  png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/insight_inves_gret_2mode_in_uttrakhand/",ind_dnsty_resul,".png",sep = ""))
  
  
  op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
  #plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
  plot(density_site, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
  par(new=T)
  plot(density_site_silulated, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 22,lty=1, lwd=3,axes = F,col="red",points.bg="red")
  arrows(0,0,x_bar,y_bar, lwd = 3, lty = 1)
  par(new=T)
  arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "red")
  par(new=T)
  axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
                labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
                modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
  par(new=F)
  #legend("topright",legend = c("Dates of peak RF","Simulated dates from UD","Density of peak RF", "Density of simulated dates"),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","black","red"))
  dev.off()
  
  
  
}




## Bimodalaty search in tamil nadu__________________________________
a<- which(all_lat_serially_IMD < 12 & all_lat_serially_IMD > 10)
b<- which(all_long_serially_IMD < 80 & all_long_serially_IMD > 79)

c<- intersect(a,b)


for (i in 1:length(c)){
  set.seed(100)
  simulated_data <- rvonmises(n=nrow(peak_rf_dates), mu=circular(pi), kappa=0)
  
  ind_dnsty_resul<- c[i]
  data<- dates_of_peak_rainfall[,ind_dnsty_resul]
  
  
  data_na_rm_simulated<- simulated_data[!is.na(simulated_data)]
  data_na_rm<- data[!is.na(data)]
  
  #mle_estimates<- mle.vonmises(data_na_rm)
  mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
  mle_estimates_simulated<- mle.vonmises(data_na_rm_simulated)
  
  
  
  #density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
  density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
  density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated$kappa)
  
  
  cos_data<- cos(data_na_rm)
  sin_data<- sin(data_na_rm)
  x_bar<- (sum(cos_data))/length(data_na_rm)
  y_bar<- (sum(sin_data))/length(data_na_rm)
  #plot(as.circular(rho.circular(data_na_rm)))
  
  
  cos_data_simulation<- cos(data_na_rm_simulated)
  sin_data_simulation<- sin(data_na_rm_simulated)
  x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
  y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)
  
  
  ####PLots_____________________________________________________________________________________________
  ##Density and resultant length
  
  
  #png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/bimodal_density_resultant_length_",all_intersec[i],".png",sep = ""))
  #png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/simulated_bimodal_data_density_leg.png",sep = ""))
  png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/insight_inves_bimodal_in_tamil_nadu/",ind_dnsty_resul,".png",sep = ""))
  
  
  op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
  #plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
  plot(density_site, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
  par(new=T)
  plot(density_site_silulated, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 22,lty=1, lwd=3,axes = F,col="red",points.bg="red")
  arrows(0,0,x_bar,y_bar, lwd = 3, lty = 1)
  par(new=T)
  arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "red")
  par(new=T)
  axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
                labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
                modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
  par(new=F)
  #legend("topright",legend = c("Dates of peak RF","Simulated dates from UD","Density of peak RF", "Density of simulated dates"),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","black","red"))
  dev.off()
  
  
  
}

























## Finding actual number of modes__________________________________________ in the whole region
# _____________________________________________________________________________________________

set.seed(100)
simulated_data <- rvonmises(n=nrow(peak_rf_dates), mu=circular(pi), kappa=0)
data_na_rm_simulated<- simulated_data[!is.na(simulated_data)]
mle_estimates_simulated<- mle.vonmises(data_na_rm_simulated)
density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated$kappa)

sim_den_calcu<- density_site_silulated$y


number_of_modes_peak_RF<- matrix(nrow = 1, ncol = ncol(dates_of_peak_rainfall))
for (i in 1:ncol(dates_of_peak_rainfall)){
  
  temp_datesrf<- dates_of_peak_rainfall[,i]
  data_na_rm<- temp_datesrf[!is.na(temp_datesrf)]
  
  mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
  density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
  
  den_calcu<- density_site$y
  
  diff_den_calc<- den_calcu - sim_den_calcu
  
  ind_positive<- which(diff_den_calc > 0)
  
  strt_pnt_positive_diff<- matrix(ind_positive[1],nrow = 1, ncol = 1)
  end_pnt_positive_diff<- matrix(nrow = 1, ncol = 1)
  for (p in 1:(length(ind_positive) - 1)){
    dif_betn_ind<- abs(ind_positive[p+1] - ind_positive[p])
    
    if (dif_betn_ind > 1){
      end_pnt_positive_diff<- rbind(end_pnt_positive_diff,ind_positive[p])
      strt_pnt_positive_diff<- rbind(strt_pnt_positive_diff, ind_positive[p+1])
    }
    
    
    
  }
  end_pnt_positive_diff<- rbind(end_pnt_positive_diff,ind_positive[length(ind_positive)])
  
  
  strt_pnt<- strt_pnt_positive_diff[,1]
  end_pnt<- end_pnt_positive_diff[-1,]
  
  
  
  if (diff_den_calc[1] > 0 & diff_den_calc[length(diff_den_calc)] > 0){
    
    strt_pnt<- strt_pnt[-1]
    end_pnt<- end_pnt[-length(end_pnt)]
    
  }
  
  
  number_of_modes_peak_RF[1,i]<- length(strt_pnt)
  
  cat(i)
}





## Plot for number of modality ____________________




#rep_index<- c(970,1051,2271)


png(filename="D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/number_modality.png")
map(india_watershed_di,xlim=c(67.5,98), ylim=c(6,38),lwd=2)
#map(india_watershed_di,xlim=c(min_long,max_long), ylim=c(min_lat,max_lat),lwd=2)
title(main="Modality of RF Seasonality")
map.axes(cex.axis=1.2,lwd=3,at=c(long_at,lat_at),
         labels=c('70°E','75°E','80°E','85°E','90°E','95°E','10°N','15°N','20°N','25°N',"30°N",'35°N'))
north(xy=cbind(95,35),type=2)

for (i in 1:length(all_lat_serially_IMD)){
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  
  if (number_of_modes_peak_RF[1,i] == 1){
    col_pnt<- "seashell4"
  }else if (number_of_modes_peak_RF[1,i] == 2){
    col_pnt<- "yellow3"
  }else{
    col_pnt<- "darkgreen"
  }
  points(temp_long,temp_lat,col = col_pnt, cex = 1,pch=16)
  
  
  

  
  
}

points(82.5,13,col = "seashell4", cex = 1.5,pch=15)
points(82.5,11,col = "yellow3", cex = 1.5,pch=15)
points(82.5,9,col = "darkgreen", cex = 1.5,pch=15)

text(86.2, 13.2,"Unimodal",cex=1)
text(86.2, 11.2,"Bimodal",cex=1)
text(86.2, 9.2,">2 Modes",cex=1)



#points ______marking the unimodal and bimodal
#if (i == 970){
  i<- 970
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  points(temp_long,temp_lat,col = "red", cex = 2,pch=3, lwd=2)
  points(87,19,col = "red", cex = 2,pch=3, lwd=2)
  text(93, 19.2,paste(temp_long,"°E;",temp_lat,"°N",sep = ""),cex=1)
#}else if (i == 1051){
  i<- 1051
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  points(temp_long,temp_lat,col = "red", cex = 2,pch=4, lwd=2)
  points(87,17.5,col = "red", cex = 2,pch=4,lwd=2)
  text(93, 17.3,paste(temp_long,"°E;",temp_lat,"°N",sep = ""),cex=1)
#}else if (i ==2271){
  i<- 2271
  temp_lat<- all_lat_serially_IMD[i]
  temp_long<- all_long_serially_IMD[i]
  points(temp_long,temp_lat,col = "red", cex = 2,pch=8, lwd=2)
  points(87,15.5,col = "red", cex = 2,pch=8, lwd=2)
  text(93, 15.3,paste(temp_long,"°E;",temp_lat,"°N",sep = ""),cex=1)
#}


dev.off()














## ____________________________________________ Insight investigation

a<- which(all_long_serially_IMD > 77.5 & all_long_serially_IMD < 81)
b<- which(all_lat_serially_IMD > 30 & all_lat_serially_IMD < 32)
intersect(a,b)

set.seed(100)
simulated_data <- rvonmises(n=nrow(peak_rf_dates), mu=circular(pi), kappa=0)


all_intersec<- intersect(a,b)
#resultant_date_rf[1,2271]

for (i in 1:length(all_intersec)){ ###____________________To be activated for search

#ind_dnsty_resul<- 2271
ind_dnsty_resul<- all_intersec[i]

data<- dates_of_peak_rainfall[,ind_dnsty_resul]


data_na_rm_simulated<- simulated_data[!is.na(simulated_data)]
data_na_rm<- data[!is.na(data)]

#mle_estimates<- mle.vonmises(data_na_rm)
mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
mle_estimates_simulated<- mle.vonmises(data_na_rm_simulated)



#density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated$kappa)


cos_data<- cos(data_na_rm)
sin_data<- sin(data_na_rm)
x_bar<- (sum(cos_data))/length(data_na_rm)
y_bar<- (sum(sin_data))/length(data_na_rm)
#plot(as.circular(rho.circular(data_na_rm)))


cos_data_simulation<- cos(data_na_rm_simulated)
sin_data_simulation<- sin(data_na_rm_simulated)
x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)


####PLots_____________________________________________________________________________________________
##Density and resultant length


png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/insight_investigation/bimodal_density_resultant_length_",all_intersec[i],".png",sep = ""))
#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/simulated_bimodal_data_density_leg.png",sep = ""))
#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/insight_investigation/simulated_bimodal_data_density.png",sep = ""))


op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
#plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
plot(density_site, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
par(new=T)
plot(density_site_silulated, main = paste("Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=F, xlim=c(-1,1.2),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 22,lty=1, lwd=3,axes = F,col="red",points.bg="red")
arrows(0,0,x_bar,y_bar, lwd = 3, lty = 1)
par(new=T)
arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "red")
par(new=T)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)
#legend("topright",legend = c("Dates of peak RF","Simulated dates from UD","Density of peak RF", "Density of simulated dates"),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","black","red"))
dev.off()



cat(i)
}











## Plot for the site with significant diff in mean and resultant data ___ 

#_____________________________________________________ Common with Watson test
a<- which(diff_resultant_length_two_sample_test > 0.5)
#two_sample_test_watson[1,a]

b<- which(diff_mean_two_sample_test > 0.5 & diff_mean_two_sample_test < 1.1)
#two_sample_test_watson[1,b]

c<- which(two_sample_test_watson == 2)

intersect_change_mean_and_resultant<- intersect(a,b)
inter_betn_smal_mean_with_5per_sig<- intersect(b,c)

ind_max_diff_in_rl<- which(diff_resultant_length_two_sample_test == max(diff_resultant_length_two_sample_test))
ind_max_dif_in_mean_date<- which(diff_mean_two_sample_test == max(diff_mean_two_sample_test))


# I choose ind #4819 to demonstrate the difference of resultant length
# And choose ind #1051 to demonstrate the difference of resultant length

## Plot for the difference in the mean____________________________________________________--
#i<- 1051 # Good example for bimodal, Also bimodality is significant here
#i<- 2717 # or 2821 Good example for uniimodal

#i<- 1051
i<- 2821


#for (i in 1:length(inter_betn_smal_mean_with_5per_sig)){ #activate this for the search

ind_dnsty_resul<- i # deactivate this for the search
#ind_dnsty_resul<- inter_betn_smal_mean_with_5per_sig[i] #activate this for the search
  
first_data<- dates_of_peak_rainfall[1:half_rec,ind_dnsty_resul]
sec_data<- dates_of_peak_rainfall[(half_rec+1):nrow(dates_of_peak_rainfall),ind_dnsty_resul]

data_na_rm<- first_data[!is.na(first_data)]
data_na_rm_simulated<- sec_data[!is.na(sec_data)]





fst_md<- mean.circular(first_data)
if (fst_md < 0){
  fst_md<- (2*pi) + fst_md
}

fst_md_J<- round(((fst_md/(2*pi))*365),0)
fst_rl<- round(rho.circular(first_data),2)

sec_md<- mean.circular(sec_data)
if (sec_md < 0){
  sec_md<- (2*pi) + sec_md
}

sec_md_J<- round(((sec_md/(2*pi))*365),0)
sec_rl<- round(rho.circular(sec_data),2)





#mle_estimates<- mle.vonmises(data_na_rm)
mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
mle_estimates_simulated<- bw.cv.ml.circular(data_na_rm_simulated,kernel = c("vonmises"))



#density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated)


cos_data<- cos(data_na_rm)
sin_data<- sin(data_na_rm)
x_bar<- (sum(cos_data))/length(data_na_rm)
y_bar<- (sum(sin_data))/length(data_na_rm)
#plot(as.circular(rho.circular(data_na_rm)))


cos_data_simulation<- cos(data_na_rm_simulated)
sin_data_simulation<- sin(data_na_rm_simulated)
x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)




#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/watson_two_sample/change_in_MD_density_le_",ind_dnsty_resul,".png",sep = ""))
png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/watson_two_sample/change_in_MD_density_",ind_dnsty_resul,".png",sep = ""))

op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
#plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
plot(density_site, main = paste("WT=5%,","E=",all_long_serially_IMD[ind_dnsty_resul],",N=",all_lat_serially_IMD[ind_dnsty_resul],
                                #",DM=",round(diff_mean_two_sample_test[1,ind_dnsty_resul],3),",DRL=",round(diff_resultant_length_two_sample_test[1,ind_dnsty_resul],3),
                                ",D_1=",fst_md_J, ",D_2=",sec_md_J,",R_1=",fst_rl, ",R_2=",sec_rl,
                                sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.5,1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F,col="darkslategray3")
par(new=T)
plot(density_site_silulated, main = paste("WT=5%,","E=",all_long_serially_IMD[ind_dnsty_resul],",N=",all_lat_serially_IMD[ind_dnsty_resul],
                                          #",DM=",round(diff_mean_two_sample_test[1,ind_dnsty_resul],3),",DRL=",round(diff_resultant_length_two_sample_test[1,ind_dnsty_resul],3),
                                          ",D_1=",fst_md_J, ",D_2=",sec_md_J,",R_1=",fst_rl, ",R_2=",sec_rl,
                                          sep = ""),ticks = F,points.plot=T, xlim=c(-1,1.2),ylim = c(-1.5,1), sub="", xlab="", ylab="", points.pch = 22,lty=1, lwd=3,axes = F,col="darkred")

arrows(0,0,x_bar,y_bar, lwd = 3, lty = 1,col="darkslategray3")
par(new=T)
arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "darkred")
par(new=T)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)
#legend("topright",legend = c("Dates of PRF 1901-1962","Dates of PRF 1963-2023","Density of PRF 1901-1962", "Density of PRF 1963-2023"),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","darkslategray3","darkred"))
dev.off()

#cat(i) #activate this for the search
#} #activate this for the search












## Plot for the difference in the resultant length ____________________________________________________--
a<- which(diff_resultant_length_two_sample_test > 0.5 & diff_resultant_length_two_sample_test < 0.75)
#two_sample_test_watson[1,a]

b<- which(diff_mean_two_sample_test > 0.5 & diff_mean_two_sample_test < 1.1)
#two_sample_test_watson[1,b]

c<- which(two_sample_test_watson == 2)


inter_betn_smal_rl_with_5per_sig<- intersect(a,c)

#i<- 4819 # Not appropriate figure
i<- 1681

#for (i in 1:length(inter_betn_smal_rl_with_5per_sig)){ #activate this for the search

ind_dnsty_resul<- i # deactivate this for the search
#ind_dnsty_resul<- inter_betn_smal_rl_with_5per_sig[i] #activate this for the search


first_data<- dates_of_peak_rainfall[1:half_rec,ind_dnsty_resul]
sec_data<- dates_of_peak_rainfall[(half_rec+1):nrow(dates_of_peak_rainfall),ind_dnsty_resul]

data_na_rm<- first_data[!is.na(first_data)]
data_na_rm_simulated<- sec_data[!is.na(sec_data)]



fst_md<- mean.circular(first_data)
if (fst_md < 0){
  fst_md<- (2*pi) + fst_md
}

fst_md_J<- round(((fst_md/(2*pi))*365),0)
fst_rl<- round(rho.circular(first_data),2)

sec_md<- mean.circular(sec_data)
if (sec_md < 0){
  sec_md<- (2*pi) + sec_md
}

sec_md_J<- round(((sec_md/(2*pi))*365),0)
sec_rl<- round(rho.circular(sec_data),2)




#mle_estimates<- mle.vonmises(data_na_rm)
mle_estimates<- bw.cv.ml.circular(data_na_rm,kernel = c("vonmises"))
mle_estimates_simulated<- bw.cv.ml.circular(data_na_rm_simulated,kernel = c("vonmises"))



#density_site<- density(as.circular(data_na_rm),bw=mle_estimates$kappa)
density_site<- density(as.circular(data_na_rm),bw=mle_estimates)
density_site_silulated<- density(as.circular(data_na_rm_simulated),bw=mle_estimates_simulated)


cos_data<- cos(data_na_rm)
sin_data<- sin(data_na_rm)
x_bar<- (sum(cos_data))/length(data_na_rm)
y_bar<- (sum(sin_data))/length(data_na_rm)
#plot(as.circular(rho.circular(data_na_rm)))


cos_data_simulation<- cos(data_na_rm_simulated)
sin_data_simulation<- sin(data_na_rm_simulated)
x_bar_simulation<- (sum(cos_data_simulation))/length(data_na_rm_simulated)
y_bar_simulation<- (sum(sin_data_simulation))/length(data_na_rm_simulated)




#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/watson_two_sample/change_in_MD_density_le_",ind_dnsty_resul,".png",sep = ""))
png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/circular_plot/watson_two_sample/change_in_RL_density_",ind_dnsty_resul,".png",sep = ""))

op <- par(mfrow=c(1,1),mar = c(5,0,1,3) + 0.4,oma=c(0, 0, 4, 0))
#plot(density_site, main = paste("ind=",all_intersec[i],", Longitude = ",all_long_serially_IMD[ind_dnsty_resul],", Latitude = ",all_lat_serially_IMD[ind_dnsty_resul],sep = ""),ticks = F,points.plot=T, xlim=c(-1.2,1),ylim = c(-1,1.5), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F)
plot(density_site, main = paste("WT=5%,","E=",all_long_serially_IMD[ind_dnsty_resul],",N=",all_lat_serially_IMD[ind_dnsty_resul],
                                #",DM=",round(diff_mean_two_sample_test[1,ind_dnsty_resul],3),",DRL=",round(diff_resultant_length_two_sample_test[1,ind_dnsty_resul],3),
                                ",D_1=",fst_md_J, ",D_2=",sec_md_J,",R_1=",fst_rl, ",R_2=",sec_rl,
                                sep = ""),ticks = F,points.plot=T, xlim=c(-1.1,1.1),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 16,lty=1, lwd=3,axes = F,col="darkslategray3")
par(new=T)
plot(density_site_silulated, main = paste("WT=5%,","E=",all_long_serially_IMD[ind_dnsty_resul],",N=",all_lat_serially_IMD[ind_dnsty_resul],
                                          #",DM=",round(diff_mean_two_sample_test[1,ind_dnsty_resul],3),",DRL=",round(diff_resultant_length_two_sample_test[1,ind_dnsty_resul],3),
                                          ",D_1=",fst_md_J, ",D_2=",sec_md_J,",R_1=",fst_rl, ",R_2=",sec_rl,
                                          sep = ""),ticks = F,points.plot=T, xlim=c(-1.1,1.1),ylim = c(-1.4,1.1), sub="", xlab="", ylab="", points.pch = 22,lty=1, lwd=3,axes = F,col="darkred")

arrows(0,0,x_bar,y_bar, lwd = 3, lty = 1,col="darkslategray3")
par(new=T)
arrows(0,0,x_bar_simulation,y_bar_simulation, lwd = 3, lty = 1,col = "darkred")
par(new=T)
axis.circular(at=c(0,((64*pi)/365),((120*pi)/365),((182*pi)/365),((242*pi)/365),((304*pi)/365),((364*pi)/365),((426*pi)/365),((488*pi)/365),((548*pi)/365),((610*pi)/365),((670*pi)/365)),
              labels=c('1-Jan','1-Feb','1-Mar','1-Apr','1-May','1-Jun','1-Jul','1-Aug','1-Sep','1-Oct','1-Nov','1-Dec'), units = NULL, template=NULL,
              modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE,tcl=0.09,cex = 0.9, lwd = 3)
par(new=F)
#legend("topright",legend = c("Dates of PRF 1901-1962","Dates of PRF 1963-2023","Density of PRF 1901-1962", "Density of PRF 1963-2023"),lty = c(0,0,1,1),lwd = c(0,0,3,3),pch = c(16,22,NA,NA),col = c("black","black","darkslategray3","darkred"))
dev.off()



#cat(i) #activate this for the search
#} #activate this for the search









## Plot for nonstationarity____________________ for a common site in JWM and Mardia
a<- which(p_val_JWM_test <= 0.05)
b<- which(p_val_Mardia_test <= 0.05)
d<- which(p_value_mk_test_rf <= 0.05)
c<- intersect(a,b)

interset_JWM_N_M_with_MK<- intersect(c,d)

#i<- c(2724,2725,3150,4069,4077,4864,4935,4936,4937,4938) #all of these works, There are few more but not mentioned

#i<- 4077 #good for MK test

#for (j in 1:length(interset_JWM_N_M_with_MK)){ #activate for a loop
i<- 2724
#i<- interset_JWM_N_M_with_MK[j]#activate for a loop

dates_rf<- dates_of_peak_rainfall[,i]
time_rf<- years_of_record[,1]

val_31apr<- (((31+28+31+30)/365)*2*pi)
ind_less_1may<- which(dates_rf <= val_31apr)
dates_rf[ind_less_1may]<- dates_rf[ind_less_1may] + (2*pi)

points_in_Y_axis<- ((c(start_month[4:12],(365+start_month[1:4])))/365)*2*pi

lm_model<- lm(dates_rf ~ time_rf)
b_0<- round(lm_model$coefficients[1],3)
b_1<- round(lm_model$coefficients[2],3)
shift_in_seasonality_in_day<- ((lm_model$coefficients[2])/(2*pi))*365
shift_in_seasonality_in_hr<- round(((lm_model$coefficients[2])/(2*pi))*365*24,0)

start<- time_rf[1]
en<- time_rf[length(time_rf)]
low_lim<- (((31+28+31+30)/365)*2*pi)
upp_lim<- (2*pi) + low_lim
png(paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/trend_plot_in_JWM_&_M_with_MK",i,".png",sep = ""))
op <- par(mar = c(7,6,4,2) + 0.4)
plot(time_rf, dates_rf,main =paste("@5%JWM&M&MK",",B0=",b_0,",B1=",b_1,",shift=",shift_in_seasonality_in_hr,"hr/yr",sep=""), type = "p", lwd=3, xlim = c(start,en), ylim = c(low_lim,upp_lim), xlab =  "Year", ylab = "", col="blue",cex.lab=1.5,cex.axis=1.5,cex.main=1.5,bty="n",xaxt="n",yaxt="n",pch=16)
axis(1, lwd = 2,cex.lab=1.5,cex.axis=1.5,at = seq(start, en+3,by=25))
axis(2, las=1, lwd = 2,cex.lab=1.5,cex.axis=1.5, at = c(points_in_Y_axis),labels=c("1-Apr","1-May","1-Jun","1-Jul","1-Aug","1-Sep","1-Oct","1-Nov","1-Dec","1-Jan","1-Feb","1-Mar","1-Apr"))
par(new=T)
abline(lm(dates_rf ~ time_rf), lwd=2)
legend("topright", lty = c(0, 1), lwd = c(NA,2), pch = c(16,NA), legend = c("Dates of peak RF","Fitted trend for peak RF"), col = c("blue", "Black"), text.font = 2, cex = 1.1)
par(new=F)
dev.off()



# For the loop search.... next lines
#png(filename=paste("D:/Research_IIT_Roorkee/Seasonality_india/pr/Plots/trend_JWM_and_Mardia",i,".png",sep=""))
#plot(time_rf,dates_rf)
#dev.off()

#cat(j)#activate for a loop
#}#activate for a loop







