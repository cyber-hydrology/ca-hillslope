#CA_cal <- function(flowdir, time_step, time_interval, dem, hydro_data, SoilText_Mat,
#                         resolution_cell,LAI, P_max, i_init, i_final, i_decay, Evptr,
#                         manning_coeff_sur, K_Sur, Depth_Sur, E_Sat_Sur_cons_val) {
setwd('C:/Users/82105/Documents/cellularAutomata/CA_model_hillslope_220612')
library(raster);library(tictoc)
source('Neighbor_mat.R');source('STRG_surflow2.R');source('Transition_fun2.R')
##################################################################################################################################
# Cellular Automata Setting ######################################################################################################
##################################################################################################################################
##### Time parameter setting for CA calculation
#flowdir <- c("D4") # D8, D4, 4+4N
#time_step <- c(30) #seconds
#time_interval <- c(2*3600) #seconds
flowdir <- c("D4") # D8, D4, 4+4N
time_step <- c(30) #seconds
time_interval <- c(2*3600) #seconds
hour_to_sec <- c(3600)

##### Matrix setting for CA calculation ##########################################################################################
#dem <- read.table('depth_revised.txt', head=T)
#resolution_cell <- c(500) #[mm]; 
dem <- read.table('depth_revised.txt', head=T)
resolution_cell <- c(500) #[mm]; 
m_to_mm <- c(1000)
cm_to_mm <- c(10)

dem$surface <- dem$surface-min(dem$surface);
dem$surface <- dem$surface*m_to_mm
num_col <- (max(dem$X)-min(dem$X))*2+1
num_row <- (max(dem$Y)-min(dem$Y))*2+1
sur_dem <- matrix(0, num_row, num_col); # [m]
rownames(sur_dem) <- seq(from=max(dem$Y), to=min(dem$Y), c(-resolution_cell/m_to_mm)); 
colnames(sur_dem) <- seq(from=min(dem$X), to=max(dem$X), c(resolution_cell/m_to_mm))
sur_waterdepth_dem <- matrix(0, num_row, num_col) # [m]
rownames(sur_waterdepth_dem) <- seq(from=max(dem$Y), to=min(dem$Y), c(-resolution_cell/m_to_mm)); 
colnames(sur_waterdepth_dem) <- seq(from=min(dem$X), to=max(dem$X), c(resolution_cell/m_to_mm))
sur_cellheight_dem <- matrix(0, num_row, num_col) # [m]
rownames(sur_cellheight_dem) <- seq(from=max(dem$Y), to=min(dem$Y), c(-resolution_cell/m_to_mm)); 
colnames(sur_cellheight_dem) <- seq(from=min(dem$X), to=max(dem$X), c(resolution_cell/m_to_mm))
sur_updatecell_dem <- matrix(0, num_row, num_col) # [m]
rownames(sur_updatecell_dem) <- seq(from=max(dem$Y), to=min(dem$Y), c(-resolution_cell/m_to_mm)); 
colnames(sur_updatecell_dem) <- seq(from=min(dem$X), to=max(dem$X), c(resolution_cell/m_to_mm))
sur_storage_dem    <- matrix(0, num_row, num_col) # [m]
rownames(sur_storage_dem) <- seq(from=max(dem$Y), to=min(dem$Y), c(-resolution_cell/m_to_mm)); 
colnames(sur_storage_dem) <- seq(from=min(dem$X), to=max(dem$X), c(resolution_cell/m_to_mm))
theta_Res_Sur_dem    <- matrix(0, num_row, num_col) # [m]
rownames(theta_Res_Sur_dem) <- seq(from=max(dem$Y), to=min(dem$Y), c(-resolution_cell/m_to_mm)); 
colnames(theta_Res_Sur_dem) <- seq(from=min(dem$X), to=max(dem$X), c(resolution_cell/m_to_mm))
theta_Sat_Sur_dem    <- matrix(0, num_row, num_col) # [m]
rownames(theta_Sat_Sur_dem) <- seq(from=max(dem$Y), to=min(dem$Y), c(-resolution_cell/m_to_mm)); 
colnames(theta_Sat_Sur_dem) <- seq(from=min(dem$X), to=max(dem$X), c(resolution_cell/m_to_mm))

##### Parameters in Waterbalance setting #########################################################################################
#Date_data <- c(120928)
#hydro_data <- read.csv(paste0('runoffEvent_',Date_data,'.csv'))
#manning_coeff_sur <- c(0.035);
#LAI <- c(3.76)#Leaf Area Index
#P_max <- 0.935+0.498*LAI+0.00575*LAI^2#P_max[mm]
#i_init  <- c(2.30)
#i_final <- c(0.49)
#i_decay <- c(1.58)#i_init[mm/hr], i_final[mm/hr]
#Evptr <- c(0)
Date_data <- c(120928)
hydro_data <- read.csv(paste0('runoffEvent_',Date_data,'.csv'))
manning_coeff_sur <- c(0.035);
LAI <- c(3.76)#Leaf Area Index
P_max <- 0.935+0.498*LAI+0.00575*LAI^2#P_max[mm]
i_init  <- c(2.30)
i_final <- c(0.49)
# 220530 c(1.58)
i_decay =  c(0.063) # K, i_init[mm/hr], i_final[mm/hr]
Evptr <- c(0)

R_cum <- c(0); P_pre <- c(0) # [mm]
Intercep_Total <- c(0); Infiltra_Total <- c(0); Runoff_Total <- c(0)

##### Setting for soil texture characteristics ###################################################################################
#SoilText_Mat <- read.csv("GL_SoilText.csv")
SoilText_Mat <- read.csv("GL_SoilText.csv")

SoilText_Mat <- cbind(SoilText_Mat, t(hydro_data[1,grep("10", colnames(hydro_data))]),
                      t(hydro_data[1,grep("30", colnames(hydro_data))]))
dem_soiltext <- matrix(NA, nrow(dem), 9)
colnames(dem_soiltext) <- c(colnames(SoilText_Mat)[2:8],"ASM_10","ASM_30")
dem_soiltext[,1] <- dem[,1]; dem_soiltext[,2] <- dem[,2]
order_dist <- sapply(1:nrow(dem), 
                     function(i) which.min(sapply(1:nrow(SoilText_Mat),
                                                  function(j) sqrt((dem[i,1]-SoilText_Mat[j,2])^2 + (dem[i,2]-SoilText_Mat[j,3])^2))))
for(k in 3:9){
  dem_soiltext[,k] <- sapply(1:nrow(dem), function(i) SoilText_Mat[order_dist[i],k+1])
}
#dem_soiltext[,8] <- rep(mean(SoilText_Mat[,9]), nrow(dem)) # Average of surface sm setting
#dem_soiltext[,9] <- rep(mean(SoilText_Mat[,10]), nrow(dem)) # Average of subsurface sm setting
dem_soiltext <- as.data.frame(dem_soiltext)

##### Parameters in STRG_surflow  ################################################################################################
#K_Sur=c(0.85);  
#Depth_Sur=c(100); 
#E_Sat_Sur_cons_val <- c(60) # % of effective saturation constant

K_Sur=c(0.30);  
Depth_Sur=c(100); 
E_Sat_Sur_cons_val <- c(10) # % of effective saturation constant

E_Sat_Sur=c(0);
E_Sat_Sur_cons = rep(E_Sat_Sur_cons_val, num_row); 
# dir.create(paste0(getwd(),"/Plot/Plots_",flowdir,time_step,"sec_0",K_Sur*100,'K',E_Sat_Sur_cons_val,'%E_Sat_hillslope_test'))
# setwd(paste0(getwd(),"/Plot/Plots_",flowdir,time_step,"sec_0",K_Sur*100,'K',E_Sat_Sur_cons_val,'%E_Sat_hillslop_test'))
# unlink(getwd(), recursive = TRUE)
##### Setting for initial conditions of soil texture, water balance ##############################################################
k <- c(1)
for(j in num_row:1){
  for(i in 1:num_col){
    sur_dem[j,i] <- dem$surface[k]
    sur_storage_dem[j,i] <- dem_soiltext$ASM_10[k]*Depth_Sur/cm_to_mm
    theta_Sat_Sur_dem[j,i] <- dem_soiltext$stsur[k]
    theta_Res_Sur_dem[j,i] <- dem_soiltext$rtsur[k]
    k <- k+1
  }}

x11()
layout(matrix(c(1,2), 1, 2, byrow=T))
rast <- raster(sur_dem)
extent(rast) <- c(c(min(dem$X)+0.5),max(dem$X),c(min(dem$Y)+0.5),max(dem$Y))
projection(rast) <- CRS("+proj=longlat +datum=WGS84")
plot(rast, col=terrain.colors(30), main="sur_dem")
rast <- raster(sur_storage_dem)
extent(rast) <- c(c(min(dem$X)+0.5),max(dem$X),c(min(dem$Y)+0.5),max(dem$Y))
projection(rast) <- CRS("+proj=longlat +datum=WGS84")
plot(rast, col=terrain.colors(30), main="sur_soil water storage")

##### Setting for boundary conditions including discharge, left, right sides #####################################################
stream <- c(3,-2)
y_str <- which(rownames(sur_updatecell_dem)==stream[2])
x_str <- which(colnames(sur_updatecell_dem)==stream[1])
#sur_storage_dem[y_str,x_str] <- hydro_data$WL[1]
#sur_storage_dem[y_str,x_str] <- c(0)

##### Setting for results matrix of runoff time series data ######################################################################
asm_hydro_data <- hydro_data[1,]
hydro_data <- hydro_data[-1,]
runoff_timeseries <- matrix(NA, 6, nrow(hydro_data)); 
rownames(runoff_timeseries) <- c('order','obs_runoff[m]','obs_runoff[m/s]','cal_sur_runoff[m]',
                                 'cal_tot_runoff[m]','cal_tot_runoff[m/s]')
runoff_timeseries[1,] <- c(1:nrow(hydro_data))
runoff_timeseries[2,] <- hydro_data$WL # [m]
runoff_timeseries[3,] <- (3*10^(-5))*exp(36.139*hydro_data$WL) # [m3/s]
runoff_timeseries[4:5,1] <- runoff_timeseries[2,1]
runoff_timeseries[6,1] <- runoff_timeseries[3,1]

##### Setting for results matrix of soil moisture time series data ###############################################################
sur_sm_timeseries <- matrix(NA, nrow(SoilText_Mat)*2, nrow(hydro_data))
rownames(sur_sm_timeseries) <- c(rownames(SoilText_Mat), paste0("pred_",rownames(SoilText_Mat)))
sur_sm_timeseries[1:nrow(SoilText_Mat),] <- t(hydro_data[,grep("10", colnames(hydro_data))])
sur_sm_timeseries[22:42,1] <- sur_sm_timeseries[1:21,1]
##################################################################################################################################
# Cellular Automata Calculation ##################################################################################################
##################################################################################################################################
for(time in 1:c(nrow(hydro_data)*time_interval/time_step)){
  #for(time in 1:10800){
  #tic()
  #time <- c(2)
  ################################################################################################################################
  ##### Hydrological components ##################################################################################################
  ################################################################################################################################
  
  ##### Rainfall [mm] 
  R_t <- (hydro_data$Rain[ceiling(time*time_step/time_interval)]/time_interval)*(time_step)
  
  ##### Evapotranspiration [mm]
  ET_t <- Evptr*(time_step/hour_to_sec)
  
  ##### Interception [mm/hr]
  R_cum <- R_cum+R_t
  P_cum <- P_max*(1-exp(-0.046*LAI*R_cum/P_max))
  P_t <- P_cum - P_pre
  P_pre <- P_cum
  Intercep_Total <- Intercep_Total + P_t
  
  ##### Infiltration [mm]
  I_t <- (i_final + (i_init-i_final)*exp(-i_decay*time_step/hour_to_sec))*(time_step/hour_to_sec)
  Infiltra_Sum <- c(0)
  
  ##### Effective rainfall
  R_eff <- R_t - P_t - ET_t
  ################################################################################################################################
  ##### Set-up for DEM matrix ####################################################################################################
  ################################################################################################################################
  for(y in 1:num_row){
  for(x in 1:num_col){
      
      R_Sur <- sur_storage_dem[y,x] + sur_waterdepth_dem[y,x]
      STRG_mat <- STRG_sur2(R_eff, I_t, K_Sur, R_Sur, 
                            E_Sat_Sur_cons[y], Depth_Sur, 
                            theta_Res_Sur_dem[y,x], theta_Sat_Sur_dem[y,x])
      sur_storage_dem[y,x] <- STRG_mat[1,1]; # soil moisture
      sur_waterdepth_dem[y,x] <- STRG_mat[1,2]; # water depth for runoff
      Infiltra_Sum <- Infiltra_Sum + STRG_mat[1,3] # infiltration 
      
    }}
  
  sur_waterdepth_dem[nrow(sur_dem),] <- c(0) # discharge at bottom area
  sur_waterdepth_dem[1,] <- c(0) # discharge at top area
  sur_waterdepth_dem[,1] <- c(0) # discharge at left area
  sur_waterdepth_dem[,ncol(sur_dem)] <- c(0) # discharge at right area
  
  Infiltra_Total <- Infiltra_Total + (Infiltra_Sum/(num_row * num_col)) # Infiltration update
  sur_cellheight_dem <- sur_dem + sur_storage_dem + sur_waterdepth_dem # Cellheight update
  sur_updatecell_dem <- matrix(0, num_row, num_col) # Updatecell update
  rownames(sur_updatecell_dem) <- as.numeric(rownames(sur_dem))
  colnames(sur_updatecell_dem) <- as.numeric(colnames(sur_dem))
  
  #x<-c(50); y <- c(55)
  
  ################################################################################################################################
  ##### Surface flow routing algorithm ###########################################################################################
  ################################################################################################################################
  for(y in 1:num_row){
  for(x in 1:num_col){
      
      if(sur_waterdepth_dem[y,x] <= c(0)){

        next
      } 
      sur_neighbor_cell <- matrix(NA, 6, 9); 
      colnames(sur_neighbor_cell) <- c('upleft','up','upright','left','point','right','doleft','do','doright')
      sur_neighbor_cell <- Neighbor_mat(x,y,num_col,num_row,sur_neighbor_cell,
                                        sur_dem,sur_cellheight_dem, sur_waterdepth_dem)
      sur_point_cell <- as.matrix(sur_neighbor_cell[,5]); sur_neighbor_cell <- sur_neighbor_cell[,-5]
      sur_updatecell_dem <- Transition_fun2(flowdir,manning_coeff_sur,resolution_cell,
                                            time_step,sur_neighbor_cell,sur_point_cell,sur_updatecell_dem)
    }}
  
  ################################################################################################################################
  ##### Export PNG file for water depth ##########################################################################################
  ################################################################################################################################
  #x11()
  #rast <- raster(sur_updatecell_dem)
  #extent(rast) <- c(c(min(dem$X)+0.5),max(dem$X),c(min(dem$Y)+0.5),max(dem$Y))
  #projection(rast) <- CRS("+proj=longlat +datum=WGS84")
  #plot(rast, col=terrain.colors(30), main="sur_updatecell_dem")
  # png(paste0(time_step,"sec_",time,"_plot.png"),width=4000,height=2000,res=200)
  # layout(matrix(c(1,2,3,4), 2, 2, byrow=F))
  # rast <- raster(sur_dem)
  # extent(rast) <- c(c(min(dem$X)+0.5),max(dem$X),c(min(dem$Y)+0.5),max(dem$Y))
  # projection(rast) <- CRS("+proj=longlat +datum=WGS84")
  # plot(rast, col=terrain.colors(30), main="sur_dem")
  # rast <- raster(sur_storage_dem)
  # extent(rast) <- c(c(min(dem$X)+0.5),max(dem$X),c(min(dem$Y)+0.5),max(dem$Y))
  # projection(rast) <- CRS("+proj=longlat +datum=WGS84")
  # plot(rast, col=terrain.colors(30), main="sur_Strg")
  # rast <- raster(sur_waterdepth_dem)
  # extent(rast) <- c(c(min(dem$X)+0.5),max(dem$X),c(min(dem$Y)+0.5),max(dem$Y))
  # projection(rast) <- CRS("+proj=longlat +datum=WGS84")
  # plot(rast, col=terrain.colors(30), main="sur_waterdepth_dem")
  # rast <- raster(sur_updatecell_dem)
  # extent(rast) <- c(c(min(dem$X)+0.5),max(dem$X),c(min(dem$Y)+0.5),max(dem$Y))
  # projection(rast) <- CRS("+proj=longlat +datum=WGS84")
  # plot(rast, col=terrain.colors(30), main="sur_updatecell_dem")
  # dev.off()
  
  ################################################################################################################################
  ##### Undating waterdepth and check error ######################################################################################
  ################################################################################################################################
  sur_waterdepth_dem <- sur_updatecell_dem + sur_waterdepth_dem; 
  if(c(min(sur_waterdepth_dem)+0.0001) < c(0) ){
    print("there is negative values in waterdepth")
    break
  }
  Runoff_Total <- 
    Runoff_Total+(3*10^(-5))*exp(36.139*c(sur_waterdepth_dem[y_str,x_str]/m_to_mm))*time_step # [m3]
  
  ################################################################################################################################
  ##### Write water balance values ###############################################################################################
  ################################################################################################################################
  if((time*time_step/time_interval) - floor(time*time_step/time_interval) == 0){
    
    ##### updating runoff time series 
    runoff_timeseries[4,c(floor(time*time_step/time_interval))] <- c(sur_waterdepth_dem[y_str,x_str])/m_to_mm # [m]
    runoff_timeseries[5,c(floor(time*time_step/time_interval))] <- c(sur_waterdepth_dem[y_str,x_str])/m_to_mm # [mm]
    runoff_timeseries[6,c(floor(time*time_step/time_interval))] <- 
      c(3*10^(-5))*exp(36.139*runoff_timeseries[4,c(floor(time*time_step/time_interval))]) # [m3/s]
    
    ##### updating soil moisture time series
    # 2022-05-24 287 row edit: sur_waterdepth_dem removed
    for(t in 1:c(nrow(sur_sm_timeseries)/2)){
      sm_y_point <- which.min(abs(SoilText_Mat$y[t] - as.numeric(rownames(sur_cellheight_dem))))
      sm_x_point <- which.min(abs(SoilText_Mat$x[t] - as.numeric(colnames(sur_cellheight_dem))))
      sur_sm_timeseries[21+t, c(floor(time*time_step/time_interval))] <- sur_storage_dem[sm_y_point,sm_x_point] # [mm]
    }
  }
  
  print(time_step * time)
  if(R_t == c(0) && max(sur_updatecell_dem) == c(0) && max(sur_waterdepth_dem) == c(0) ){
    break

  }
  
 }

Rain_space <- sum((hydro_data$Rain/m_to_mm)*ncol(sur_dem)*(resolution_cell/m_to_mm)*nrow(sur_dem)*(resolution_cell/m_to_mm)) # [m3]
Intercep_space  <- 
  sum((Intercep_Total/m_to_mm)*ncol(sur_dem)*(resolution_cell/m_to_mm)*nrow(sur_dem)*(resolution_cell/m_to_mm)) # [m3]
Infiltra_space  <- 
  sum((Infiltra_Total/m_to_mm)*ncol(sur_dem)*(resolution_cell/m_to_mm)*nrow(sur_dem)*(resolution_cell/m_to_mm)) # [m3]
Runoff_space_cal  <- Runoff_Total # [m3]
runoff_diff <- runoff_timeseries[2,] - runoff_timeseries[2,1]
runoff_flowrate <- (3*10^(-5))*exp(36.139*runoff_diff)
Runoff_space_obs <- sum(runoff_flowrate[runoff_flowrate>0])*length(runoff_diff[runoff_diff>0])*time_interval # [m3]

##################################################################################################################################
# Cellular Automata Outputs ######################################################################################################
##################################################################################################################################
setwd('C:/Users/82105/Documents/cellularAutomata/CA_model_hillslope_220612')
water_balance_mat <- matrix(NA, 1, 5)
colnames(water_balance_mat) <- c('Rainfall[m3]','Interception[m3]','Infiltration[m3]',
                                 'Observed runoff[m3]','Calculated runoff[m3]')
water_balance_mat[1,] <- c(round(Rain_space,2),round(Intercep_space,2),round(Infiltra_space,2),
                           round(Runoff_space_obs,2),round(Runoff_space_cal,2))

sur_sm_timeseries <- cbind(rbind(t(asm_hydro_data[1,grep("10", colnames(hydro_data))]),
                                 t(asm_hydro_data[1,grep("10", colnames(hydro_data))])),sur_sm_timeseries)

runoff_timeseries <- cbind(c(0,asm_hydro_data$WL,c(3*10^(-5))*exp(36.139*asm_hydro_data$WL),asm_hydro_data$WL,
                             asm_hydro_data$WL,c(3*10^(-5))*exp(36.139*asm_hydro_data$WL)),runoff_timeseries)

write.csv(water_balance_mat,
          paste0(getwd(),'/Results/WaterBudget_',flowdir,time_step,"sec_0",K_Sur*100,'K',E_Sat_Sur_cons_val,'%E_Sat_test12.csv'))
write.csv(runoff_timeseries,
          paste0(getwd(),'/Results/runoff_timeseries_',flowdir,time_step,"sec_0",K_Sur*100,'K',E_Sat_Sur_cons_val,'%E_Sat_test12.csv'))
write.csv(sur_sm_timeseries,
          paste0(getwd(),'/Results/sur_sm_timeseries_',flowdir,time_step,"sec_0",K_Sur*100,'K',E_Sat_Sur_cons_val,'%E_Sat_test12.csv'))
CA_results <-
  list("Water_balance_mat"=water_balance_mat,"Runoff_timeseries"=runoff_timeseries,"SM_timeseries"=sur_sm_timeseries)
toc()
return(CA_results)



