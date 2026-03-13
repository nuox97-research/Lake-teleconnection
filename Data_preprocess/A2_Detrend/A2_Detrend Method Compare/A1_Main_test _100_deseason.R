### 
rm(list = ls())

library(dplyr)
library(data.table)  # fread
library(parallel)
library(rEDM)
library(KFAS)
library(forecast)

# Set working directory
setwd("D:/Lake teleconnection code/Data_preprocess/A2_Detrend/A2_Detrend Method Compare")

X <- NULL
Y <- NULL

X[1] = 0.8
Y[1] = 0.8

for (i in 2:3000){
  X[i] <- X[i-1]*(3.8-3.8*X[i-1]-0.02*Y[i-1])
  Y[i] <- Y[i-1]*(3.5-3.5*Y[i-1]-0.1*X[i-1])
}
# plot(X,type = 'l',col='red')
# lines(Y,type = 'l',col = 'blue')

## CCM
CCM_cal = function(X,Y) {
  # Ebest
  BestE <- NULL
  Data = data.frame(X = X,Y = Y)
  
  for (i in c(1:ncol(Data))){
    Sim_r <- simplex(Data[,i],lib=c(1,floor(nrow(Data)/2)),pred=c(floor(nrow(Data)/2)+1,nrow(Data)),E=c(2:10))
    BestE <- rbind(BestE, Sim_r[which.min(Sim_r$mae),"E"][1])
  }
  
  # time lag
  indmat <- expand.grid(1:2, 1:2) 
  indmat <- indmat[indmat$Var1 != indmat$Var2, ]
  colnames(indmat)=c('Effect','Cause')
  
  lag_all <- NULL
  for (i in c(1:nrow(indmat))){
    lib = Data[,indmat[i,1]] %>% c()
    target =  Data[,indmat[i,2]] %>% c()
    Ebest = as.numeric(BestE[indmat[i,1]])
    
    lag_rho <- NULL
    block_temp = data.frame(lib,target)  
    colnames(block_temp) = c("lib","target")    
    test_df = data.frame(time = 1:nrow(block_temp),block_temp)
    colnames(test_df) = c("time","lib","target")
    
    cmap = CCM(dataFrame = test_df, E = Ebest, Tp = 0, columns = "lib",
               target = "target", libSizes = NROW(block_temp), random = FALSE, includeData = TRUE)
    rho <- cmap$CCM1_PredictStat$rho
    
    returnrow <- cbind(indmat[i,],rho)
    lag_all <- rbind(lag_all,returnrow)
  }
  
  return(lag_all)
}

CCM_perf = CCM_cal(X,Y)

####################################################################
###################### Simple_sin CCM test #########################
####################################################################
Dist_num = 6
indmat <- do.call(rbind, replicate(50, expand.grid(1:Dist_num, 1:Dist_num), simplify = FALSE))
indmat <- indmat[order(indmat[, 1], indmat[, 2]), ]
rownames(indmat) = c(1:nrow(indmat))
colnames(indmat) = c('Way_x','Way_y')

{
  Disturbance <- function(Xlength,Way){
    time = 1:Xlength
    # simple sinmodel
    if (Way == 1){
      Seasonseries <- 0.5*sin(2*pi*time/365.25)
    }
    # multiple sinmodel
    if (Way == 2){
      Seasonseries <- 0.5*sin(2*pi*time/365.25) + 0.4*sin(2*pi*time/(365.25/2))
    }
    # varying_amp
    if (Way == 3){
      base_sin = 0.5*sin(2*pi*time/365.25)
      Seasonseries <- base_sin
      for (i in c(1:ceiling(Xlength/365.25))) {
        i_0 = (i-1)*floor(365.25)
        i_1 = i * floor(365.25)
        if (i_1 > Xlength-1){
          i_1 = Xlength
        }
        Seasonseries[i_0:i_1] = base_sin[i_0:i_1] + base_sin[i_0:i_1] * (sample(-5:5, 1) / 50) #20
      }
      # plot(Seasonseries)
    }
    # varying_timing
    if (Way == 4){
      base_sin = 0.5*sin(2*pi*time/365.25)
      shifted_sin = 0.5*sin(2*pi*time/365.25 - pi/4)
      Seasonseries <- base_sin
      for (i in c(1:ceiling(Xlength/365.25))) {
        i_0 = (i-1)*floor(365.25)
        i_1 = i * floor(365.25)
        if (i_1 > Xlength-1){
          i_1 = Xlength
        }
        Seasonseries[i_0:i_1] = base_sin[i_0:i_1] + shifted_sin[i_0:i_1] * (sample(-5:5, 1) / 50) #20
      }
      # plot(Seasonseries)
    }
    # varying_amp_timing
    if (Way == 5){
      base_sin =  0.5*sin(2*pi*time/365.25)
      shifted_sin = 0.5*sin(2*pi*time/365.25 - pi/4)
      Seasonseries <- base_sin
      for (i in c(1:ceiling(Xlength/365.25))) {
        i_0 = (i-1)*floor(365.25)
        i_1 = i * floor(365.25)
        if (i_1 > Xlength-1){
          i_1 = Xlength
        }
        Seasonseries[i_0:i_1] = base_sin[i_0:i_1] + shifted_sin[i_0:i_1] * (sample(-5:5, 1) / 50) + base_sin[i_0:i_1] * (sample(-5:5, 1) / 50) #20
      }
      # plot(Seasonseries)
    }
    # extra_noise
    if (Way == 6){
      base_sin =  0.5*sin(2*pi*time/365.25)
      Seasonseries <- base_sin
      for (i in c(1:ceiling(Xlength/365.25))) {
        i_0 = (i-1)*floor(365.25)
        i_1 = i * floor(365.25)
        if (i_1 > Xlength-1){
          i_1 = Xlength
        }
        Seasonseries[i_0:i_1] = base_sin[i_0:i_1] +  rnorm(n = length(base_sin[i_0:i_1])) * (sample(-5:5, 1) / 100) #20
      }
      # plot(Seasonseries)
    }
    return(Seasonseries)
  }
  
  DetrendCompare <- function(libtar_idx){

    Direction = 1 # 1: X -> Y
    Direction2 = 2 # 2: Y -> X
    
    CCM_stl_record <- NULL
    # CCM_ssm_record <- NULL
    CCM_harm_record <- NULL
    CCM_mean_record <- NULL
    
    X_dist = X + Disturbance(length(X),libtar_idx[[1]])
    Y_dist = Y + Disturbance(length(X),libtar_idx[[2]])
    # CCM_Dist = CCM_cal(X_dist,Y_dist)
    
    # stl
    X_dist_destl <- stl(ts(X_dist,frequency = 365), s.window="periodic", robust=TRUE)
    X_dist_stl = X_dist - X_dist_destl$time.series[, 'seasonal']
    # plot(X_dist,type = 'l')
    # lines(1:3000,X_dist_stl,type = 'l',col = 'red')
    Y_dist_destl <- stl(ts(Y_dist,frequency = 365), s.window="periodic", robust=TRUE) 
    Y_dist_stl = Y_dist - Y_dist_destl$time.series[, 'seasonal']
    # Y_dist_stl_seas = Y_dist_destl$time.series[, 'seasonal']
    # plot(Y_dist,type = 'l')
    # lines(1:3000,Y_dist_stl,type = 'l',col = 'red')
    # lines(1:3000,Y_dist_stl_seas,type = 'l',col = 'blue')
    CCM_stl = CCM_cal(X_dist_stl,Y_dist_stl)
    
    result_stl = cbind(CCM_stl$rho[Direction],CCM_stl$rho[Direction2])
    CCM_stl_record <- rbind(CCM_stl_record,result_stl)
    
    # 3harmonic
    n = length(X_dist)
    order = 3
    xt = (1:n)/365
    x_rad = xt * 2 * pi
    nr_indep = order*2 + 1
    Indep_X = matrix(0,length(X_dist), nr_indep)
    Indep_X[,1] = x_rad;
    i = 2;
    for (freq in c(1:3)){
      cos_freq = cos(x_rad * freq);
      sin_freq = sin(x_rad * freq);
      Indep_X[,i] = cos_freq;
      i = i + 1;
      Indep_X[,i] = sin_freq;
      i = i + 1;
    }
    myfit <- lm(X_dist~Indep_X)
    Season = myfit$fitted.values
    Deall_X = X_dist - Season
    # plot(X_dist,type = 'l')
    # lines(1:3000,Deall_X,type = 'l',col = 'red')
    # lines(Season,type = 'l',col = 'blue')
    
    n = length(Y_dist)
    order = 3
    xt = (1:n)/365
    x_rad = xt * 2 * pi
    nr_indep = order*2 + 1
    Indep_Y = matrix(0,length(Y_dist), nr_indep)
    Indep_Y[,1] = x_rad;
    i = 2;
    for (freq in c(1:3)){
      cos_freq = cos(x_rad * freq);
      sin_freq = sin(x_rad * freq);
      Indep_Y[,i] = cos_freq;
      i = i + 1;
      Indep_Y[,i] = sin_freq;
      i = i + 1;
    }
    myfit <- lm(Y_dist~Indep_Y)
    Season = myfit$fitted.values
    Deall_Y = Y_dist - Season
    
    CCM_harm = CCM_cal(Deall_X,Deall_Y)
    result_harm = cbind(CCM_harm$rho[Direction],CCM_harm$rho[Direction2])
    CCM_harm_record <- rbind(CCM_harm_record,result_harm)
    
    # Mean
    tsx = ts(X_dist,frequency = 365)
    tstime = seq.Date(from = as.Date("2021-01-01"),by = "day",length.out = 3000)
    Data_x = data.frame(Date = tstime,data = tsx,month = month(tstime))
    Spcmonth_nomean <- NULL
    for (i in c(1:12)){
      Spcmonth = Data_x[which(Data_x$month == i),]
      meanval = mean(Spcmonth$data)
      Spcmonth$data = Spcmonth$data-meanval
      Spcmonth_nomean = rbind(Spcmonth_nomean,Spcmonth)
    }
    Spcmonth_nomean = Spcmonth_nomean[order(Spcmonth_nomean$Date),]
    Detrendall_X = Spcmonth_nomean$data
    # plot(X_dist,type = 'l')
    # lines(1:3000,Detrendall_X,type = 'l',col = 'red')
    # LFit_X = lm(Spcmonth_nomean$data~Spcmonth_nomean$Date)
    # Year_X = LFit_X$fitted.values
    # Detrendall_X = Spcmonth_nomean$data - Year_X
    
    tsy = ts(Y_dist,frequency = 365)
    tstime = seq.Date(from = as.Date("2021-01-01"),by = "day",length.out = 3000)
    Data_y = data.frame(Date = tstime,data = tsy,month = month(tstime))
    Spcmonth_nomean <- NULL
    for (i in c(1:12)){
      Spcmonth = Data_y[which(Data_y$month == i),]
      meanval = mean(Spcmonth$data)
      Spcmonth$data = Spcmonth$data-meanval
      Spcmonth_nomean = rbind(Spcmonth_nomean,Spcmonth)
    }
    Spcmonth_nomean = Spcmonth_nomean[order(Spcmonth_nomean$Date),]
    Detrendall_Y = Spcmonth_nomean$data
    # LFit_Y = lm(Spcmonth_nomean$data~Spcmonth_nomean$Date)
    # Year_Y = LFit_Y$fitted.values
    # Detrendall_Y = Spcmonth_nomean$data - Year_Y
    
    CCM_mean = CCM_cal(Detrendall_X,Detrendall_Y)
    result_mean = cbind(CCM_mean$rho[Direction],CCM_mean$rho[Direction2])
    CCM_mean_record <- rbind(CCM_mean_record,result_mean)
    
    X2Y <- cbind(libtar_idx[[1]],libtar_idx[[2]],CCM_stl_record,CCM_harm_record,CCM_mean_record) %>% data.frame()

    return(X2Y)
  }
  
  multifun = function(id){
    library(rEDM)
    library(dplyr)  # %>%
    library(data.table)
    library(KFAS)
    library(forecast)
    
    df = Index_list[[id]]
    result_df = apply(X = df, MARGIN=1, FUN = DetrendCompare) %>% unlist() %>% as.character() %>% t() %>% data.frame()
    fwrite(result_df,
           file = paste0(result_path, 'AMethods_',id,'.csv'), row.names = F)
  }
  
  cores = 2
  cl = makeCluster(cores)
  result_path = paste0('D:\\Lake teleconnection code\\Data_preprocess\\A2_Detrend\\A2_Detrend Method Compare\\')

  n = cores*10
  nr = nrow(indmat)
  # Split the mission for all cores
  Index_list = split(indmat, rep(1:n, each=ceiling(nr/n), length.out=nr))
  
  print(Sys.time())
  clusterExport(cl,varlist = c('result_path', 'Index_list', 'CCM_cal','Disturbance','DetrendCompare','multifun','X','Y'))
  parLapply(cl = cl, X = 1:n, fun = multifun)
  stopCluster(cl)  # Release CPU
  clusterEvalQ
  print(Sys.time())
  
  ## rind multi-result
  file_list = list.files(result_path, pattern = paste0('^AMethods_.*\\.csv'), full.names = T)
  id<- sub("^.*AMethods_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
  file_list = file_list[order(id)]
  Total_result = lapply(file_list, FUN = function(p){
    fread(p) %>% data.frame()}) %>% do.call(rbind, .)
  fwrite(Total_result,paste0(result_path, 'Methods_Compare.csv'),row.names = F)
  file.remove(file_list)
}




