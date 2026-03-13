### Full Code
### Lake CCM
rm(list = ls())
library(dplyr) 
library(rEDM) # Empirical dynamical modeling for CCM analysis
library(Kendall) # Kendall's tau test for the convergence of CCM
library(data.table)  # fread
library(parallel)


######################################################################################
##################### Loading time series data for Lake ##############################
######################################################################################
# Set working directory
setwd("D:/TSI 02-12/CCM/Moving Window/Moving window CCM deseason")

folder_path <- 'D:/TSI 02-12/CCM/Moving Window/Matlab_datacompile/Harm_detrendall'
Date_file_names <- list.files(folder_path, pattern = paste0('^Date_.*\\.csv'), full.names = T)
id<- sub("^.*Date_MW(.*)\\.csv$", "\\1", Date_file_names) %>% as.numeric()
Date_file_names <- Date_file_names[order(id)]

TSI_file_names <- list.files(folder_path, pattern = paste0('^TSI_.*\\.csv'), full.names = T)
TSI_file_names <- TSI_file_names[order(id)]

Num_MW <- length(Date_file_names) %>% as.numeric()

#####
{
  for (k in 1:Num_MW){
    cat('Process:',k,'/',Num_MW,'\n')
    
    Worldlake <- fread(TSI_file_names[k]) # Loading dataset with raw variables
    colnames(Worldlake) = c('LakeID','Lat','Lon',as.character(c(1:(ncol(Worldlake)-3)))) # rename the variables
    nlakes = nrow(Worldlake) # time series length
    nyear = ncol(Worldlake)-3
    
    ##################### calculate the best E for each residual time series ##############
    {
      # define function
      bestE_fun = function(row){
        FUI_value = t(row) %>% c() %>% tail(., -3)
        FUI_value = data.frame(TSI = FUI_value)
        
        Sim_r <- simplex(FUI_value[,1],lib=c(1,floor(nrow(FUI_value)/2)),pred=c(floor(nrow(FUI_value)/2)+1,nrow(FUI_value)),E=c(2:11))
        Ebest <-Sim_r[which.min(Sim_r$mae),"E"][1]
        return_row = c(lon = row[[3]], lat = row[[2]], Ebest = Ebest)
        return(return_row)
      }
      
      multifun = function(id){
        library(rEDM)
        library(dplyr)  # %>%
        library(data.table)
        
        df = df_list[[id]]
        result_df = apply(X = df, MARGIN=1, FUN = bestE_fun) %>% t() %>% data.frame()
        fwrite(result_df,
               file = paste0(result_path, 'Ebest_',id,'.csv'), row.names = F)
      }
      
      ### CCM analysis for all causality links
      ### Parallel operation
      cores = 2
      cl = makeCluster(cores)
      result_path = paste0('D:\\TSI 02-12\\CCM\\Moving Window\\Moving window CCM deseason\\')
      
      n = cores*10
      df_list = split(Worldlake, rep(1:n,each=ceiling(nlakes/n),length.out=nlakes)) # the nrow must be larger than the number of groups
      
      print(Sys.time())
      clusterExport(cl,varlist = c('result_path', 'df_list', 'bestE_fun','multifun'))
      parLapply(cl = cl, X = 1:n, fun = multifun)
      stopCluster(cl)  # Release CPU
      clusterEvalQ
      print(Sys.time())
      
      # rbind multi-progress results
      file_list = list.files(result_path, pattern = paste0('^Ebest_.*\\.csv'), full.names = T)
      id<- sub("^.*Ebest_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
      file_list = file_list[order(id)]
      Ebest_df = lapply(file_list, FUN = function(p){
        read.csv(p) %>% data.frame()}) %>% do.call(rbind, .)
      # write.csv(Total_result ,paste0(result_path, 'Total_Ebest.csv'),row.names = F)
      file.remove(file_list)
    }
    
    ################################ Test the time delay #################################
    {
      laf = function(x,y,lagf){
        n <- NROW(x)
        x.t <- x
        y.t <- y
        if(lagf <= 0){x.t = x.t[(1-lagf):n];y.t = y.t[1:(n+lagf)]} 
        if(lagf > 0){x.t = x.t[1:(n-lagf)];y.t = y.t[(1+lagf):n]}  
        return(cbind(x.t,y.t))
      }
      
      lagcal_fun <- function(libtar_idx){
        lib = t(Worldlake[libtar_idx[[1]], ]) %>% c() %>% tail(., -3)
        target =  t(Worldlake[libtar_idx[[2]], ]) %>% c() %>% tail(., -3)
        Ebest = as.numeric(Ebest_df[libtar_idx[[1]], 3])
        
        lag_rho <- NULL
        for (lags in -9:0){
          block_temp = data.frame(laf(lib,target,lagf=lags))  
          colnames(block_temp) = c("lib","target")    
          test_df = data.frame(time = 1:nrow(block_temp),block_temp)
          colnames(test_df) = c("time","lib","target")
          
          cmap = CCM(dataFrame = test_df, E = Ebest, Tp = 0, columns = "lib",
                     target = "target", libSizes = NROW(block_temp), random = FALSE, includeData = TRUE)
          rho <- cmap$CCM1_PredictStat$rho
          
          lag_rho <- cbind(lag_rho, rho)
        }
        return(cbind(lib = Worldlake[libtar_idx[[1]],'LakeID'] %>% as.numeric(), target = Worldlake[libtar_idx[[2]],'LakeID'] %>% as.numeric(), lag_rho))
      }
      
      multifun = function(id){
        library(rEDM)
        library(dplyr)  # %>%
        library(data.table)
        
        df = Index_list[[id]]
        result_df = apply(X = df, MARGIN=1, FUN = lagcal_fun) %>% t() %>% data.frame()
        fwrite(result_df,
               file = paste0(result_path, 'Lags_',id,'.csv'), row.names = F)
      }
      
      ### CCM analysis for all causality links
      ### Parallel operation
      # Pre-set all the candidate causal links
      indmat <- expand.grid(1:nlakes, 1:nlakes) 
      indmat <- indmat[indmat$Var1 != indmat$Var2, ]
      colnames(indmat)=c('Effect','Cause')
      
      # Preparation for the parallel loop
      cores = 30
      cl = makeCluster(cores)
      result_path = paste0('D:\\TSI 02-12\\CCM\\Moving Window\\Moving window CCM deseason\\')
      
      n = cores*10
      nr = nrow(indmat)
      # Split the mission for all cores
      Index_list = split(indmat, rep(1:n, each=ceiling(nr/n), length.out=nr))
      
      print(Sys.time())
      clusterExport(cl,varlist = c('result_path', 'Index_list', 'laf', 'lagcal_fun','multifun','Ebest_df','Worldlake'))
      parLapply(cl = cl, X = 1:n, fun = multifun)
      stopCluster(cl)  # Release CPU
      clusterEvalQ
      print(Sys.time())
      
      ## rind multi-result
      file_list = list.files(result_path, pattern = paste0('^Lags_.*\\.csv'), full.names = T)
      id<- sub("^.*Lags_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
      file_list = file_list[order(id)]
      Lags_Record = lapply(file_list, FUN = function(p){
        fread(p) %>% data.frame()}) %>% do.call(rbind, .)
      colnames(Lags_Record) = c('lib','target','rho','rho','rho','rho','rho','rho','rho','rho','rho','rho')
      head(Lags_Record)
      # fwrite(Total_result,paste0(result_path, 'Lags_Record.csv'),row.names = F)
      file.remove(file_list)
      
    }
    
    ################################# find maximum lag ##############################
    {
      # define function
      Maxlag_fun = function(row){
        a <- row[3:length(row)]
        a <- as.numeric(t(a))
        a[is.na(a)] <- 0
        a[which(a<0)] <- 0
        optlag <- which.max(a)
        maxrho <- a[which.max(a)]
        
        out <- cbind(row[1],row[2],optlag-10,maxrho)
        return(out)
      }
      
      multifun = function(id){
        library(dplyr)  # %>%
        library(data.table)
        
        df = df_list[[id]]
        result_df = apply(X = df, MARGIN=1, FUN = Maxlag_fun) %>% t() %>% data.frame()
        fwrite(result_df,
               file = paste0(result_path, 'Maxlag_',id,'.csv'), row.names = F)
      }
      
      ### CCM analysis for all causality links
      ### Parallel operation
      cores = 2
      cl = makeCluster(cores)
      result_path = paste0('D:\\TSI 02-12\\CCM\\Moving Window\\Moving window CCM deseason\\')
      
      n = cores*10
      df_list = split(Lags_Record, rep(1:n,each=ceiling(nlakes/n),length.out=nlakes)) # the nrow must be larger than the number of groups
      
      print(Sys.time())
      clusterExport(cl,varlist = c('result_path', 'df_list', 'Maxlag_fun','multifun'))
      parLapply(cl = cl, X = 1:n, fun = multifun)
      stopCluster(cl)  # Release CPU
      # clusterEvalQ
      print(Sys.time())
      
      # rbind multi-progress results
      file_list = list.files(result_path, pattern = paste0('^Maxlag_.*\\.csv'), full.names = T)
      id<- sub("^.*Maxlag_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
      file_list = file_list[order(id)]
      Optimal_lag = lapply(file_list, FUN = function(p){
        fread(p) %>% data.frame()}) %>% do.call(rbind, .)
      colnames(Optimal_lag) = c('lib','target','Optlag','Maxrho')
      # fwrite(Total_result ,paste0(result_path, 'Optimal_lags.csv'),row.names = F)
      file.remove(file_list)
    }
    
    ################## Kendall's tau test and Student's t-test significance #################
    {
      laf = function(x,y,lagf){
        n <- NROW(x)
        x.t <- x
        y.t <- y
        if(lagf <= 0){x.t = x.t[(1-lagf):n];y.t = y.t[1:(n+lagf)]} 
        if(lagf > 0){x.t = x.t[1:(n-lagf)];y.t = y.t[(1+lagf):n]}  
        return(cbind(x.t,y.t))
      }
      
      Sig_fun <- function(libtar_idx){
        library(Kendall)
        libnum = libtar_idx[[1]]
        tarnum = libtar_idx[[2]]
        libloc = which(Worldlake$LakeID == libnum)
        tarloc = which(Worldlake$LakeID == tarnum)
        
        lib = t(Worldlake[libloc, ]) %>% c() %>% tail(., -3)
        target =  t(Worldlake[tarloc, ]) %>% c() %>% tail(., -3)
        Ebest = as.numeric(Ebest_df[libloc, 3])
        
        k <- which(Optimal_lag$lib == libnum & Optimal_lag$target == tarnum)
        optlag = Optimal_lag[k,]$Optlag
        
        block_temp = data.frame(laf(lib,target,lagf=optlag))  
        colnames(block_temp) = c("lib","target") 
        
        libsizes <- seq(11, nrow(block_temp), by = 5)
        out.temp <- ccm(block = block_temp,E = Ebest,lib_column="lib",target_column="target", 
                        lib_sizes = libsizes,random_libs = FALSE,stats_only = T, )
        
        out.temp$`lib:target`[is.na(out.temp$`lib:target`)]<-0                
        out.temp$`lib:target`[which(!is.finite(out.temp$`lib:target`))] <- 0  
        
        out.temp_means <- out.temp$`lib:target`
        difference<-out.temp_means[length(out.temp_means)]-out.temp_means[1]
        
        # Kendall's tau test and Student's t-test significance
        kend <- MannKendall(out.temp$`lib:target`);  # Kendall's tau test for monototic increase,both (1) Kendall's tau and (2) terminal rho are significantly larger than zero
        
        # Fisher's Z test
        # method 1
        np<-nrow(block_temp);
        crirho <- qt(0.95,np-1)/(np-2+qt(0.95,np-1)^2);
        # method 2
        rho.Lmax = out.temp$`lib:target`[which.max(libsizes)]
        rho.Lmin = out.temp$`lib:target`[1]
        N = tail(out.temp$LibSize,1) 
        z=abs(0.5*(log((1+rho.Lmax)/(1-rho.Lmax))-log((1+rho.Lmin)/(1-rho.Lmin)))*(2/(N-3))^-0.5)
        z.p=(1-pnorm(z))
        
        significance_monototic_1 <- (kend$tau[1]>0)*(kend$sl[1]<0.05)*(out.temp_means[length(out.temp_means)]>crirho) # ccm.sig records the significance of each CCM
        significance_monototic_2 <- (kend$sl[1]<0.05)*(kend$tau[1]>0)*(z.p<0.05)
        
        rho_mono<-out.temp_means[length(out.temp_means)]
        P_value_mono<-kend$sl[1]
        
        return(cbind(lib = Worldlake[libloc,'LakeID'] %>% as.numeric(), target = Worldlake[tarloc,'LakeID'] %>% as.numeric(), optlag, difference, P_value_mono, rho_mono, significance_monototic_1, significance_monototic_2))
      }
      
      multifun = function(id){
        library(rEDM)
        library(dplyr)  # %>%
        library(data.table)
        
        df = Index_list[[id]]
        result_df = apply(X = df, MARGIN=1, FUN = Sig_fun) %>% t() %>% data.frame()
        fwrite(result_df,
               file = paste0(result_path, 'SigLags_',id,'.csv'), row.names = F)
      }
      
      ### CCM analysis for all causality links
      ### Parallel operation
      Optimal_lag <- Optimal_lag[order(Optimal_lag$target),] 
      Optimal_lag <- subset(Optimal_lag,Optimal_lag$Maxrho>0)
      
      # Pre-set all the candidate causal links
      indmat <- Optimal_lag[,c(1,2)]
      colnames(indmat)=c('Effect','Cause')
      
      # Preparation for the parallel loop
      cores = 20
      cl = makeCluster(cores)
      result_path = paste0('D:\\TSI 02-12\\CCM\\Moving Window\\Moving window CCM deseason\\')
      
      n = cores*10
      nr = nrow(indmat)
      # Split the mission for all cores
      Index_list = split(indmat, rep(1:n, each=ceiling(nr/n), length.out=nr))
      
      print(Sys.time())
      clusterExport(cl,varlist = c('result_path', 'Index_list', 'laf', 'Sig_fun','multifun','Ebest_df','Optimal_lag','Worldlake'))
      parLapply(cl = cl, X = 1:n, fun = multifun)
      stopCluster(cl)  # Release CPU
      clusterEvalQ
      print(Sys.time())
      
      ## rind multi-result
      file_list = list.files(result_path, pattern = paste0('^SigLags_.*\\.csv'), full.names = T)
      id<- sub("^.*SigLags_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
      file_list = file_list[order(id)]
      SigLags = lapply(file_list, FUN = function(p){
        fread(p) %>% data.frame()}) %>% do.call(rbind, .)
      # fwrite(Total_result ,paste0(result_path, 'SigLags_Record.csv'),row.names = F)
      file.remove(file_list)
      
      colnames(SigLags) <- c('lib', 'target', 'optlag', 'difference', 'P_value_mono', 'rho_mono', 'significance_monototic_ttest', 'significance_monototic_Fztest')
      SigLags$significance <- ifelse(SigLags$significance_monototic_Fztest==1,1,0) 
      # head(SigLags)
      
      Maxlag_rhodifer <- cbind(SigLags,Optimal_lag[,4])
      colnames(Maxlag_rhodifer) = c('lib','target','optlag','difference','P_value_mono','rho_mono','significance_monototic_ttest','significance_monototic_Fztest','Significance','Optlag_rho')
      # head(Maxlag_rhodifer)
      # fwrite(Maxlag_rhodifer,"Maxlag_rhodifer.csv")
    }
    
    #########################################################################################
    ##################################### Null surrogates ###################################
    #########################################################################################
    
    ############################### generating surrogates #############################
    {
      Date <- fread(Date_file_names[k])
      
      yearday_anom <- function(t,x){
        # t: date formatted with POSIXt
        # x: time-series values to compute seasonal mean and anomaly
        doy <- as.numeric(strftime(t, format = "%j"))
        I_use <- which(!is.na(x))
        # create time indices to use for smoothing, replicating data to "wrap around"
        doy_sm <- rep(doy[I_use],3) + rep(c(-366,0,366),each=length(I_use))
        x_sm <- rep(x[I_use],3)
        xsp <- smooth.spline(doy_sm, y = x_sm, w = NULL, spar = 0.8, cv = NA,
                             all.knots = TRUE,keep.data = TRUE, df.offset = 0)
        xbars <- data.frame(t=t,doy=doy) %>%
          left_join(data.frame(doy=xsp$x,xbar=xsp$y),by='doy') 
        xbar <- xbars$xbar 
        
        out = data.frame(t=t,mean=xbar,anomaly=(x - xbar))
        names(out) <- c('date','mean','anomaly')
        return(out)
      }
      
      Surrogate_TSI<-NULL
      for (i in 1:dim(Worldlake)[1]) {
        
        df.in <- data.frame(Date, t(Worldlake[i,]) %>% tail(-3,))
        colnames(df.in) = c('date','TSI')
        
        set.seed(1000) 
        
        Mean_anom <- yearday_anom(df.in$date,df.in$TSI)
        TSI.mean <- Mean_anom$mean
        TSI.anom <- Mean_anom$anomaly
        
        TSI.sample <- do.call(cbind,
                              lapply(1:100, function(ii) {
                                I_na <- is.na(TSI.anom)
                                out <- TSI.mean
                                out[I_na] <- NA
                                out[!I_na] <- out[!I_na] + sample(TSI.anom[!I_na],sum(!I_na),replace = FALSE)
                                return(out)
                              }))
        
        Surrogate_TSI<-cbind(Surrogate_TSI,TSI.sample)
      }
    }
    
    ################################ rho of surrogates ################################
    {
      laf <- function(x,y,lagf){
        n <- NROW(x)
        x.t <- x
        y.t <- y
        if(lagf <= 0){x.t = x.t[(1-lagf):n];y.t = y.t[1:(n+lagf)]} 
        if(lagf > 0){x.t = x.t[1:(n-lagf)];y.t = y.t[(1+lagf):n]}  
        return(cbind(x.t,y.t))
      }
      
      Surrogate_fun <- function(libtar_idx){
        libnum = libtar_idx[[1]]
        tarnum = libtar_idx[[2]]
        libloc = which(Worldlake$LakeID == libnum)
        tarloc = which(Worldlake$LakeID == tarnum)
        
        Ebest <- as.numeric(Ebest_df[libloc, 3])
        max_range<-c(1:100)+(tarloc-1)*100
        lib <- t(Worldlake[libloc, ]) %>% c() %>% tail(., -3)
        
        lag_rho<-NULL
        for (i_surr in  max_range){
          dat_temp <- data.frame(lib = lib, target = Surrogate_TSI[,i_surr])
          
          lag_rho_surroga <- do.call(cbind,
                                     lapply(0, function(ii) {
                                       block_temp = data.frame(laf(dat_temp[,1],dat_temp[,2],lagf=ii))  
                                       colnames(block_temp) = c("lib","target") 
                                       lib_ccm <- c(1,NROW(block_temp))
                                       pred_ccm <- c(1,NROW(block_temp))
                                       
                                       df.out.ccm <- ccm(block=block_temp, E=Ebest, lib=lib_ccm, pred =pred_ccm,       
                                                         lib_column="lib", target_column="target", lib_sizes=NROW(block_temp),    
                                                         exclusion_radius=0, random_libs = FALSE, num_sample=1,)
                                       return(df.out.ccm$`lib:target`)
                                     }))
          
          lag_rho<-cbind(lag_rho,lag_rho_surroga)
        }
        return(cbind(libtar_idx[[1]],libtar_idx[[2]],lag_rho))
      }
      
      multifun <- function(id){
        library(rEDM)
        library(dplyr)  # %>%
        library(data.table)
        
        df = Index_list[[id]]
        result_df = apply(X = df, MARGIN=1, FUN=Surrogate_fun) %>% t() %>% data.frame()
        fwrite(result_df,
               file = paste0(result_path, 'ASurrogate_',id,'.csv'), row.names = F)
      }
      
      ### CCM analysis for all causality links
      ### Parallel operation
      Maxlag_rhodifer<-subset(Maxlag_rhodifer,Significance==1)  
      # dim(Maxlag_rhodifer) 
      
      # Pre-set all the candidate causal links
      indmat <- cbind(Maxlag_rhodifer$lib, Maxlag_rhodifer$target) %>% data.frame()
      colnames(indmat)=c('Effect','Cause')
      
      # Preparation for the parallel loop
      cores = 20
      cl = makeCluster(cores)
      result_path = paste0('D:\\TSI 02-12\\CCM\\Moving Window\\Moving window CCM deseason\\')
      
      n = cores*10
      nr = nrow(indmat)
      # Split the mission for all cores
      Index_list = split(indmat, rep(1:n, each=ceiling(nr/n), length.out=nr))
      
      print(Sys.time())
      clusterExport(cl,varlist = c('result_path', 'Index_list', 'laf', 'Surrogate_fun', 'multifun', 'Ebest_df', 'Worldlake', 'Surrogate_TSI'))
      parLapply(cl = cl, X = 1:n, fun = multifun)
      stopCluster(cl)  # Release CPU
      clusterEvalQ
      print(Sys.time())
      
      ## rind multi-result
      file_list = list.files(result_path, pattern = paste0('^ASurrogate_.*\\.csv'), full.names = T)
      id<- sub("^.*ASurrogate_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
      file_list = file_list[order(id)]
      rho_surrodata1 = lapply(file_list, FUN = function(p){ 
        fread(p) %>% data.frame()}) %>% do.call(rbind, .)
      colnames(rho_surrodata1) = c('lib','target',rep("out.temp", times=100))
      # fwrite(Total_result ,paste0(result_path, 'NullSurrogate_rho.csv'),row.names = F)
      file.remove(file_list)
      
    }
    
    ######################### Test the significant of surrogate #######################
    {
      rho_surrodata<-rho_surrodata1[1:dim(rho_surrodata1)[1],3:dim(rho_surrodata1)[2]]
      # head(rho_surrodata)
      # dim(rho_surrodata) # 
      
      Surrosig_fun <- function(nrow){
        aa <- rho_surrodata[nrow,] %>% data.frame()
        aa[which(aa<0)] <- 0
        surrorho <- Maxlag_rhodifer[nrow,]
        bb<-(sum(surrorho[[10]] < aa))/(length(aa))
        return(bb) 
      }
      
      multifun = function(id){
        library(dplyr)  # %>%
        library(data.table)
        
        df = df_list[[id]]
        result_df = lapply(X = df, FUN = Surrosig_fun) %>% do.call(rbind, .) %>% data.frame()
        fwrite(result_df,
               file = paste0(result_path, 'Surrosig_',id,'.csv'))
      }
      
      
      ### CCM analysis for all causality links
      ### Parallel operation
      cores = 2
      cl = makeCluster(cores)
      result_path   = paste0('D:\\TSI 02-12\\CCM\\Moving Window\\Moving window CCM deseason\\')
      
      n = cores*10
      ididx <- c(1:dim(Maxlag_rhodifer)[1]) %>% data.frame()
      nlen <- dim(ididx)[1]
      df_list = split(c(1:nlen), rep(1:n,each=ceiling(nlen/n),length.out=nlen)) # the nrow must be larger than the number of groups
      
      print(Sys.time())
      clusterExport(cl,varlist = c('result_path', 'df_list', 'Surrosig_fun','multifun','rho_surrodata','Maxlag_rhodifer'))
      parLapply(cl = cl, X = 1:n, fun = multifun)
      stopCluster(cl)  # Release CPU
      # clusterEvalQ
      print(Sys.time())
      
      # rbind multi-progress results
      file_list = list.files(result_path, pattern = paste0('^Surrosig_.*\\.csv'), full.names = T)
      id<- sub("^.*Surrosig_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
      file_list = file_list[order(id)]
      interaction_significance = lapply(file_list, FUN = function(p){
        fread(p) %>% data.frame()}) %>% do.call(rbind, .)
      colnames(interaction_significance) <- c('Surrosig')
      # fwrite(Total_result ,paste0(result_path, 'Surrosig_total.csv'),row.names = F)
      file.remove(file_list)
      
      Maxlag_rhodifer <- cbind(Maxlag_rhodifer,interaction_significance)
      Maxlag_rhodifer_sig <- subset(Maxlag_rhodifer, Surrosig<0.05)
      # dim(Maxlag_rhodifer_sig)
      # head(Maxlag_rhodifer_sig)
      # fwrite(Maxlag_rhodifer_sig,'Maxlag_rhodifer_sig.csv')
    }
    
    ########################### transfer CCM result to matrix #########################
    {
      CCMrho <- Maxlag_rhodifer_sig
      head(CCMrho)

      CCMresult <- data.frame(CCMrho[,c(1,2,3,10)])
      head(CCMresult)
      
      library(geosphere)
      Dist <- NULL
      for (i in 1:nrow(CCMresult)){
        libnum = CCMresult$lib[i]
        tarnum = CCMresult$target[i]
        libloc = which(Worldlake$LakeID == libnum)
        tarloc = which(Worldlake$LakeID == tarnum)
        OLoc = cbind(Worldlake[libloc,"Lon"],Worldlake[libloc,"Lat"])
        DLoc = cbind(Worldlake[tarloc,"Lon"],Worldlake[tarloc,"Lat"])
        Dist = rbind(Dist,distGeo(OLoc,DLoc)/1000)
      }
      
      CCMresult <- cbind(CCMresult,Dist)
      head(CCMresult)
      
      fwrite(CCMresult,paste0('CCM_MW_',substring(TSI_file_names[k],75)))
      
      ##########
      # CCMrho_matrix = matrix(0,nlakes,nlakes)
      # for (i in c(1:nrow(CCMresult))){
      #   lib = CCMresult[i,1]
      #   target = CCMresult[i,2]
      #   libloc = which(Worldlake$LakeID == lib)
      #   tarloc = which(Worldlake$LakeID == target)
      #   CCMrho_matrix[tarloc,libloc] = CCMresult[i,4]
      # }
      # 
      # CCMrho_matrix <- data.frame(CCMrho_matrix)
      # 
      # fwrite(CCMrho_matrix,'CCM_rho_sig.csv')
    }
    
  }
}





















   















