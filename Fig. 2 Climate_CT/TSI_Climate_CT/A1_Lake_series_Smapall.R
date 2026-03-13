##### Here, we take data from 2002 to 2012 as an exemple to calculate CT-Climate-TSI relations
#### The computational procedures remain consistent for TSI (2016-2024) 
#### Alongside this workflow for processing 2002-2012 data, 
#### we directly provide finalized results for the TSI (2016-2024) in this folder to support downstream simulations


### Full Code
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
setwd("D:/Lake teleconnection code/Fig. 2 Climate_CT/TSI_Climate_CT")
source('Function_Needed.R')

#####
##### Time series
print(Sys.time())
TSI <- fread("D:/Lake teleconnection code/Data_preprocess/A2_Detrend/De_Seasonal_FilterFilled.csv") %>% data.frame()

nlakes = nrow(TSI) # time series length
nyear = ncol(TSI)-3

# climate
Temp2m <- fread("De_Seasonal_Temp.csv") %>% data.frame() 
Preciptn <- fread("De_Seasonal_Preciptn.csv") %>% data.frame() 
SLH <- fread("De_Seasonal_SLH.csv") %>% data.frame() 
SM <- fread("De_Seasonal_SoilMoisture.csv") %>% data.frame() 
SSH <- fread("De_Seasonal_SSH.csv") %>% data.frame() 
SurPressure <- fread("De_Seasonal_SurfacePressure.csv") %>% data.frame() 
Temp2mmax <- fread("De_Seasonal_Tempmax.csv") %>% data.frame() 
SolarRadiation <- fread("De_Seasonal_SolarRadiation.csv") %>% data.frame()  
Windspeed <- fread("De_Seasonal_Windspeed.csv") %>% data.frame()  

# Climate to Lake
Climate_Lake_CCM <- fread('Sig_CCM_result.csv')
Climate_Lake_Smap <- fread('Smap_results.csv')

# CT
CT <- fread("CT_De_Seasonal.csv") %>% data.frame()  # Loading dataset with raw variables
colnames(CT) <- c('AO','PDO','PNA','WP','EPNP','ENSO','AMO','TNA','TSA','NAO','EA','CAR','IOD','AAO')

CT_CT <- fread('CCM_rho_sig.csv') %>% as.matrix()

# potentially direct linked
Lake_Lake <- fread('CCM_Sync_harm.csv')

######################################################################################
####################### Test the time delay & Significance ###########################
######################################################################################
{
  Sig_fun <- function(nrow){
    row = Lake_Lake[nrow,]
    libno = as.numeric(row[,1])
    tarno = as.numeric(row[,2])
    
    {# start
      lib.climate.ccmresult <- Climate_Lake_CCM[which(Climate_Lake_CCM$LakeNo == libno),]
      tar.climate.ccmresult <- Climate_Lake_CCM[which(Climate_Lake_CCM$LakeNo == tarno),]
      lib.climate.ccmresult.loc <- (which(lib.climate.ccmresult>0)-1) %>% tail(-1,)
      tar.climate.ccmresult.loc <- (which(tar.climate.ccmresult>0)-1) %>% tail(-1,)
      
      Lake_smap <- NULL
      
      if (length(lib.climate.ccmresult.loc)>0 & length(tar.climate.ccmresult.loc)>0){ # Lakes are influenced by CFs
        Temp.libloc <- which(Temp2m$LakeNo == libno)
        SLH.libloc <- which(SLH$LakeNo == libno)
        lib.climate.block =  Generate_Climate_Block(Temp.libloc,SLH.libloc,Temp2m,Preciptn,SLH,SM,SSH,SurPressure,Temp2mmax,SolarRadiation,Windspeed)
        
        Temp.tarloc <- which(Temp2m$LakeNo == tarno)
        SLH.tarloc <- which(SLH$LakeNo == tarno)
        tar.climate.block =  Generate_Climate_Block(Temp.tarloc,SLH.tarloc,Temp2m,Preciptn,SLH,SM,SSH,SurPressure,Temp2mmax,SolarRadiation,Windspeed)
        
        climatesmap_lib = Climate_Lake_Smap[which(Climate_Lake_Smap$LakeNo == libno),]
        climatesmap_tar = Climate_Lake_Smap[which(Climate_Lake_Smap$LakeNo == tarno),]
        
        lib_cf_block = lib.climate.block[,lib.climate.ccmresult.loc,drop = FALSE]
        tar_cf_block = tar.climate.block[,tar.climate.ccmresult.loc,drop = FALSE]
        
        for (libCF in c(1:length(lib.climate.ccmresult.loc))){
          libcli_num = lib.climate.ccmresult.loc[libCF]
          for (tarCF in c(1:length(tar.climate.ccmresult.loc))){
            tarcli_num = tar.climate.ccmresult.loc[tarCF]
            lib_cf_strength = climatesmap_lib[[as.numeric(libcli_num)+3]]
            tar_cf_strength = climatesmap_tar[[as.numeric(tarcli_num)+3]]
            
            block_climate = data.frame(lib_cf_block[,libCF,drop=FALSE], tar_cf_block[,tarCF,drop=FALSE])
            ifcli_cli_linked = Climate.Climate.CCMcheck(block_climate,Date)
            
            if (any(as.logical(ifcli_cli_linked))){ # climate factors are linked
              lib_ct_climate = Climate.CT.CCMcheck(block_climate[,1,drop=FALSE],CT,Date)
              tar_ct_climate = Climate.CT.CCMcheck(block_climate[,2,drop=FALSE],CT,Date)
              
              lib_cli_ct_loc = which(lib_ct_climate > 0)
              tar_cli_ct_loc = which(tar_ct_climate > 0)
              lib_smap_ct_strength = CT.Climate.SmapFunc(lib_ct_climate, lib_cf_block[,libCF,drop=FALSE], CT, alllength=320)
              tar_smap_ct_strength = CT.Climate.SmapFunc(tar_ct_climate, tar_cf_block[,tarCF,drop=FALSE], CT, alllength=320)
              
              if (sum(!is.na(lib_smap_ct_strength))>0 & sum(!is.na(tar_smap_ct_strength))>0){ # both climate are influenced by CT
                for (m in lib_cli_ct_loc){
                  for (n in tar_cli_ct_loc){
                    if (m == n){
                      ct_relation1 = 1
                      ct_relation2 = 1
                    }else{
                      ct_relation1 = CT_CT[m, n]
                      ct_relation2 = CT_CT[n, m]
                    }
                    if (any(as.logical(ct_relation1)) | any(as.logical(ct_relation2))){
                      if (m == n){
                        if (libcli_num == tarcli_num){
                          result_row = c(libno, tarno, libcli_num, tarcli_num, m, n,
                                         lib_cf_strength, tar_cf_strength, lib_smap_ct_strength[m], tar_smap_ct_strength[n], 1)
                        }else{
                          result_row = c(libno, tarno, libcli_num, tarcli_num, m, n,
                                         lib_cf_strength, tar_cf_strength, lib_smap_ct_strength[m], tar_smap_ct_strength[n], 3)
                        }
                      }else{
                        if (libcli_num == tarcli_num){
                          result_row = c(libno, tarno, libcli_num, tarcli_num, m, n,
                                         lib_cf_strength, tar_cf_strength, lib_smap_ct_strength[m], tar_smap_ct_strength[n], 2)
                        }else{
                          result_row = c(libno, tarno, libcli_num, tarcli_num, m, n,
                                         lib_cf_strength, tar_cf_strength, lib_smap_ct_strength[m], tar_smap_ct_strength[n], 4)
                        }
                      }
                    }else{ # 9 means CFs are linked, they are influenced by CTs, but CTs are not linked
                      result_row = c(libno, tarno, libcli_num, tarcli_num, m, n,
                                     lib_cf_strength, tar_cf_strength, lib_smap_ct_strength[m], tar_smap_ct_strength[n], 9)
                    }
                    Lake_smap = rbind(Lake_smap,result_row)
                  }
                }
              }else{
                # 10 means lakes are influenced by CFs, but at least one CF is not influenced by CTs
                result_row = c(libno, tarno, libcli_num, tarcli_num, NaN, NaN,
                               lib_cf_strength, tar_cf_strength, NaN, NaN, 10)
                Lake_smap = rbind(Lake_smap,result_row)
              }
              
            }else{ # climate factors are not linked
              lib_ct_climate = Climate.CT.CCMcheck(block_climate[,1,drop=FALSE],CT,Date)
              tar_ct_climate = Climate.CT.CCMcheck(block_climate[,2,drop=FALSE],CT,Date)
              
              lib_cli_ct_loc = which(lib_ct_climate > 0)
              tar_cli_ct_loc = which(tar_ct_climate > 0)
              lib_smap_ct_strength = CT.Climate.SmapFunc(lib_ct_climate, lib_cf_block[,libCF,drop=FALSE], CT, alllength=320)
              tar_smap_ct_strength = CT.Climate.SmapFunc(tar_ct_climate, tar_cf_block[,tarCF,drop=FALSE], CT, alllength=320)
              
              if (sum(!is.na(lib_smap_ct_strength))>0 & sum(!is.na(tar_smap_ct_strength))>0){ # both climate are influenced by CT
                for (m in lib_cli_ct_loc){
                  for (n in tar_cli_ct_loc){
                    if (m == n){
                      ct_relation1 = 1
                      ct_relation2 = 1
                    }else{
                      ct_relation1 = CT_CT[m, n]
                      ct_relation2 = CT_CT[n, m]
                    }
                    # Check if there is any relationship
                    if (any(as.logical(ct_relation1)) | any(as.logical(ct_relation2))){
                      if (m==n){
                        if (libcli_num == tarcli_num){
                          result_row = c(libno, tarno, libcli_num, tarcli_num, m, n,
                                         lib_cf_strength, tar_cf_strength,lib_smap_ct_strength[m], tar_smap_ct_strength[n], 5)
                        }else{
                          result_row = c(libno, tarno, libcli_num, tarcli_num, m, n,
                                         lib_cf_strength, tar_cf_strength,lib_smap_ct_strength[m], tar_smap_ct_strength[n], 7)
                        }
                      }# if m==n
                      else{
                        if (libcli_num == tarcli_num){
                          result_row = c(libno, tarno, libcli_num, tarcli_num, m, n,
                                         lib_cf_strength, tar_cf_strength,lib_smap_ct_strength[m], tar_smap_ct_strength[n], 6)
                        }else{
                          result_row = c(libno, tarno, libcli_num, tarcli_num, m, n,
                                         lib_cf_strength, tar_cf_strength,lib_smap_ct_strength[m], tar_smap_ct_strength[n], 8)
                        }
                      }
                      Lake_smap = rbind(Lake_smap,result_row) 
                    }else{ # 11 means CFs are not linked, they are influenced by CTs, but CTs are not linked
                      result_row = c(libno, tarno, libcli_num, tarcli_num, m, n,
                                     lib_cf_strength, tar_cf_strength,lib_smap_ct_strength[m], tar_smap_ct_strength[n], 11)
                      Lake_smap = rbind(Lake_smap,result_row) 
                    }
                  } # n
                } # m
              }else{
                # 12 means lakes are influenced by CFs, CFs are not linked, but at least one CF is not influenced by CTs
                result_row = c(libno, tarno, libcli_num, tarcli_num, NaN, NaN,
                               lib_cf_strength, tar_cf_strength, NaN, NaN, 12)
                Lake_smap = rbind(Lake_smap,result_row)
              }
              
            } # else
          } # tarCF
        } # libCF
      }else{
        # 13 means at least one lake is not influenced by CFs, not to mention CTs
        result_row = c(libno, tarno, NaN, NaN, NaN, NaN,
                       NaN, NaN, NaN, NaN, 13)
        Lake_smap = rbind(Lake_smap,result_row)
      }
    } # end
    sorted_lake_smap = Lake_smap[order(Lake_smap[, ncol(Lake_smap)]), ]
    return(sorted_lake_smap)
  }
  
  multifun = function(id){
    library(rEDM)
    library(dplyr)  # %>%
    library(data.table)
    
    df = df_list[[id]]
    result_df = lapply(X = df, FUN = Sig_fun) %>% do.call(rbind, .) %>% data.frame()
    fwrite(result_df,
           file = paste0(result_path, 'Smap_',id,'.csv'), row.names = F)
  }
  
  ### CCM analysis for all causality links
  ### Parallel operation
  Date = fread("D:/Lake teleconnection code/Data_preprocess/Date.csv")
  
  # Preparation for the parallel loop
  cores = 36
  cl = makeCluster(cores)
  result_path = paste0('D:\\Lake teleconnection code\\Fig. 2 Climate_CT\\TSI_Climate_CT\\')
  n = cores*10
  nr = nrow(Lake_Lake)
  
  # Pre-set all the candidate causal links
  Index_list = split(Lake_Lake, rep(1:n,each=ceiling(nr/n),length.out=nr)) # the nrow must be larger than the number of groups
  df_list = split(1:nr, rep(1:n,each=ceiling(nr/n),length.out=nr))
  
  print(Sys.time())
  clusterExport(cl,varlist = c('result_path', 'df_list', 'Sig_fun','multifun', 'laf', 'yearday_anom','Embed2','CT.Climate.SmapFunc',
                               'Generate_Climate_Block', 'Lag_fun','Kendall_fun','Surrosig_fun','Climate.Climate.CCMcheck','Climate.CT.CCMcheck','Classifi.cli',
                               'Date','TSI', 'Temp2m','Preciptn','SLH','SM','SSH','SurPressure','Temp2mmax','SolarRadiation','Windspeed',
                               'Climate_Lake_CCM','Climate_Lake_Smap','CT','CT_CT','Lake_Lake'))
  parLapply(cl = cl, X = 1:length(df_list), fun = multifun)
  stopCluster(cl)  # Release CPU
  clusterEvalQ
  print(Sys.time())
  
  ## rind multi-result
  file_list = list.files(result_path, pattern = paste0('^Smap_.*\\.csv'), full.names = T)
  id<- sub("^.*Smap_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
  file_list = file_list[order(id)]
  Total_result = lapply(file_list, FUN = function(p){
    fread(p) %>% data.frame()}) %>% do.call(rbind, .)
  fwrite(Total_result,paste0(result_path, 'Lake_Climate_CT_Smap_0212'),row.names = F)
  file.remove(file_list)
}
