library(metafor)
library(dplyr)
################################################
#                                              #
#          Supplementary Table 17              #
#                                              #
################################################
# ==================== Define general function ====================
fit_models_and_get_results <- function(data_subset, group_name) {
  results <- data.frame()
  if (nrow(data_subset) < 12) {
    return(results)
  }
  data_subset$Latitude_abs <- abs(data_subset$Latitude)
  
  # model 1
  try({
    rmod.linear <- rma.mv(yi, vi, 
                          mods = ~ Latitude_abs,
                          method = "REML",
                          random = list(~1 | Paper_id / Observation_id), 
                          data = data_subset)
    
    results <- rbind(results, data.frame(
      Group = group_name,
      Model = "Linear",
      Formula = "yi ~ Latitude_abs",
      Intercept = coef(rmod.linear)[1],
      Slope1 = coef(rmod.linear)[2],
      Slope2 = NA,
      AIC = AIC(rmod.linear),
      AICc = AICc(rmod.linear),
      stringsAsFactors = FALSE
    ))
  }, silent = TRUE)
  
  # model 2
  try({
    data_quad <- data_subset
    data_quad$Latitude_abs_sq <- data_quad$Latitude_abs^2
    
    rmod.quadratic <- rma.mv(yi, vi, 
                             mods = ~ Latitude_abs + Latitude_abs_sq,
                             method = "REML",
                             random = list(~1 | Paper_id / Observation_id), 
                             data = data_quad)
    
    results <- rbind(results, data.frame(
      Group = group_name,
      Model = "Quadratic (raw)",
      Formula = "yi ~ Latitude_abs + Latitude_abs^2",
      Intercept = coef(rmod.quadratic)[1],
      Slope1 = coef(rmod.quadratic)[2],
      Slope2 = coef(rmod.quadratic)[3],
      AIC = AIC(rmod.quadratic),
      AICc = AICc(rmod.quadratic),
      stringsAsFactors = FALSE
    ))
  }, silent = TRUE)
  
  # model 3
  try({
    data_log <- data_subset
    data_log$logdist <- log(data_log$Latitude_abs)
    data_log$logdist_sq <- data_log$logdist^2
    
    rmod.log_quadratic <- rma.mv(yi, vi, 
                                 mods = ~ logdist_sq + logdist,
                                 method = "REML",
                                 random = list(~1 | Paper_id / Observation_id), 
                                 data = data_log)
    
    results <- rbind(results, data.frame(
      Group = group_name,
      Model = "Log-quadratic",
      Formula = "yi ~ log(Latitude_abs) + log(Latitude_abs)^2",
      Intercept = coef(rmod.log_quadratic)[1],
      Slope1 = coef(rmod.log_quadratic)[3],  # logdist
      Slope2 = coef(rmod.log_quadratic)[2],  # logdist_sq
      AIC = AIC(rmod.log_quadratic),
      AICc = AICc(rmod.log_quadratic),
      stringsAsFactors = FALSE
    ))
  }, silent = TRUE)
  
  return(results)
}


data <- read.csv("Supplementary Data 1.csv", header = TRUE)
# ==================== Taxa ====================
taxa_list <- c("Algae+Protist", "Radiata", "Lophotrochozoa", "Ecdysozoa", "Deuterostomia")

all_taxa_results <- data.frame()
for (taxa_name in taxa_list) {
  cat("  分析 Taxa:", taxa_name, "\n")
  data_subset <- subset(data, Taxa == taxa_name)
  taxa_results <- fit_models_and_get_results(data_subset, paste("Taxa:", taxa_name))
  all_taxa_results <- rbind(all_taxa_results, taxa_results)
}

write.csv(all_taxa_results, "Model_Comparison_All_Taxa.csv", row.names = FALSE)


# ==================== A-Taxa ====================

data_stressorA <- subset(data, Stressor == "A")
stressorA_taxa_results <- data.frame()
for (taxa_name in taxa_list) {
  cat("  分析 Taxa:", taxa_name, "(Stressor A)\n")
  data_subset <- subset(data_stressorA, Taxa == taxa_name)
  taxa_results <- fit_models_and_get_results(data_subset, paste("Taxa:", taxa_name, "(Stressor A)"))
  stressorA_taxa_results <- rbind(stressorA_taxa_results, taxa_results)
}
write.csv(stressorA_taxa_results, "Model_Comparison_StressorA_All_Taxa.csv", row.names = FALSE)


# ==================== W-Taxa ====================
data_stressorW <- subset(data, Stressor == "W")

stressorW_taxa_results <- data.frame()
for (taxa_name in taxa_list) {
  cat("  分析 Taxa:", taxa_name, "(Stressor W)\n")
  data_subset <- subset(data_stressorW, Taxa == taxa_name)
  taxa_results <- fit_models_and_get_results(data_subset, paste("Taxa:", taxa_name, "(Stressor W)"))
  stressorW_taxa_results <- rbind(stressorW_taxa_results, taxa_results)
}
write.csv(stressorW_taxa_results, "Model_Comparison_StressorW_All_Taxa.csv", row.names = FALSE)


# ==================== Stressor+ Overall ====================
stressor_groups <- list(
  "W" = subset(data, Stressor == "W"),
  "A" = subset(data, Stressor == "A"),
  "A-W" = subset(data, Stressor == "A-W"),
  "H" = subset(data, Stressor == "H"),
  "Overall" = data
)

all_stressor_results <- data.frame()
for (stressor_name in names(stressor_groups)) {
  cat("  分析 Stressor:", stressor_name, "\n")
  data_subset <- stressor_groups[[stressor_name]]
  stressor_results <- fit_models_and_get_results(data_subset, paste("Stressor:", stressor_name))
  all_stressor_results <- rbind(all_stressor_results, stressor_results)
}
write.csv(all_stressor_results, "Model_Comparison_All_Stressors.csv", row.names = FALSE)





################################################
#                                              #
#        Fig. 3_Supplementary Table 4          #
#                                              #
################################################
################# Fig. 3a_All taxa #################
######### Overall #########
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix

#                estimate      se     zval    pval    ci.lb    ci.ub    
# intrcpt        -0.1946  0.0801  -2.4306  0.0151  -0.3515  -0.0377  * 
# Latitude_abs    0.0034  0.0021   1.5751  0.1152  -0.0008   0.0076   

######## Acidification ########
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data <- subset(data, Stressor == "A")
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data$Latitude_abs_sq <- data$Latitude_abs^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs + Latitude_abs_sq,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                    estimate      se     zval    pval    ci.lb    ci.ub    
# intrcpt            0.4957  0.2254   2.1987  0.0279   0.0538   0.9375  * 
# Latitude_abs      -0.0247  0.0120  -2.0618  0.0392  -0.0481  -0.0012  * 
# Latitude_abs_sq    0.0003  0.0002   1.7882  0.0737  -0.0000   0.0006  .    

######### Warming ########
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data <- subset(data, Stressor == "W")
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                estimate      se     zval    pval    ci.lb    ci.ub      
# intrcpt        -0.5828  0.0948  -6.1453  <.0001  -0.7687  -0.3969  *** 
# Latitude_abs    0.0129  0.0027   4.8858  <.0001   0.0078   0.0181  *** 

######### Acidification-Warming ########
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Stressor == "A-W")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#               estimate      se     zval    pval     ci.lb    ci.ub    
# intrcpt      21.2112  9.4337   2.2484  0.0245    2.7214  39.7010  * 
# logdist_sq    2.1035  0.8785   2.3944  0.0166    0.3816   3.8254  * 
# logdist     -13.5354  5.7916  -2.3371  0.0194  -24.8868  -2.1840  *   

######### Oxygen #########
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data, Stressor == "H")
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                 estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt         0.6424  0.6123   1.0491  0.2941  -0.5577  1.8425    
# Latitude_abs   -0.0214  0.0134  -1.5972  0.1102  -0.0477  0.0049    

################# Fig. 3b #################
################ Overall ################
# Algae+Protist
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Algae+Protist")
nrow(data)
data$Latitude_abs_sq <- data$Latitude_abs^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs + Latitude_abs_sq,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                   estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt            0.4963  0.3945   1.2582  0.2083  -0.2769  1.2695    
# Latitude_abs      -0.0198  0.0203  -0.9760  0.3291  -0.0595  0.0199    
# Latitude_abs_sq    0.0003  0.0002   1.1857  0.2358  -0.0002  0.0007  

# Radiata
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Radiata")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                 estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt        -0.3295  0.2221  -1.4839  0.1378  -0.7648  0.1057    
# Latitude_abs    0.0028  0.0086   0.3293  0.7419  -0.0141  0.0198   

# Lophotrochozoa
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Lophotrochozoa")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#               estimate      se     zval    pval     ci.lb    ci.ub    
# intrcpt       7.7738  5.6219   1.3828  0.1667   -3.2450  18.7926    
# logdist_sq    0.5968  0.4240   1.4073  0.1593   -0.2344   1.4279    
# logdist      -4.3334  3.0884  -1.4031  0.1606  -10.3866   1.7198    

# Ecdysozoa
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Ecdysozoa")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#             estimate      se     zval    pval     ci.lb    ci.ub     
# intrcpt      -6.2432  2.1277  -2.9342  0.0033  -10.4134  -2.0730  ** 
# logdist_sq   -0.5876  0.2070  -2.8382  0.0045   -0.9934  -0.1818  ** 
# logdist       3.8152  1.3369   2.8538  0.0043    1.1950   6.4355  ** 


# Deuterostomia
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Deuterostomia")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#               estimate      se     zval    pval     ci.lb    ci.ub    
# intrcpt      -6.9244  3.3521  -2.0657  0.0389  -13.4945  -0.3544  * 
# logdist_sq   -0.6239  0.3119  -2.0005  0.0454   -1.2351  -0.0126  * 
# logdist       4.1387  2.0566   2.0123  0.0442    0.1077   8.1696  * 


######### Acidification #######
# Algae+Protist+Latitude
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Algae+Protist")
data <- subset(data, Stressor == "A")
nrow(data)
data$Latitude_abs_sq <- data$Latitude_abs^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs + Latitude_abs_sq,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt            0.6581  0.4434   1.4842  0.1378  -0.2110  1.5273    
# Latitude_abs      -0.0255  0.0237  -1.0747  0.2825  -0.0719  0.0210    
# Latitude_abs_sq    0.0003  0.0003   1.0522  0.2927  -0.0003  0.0009    



# Radiata+Latitude
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Radiata")
data <- subset(data, Stressor == "A")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                 estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt         0.1185  0.1846   0.6423  0.5207  -0.2432  0.4803    
# Latitude_abs   -0.0104  0.0059  -1.7457  0.0809  -0.0220  0.0013  . 

# Lophotrochozoa+Latitude
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Lophotrochozoa")
data <- subset(data, Stressor == "A")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#              estimate       se     zval    pval     ci.lb    ci.ub    
# intrcpt      15.4166  11.3364   1.3599  0.1739   -6.8023  37.6355    
# logdist_sq    1.2083   0.8947   1.3505  0.1768   -0.5453   2.9618    
# logdist      -8.6561   6.3739  -1.3580  0.1745  -21.1488   3.8367     


# Ecdysozoa+Latitude
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Ecdysozoa")
data <- subset(data, Stressor == "A")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#               estimate      se     zval    pval    ci.lb    ci.ub    
# intrcpt      -5.0412  2.1106  -2.3885  0.0169  -9.1778  -0.9045  * 
# logdist_sq   -0.5171  0.2080  -2.4857  0.0129  -0.9248  -0.1094  * 
# logdist       3.2709  1.3274   2.4642  0.0137   0.6693   5.8725  * 

# Deuterostomia+Latitude
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Deuterostomia")
data <- subset(data, Stressor == "A")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#               estimate       se     zval    pval     ci.lb     ci.ub    
# intrcpt      58.1049  23.6127   2.4607  0.0139   11.8249  104.3850  * 
# logdist_sq    3.6786   1.6626   2.2125  0.0269    0.4199    6.9372  * 
# logdist     -29.3665  12.5319  -2.3433  0.0191  -53.9286   -4.8044  * 

######## Warming ########
# Algae+Protist+Latitude
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Algae+Protist")
data <- subset(data, Stressor == "W")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt        -0.3222  0.2571  -1.2535  0.2100  -0.8260  0.1816    
# Latitude_abs    0.0136  0.0057   2.4038  0.0162   0.0025  0.0247  * 

# Radiata+Latitude
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Radiata")
data <- subset(data, Stressor == "W")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                 estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt        -0.9408  0.7479  -1.2578  0.2084  -2.4068  0.5251    
# Latitude_abs    0.0301  0.0366   0.8236  0.4102  -0.0416  0.1019    


# Lophotrochozoa+Latitude
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Lophotrochozoa")
data <- subset(data, Stressor == "W")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt        -0.4111  0.2780  -1.4790  0.1392  -0.9559  0.1337    
# Latitude_abs    0.0090  0.0065   1.3861  0.1657  -0.0037  0.0216    


# Ecdysozoa+Latitude
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Ecdysozoa")
data <- subset(data, Stressor == "W")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#              estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt      -3.6095  2.7954  -1.2912  0.1966  -9.0884  1.8694    
# logdist_sq   -0.2686  0.2848  -0.9431  0.3456  -0.8267  0.2896    
# logdist       1.9135  1.8005   1.0628  0.2879  -1.6154  5.4424  


# Deuterostomia+Latitude
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Deuterostomia")
data <- subset(data, Stressor == "W")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#              estimate      se     zval    pval     ci.lb    ci.ub     
# intrcpt     -11.6397  3.6992  -3.1466  0.0017  -18.8900  -4.3894  ** 
# logdist_sq   -1.1898  0.3951  -3.0116  0.0026   -1.9642  -0.4155  ** 
# logdist       7.4526  2.4445   3.0487  0.0023    2.6615  12.2436  ** 











