
library(metafor)    
library(multcomp)  
library(ape)        
library(phytools)   
library(dplyr)     
library(Matrix)        
library(parallel) 







# Calculate Blomberg's K
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$species <- tolower(gsub(" ", "_", data$species))
tree<-read.tree("phylogeny.tre")
vcov <- ape::vcv(tree)
species_avg_effect <- aggregate(yi ~ species, data = data, FUN = mean)
effect_sizes <- species_avg_effect$yi
names(effect_sizes) <- species_avg_effect$species

result_k <- phylosig(tree, effect_sizes, method = "K", test = TRUE, nsim = 10000)
print(result_k )
# Phylogenetic signal K : 0.0390648 
# P-value (based on 10000 randomizations) : 0.3481 



################################################
#                                              #
#        Fig. 2.d_Supplementary Table 1.       #
#                                              #
################################################
############## overall #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
# estimate      se     zval    pval    ci.lb    ci.ub     
# -0.0767  0.0283  -2.7125  0.0067  -0.1321  -0.0213  ** 



############## Stressor #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Stressor<-factor(data$Stressor)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
#               estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA     -0.0002  0.0347  -0.0058  0.9954  -0.0682   0.0678      
# StressorA-W   -0.2458  0.1176  -2.0900  0.0366  -0.4764  -0.0153    * 
# StressorH     -0.3229  0.1143  -2.8247  0.0047  -0.5470  -0.0989   ** 
# StressorW     -0.1480  0.0434  -3.4054  0.0007  -0.2331  -0.0628  *** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("A"=1,"A-W"=1,"H"=1,"W"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#               Estimate Std. Error z value Pr(>|z|)   
# A-W - A == 0 -0.24562    0.11909  -2.062  0.03916 * 
# H - A == 0   -0.32275    0.11948  -2.701  0.00691 **
# W - A == 0   -0.14776    0.05397  -2.738  0.00619 **
# H - A-W == 0 -0.07713    0.16403  -0.470  0.63821   
# W - A-W == 0  0.09786    0.12123   0.807  0.41952   
# W - H == 0    0.17499    0.12231   1.431  0.15250   



############## Taxa ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Taxa<-factor(data$Taxa)
rep_data <- as.data.frame(table(data$Taxa))
colnames(rep_data) <- c("Taxa","rep")
data <- merge(data, rep_data, by = c("Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                       estimate      se     zval    pval    ci.lb    ci.ub     
# TaxaAlgae+Protist     0.1826  0.0582   3.1372  0.0017   0.0685   0.2967  ** 
# TaxaDeuterostomia    -0.2047  0.0819  -2.4988  0.0125  -0.3652  -0.0441   * 
# TaxaEcdysozoa        -0.1937  0.0606  -3.1951  0.0014  -0.3126  -0.0749  ** 
# TaxaLophotrochozoa   -0.0653  0.0391  -1.6703  0.0949  -0.1419   0.0113   . 
# TaxaRadiata          -0.2512  0.0813  -3.0907  0.0020  -0.4105  -0.0919  ** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Algae+Protist"=1,"Deuterostomia"=1,"Ecdysozoa"=1,"Lophotrochozoa"=1,"Radiata"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                     Estimate Std. Error z value Pr(>|z|)    
# Deuterostomia - Algae+Protist == 0  -0.38727    0.10048  -3.854 0.000116 ***
# Ecdysozoa - Algae+Protist == 0      -0.37631    0.08402  -4.479 7.50e-06 ***
# Lophotrochozoa - Algae+Protist == 0 -0.24787    0.06868  -3.609 0.000307 ***
# Radiata - Algae+Protist == 0        -0.43375    0.09996  -4.339 1.43e-05 ***
# Ecdysozoa - Deuterostomia == 0       0.01096    0.10173   0.108 0.914221    
# Lophotrochozoa - Deuterostomia == 0  0.13940    0.09068   1.537 0.124210    
# Radiata - Deuterostomia == 0        -0.04648    0.11539  -0.403 0.687054    
# Lophotrochozoa - Ecdysozoa == 0      0.12845    0.07161   1.794 0.072866 .  
# Radiata - Ecdysozoa == 0            -0.05744    0.10139  -0.567 0.571036    
# Radiata - Lophotrochozoa == 0       -0.18589    0.09018  -2.061 0.039268 *  



############## Habitat ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Habitat<-factor(data$Habitat)
rep_data <- as.data.frame(table(data$Habitat))
colnames(rep_data) <- c("Habitat","rep")
data <- merge(data, rep_data, by = c("Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                    estimate      se     zval    pval    ci.lb    ci.ub     
# HabitatPolar       0.3924  0.1906   2.0590  0.0395   0.0189   0.7659   * 
# HabitatTemperate   -0.0498  0.0487  -1.0210  0.3073  -0.1453   0.0458     
# HabitatTropical    -0.1041  0.0341  -3.0527  0.0023  -0.1710  -0.0373  ** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Polar"=1,"Temperate"=1,"Tropical"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                           Estimate Std. Error z value Pr(>|z|)  
# Temperate - Polar == 0   -0.44215    0.19670  -2.248   0.0246 *
# Tropical - Polar == 0    -0.49651    0.19360  -2.565   0.0103 *
# Tropical - Temperate == 0 -0.05436    0.05947  -0.914   0.3607  


############## Trait ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Trait<-factor(data$Trait)
rep_data <- as.data.frame(table(data$Trait))
colnames(rep_data) <- c("Trait","rep")
data <- merge(data, rep_data, by = c("Trait"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trait - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                     estimate      se     zval    pval    ci.lb    ci.ub    
# TraitAbundance       0.0997  0.1152   0.8655  0.3868  -0.1261   0.3254    
# TraitAccumulation   -0.1013  0.0422  -2.4028  0.0163  -0.1840  -0.0187  * 
# TraitBehaviour      -0.1357  0.1355  -1.0014  0.3167  -0.4013   0.1299    
# TraitGrowth         -0.0502  0.0555  -0.9041  0.3659  -0.1589   0.0586    
# TraitMetabolism     -0.0632  0.0432  -1.4635  0.1433  -0.1479   0.0214    
# TraitReproduction   -0.1502  0.1014  -1.4812  0.1386  -0.3489   0.0485    
# TraitSurvival       -0.1177  0.0884  -1.3307  0.1833  -0.2909   0.0556    
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Abundance"=1,"Accumulation"=1,"Behaviour"=1,
                                     "Growth"=1,"Metabolism"=1,"Reproduction"=1,"Survival"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                  Estimate Std. Error z value Pr(>|z|)  
# Accumulation - Abundance == 0    -0.20102    0.12141  -1.656   0.0978 .
# Behaviour - Abundance == 0       -0.23539    0.17705  -1.330   0.1837  
# Growth - Abundance == 0          -0.14984    0.12374  -1.211   0.2259  
# Metabolism - Abundance == 0      -0.16291    0.12005  -1.357   0.1748  
# Reproduction - Abundance == 0    -0.24986    0.15210  -1.643   0.1004  
# Survival - Abundance == 0        -0.21734    0.14274  -1.523   0.1279  
# Behaviour - Accumulation == 0    -0.03437    0.13970  -0.246   0.8057  
# Growth - Accumulation == 0        0.05118    0.06536   0.783   0.4336  
# Metabolism - Accumulation == 0    0.03811    0.05541   0.688   0.4915  
# Reproduction - Accumulation == 0 -0.04884    0.10667  -0.458   0.6471  
# Survival - Accumulation == 0     -0.01632    0.09491  -0.172   0.8635  
# Growth - Behaviour == 0           0.08555    0.14243   0.601   0.5481  
# Metabolism - Behaviour == 0       0.07248    0.14089   0.514   0.6069  
# Reproduction - Behaviour == 0    -0.01447    0.16697  -0.087   0.9309  
# Survival - Behaviour == 0         0.01805    0.15782   0.114   0.9090  
# Metabolism - Growth == 0         -0.01307    0.06514  -0.201   0.8410  
# Reproduction - Growth == 0       -0.10002    0.10756  -0.930   0.3524  
# Survival - Growth == 0           -0.06750    0.10005  -0.675   0.4999  
# Reproduction - Metabolism == 0   -0.08695    0.10794  -0.806   0.4205  
# Survival - Metabolism == 0       -0.05443    0.09634  -0.565   0.5721  
# Survival - Reproduction == 0      0.03252    0.12771   0.255   0.7990  


############## Metal ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Metal<-factor(data$Metal)
rep_data <- as.data.frame(table(data$Metal))
colnames(rep_data) <- c("Metal","rep")
data <- merge(data, rep_data, by = c("Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Metal - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                  estimate      se     zval    pval    ci.lb    ci.ub     
# MetalArsenic      0.0938  0.1572   0.5970  0.5505  -0.2142   0.4018     
# MetalCadmium     -0.0537  0.0552  -0.9734  0.3303  -0.1619   0.0545     
# MetalCobalt      -0.0469  0.1214  -0.3859  0.6996  -0.2849   0.1912     
# MetalCopper      -0.1432  0.0458  -3.1252  0.0018  -0.2330  -0.0534  ** 
# MetalIron         0.2133  0.1210   1.7634  0.0778  -0.0238   0.4505   . 
# MetalLead        -0.3522  0.1086  -3.2417  0.0012  -0.5651  -0.1393  ** 
# MetalManganese   -0.1030  0.1413  -0.7290  0.4660  -0.3798   0.1739     
# MetalMercury     -0.0401  0.0759  -0.5287  0.5970  -0.1889   0.1086     
# MetalNickel      -0.0085  0.1780  -0.0478  0.9619  -0.3573   0.3403     
# MetalSilver      -0.1249  0.1190  -1.0495  0.2939  -0.3581   0.1083     
# MetalZinc         0.0328  0.0889   0.3695  0.7118  -0.1414   0.2071     
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Arsenic "=1,"Cadmium"=1,"Cobalt"=1,
                                     "Copper"=1,"Iron"=1,"Lead "=1,"Manganese"=1,
                                     "Mercury"=1,"Nickel"=1,"Silver"=1,"Zinc"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                           Estimate Std. Error z value Pr(>|z|)    
# Cadmium - Arsenic  == 0   -0.147550   0.165115  -0.894 0.371526    
# Cobalt - Arsenic  == 0    -0.140674   0.198492  -0.709 0.478501    
# Copper - Arsenic  == 0    -0.236996   0.163450  -1.450 0.147070    
# Iron - Arsenic  == 0       0.119523   0.198323   0.603 0.546729    
# Lead  - Arsenic  == 0     -0.446004   0.190214  -2.345 0.019040 *  
# Manganese - Arsenic  == 0 -0.196793   0.211171  -0.932 0.351383    
# Mercury - Arsenic  == 0   -0.133939   0.174103  -0.769 0.441711    
# Nickel - Arsenic  == 0    -0.102321   0.237402  -0.431 0.666466    
# Silver - Arsenic  == 0    -0.218686   0.196721  -1.112 0.266286    
# Zinc - Arsenic  == 0      -0.060965   0.179941  -0.339 0.734757    
# Cobalt - Cadmium == 0      0.006875   0.131994   0.052 0.958459    
# Copper - Cadmium == 0     -0.089446   0.070219  -1.274 0.202730    
# Iron - Cadmium == 0        0.267072   0.132968   2.009 0.044585 *  
# Lead  - Cadmium == 0      -0.298454   0.121001  -2.467 0.013643 *  
# Manganese - Cadmium == 0  -0.049243   0.150229  -0.328 0.743072    
# Mercury - Cadmium == 0     0.013611   0.093505   0.146 0.884267    
# Nickel - Cadmium == 0      0.045228   0.186304   0.243 0.808186    
# Silver - Cadmium == 0     -0.071136   0.124466  -0.572 0.567642    
# Zinc - Cadmium == 0        0.086585   0.103390   0.837 0.402337    
# Copper - Cobalt == 0      -0.096322   0.129533  -0.744 0.457115    
# Iron - Cobalt == 0         0.260197   0.171028   1.521 0.128166    
# Lead  - Cobalt == 0       -0.305329   0.162840  -1.875 0.060790 .  
# Manganese - Cobalt == 0   -0.056118   0.152039  -0.369 0.712050    
# Mercury - Cobalt == 0      0.006736   0.142440   0.047 0.962285    
# Nickel - Cobalt == 0       0.038353   0.215378   0.178 0.858665    
# Silver - Cobalt == 0      -0.078011   0.165537  -0.471 0.637454    
# Zinc - Cobalt == 0         0.079709   0.146905   0.543 0.587410    
# Iron - Copper == 0         0.356519   0.129352   2.756 0.005848 ** 
# Lead  - Copper == 0       -0.209008   0.116210  -1.799 0.072093 .  
# Manganese - Copper == 0    0.040203   0.148201   0.271 0.786180    
# Mercury - Copper == 0      0.103057   0.088512   1.164 0.244288    
# Nickel - Copper == 0       0.134675   0.183723   0.733 0.463538    
# Silver - Copper == 0       0.018310   0.126298   0.145 0.884729    
# Zinc - Copper == 0         0.176031   0.097413   1.807 0.070752 .  
# Lead  - Iron == 0         -0.565526   0.162597  -3.478 0.000505 ***
# Manganese - Iron == 0     -0.316316   0.185824  -1.702 0.088711 .  
# Mercury - Iron == 0       -0.253462   0.142808  -1.775 0.075923 .  
# Nickel - Iron == 0        -0.221844   0.215180  -1.031 0.302555    
# Silver - Iron == 0        -0.338208   0.169644  -1.994 0.046192 *  
# Zinc - Iron == 0          -0.180488   0.149795  -1.205 0.228244    
# Manganese - Lead  == 0     0.249211   0.178080   1.399 0.161683    
# Mercury - Lead  == 0       0.312065   0.132228   2.360 0.018273 *  
# Nickel - Lead  == 0        0.343682   0.208482   1.649 0.099250 .  
# Silver - Lead  == 0        0.227318   0.160818   1.414 0.157504    
# Zinc - Lead  == 0          0.385039   0.139382   2.762 0.005736 ** 
# Mercury - Manganese == 0   0.062854   0.159242   0.395 0.693060    
# Nickel - Manganese == 0    0.094472   0.227119   0.416 0.677442    
# Silver - Manganese == 0   -0.021893   0.179500  -0.122 0.902926    
# Zinc - Manganese == 0      0.135828   0.162550   0.836 0.403377    
# Nickel - Mercury == 0      0.031618   0.193450   0.163 0.870171    
# Silver - Mercury == 0     -0.084747   0.140094  -0.605 0.545226    
# Zinc - Mercury == 0        0.072974   0.116004   0.629 0.529309    
# Silver - Nickel == 0      -0.116364   0.213999  -0.544 0.586606    
# Zinc - Nickel == 0         0.041356   0.197861   0.209 0.834435    
# Zinc - Silver == 0         0.157721   0.144856   1.089 0.276237    


################################################
#                                              #
#                   Fig. 4                     #
#                                              #
################################################
############## Stressor:Taxa_Supplementary Table 5 #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Taxa<-factor(data$Taxa)
data$Stressor<-factor(data$Stressor)
rep_data <- as.data.frame(table(data$Stressor,data$Taxa))
colnames(rep_data) <- c("Stressor","Taxa","rep")
data <- merge(data, rep_data, by = c("Stressor","Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                   estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:TaxaAlgae+Protist       0.1774  0.0602   2.9475  0.0032   0.0594   0.2953   ** 
# StressorA-W:TaxaAlgae+Protist    -0.1083  0.2762  -0.3922  0.6949  -0.6497   0.4330      
# StressorW:TaxaAlgae+Protist       0.2462  0.1296   1.8993  0.0575  -0.0079   0.5003    . 
# StressorA:TaxaDeuterostomia      -0.1143  0.1011  -1.1305  0.2583  -0.3126   0.0839      
# StressorA-W:TaxaDeuterostomia    -0.1913  0.2043  -0.9362  0.3492  -0.5916   0.2091      
# StressorW:TaxaDeuterostomia      -0.3381  0.1185  -2.8524  0.0043  -0.5705  -0.1058   ** 
# StressorA:TaxaEcdysozoa           0.0555  0.0874   0.6352  0.5253  -0.1158   0.2269      
# StressorH:TaxaEcdysozoa          -0.4381  0.1861  -2.3540  0.0186  -0.8029  -0.0733    * 
# StressorW:TaxaEcdysozoa          -0.3307  0.0817  -4.0501  <.0001  -0.4908  -0.1707  *** 
# StressorA:TaxaLophotrochozoa     -0.0530  0.0464  -1.1412  0.2538  -0.1439   0.0380      
# StressorA-W:TaxaLophotrochozoa   -0.2165  0.1791  -1.2092  0.2266  -0.5675   0.1344      
# StressorH:TaxaLophotrochozoa     -0.2951  0.1404  -2.1019  0.0356  -0.5703  -0.0199    * 
# StressorW:TaxaLophotrochozoa     -0.0344  0.0620  -0.5545  0.5792  -0.1559   0.0871      
# StressorA:TaxaRadiata            -0.1731  0.1188  -1.4572  0.1451  -0.4060   0.0597      
# StressorW:TaxaRadiata            -0.3117  0.1028  -3.0329  0.0024  -0.5131  -0.1103   ** 




############## Stressor:Habitat_Supplementary Table 6 #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Habitat<-factor(data$Habitat)
data$Stressor<-factor(data$Stressor)
rep_data <- as.data.frame(table(data$Stressor,data$Habitat))
colnames(rep_data) <- c("Stressor","Habitat","rep")
data <- merge(data, rep_data, by = c("Stressor","Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:HabitatPolar         0.5137  0.1983   2.5899  0.0096   0.1249   0.9025   ** 
# StressorA:HabitatTemperate     -0.0325  0.0576  -0.5654  0.5718  -0.1453   0.0803      
# StressorA-W:HabitatTemperate   -0.0446  0.2130  -0.2095  0.8341  -0.4622   0.3729      
# StressorH:HabitatTemperate     -0.4022  0.1148  -3.5027  0.0005  -0.6273  -0.1771  *** 
# StressorW:HabitatTemperate      0.0915  0.0736   1.2429  0.2139  -0.0528   0.2357      
# StressorA:HabitatTropical       0.0217  0.0381   0.5693  0.5692  -0.0530   0.0964      
# StressorA-W:HabitatTropical    -0.3280  0.1359  -2.4129  0.0158  -0.5944  -0.0616    * 
# StressorH:HabitatTropical       0.2278  0.2993   0.7611  0.4466  -0.3588   0.8144      
# StressorW:HabitatTropical      -0.3018  0.0495  -6.1025  <.0001  -0.3987  -0.2049  *** 

############## Stressor:Habitat:Taxa_Supplementary Table 7 #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Habitat<-factor(data$Habitat)
data$Stressor<-factor(data$Stressor)
data$Taxa<-factor(data$Taxa)
rep_data <- as.data.frame(table(data$Stressor,data$Habitat,data$Taxa))
colnames(rep_data) <- c("Stressor","Habitat","Taxa","rep")
data <- merge(data, rep_data, by = c("Stressor","Habitat","Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Habitat:Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                  estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:HabitatPolar:TaxaAlgae+Protist          0.7031  0.2649   2.6545  0.0079   0.1840   1.2222   ** 
# StressorA:HabitatTemperate:TaxaAlgae+Protist       0.1989  0.1203   1.6532  0.0983  -0.0369   0.4346    . 
# StressorW:HabitatTemperate:TaxaAlgae+Protist       0.2731  0.2057   1.3276  0.1843  -0.1301   0.6762      
# StressorA:HabitatTropical:TaxaAlgae+Protist        0.1706  0.0679   2.5114  0.0120   0.0375   0.3038    * 
# StressorA-W:HabitatTropical:TaxaAlgae+Protist     -0.1827  0.2742  -0.6664  0.5052  -0.7201   0.3547      
# StressorW:HabitatTropical:TaxaAlgae+Protist       -0.0399  0.1979  -0.2019  0.8400  -0.4277   0.3478      
# StressorA:HabitatTemperate:TaxaDeuterostomia      -0.3052  0.1384  -2.2062  0.0274  -0.5764  -0.0341    * 
# StressorA:HabitatTropical:TaxaDeuterostomia        0.0791  0.1379   0.5738  0.5661  -0.1911   0.3493      
# StressorA-W:HabitatTropical:TaxaDeuterostomia     -0.1213  0.2030  -0.5978  0.5500  -0.5192   0.2765      
# StressorW:HabitatTropical:TaxaDeuterostomia       -0.3150  0.1163  -2.7095  0.0067  -0.5428  -0.0871   ** 
# StressorA:HabitatTemperate:TaxaEcdysozoa          -0.1835  0.2350  -0.7809  0.4348  -0.6442   0.2771      
# StressorH:HabitatTemperate:TaxaEcdysozoa          -0.4306  0.1812  -2.3763  0.0175  -0.7858  -0.0754    * 
# StressorW:HabitatTemperate:TaxaEcdysozoa          -0.1859  0.2105  -0.8831  0.3772  -0.5984   0.2266      
# StressorA:HabitatTropical:TaxaEcdysozoa            0.0949  0.0897   1.0577  0.2902  -0.0810   0.2708      
# StressorW:HabitatTropical:TaxaEcdysozoa           -0.3498  0.0847  -4.1323  <.0001  -0.5157  -0.1839  *** 
# StressorW:HabitatPolar:TaxaLophotrochozoa         0.3080  0.2761   1.1154  0.2647  -0.2332   0.8491      
# StressorA:HabitatTemperate:TaxaLophotrochozoa     -0.0263  0.0728  -0.3618  0.7175  -0.1690   0.1164      
# StressorA-W:HabitatTemperate:TaxaLophotrochozoa   -0.0410  0.2111  -0.1944  0.8459  -0.4549   0.3728      
# StressorH:HabitatTemperate:TaxaLophotrochozoa     -0.4451  0.1550  -2.8721  0.0041  -0.7488  -0.1414   ** 
# StressorW:HabitatTemperate:TaxaLophotrochozoa      0.0909  0.0789   1.1511  0.2497  -0.0639   0.2456      
# StressorA:HabitatTropical:TaxaLophotrochozoa      -0.0632  0.0566  -1.1158  0.2645  -0.1742   0.0478      
# StressorA-W:HabitatTropical:TaxaLophotrochozoa    -0.5956  0.3300  -1.8050  0.0711  -1.2423   0.0511    . 
# StressorH:HabitatTropical:TaxaLophotrochozoa       0.2289  0.2931   0.7809  0.4349  -0.3456   0.8034      
# StressorW:HabitatTropical:TaxaLophotrochozoa      -0.2698  0.1003  -2.6909  0.0071  -0.4663  -0.0733   ** 
# StressorA:HabitatTropical:TaxaRadiata             -0.1720  0.1159  -1.4838  0.1379  -0.3993   0.0552      
# StressorW:HabitatTropical:TaxaRadiata             -0.3136  0.1001  -3.1336  0.0017  -0.5098  -0.1175   ** 

################################################
#                                              #
#                    Fig. 5                    #
#                                              #
################################################
####### Trophic.level #########
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Trophic.level<-factor(data$Trophic.level)
rep_data <- as.data.frame(table(data$Trophic.level))
colnames(rep_data) <- c("Trophic.level","rep")
data <- merge(data, rep_data, by = c("Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                  estimate      se     zval    pval    ci.lb    ci.ub      
# Trophic.levelHerbivore          -0.0966  0.0358  -2.6973  0.0070  -0.1667  -0.0264   ** 
# Trophic.levelPredator           -0.1799  0.0508  -3.5399  0.0004  -0.2795  -0.0803  *** 
# Trophic.levelPrimary producer    0.1176  0.0567   2.0749  0.0380   0.0065   0.2287    * 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                   Estimate Std. Error z value Pr(>|z|)    
# Predator - Herbivore == 0         -0.08333    0.06158  -1.353  0.17596    
# Primary producer - Herbivore == 0  0.21419    0.06583   3.253  0.00114 ** 
# Primary producer - Predator == 0   0.29752    0.07468   3.984 6.78e-05 ***

######### Stressor:Trophic.level_Supplementary Table 9 #######
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Trophic.level<-factor(data$Trophic.level)
data$Stressor<-factor(data$Stressor)
rep_data <- as.data.frame(table(data$Stressor,data$Trophic.level))
colnames(rep_data) <- c("Stressor","Trophic.level","rep")
data <- merge(data, rep_data, by = c("Stressor","Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trophic.level- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                            estimate      se     zval    pval    ci.lb    ci.ub     
# StressorA:Trophic.levelHerbivore          -0.0136  0.0475  -0.2870  0.7741  -0.1068   0.0795     
# StressorA-W:Trophic.levelHerbivore        -0.2567  0.1888  -1.3595  0.1740  -0.6267   0.1134     
# StressorH:Trophic.levelHerbivore          -0.2578  0.1407  -1.8328  0.0668  -0.5335   0.0179   . 
# StressorW:Trophic.levelHerbivore          -0.1764  0.0552  -3.1933  0.0014  -0.2847  -0.0681  ** 
# StressorA:Trophic.levelPredator           -0.1081  0.0668  -1.6195  0.1053  -0.2390   0.0227     
# StressorA-W:Trophic.levelPredator         -0.2968  0.1507  -1.9693  0.0489  -0.5923  -0.0014   * 
# StressorH:Trophic.levelPredator           -0.4360  0.1864  -2.3386  0.0194  -0.8014  -0.0706   * 
# StressorW:Trophic.levelPredator           -0.2081  0.0730  -2.8515  0.0044  -0.3512  -0.0651  ** 
# StressorA:Trophic.levelPrimary producer    0.1233  0.0627   1.9678  0.0491   0.0005   0.2462   * 
# StressorW:Trophic.levelPrimary producer    0.1277  0.1233   1.0350  0.3007  -0.1141   0.3694     


######### Trophic.level:Metal_Supplementary Table 10 #########
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Metal<-factor(data$Metal)
data$Trophic.level<-factor(data$Trophic.level)
rep_data <- as.data.frame(table(data$Trophic.level,data$Metal))
colnames(rep_data) <- c("Trophic.level","Metal","rep")
data <- merge(data, rep_data, by = c("Trophic.level","Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level:Metal- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                              estimate      se     zval    pval    ci.lb    ci.ub      
# Trophic.levelHerbivore:MetalArsenic           0.1827  0.1838   0.9943  0.3201  -0.1774   0.5429      
# Trophic.levelPrimary producer:MetalArsenic    0.2048  0.1914   1.0701  0.2846  -0.1703   0.5799      
# Trophic.levelHerbivore:MetalCadmium          -0.1478  0.0641  -2.3043  0.0212  -0.2735  -0.0221    * 
# Trophic.levelPredator:MetalCadmium           -0.0421  0.1253  -0.3362  0.7367  -0.2877   0.2034      
# Trophic.levelPrimary producer:MetalCadmium    0.3744  0.1421   2.6349  0.0084   0.0959   0.6529   ** 
# Trophic.levelHerbivore:MetalCobalt           -0.0423  0.1957  -0.2161  0.8289  -0.4258   0.3412      
# Trophic.levelPredator:MetalCobalt            -0.1027  0.1613  -0.6365  0.5244  -0.4189   0.2135      
# Trophic.levelPrimary producer:MetalCobalt     0.1614  0.3117   0.5179  0.6046  -0.4495   0.7723      
# Trophic.levelHerbivore:MetalCopper           -0.1719  0.0598  -2.8738  0.0041  -0.2891  -0.0546   ** 
# Trophic.levelPredator:MetalCopper            -0.2126  0.0840  -2.5324  0.0113  -0.3772  -0.0481    * 
# Trophic.levelPrimary producer:MetalCopper    -0.0218  0.0888  -0.2450  0.8064  -0.1958   0.1523      
# Trophic.levelPrimary producer:MetalIron       0.2064  0.1142   1.8079  0.0706  -0.0174   0.4302    . 
# Trophic.levelHerbivore:MetalLead             -0.1690  0.1164  -1.4522  0.1464  -0.3971   0.0591      
# Trophic.levelPredator:MetalLead              -1.1289  0.2437  -4.6320  <.0001  -1.6066  -0.6512  *** 
# Trophic.levelHerbivore:MetalManganese        -0.0084  0.2862  -0.0294  0.9766  -0.5693   0.5525      
# Trophic.levelPredator:MetalManganese         -0.1581  0.1638  -0.9654  0.3343  -0.4792   0.1629      
# Trophic.levelHerbivore:MetalMercury           0.0415  0.0909   0.4567  0.6479  -0.1366   0.2196      
# Trophic.levelPredator:MetalMercury           -0.1616  0.1154  -1.4000  0.1615  -0.3879   0.0646      
# Trophic.levelHerbivore:MetalNickel           -0.0060  0.1766  -0.0342  0.9728  -0.3522   0.3401      
# Trophic.levelHerbivore:MetalSilver           -0.1601  0.1443  -1.1092  0.2673  -0.4429   0.1228      
# Trophic.levelPredator:MetalSilver            -0.1366  0.2025  -0.6746  0.4999  -0.5335   0.2603      
# Trophic.levelHerbivore:MetalZinc              0.0712  0.1435   0.4959  0.6200  -0.2101   0.3524      
# Trophic.levelPredator:MetalZinc              -0.0907  0.1396  -0.6499  0.5157  -0.3643   0.1829      
# Trophic.levelPrimary producer:MetalZinc       0.1837  0.1744   1.0533  0.2922  -0.1581   0.5254      


############## Habitat:Trophic.level_Supplementary Table 11 #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Trophic.level<-factor(data$Trophic.level)
data$Habitat<-factor(data$Habitat)
rep_data <- as.data.frame(table(data$Habitat,data$Trophic.level))
colnames(rep_data) <- c("Habitat","Trophic.level","rep")
data <- merge(data, rep_data, by = c("Habitat","Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Habitat:Trophic.level- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                  estimate      se     zval    pval    ci.lb    ci.ub      
# HabitatPolar:Trophic.levelHerbivore              0.3100  0.3084   1.0053  0.3148  -0.2944   0.9145      
# HabitatTemperate:Trophic.levelHerbivore          -0.0889  0.0621  -1.4313  0.1523  -0.2107   0.0328      
# HabitatTropical:Trophic.levelHerbivore           -0.1094  0.0441  -2.4819  0.0131  -0.1958  -0.0230    * 
# HabitatTemperate:Trophic.levelPredator           -0.0957  0.0844  -1.1333  0.2571  -0.2611   0.0698      
# HabitatTropical:Trophic.levelPredator            -0.2299  0.0635  -3.6200  0.0003  -0.3544  -0.1054  *** 
# HabitatPolar:Trophic.levelPrimary producer       0.4416  0.2305   1.9159  0.0554  -0.0102   0.8934    . 
# HabitatTemperate:Trophic.levelPrimary producer    0.1905  0.1213   1.5701  0.1164  -0.0473   0.4283      
# HabitatTropical:Trophic.levelPrimary producer     0.0651  0.0667   0.9773  0.3284  -0.0655   0.1958      

################################################
#                                              #
#                    Fig. 6.c                  #
#                                              #
################################################
############## Stressor:Metal #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Metal<-factor(data$Metal)
data$Stressor<-factor(data$Stressor)
rep_data <- as.data.frame(table(data$Stressor,data$Metal))
colnames(rep_data) <- c("Stressor","Metal","rep")
data <- merge(data, rep_data, by = c("Stressor","Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Metal - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                             estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:MetalArsenic      0.3556  0.2171   1.6377  0.1015  -0.0700   0.7812      
# StressorW:MetalArsenic     -0.1309  0.2096  -0.6244  0.5324  -0.5418   0.2800      
# StressorA:MetalCadmium      0.0250  0.0666   0.3763  0.7067  -0.1054   0.1555      
# StressorA-W:MetalCadmium   -0.2179  0.2533  -0.8603  0.3896  -0.7143   0.2785      
# StressorH:MetalCadmium     -0.3532  0.1321  -2.6732  0.0075  -0.6122  -0.0942   ** 
# StressorW:MetalCadmium     -0.0304  0.1187  -0.2561  0.7978  -0.2631   0.2023      
# StressorA:MetalCobalt      -0.0769  0.1391  -0.5524  0.5807  -0.3495   0.1958      
# StressorW:MetalCobalt       0.2714  0.1975   1.3743  0.1694  -0.1156   0.6584      
# StressorA:MetalCopper      -0.1006  0.0574  -1.7526  0.0797  -0.2130   0.0119    . 
# StressorA-W:MetalCopper    -0.2851  0.2811  -1.0141  0.3105  -0.8361   0.2659      
# StressorW:MetalCopper      -0.2021  0.0691  -2.9241  0.0035  -0.3375  -0.0666   ** 
# StressorA:MetalIron         0.1436  0.1436   1.0004  0.3171  -0.1378   0.4250      
# StressorW:MetalIron         0.3362  0.1994   1.6863  0.0917  -0.0546   0.7271    . 
# StressorA:MetalLead        -0.1290  0.1787  -0.7218  0.4704  -0.4793   0.2213      
# StressorW:MetalLead        -0.4647  0.1316  -3.5323  0.0004  -0.7225  -0.2068  *** 
# StressorA:MetalManganese    0.1376  0.1814   0.7584  0.4482  -0.2180   0.4932      
# StressorW:MetalManganese    0.0872  0.2296   0.3797  0.7042  -0.3628   0.5372      
# StressorA:MetalMercury      0.1594  0.1009   1.5802  0.1141  -0.0383   0.3572      
# StressorA-W:MetalMercury   -0.1962  0.1801  -1.0892  0.2761  -0.5491   0.1568      
# StressorW:MetalMercury     -0.2004  0.0986  -2.0316  0.0422  -0.3937  -0.0071    * 
# StressorA:MetalNickel       0.1973  0.2481   0.7952  0.4265  -0.2890   0.6837      
# StressorW:MetalNickel      -0.1964  0.2337  -0.8402  0.4008  -0.6545   0.2617      
# StressorA:MetalSilver      -0.0980  0.1326  -0.7389  0.4600  -0.3578   0.1619      
# StressorW:MetalSilver      -0.0919  0.2534  -0.3625  0.7170  -0.5885   0.4048      
# StressorA:MetalZinc         0.0223  0.0962   0.2318  0.8167  -0.1662   0.2108      
# StressorW:MetalZinc         0.0067  0.2561   0.0262  0.9791  -0.4953   0.5087      



################################################
#                                              #
#            Extended Data Fig. 2              #
#                                              #
################################################
######### Stressor:Trait_Supplementary Table 2 #########
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Trait<-factor(data$Trait)
data$Stressor<-factor(data$Stressor)
rep_data <- as.data.frame(table(data$Stressor,data$Trait))
colnames(rep_data) <- c("Stressor","Trait","rep")
data <- merge(data, rep_data, by = c("Stressor","Trait"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trait - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:TraitAbundance         0.2255  0.1713   1.3167  0.1879  -0.1102   0.5612      
# StressorW:TraitAbundance         0.0081  0.1531   0.0529  0.9578  -0.2920   0.3083      
# StressorA:TraitAccumulation     -0.0073  0.0485  -0.1495  0.8812  -0.1023   0.0878      
# StressorA-W:TraitAccumulation   -0.2170  0.1663  -1.3049  0.1919  -0.5428   0.1089      
# StressorH:TraitAccumulation     -0.2367  0.1417  -1.6710  0.0947  -0.5144   0.0409    . 
# StressorW:TraitAccumulation     -0.2590  0.0763  -3.3951  0.0007  -0.4085  -0.1095  *** 
# StressorA:TraitBehaviour         0.0114  0.2004   0.0569  0.9546  -0.3814   0.4042      
# StressorW:TraitBehaviour        -0.1632  0.2028  -0.8045  0.4211  -0.5607   0.2343      
# StressorA:TraitGrowth            0.0161  0.0639   0.2511  0.8017  -0.1092   0.1413      
# StressorW:TraitGrowth           -0.1275  0.1014  -1.2569  0.2088  -0.3263   0.0713      
# StressorA:TraitMetabolism       -0.0473  0.0551  -0.8571  0.3914  -0.1554   0.0608      
# StressorA-W:TraitMetabolism     -0.5495  0.2231  -2.4626  0.0138  -0.9869  -0.1122    * 
# StressorH:TraitMetabolism       -0.4201  0.2116  -1.9853  0.0471  -0.8349  -0.0054    * 
# StressorW:TraitMetabolism       -0.0155  0.0674  -0.2299  0.8182  -0.1476   0.1166      
# StressorA:TraitReproduction      0.1545  0.1370   1.1277  0.2594  -0.1141   0.4231      
# StressorW:TraitReproduction     -0.4777  0.1459  -3.2741  0.0011  -0.7636  -0.1917   ** 
# StressorA:TraitSurvival          0.0627  0.1443   0.4343  0.6641  -0.2202   0.3456      
# StressorW:TraitSurvival         -0.2617  0.1117  -2.3441  0.0191  -0.4806  -0.0429    * 

############## Stressor:Trait:Habitat_Supplementary Table 3 ########## 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Habitat<-factor(data$Habitat)
data$Stressor<-factor(data$Stressor)
data$Trait<-factor(data$Trait)
rep_data <- as.data.frame(table(data$Stressor,data$Trait,data$Habitat))
colnames(rep_data) <- c("Stressor","Trait","Habitat","rep")
data <- merge(data, rep_data, by = c("Stressor","Trait","Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trait:Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                  estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:TraitMetabolism:HabitatPolar           0.6592  0.2395   2.7522  0.0059   0.1897   1.1286   ** 
# StressorA:TraitAbundance:HabitatTemperate         0.3739  0.2068   1.8084  0.0705  -0.0313   0.7792    . 
# StressorA:TraitAccumulation:HabitatTemperate     -0.0604  0.0782  -0.7732  0.4394  -0.2136   0.0927      
# StressorA-W:TraitAccumulation:HabitatTemperate   -0.1561  0.2669  -0.5848  0.5587  -0.6792   0.3670      
# StressorH:TraitAccumulation:HabitatTemperate     -0.3778  0.1415  -2.6698  0.0076  -0.6552  -0.1005   ** 
# StressorW:TraitAccumulation:HabitatTemperate      0.1832  0.1187   1.5434  0.1227  -0.0494   0.4157      
# StressorA:TraitBehaviour:HabitatTemperate        -0.1229  0.3146  -0.3905  0.6962  -0.7395   0.4938      
# StressorA:TraitGrowth:HabitatTemperate           -0.0233  0.1408  -0.1653  0.8687  -0.2993   0.2527      
# StressorW:TraitGrowth:HabitatTemperate           -0.2985  0.1938  -1.5402  0.1235  -0.6783   0.0813      
# StressorA:TraitMetabolism:HabitatTemperate        0.0155  0.1227   0.1262  0.8996  -0.2250   0.2559      
# StressorH:TraitMetabolism:HabitatTemperate       -0.3950  0.2476  -1.5948  0.1108  -0.8803   0.0904      
# StressorW:TraitMetabolism:HabitatTemperate        0.2229  0.1129   1.9747  0.0483   0.0017   0.4442    * 
# StressorW:TraitSurvival:HabitatTemperate         -0.2067  0.2275  -0.9085  0.3636  -0.6525   0.2392      
# StressorA:TraitAbundance:HabitatTropical         -0.1112  0.2866  -0.3880  0.6980  -0.6730   0.4506      
# StressorW:TraitAbundance:HabitatTropical         -0.1359  0.1559  -0.8719  0.3833  -0.4415   0.1696      
# StressorA:TraitAccumulation:HabitatTropical       0.0439  0.0582   0.7538  0.4509  -0.0702   0.1580      
# StressorA-W:TraitAccumulation:HabitatTropical    -0.2157  0.1998  -1.0794  0.2804  -0.6074   0.1760      
# StressorW:TraitAccumulation:HabitatTropical      -0.5635  0.0967  -5.8264  <.0001  -0.7530  -0.3739  *** 
# StressorA:TraitBehaviour:HabitatTropical          0.0789  0.2468   0.3197  0.7492  -0.4047   0.5625      
# StressorW:TraitBehaviour:HabitatTropical         -0.2099  0.2098  -1.0004  0.3171  -0.6211   0.2013      
# StressorA:TraitGrowth:HabitatTropical             0.0361  0.0693   0.5210  0.6023  -0.0998   0.1721      
# StressorW:TraitGrowth:HabitatTropical            -0.0905  0.1125  -0.8040  0.4214  -0.3110   0.1301      
# StressorA:TraitMetabolism:HabitatTropical        -0.0579  0.0595  -0.9737  0.3302  -0.1745   0.0587      
# StressorA-W:TraitMetabolism:HabitatTropical      -0.5889  0.2154  -2.7340  0.0063  -1.0110  -0.1667   ** 
# StressorW:TraitMetabolism:HabitatTropical        -0.2187  0.0838  -2.6110  0.0090  -0.3829  -0.0545   ** 
# StressorA:TraitReproduction:HabitatTropical       0.2162  0.1377   1.5702  0.1164  -0.0537   0.4861      
# StressorW:TraitReproduction:HabitatTropical      -0.5136  0.1407  -3.6508  0.0003  -0.7894  -0.2379  *** 
# StressorA:TraitSurvival:HabitatTropical           0.1178  0.1406   0.8383  0.4018  -0.1576   0.3933      
# StressorW:TraitSurvival:HabitatTropical          -0.3079  0.1214  -2.5356  0.0112  -0.5459  -0.0699    * 


################################################
#                                              #
#  Extended Data Fig. 3_Supplementary Table 8  #
#                                              #
################################################
###### Trophic.level-Accumulation ######
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data <- subset(data, Trait == "Accumulation")
nrow(data)
data$Trophic.level<-factor(data$Trophic.level)
rep_data <- as.data.frame(table(data$Trophic.level))
colnames(rep_data) <- c("Trophic.level","rep")
data <- merge(data, rep_data, by = c("Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                  estimate      se     zval    pval    ci.lb    ci.ub    
# Trophic.levelHerbivore          -0.0711  0.0613  -1.1598  0.2461  -0.1913   0.0491    
# Trophic.levelPredator           -0.2149  0.0854  -2.5175  0.0118  -0.3822  -0.0476  * 
# Trophic.levelPrimary producer    0.3402  0.1515   2.2458  0.0247   0.0433   0.6372  * 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                   Estimate Std. Error z value Pr(>|z|)   
# Predator - Herbivore == 0          -0.1438     0.1051  -1.368  0.17119   
# Primary producer - Herbivore == 0   0.4113     0.1553   2.648  0.00809 **
# Primary producer - Predator == 0    0.5552     0.1739   3.192  0.00141 **

############## Trophic.level-Abundance #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data <- subset(data, Trait == "Abundance")
nrow(data)
data$Trophic.level<-factor(data$Trophic.level)
rep_data <- as.data.frame(table(data$Trophic.level))
colnames(rep_data) <- c("Trophic.level","rep")
data <- merge(data, rep_data, by = c("Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Trophic.levelHerbivore          -0.1820  0.4154  -0.4381  0.6613  -0.9961  0.6321    
# Trophic.levelPredator           -0.0873  0.3618  -0.2412  0.8094  -0.7964  0.6219    
# Trophic.levelPrimary producer    0.2546  0.2869   0.8875  0.3748  -0.3077  0.8169    

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                    Estimate Std. Error z value Pr(>|z|)
# Predator - Herbivore == 0          0.09472    0.55085   0.172    0.863
# Primary producer - Herbivore == 0  0.43659    0.50479   0.865    0.387
# Primary producer - Predator == 0   0.34187    0.46175   0.740    0.459






