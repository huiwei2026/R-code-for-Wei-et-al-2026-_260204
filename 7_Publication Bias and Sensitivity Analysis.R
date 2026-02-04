##################### Checking publication Bias ###############
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
if (any(data$vi <= 0)) {
  data$vi[data$vi <= 0] <- 1e-6
  summary(data$vi)
  data[data$vi < 0.001, ]
  data <- data[data$vi >= 0.001, ]
}
rmod<-rma(yi,vi,data=data,method= "REML", subset = vi>0.00001)
funnel(rmod)
rmod.tf<-trimfill(rmod, ma.common= FALSE, estimator="R0")
funnel(rmod.tf)
fsn(yi, vi, data=data, type="Rosenthal", alpha=.05)
# Fail-safe N Calculation Using the Rosenthal Approach
# 
# Observed Significance Level: <.0001
# Target Significance Level:   0.05
# 
# Fail-safe N: 216952


### A. Removing one study that accounts for the largest number of experiment ###
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Paper_id<-factor(data$Paper_id)
summary(data$Paper_id)
# 134     113      77      21      89     125     127     139      58      62 
# 24      20      19      15      14      13      13      13      12      11 
# 121      39      80      41     118     143      25      61      67      69 
# 11      10      10       9       9       9       8       8       8       8 
# 141      37      38      43      55      59      66      81     102     108 
# 7       6       6       6       6       6       6       6       6       6 
# 117     133     147      12      13      86     100     122     124       1 
# 6       6       6       5       5       5       5       5       5       4 
# 5      14      17      19      22      24      28      40      42      53 
# 4       4       4       4       4       4       4       4       4       4 
# 60      83      84      88     105     106     107     128     135     136 
# 4       4       4       4       4       4       4       4       4       4 
# 140     145     148       6      29      30      31      33      35      45 
# 4       4       4       3       3       3       3       3       3       3 
# 47      51      71      73      76      78      79      90      95     109 
# 3       3       3       3       3       3       3       3       3       3 
# 112     129     142     146       3       4       8       9      10      15 
# 3       3       3       3       2       2       2       2       2       2 
# 16      18      20      26      32      34      36      44      46 (Other) 
# 2       2       2       2       2       2       2       2       2      71 
############## overall #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
# estimate      se     zval    pval    ci.lb    ci.ub     
# -0.0821  0.0279  -2.9465  0.0032  -0.1367  -0.0275  ** 

############## Stressor #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
nrow(data)
data$Stressor<-factor(data$Stressor)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
# estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA     -0.0080  0.0345  -0.2307  0.8175  -0.0756   0.0597      
# StressorA-W   -0.2487  0.1197  -2.0778  0.0377  -0.4833  -0.0141    * 
# StressorH     -0.3221  0.1139  -2.8278  0.0047  -0.5453  -0.0988   ** 
# StressorW     -0.1504  0.0430  -3.5001  0.0005  -0.2346  -0.0662  *** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("A"=1,"A-W"=1,"H"=1,"W"=1), 
                                   type="Tukey"))), test=adjusted("none"))
# Estimate Std. Error z value Pr(>|z|)   
# A-W - A == 0 -0.24073    0.12146  -1.982  0.04748 * 
# H - A == 0   -0.31410    0.11901  -2.639  0.00831 **
# W - A == 0   -0.14242    0.05374  -2.650  0.00805 **
# H - A-W == 0 -0.07337    0.16522  -0.444  0.65698   
# W - A-W == 0  0.09831    0.12356   0.796  0.42622   
# W - H == 0    0.17169    0.12173   1.410  0.15842   



############## Taxa ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
# estimate      se     zval    pval    ci.lb    ci.ub     
# TaxaAlgae+Protist     0.1967  0.0622   3.1612  0.0016   0.0747   0.3186  ** 
# TaxaDeuterostomia    -0.2026  0.0821  -2.4679  0.0136  -0.3634  -0.0417   * 
# TaxaEcdysozoa        -0.1925  0.0606  -3.1749  0.0015  -0.3114  -0.0737  ** 
# TaxaLophotrochozoa   -0.0802  0.0398  -2.0129  0.0441  -0.1582  -0.0021   * 
# TaxaRadiata          -0.2517  0.0818  -3.0784  0.0021  -0.4120  -0.0915  ** 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Algae+Protist"=1,"Deuterostomia"=1,"Ecdysozoa"=1,"Lophotrochozoa"=1,"Radiata"=1), 
                                   type="Tukey"))), test=adjusted("none"))
# Estimate Std. Error z value Pr(>|z|)    
# Deuterostomia - Algae+Protist == 0  -0.39923    0.10299  -3.876 0.000106 ***
# Ecdysozoa - Algae+Protist == 0      -0.38920    0.08688  -4.480 7.47e-06 ***
# Lophotrochozoa - Algae+Protist == 0 -0.27683    0.07387  -3.748 0.000178 ***
# Radiata - Algae+Protist == 0        -0.44839    0.10275  -4.364 1.28e-05 ***
# Ecdysozoa - Deuterostomia == 0       0.01004    0.10189   0.098 0.921541    
# Lophotrochozoa - Deuterostomia == 0  0.12240    0.09115   1.343 0.179325    
# Radiata - Deuterostomia == 0        -0.04916    0.11586  -0.424 0.671362    
# Lophotrochozoa - Ecdysozoa == 0      0.11237    0.07204   1.560 0.118819    
# Radiata - Ecdysozoa == 0            -0.05919    0.10180  -0.581 0.560938    
# Radiata - Lophotrochozoa == 0       -0.17156    0.09095  -1.886 0.059251 .  



############## Habitat ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
# HabitatFrigid       0.3981  0.1863   2.1367  0.0326   0.0329   0.7632    * 
# HabitatTemperate   -0.0482  0.0473  -1.0176  0.3089  -0.1409   0.0446      
# HabitatTropical    -0.1139  0.0335  -3.4029  0.0007  -0.1795  -0.0483  *** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Frigid"=1,"Temperate"=1,"Tropical"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                           Estimate Std. Error z value Pr(>|z|)   
# Temperate - Frigid == 0   -0.44621    0.19221  -2.321  0.02026 * 
# Tropical - Frigid == 0    -0.51195    0.18927  -2.705  0.00683 **
# Tropical - Temperate == 0 -0.06574    0.05794  -1.135  0.25657   


############## Trait ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                      estimate      se     zval    pval    ci.lb    ci.ub     
# TraitAbundance       0.0975  0.1161   0.8397  0.4011  -0.1301   0.3251     
# TraitAccumulation   -0.1209  0.0430  -2.8085  0.0050  -0.2052  -0.0365  ** 
# TraitBehaviour      -0.1403  0.1373  -1.0217  0.3069  -0.4094   0.1288     
# TraitGrowth         -0.0514  0.0578  -0.8901  0.3734  -0.1646   0.0618     
# TraitMetabolism     -0.0563  0.0441  -1.2769  0.2016  -0.1428   0.0301     
# TraitReproduction   -0.1521  0.1032  -1.4745  0.1404  -0.3543   0.0501     
# TraitSurvival       -0.1232  0.0892  -1.3809  0.1673  -0.2981   0.0517    

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Abundance"=1,"Accumulation"=1,"Behaviour"=1,
                                     "Growth"=1,"Metabolism"=1,"Reproduction"=1,"Survival"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                  Estimate Std. Error z value Pr(>|z|)  
# Accumulation - Abundance == 0    -0.218378   0.122886  -1.777   0.0756 .
# Behaviour - Abundance == 0       -0.237805   0.179181  -1.327   0.1844  
# Growth - Abundance == 0          -0.148916   0.125835  -1.183   0.2366  
# Metabolism - Abundance == 0      -0.153833   0.121601  -1.265   0.2058  
# Reproduction - Abundance == 0    -0.249614   0.154202  -1.619   0.1055  
# Survival - Abundance == 0        -0.220701   0.144263  -1.530   0.1261  
# Behaviour - Accumulation == 0    -0.019427   0.141996  -0.137   0.8912  
# Growth - Accumulation == 0        0.069462   0.069129   1.005   0.3150  
# Metabolism - Accumulation == 0    0.064546   0.057898   1.115   0.2649  
# Reproduction - Accumulation == 0 -0.031236   0.109121  -0.286   0.7747  
# Survival - Accumulation == 0     -0.002323   0.096434  -0.024   0.9808  
# Growth - Behaviour == 0           0.088889   0.145237   0.612   0.5405  
# Metabolism - Behaviour == 0       0.083973   0.143138   0.587   0.5574  
# Reproduction - Behaviour == 0    -0.011809   0.169812  -0.070   0.9446  
# Survival - Behaviour == 0         0.017104   0.160300   0.107   0.9150  
# Metabolism - Growth == 0         -0.004917   0.068571  -0.072   0.9428  
# Reproduction - Growth == 0       -0.100698   0.110641  -0.910   0.3628  
# Survival - Growth == 0           -0.071785   0.102400  -0.701   0.4833  
# Reproduction - Metabolism == 0   -0.095781   0.110348  -0.868   0.3854  
# Survival - Metabolism == 0       -0.066868   0.097826  -0.684   0.4943  
# Survival - Reproduction == 0      0.028913   0.130253   0.222   0.8243  


############## Metal####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
# MetalArsenic     -0.1261  0.2163  -0.5831  0.5598  -0.5501   0.2978     
# MetalCadmium     -0.0576  0.0553  -1.0416  0.2976  -0.1659   0.0507     
# MetalCobalt      -0.0414  0.1222  -0.3391  0.7346  -0.2809   0.1980     
# MetalCopper      -0.1441  0.0458  -3.1482  0.0016  -0.2338  -0.0544  ** 
# MetalIron         0.2097  0.1201   1.7455  0.0809  -0.0258   0.4452   . 
# MetalLead        -0.3503  0.1092  -3.2075  0.0013  -0.5643  -0.1362  ** 
# MetalManganese   -0.0986  0.1430  -0.6894  0.4906  -0.3788   0.1816     
# MetalMercury     -0.0384  0.0753  -0.5106  0.6096  -0.1860   0.1091     
# MetalNickel      -0.0125  0.1756  -0.0709  0.9434  -0.3566   0.3317     
# MetalSilver      -0.1234  0.1207  -1.0224  0.3066  -0.3599   0.1131     
# MetalZinc         0.0344  0.0898   0.3825  0.7021  -0.1417   0.2104     
  
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Arsenic "=1,"Cadmium"=1,"Cobalt"=1,
                                     "Copper"=1,"Iron"=1,"Lead "=1,"Manganese"=1,
                                     "Mercury"=1,"Nickel"=1,"Silver"=1,"Zinc"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                            Estimate Std. Error z value Pr(>|z|)    
# Cadmium - Arsenic  == 0    0.068585   0.221469   0.310 0.756801    
# Cobalt - Arsenic  == 0     0.084711   0.248307   0.341 0.732987    
# Copper - Arsenic  == 0    -0.017951   0.220815  -0.081 0.935208    
# Iron - Arsenic  == 0       0.335867   0.247442   1.357 0.174668    
# Lead  - Arsenic  == 0     -0.224148   0.241212  -0.929 0.352756    
# Manganese - Arsenic  == 0  0.027587   0.259130   0.106 0.915216    
# Mercury - Arsenic  == 0    0.087699   0.228522   0.384 0.701152    
# Nickel - Arsenic  == 0     0.113685   0.278594   0.408 0.683226    
# Silver - Arsenic  == 0     0.002762   0.247214   0.011 0.991085    
# Zinc - Arsenic  == 0       0.160494   0.233425   0.688 0.491730    
# Cobalt - Cadmium == 0      0.016126   0.132882   0.121 0.903409    
# Copper - Cadmium == 0     -0.086536   0.070409  -1.229 0.219053    
# Iron - Cadmium == 0        0.267282   0.132238   2.021 0.043257 *  
# Lead  - Cadmium == 0      -0.292733   0.121624  -2.407 0.016090 *  
# Manganese - Cadmium == 0  -0.040998   0.152020  -0.270 0.787400    
# Mercury - Cadmium == 0     0.019113   0.093085   0.205 0.837312    
# Nickel - Cadmium == 0      0.045099   0.184058   0.245 0.806435    
# Silver - Cadmium == 0     -0.065823   0.126468  -0.520 0.602734    
# Zinc - Cadmium == 0        0.091908   0.104351   0.881 0.378448    
# Copper - Cobalt == 0      -0.102662   0.130264  -0.788 0.430631    
# Iron - Cobalt == 0         0.251156   0.171021   1.469 0.141950    
# Lead  - Cobalt == 0       -0.308859   0.163796  -1.886 0.059344 .  
# Manganese - Cobalt == 0   -0.057124   0.155591  -0.367 0.713513    
# Mercury - Cobalt == 0      0.002987   0.142849   0.021 0.983315    
# Nickel - Cobalt == 0       0.028973   0.213862   0.135 0.892235    
# Silver - Cobalt == 0      -0.081949   0.167707  -0.489 0.625093    
# Zinc - Cobalt == 0         0.075782   0.148446   0.511 0.609698    
# Iron - Copper == 0         0.353818   0.128561   2.752 0.005921 ** 
# Lead  - Copper == 0       -0.206197   0.116891  -1.764 0.077730 .  
# Manganese - Copper == 0    0.045538   0.149867   0.304 0.761236    
# Mercury - Copper == 0      0.105650   0.087986   1.201 0.229846    
# Nickel - Copper == 0       0.131636   0.181421   0.726 0.468095    
# Silver - Copper == 0       0.020713   0.128023   0.162 0.871469    
# Zinc - Copper == 0         0.178445   0.098447   1.813 0.069894 .  
# Lead  - Iron == 0         -0.560015   0.162361  -3.449 0.000562 ***
# Manganese - Iron == 0     -0.308280   0.186616  -1.652 0.098545 .  
# Mercury - Iron == 0       -0.248168   0.141780  -1.750 0.080053 .  
# Nickel - Iron == 0        -0.222182   0.212751  -1.044 0.296332    
# Silver - Iron == 0        -0.333105   0.170256  -1.956 0.050407 .  
# Zinc - Iron == 0          -0.175373   0.149711  -1.171 0.241433    
# Manganese - Lead  == 0     0.251735   0.179806   1.400 0.161502    
# Mercury - Lead  == 0       0.311847   0.132380   2.356 0.018488 *  
# Nickel - Lead  == 0        0.337832   0.206761   1.634 0.102274    
# Silver - Lead  == 0        0.226910   0.162504   1.396 0.162614    
# Zinc - Lead  == 0          0.384641   0.140512   2.737 0.006192 ** 
# Mercury - Manganese == 0   0.060111   0.160585   0.374 0.708161    
# Nickel - Manganese == 0    0.086097   0.226359   0.380 0.703680    
# Silver - Manganese == 0   -0.024825   0.182357  -0.136 0.891715    
# Zinc - Manganese == 0      0.132906   0.164922   0.806 0.420316    
# Nickel - Mercury == 0      0.025986   0.191030   0.136 0.891797    
# Silver - Mercury == 0     -0.084937   0.141311  -0.601 0.547800    
# Zinc - Mercury == 0        0.072795   0.116403   0.625 0.531728    
# Silver - Nickel == 0      -0.110922   0.212999  -0.521 0.602531    
# Zinc - Nickel == 0         0.046809   0.196306   0.238 0.811533    
# Zinc - Silver == 0         0.157731   0.147112   1.072 0.283637    




############## Trophic.level #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# Trophic.levelHerbivore          -0.1077  0.0361  -2.9811  0.0029  -0.1785  -0.0369   ** 
# Trophic.levelPredator           -0.1787  0.0507  -3.5244  0.0004  -0.2780  -0.0793  *** 
# Trophic.levelPrimary producer    0.1274  0.0601   2.1192  0.0341   0.0096   0.2452    * 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                    Estimate Std. Error z value Pr(>|z|)    
# Predator - Herbivore == 0         -0.07096    0.06176  -1.149 0.250532    
# Primary producer - Herbivore == 0  0.23508    0.07011   3.353 0.000799 ***
# Primary producer - Predator == 0   0.30605    0.07727   3.961 7.47e-05 ***








############## Stressor:Taxa #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
# StressorA:TaxaAlgae+Protist       0.1856  0.0649   2.8595  0.0042   0.0584   0.3127   ** 
# StressorA-W:TaxaAlgae+Protist    -0.1204  0.2794  -0.4310  0.6665  -0.6680   0.4272      
# StressorW:TaxaAlgae+Protist       0.2450  0.1295   1.8917  0.0585  -0.0088   0.4987    . 
# StressorA:TaxaDeuterostomia      -0.1118  0.1011  -1.1058  0.2688  -0.3101   0.0864      
# StressorA-W:TaxaDeuterostomia    -0.1786  0.2069  -0.8633  0.3880  -0.5840   0.2269      
# StressorW:TaxaDeuterostomia      -0.3355  0.1192  -2.8141  0.0049  -0.5692  -0.1018   ** 
# StressorA:TaxaEcdysozoa           0.0583  0.0866   0.6732  0.5008  -0.1115   0.2281      
# StressorH:TaxaEcdysozoa          -0.4313  0.1868  -2.3091  0.0209  -0.7974  -0.0652    * 
# StressorW:TaxaEcdysozoa          -0.3316  0.0811  -4.0871  <.0001  -0.4906  -0.1726  *** 
# StressorA:TaxaLophotrochozoa     -0.0733  0.0473  -1.5506  0.1210  -0.1660   0.0194      
# StressorA-W:TaxaLophotrochozoa   -0.2241  0.1821  -1.2308  0.2184  -0.5810   0.1328      
# StressorH:TaxaLophotrochozoa     -0.2977  0.1410  -2.1117  0.0347  -0.5741  -0.0214    * 
# StressorW:TaxaLophotrochozoa     -0.0361  0.0622  -0.5816  0.5609  -0.1580   0.0857      
# StressorA:TaxaRadiata            -0.1723  0.1195  -1.4416  0.1494  -0.4065   0.0619      
# StressorW:TaxaRadiata            -0.3133  0.1028  -3.0460  0.0023  -0.5149  -0.1117   ** 




############## Stressor:Habitat #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:HabitatFrigid         0.5166  0.1936   2.6688  0.0076   0.1372   0.8960   ** 
# StressorA:HabitatTemperate     -0.0309  0.0557  -0.5543  0.5794  -0.1401   0.0784      
# StressorA-W:HabitatTemperate   -0.0418  0.2171  -0.1923  0.8475  -0.4674   0.3838      
# StressorH:HabitatTemperate     -0.3991  0.1136  -3.5116  0.0004  -0.6218  -0.1763  *** 
# StressorW:HabitatTemperate      0.0852  0.0722   1.1815  0.2374  -0.0562   0.2267      
# StressorA:HabitatTropical       0.0082  0.0378   0.2173  0.8280  -0.0659   0.0824      
# StressorA-W:HabitatTropical    -0.3263  0.1373  -2.3770  0.0175  -0.5954  -0.0573    * 
# StressorH:HabitatTropical       0.2272  0.2997   0.7581  0.4484  -0.3602   0.8146      
# StressorW:HabitatTropical      -0.3030  0.0484  -6.2562  <.0001  -0.3980  -0.2081  *** 

############## Stressor:Habitat:Taxa#################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                                                   estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:HabitatFrigid:TaxaAlgae+Protist          0.7032  0.2610   2.6938  0.0071   0.1916   1.2149   ** 
# StressorA:HabitatTemperate:TaxaAlgae+Protist       0.1856  0.1176   1.5782  0.1145  -0.0449   0.4161      
# StressorW:HabitatTemperate:TaxaAlgae+Protist       0.2599  0.2016   1.2891  0.1974  -0.1353   0.6551      
# StressorA:HabitatTropical:TaxaAlgae+Protist        0.1791  0.0749   2.3913  0.0168   0.0323   0.3260    * 
# StressorA-W:HabitatTropical:TaxaAlgae+Protist     -0.1923  0.2753  -0.6987  0.4848  -0.7319   0.3472      
# StressorW:HabitatTropical:TaxaAlgae+Protist       -0.0409  0.1971  -0.2075  0.8356  -0.4272   0.3454      
# StressorA:HabitatTemperate:TaxaDeuterostomia      -0.2972  0.1361  -2.1840  0.0290  -0.5639  -0.0305    * 
# StressorA:HabitatTropical:TaxaDeuterostomia        0.0841  0.1374   0.6123  0.5403  -0.1852   0.3535      
# StressorA-W:HabitatTropical:TaxaDeuterostomia     -0.1084  0.2039  -0.5315  0.5951  -0.5081   0.2913      
# StressorW:HabitatTropical:TaxaDeuterostomia       -0.3142  0.1160  -2.7089  0.0068  -0.5415  -0.0869   ** 
# StressorA:HabitatTemperate:TaxaEcdysozoa          -0.1813  0.2350  -0.7718  0.4403  -0.6418   0.2792      
# StressorH:HabitatTemperate:TaxaEcdysozoa          -0.4204  0.1804  -2.3304  0.0198  -0.7739  -0.0668    * 
# StressorW:HabitatTemperate:TaxaEcdysozoa          -0.1963  0.2083  -0.9428  0.3458  -0.6045   0.2118      
# StressorA:HabitatTropical:TaxaEcdysozoa            0.0982  0.0871   1.1272  0.2597  -0.0725   0.2689      
# StressorW:HabitatTropical:TaxaEcdysozoa           -0.3497  0.0826  -4.2333  <.0001  -0.5117  -0.1878  *** 
# StressorW:HabitatFrigid:TaxaLophotrochozoa         0.3123  0.2724   1.1462  0.2517  -0.2217   0.8463      
# StressorA:HabitatTemperate:TaxaLophotrochozoa     -0.0230  0.0711  -0.3236  0.7463  -0.1624   0.1164      
# StressorA-W:HabitatTemperate:TaxaLophotrochozoa   -0.0376  0.2143  -0.1753  0.8609  -0.4577   0.3825      
# StressorH:HabitatTemperate:TaxaLophotrochozoa     -0.4473  0.1542  -2.9000  0.0037  -0.7496  -0.1450   ** 
# StressorW:HabitatTemperate:TaxaLophotrochozoa      0.0897  0.0785   1.1431  0.2530  -0.0641   0.2436      
# StressorA:HabitatTropical:TaxaLophotrochozoa      -0.0973  0.0576  -1.6898  0.0911  -0.2101   0.0156    . 
# StressorA-W:HabitatTropical:TaxaLophotrochozoa    -0.6119  0.3329  -1.8379  0.0661  -1.2644   0.0406    . 
# StressorH:HabitatTropical:TaxaLophotrochozoa       0.2282  0.2944   0.7751  0.4383  -0.3489   0.8053      
# StressorW:HabitatTropical:TaxaLophotrochozoa      -0.2699  0.0993  -2.7167  0.0066  -0.4646  -0.0752   ** 
# StressorA:HabitatTropical:TaxaRadiata             -0.1707  0.1159  -1.4733  0.1407  -0.3978   0.0564      
# StressorW:HabitatTropical:TaxaRadiata             -0.3163  0.0993  -3.1845  0.0015  -0.5110  -0.1216   ** 


############## Stressor:Trophic.level #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                                           estimate      se     zval    pval    ci.lb    ci.ub     
# StressorA:Trophic.levelHerbivore          -0.0288  0.0484  -0.5951  0.5518  -0.1237   0.0661     
# StressorA-W:Trophic.levelHerbivore        -0.2637  0.1919  -1.3743  0.1693  -0.6397   0.1124     
# StressorH:Trophic.levelHerbivore          -0.2611  0.1414  -1.8470  0.0647  -0.5381   0.0160   . 
# StressorW:Trophic.levelHerbivore          -0.1783  0.0552  -3.2321  0.0012  -0.2864  -0.0702  ** 
# StressorA:Trophic.levelPredator           -0.1063  0.0668  -1.5911  0.1116  -0.2372   0.0246     
# StressorA-W:Trophic.levelPredator         -0.2915  0.1532  -1.9026  0.0571  -0.5917   0.0088   . 
# StressorH:Trophic.levelPredator           -0.4281  0.1869  -2.2912  0.0220  -0.7943  -0.0619   * 
# StressorW:Trophic.levelPredator           -0.2090  0.0734  -2.8454  0.0044  -0.3529  -0.0650  ** 
# StressorA:Trophic.levelPrimary producer    0.1271  0.0674   1.8845  0.0595  -0.0051   0.2593   . 
# StressorW:Trophic.levelPrimary producer    0.1287  0.1236   1.0405  0.2981  -0.1137   0.3710     


############## Habitat:Trophic.level #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                                                  estimate      se     zval    pval    ci.lb 
# HabitatFrigid:Trophic.levelHerbivore              0.3143  0.3052   1.0300  0.3030  -0.2838 
# HabitatTemperate:Trophic.levelHerbivore          -0.0880  0.0616  -1.4268  0.1536  -0.2088 
# HabitatTropical:Trophic.levelHerbivore           -0.1254  0.0444  -2.8241  0.0047  -0.2124 
# HabitatTemperate:Trophic.levelPredator           -0.0914  0.0833  -1.0966  0.2728  -0.2546 
# HabitatTropical:Trophic.levelPredator            -0.2299  0.0632  -3.6404  0.0003  -0.3537 
# HabitatFrigid:Trophic.levelPrimary producer       0.4459  0.2292   1.9460  0.0517  -0.0032 
# HabitatTemperate:Trophic.levelPrimary producer    0.1841  0.1196   1.5399  0.1236  -0.0502 
# HabitatTropical:Trophic.levelPrimary producer     0.0721  0.0723   0.9971  0.3187  -0.0696 
#                                                    ci.ub      
# HabitatFrigid:Trophic.levelHerbivore             0.9124      
# HabitatTemperate:Trophic.levelHerbivore          0.0329      
# HabitatTropical:Trophic.levelHerbivore          -0.0384   ** 
#   HabitatTemperate:Trophic.levelPredator           0.0719      
# HabitatTropical:Trophic.levelPredator           -0.1061  *** 
#   HabitatFrigid:Trophic.levelPrimary producer      0.8951    . 
# HabitatTemperate:Trophic.levelPrimary producer   0.4185      
# HabitatTropical:Trophic.levelPrimary producer    0.2137      


############## Trophic.level:Metal #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                                               estimate      se     zval    pval    ci.lb    ci.ub      
# Trophic.levelHerbivore:MetalArsenic          -0.0713  0.2627  -0.2713  0.7862  -0.5862   0.4437      
# Trophic.levelHerbivore:MetalCadmium          -0.1498  0.0640  -2.3428  0.0191  -0.2752  -0.0245    * 
# Trophic.levelPredator:MetalCadmium           -0.0445  0.1256  -0.3542  0.7232  -0.2906   0.2017      
# Trophic.levelPrimary producer:MetalCadmium    0.3788  0.1428   2.6524  0.0080   0.0989   0.6587   ** 
# Trophic.levelHerbivore:MetalCobalt           -0.0409  0.1982  -0.2064  0.8365  -0.4295   0.3476      
# Trophic.levelPredator:MetalCobalt            -0.0902  0.1613  -0.5590  0.5762  -0.4064   0.2260      
# Trophic.levelPrimary producer:MetalCobalt     0.1706  0.3155   0.5407  0.5887  -0.4477   0.7889      
# Trophic.levelHerbivore:MetalCopper           -0.1716  0.0596  -2.8769  0.0040  -0.2884  -0.0547   ** 
# Trophic.levelPredator:MetalCopper            -0.2147  0.0844  -2.5432  0.0110  -0.3802  -0.0492    * 
# Trophic.levelPrimary producer:MetalCopper    -0.0224  0.0892  -0.2506  0.8021  -0.1971   0.1524      
# Trophic.levelPrimary producer:MetalIron       0.2015  0.1132   1.7801  0.0751  -0.0204   0.4233    . 
# Trophic.levelHerbivore:MetalLead             -0.1656  0.1171  -1.4140  0.1574  -0.3951   0.0639      
# Trophic.levelPredator:MetalLead              -1.1290  0.2462  -4.5852  <.0001  -1.6116  -0.6464  *** 
# Trophic.levelHerbivore:MetalManganese        -0.0095  0.2916  -0.0325  0.9741  -0.5811   0.5621      
# Trophic.levelPredator:MetalManganese         -0.1454  0.1642  -0.8856  0.3758  -0.4672   0.1764      
# Trophic.levelHerbivore:MetalMercury           0.0441  0.0902   0.4887  0.6251  -0.1326   0.2208      
# Trophic.levelPredator:MetalMercury           -0.1579  0.1137  -1.3887  0.1649  -0.3808   0.0650      
# Trophic.levelHerbivore:MetalNickel           -0.0125  0.1728  -0.0726  0.9422  -0.3513   0.3262      
# Trophic.levelHerbivore:MetalSilver           -0.1579  0.1460  -1.0821  0.2792  -0.4440   0.1281      
# Trophic.levelPredator:MetalSilver            -0.1267  0.2053  -0.6172  0.5371  -0.5290   0.2756      
# Trophic.levelHerbivore:MetalZinc              0.0757  0.1457   0.5196  0.6033  -0.2099   0.3614      
# Trophic.levelPredator:MetalZinc              -0.0923  0.1405  -0.6573  0.5110  -0.3677   0.1830      
# Trophic.levelPrimary producer:MetalZinc       0.1882  0.1748   1.0763  0.2818  -0.1545   0.5309      




############## Stressor:Metal ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                            estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:MetalArsenic     -0.1258  0.2120  -0.5933  0.5530  -0.5414   0.2898      
# StressorA:MetalCadmium      0.0223  0.0668   0.3344  0.7381  -0.1086   0.1533      
# StressorA-W:MetalCadmium   -0.2255  0.2583  -0.8732  0.3825  -0.7317   0.2807      
# StressorH:MetalCadmium     -0.3521  0.1324  -2.6584  0.0079  -0.6116  -0.0925   ** 
# StressorW:MetalCadmium     -0.0307  0.1198  -0.2565  0.7975  -0.2654   0.2040      
# StressorA:MetalCobalt      -0.0744  0.1407  -0.5283  0.5973  -0.3502   0.2015      
# StressorW:MetalCobalt       0.2683  0.2002   1.3397  0.1803  -0.1242   0.6607      
# StressorA:MetalCopper      -0.0985  0.0575  -1.7127  0.0868  -0.2112   0.0142    . 
# StressorA-W:MetalCopper    -0.2858  0.2862  -0.9988  0.3179  -0.8467   0.2751      
# StressorW:MetalCopper      -0.2067  0.0692  -2.9887  0.0028  -0.3423  -0.0712   ** 
# StressorA:MetalIron         0.1375  0.1430   0.9612  0.3365  -0.1428   0.4178      
# StressorW:MetalIron         0.3350  0.1977   1.6946  0.0902  -0.0525   0.7226    . 
# StressorA:MetalLead        -0.1289  0.1805  -0.7141  0.4752  -0.4826   0.2248      
# StressorW:MetalLead        -0.4596  0.1323  -3.4739  0.0005  -0.7189  -0.2003  *** 
# StressorA:MetalManganese    0.1400  0.1848   0.7574  0.4488  -0.2222   0.5021      
# StressorW:MetalManganese    0.0851  0.2338   0.3640  0.7159  -0.3731   0.5433      
# StressorA:MetalMercury      0.1597  0.1004   1.5909  0.1116  -0.0370   0.3565      
# StressorA-W:MetalMercury   -0.1977  0.1829  -1.0804  0.2799  -0.5562   0.1609      
# StressorW:MetalMercury     -0.1993  0.0992  -2.0083  0.0446  -0.3938  -0.0048    * 
# StressorA:MetalNickel       0.2016  0.2487   0.8106  0.4176  -0.2859   0.6891      
# StressorW:MetalNickel      -0.1981  0.2278  -0.8696  0.3845  -0.6447   0.2484      
# StressorA:MetalSilver      -0.0947  0.1345  -0.7043  0.4813  -0.3583   0.1689      
# StressorW:MetalSilver      -0.1035  0.2580  -0.4013  0.6882  -0.6091   0.4021      
# StressorA:MetalZinc         0.0245  0.0972   0.2519  0.8011  -0.1661   0.2150      
# StressorW:MetalZinc         0.0097  0.2603   0.0373  0.9703  -0.5004   0.5198    




############## Stressor:Trait ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                                  estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:TraitAbundance         0.2163  0.1722   1.2561  0.2091  -0.1212   0.5537      
# StressorW:TraitAbundance         0.0065  0.1540   0.0419  0.9665  -0.2954   0.3084      
# StressorA:TraitAccumulation     -0.0328  0.0499  -0.6570  0.5112  -0.1305   0.0650      
# StressorA-W:TraitAccumulation   -0.2142  0.1700  -1.2595  0.2079  -0.5474   0.1191      
# StressorH:TraitAccumulation     -0.2402  0.1402  -1.7133  0.0867  -0.5150   0.0346    . 
# StressorW:TraitAccumulation     -0.2572  0.0766  -3.3590  0.0008  -0.4073  -0.1071  *** 
# StressorA:TraitBehaviour         0.0078  0.2041   0.0380  0.9697  -0.3922   0.4077      
# StressorW:TraitBehaviour        -0.1660  0.2059  -0.8063  0.4201  -0.5695   0.2375      
# StressorA:TraitGrowth            0.0175  0.0673   0.2597  0.7951  -0.1144   0.1494      
# StressorW:TraitGrowth           -0.1403  0.1022  -1.3732  0.1697  -0.3405   0.0599      
# StressorA:TraitMetabolism       -0.0409  0.0571  -0.7170  0.4734  -0.1529   0.0710      
# StressorA-W:TraitMetabolism     -0.5579  0.2281  -2.4464  0.0144  -1.0049  -0.1109    * 
# StressorH:TraitMetabolism       -0.4164  0.2134  -1.9507  0.0511  -0.8347   0.0020    . 
# StressorW:TraitMetabolism       -0.0139  0.0677  -0.2060  0.8368  -0.1465   0.1187      
# StressorA:TraitReproduction      0.1608  0.1395   1.1534  0.2488  -0.1125   0.4342      
# StressorW:TraitReproduction     -0.4882  0.1483  -3.2918  0.0010  -0.7788  -0.1975  *** 
# StressorA:TraitSurvival          0.0619  0.1469   0.4214  0.6735  -0.2260   0.3498      
# StressorW:TraitSurvival         -0.2654  0.1115  -2.3789  0.0174  -0.4840  -0.0467    * 


############## Stressor:Trait:Habitat  #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                                                   estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:TraitMetabolism:HabitatFrigid           0.6633  0.2373   2.7954  0.0052   0.1982   1.1283   ** 
# StressorA:TraitAbundance:HabitatTemperate         0.3551  0.2062   1.7222  0.0850  -0.0490   0.7591    . 
# StressorA:TraitAccumulation:HabitatTemperate     -0.0566  0.0757  -0.7485  0.4542  -0.2049   0.0917      
# StressorA-W:TraitAccumulation:HabitatTemperate   -0.1564  0.2734  -0.5721  0.5672  -0.6923   0.3794      
# StressorH:TraitAccumulation:HabitatTemperate     -0.3725  0.1394  -2.6719  0.0075  -0.6458  -0.0993   ** 
# StressorW:TraitAccumulation:HabitatTemperate      0.1775  0.1188   1.4940  0.1352  -0.0554   0.4104      
# StressorA:TraitBehaviour:HabitatTemperate        -0.1154  0.3166  -0.3647  0.7154  -0.7359   0.5050      
# StressorA:TraitGrowth:HabitatTemperate           -0.0317  0.1406  -0.2258  0.8214  -0.3073   0.2438      
# StressorW:TraitGrowth:HabitatTemperate           -0.3005  0.1959  -1.5339  0.1251  -0.6844   0.0835      
# StressorA:TraitMetabolism:HabitatTemperate        0.0155  0.1226   0.1263  0.8995  -0.2248   0.2558      
# StressorH:TraitMetabolism:HabitatTemperate       -0.4023  0.2493  -1.6139  0.1065  -0.8908   0.0862      
# StressorW:TraitMetabolism:HabitatTemperate        0.2219  0.1127   1.9693  0.0489   0.0010   0.4428    * 
# StressorW:TraitSurvival:HabitatTemperate         -0.2092  0.2142  -0.9770  0.3286  -0.6290   0.2105      
# StressorA:TraitAbundance:HabitatTropical         -0.1041  0.2868  -0.3631  0.7165  -0.6663   0.4580      
# StressorW:TraitAbundance:HabitatTropical         -0.1298  0.1563  -0.8301  0.4065  -0.4362   0.1767      
# StressorA:TraitAccumulation:HabitatTropical       0.0040  0.0616   0.0650  0.9481  -0.1166   0.1247      
# StressorA-W:TraitAccumulation:HabitatTropical    -0.2069  0.2043  -1.0125  0.3113  -0.6074   0.1936      
# StressorW:TraitAccumulation:HabitatTropical      -0.5518  0.0969  -5.6921  <.0001  -0.7418  -0.3618  *** 
# StressorA:TraitBehaviour:HabitatTropical          0.0695  0.2528   0.2749  0.7834  -0.4259   0.5649      
# StressorW:TraitBehaviour:HabitatTropical         -0.2089  0.2129  -0.9811  0.3265  -0.6261   0.2084      
# StressorA:TraitGrowth:HabitatTropical             0.0374  0.0737   0.5077  0.6116  -0.1071   0.1820      
# StressorW:TraitGrowth:HabitatTropical            -0.1056  0.1130  -0.9338  0.3504  -0.3271   0.1160      
# StressorA:TraitMetabolism:HabitatTropical        -0.0567  0.0619  -0.9158  0.3598  -0.1780   0.0646      
# StressorA-W:TraitMetabolism:HabitatTropical      -0.5931  0.2202  -2.6933  0.0071  -1.0248  -0.1615   ** 
# StressorW:TraitMetabolism:HabitatTropical        -0.2227  0.0843  -2.6409  0.0083  -0.3880  -0.0574   ** 
# StressorA:TraitReproduction:HabitatTropical       0.2201  0.1401   1.5712  0.1161  -0.0545   0.4947      
# StressorW:TraitReproduction:HabitatTropical      -0.5193  0.1428  -3.6357  0.0003  -0.7993  -0.2394  *** 
# StressorA:TraitSurvival:HabitatTropical           0.1148  0.1431   0.8025  0.4223  -0.1656   0.3952      
# StressorW:TraitSurvival:HabitatTropical          -0.3076  0.1231  -2.4982  0.0125  -0.5490  -0.0663    * 


############## Trophic.level-Accumulation #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                                 estimate      se     zval    pval    ci.lb    ci.ub    
# Trophic.levelHerbivore          -0.1021  0.0624  -1.6351  0.1020  -0.2245   0.0203    
# Trophic.levelPredator           -0.2097  0.0838  -2.5029  0.0123  -0.3739  -0.0455  * 
# Trophic.levelPrimary producer    0.4958  0.2016   2.4594  0.0139   0.1007   0.8909  * 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
# Estimate Std. Error z value Pr(>|z|)   
# Predator - Herbivore == 0          -0.1076     0.1045  -1.030  0.30302   
# Primary producer - Herbivore == 0   0.5979     0.2110   2.833  0.00461 **
# Primary producer - Predator == 0    0.7055     0.2183   3.232  0.00123 **


############## Trophic.level-Abundance #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data <- subset(data, Trait == "Abundance")
nrow(data)
data<-subset (data,data$Paper_id!=134)
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
#                                     Estimate Std. Error z value Pr(>|z|)
# Predator - Herbivore == 0          0.09472    0.55085   0.172    0.863
# Primary producer - Herbivore == 0  0.43659    0.50479   0.865    0.387
# Primary producer - Predator == 0   0.34187    0.46175   0.740    0.459











### B. Removing one trait of experiment ###

############## overall #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
# estimate      se     zval    pval    ci.lb    ci.ub     
# -0.0852  0.0284  -3.0051  0.0027  -0.1408  -0.0296  ** 

rm(list = ls())



############## Stressor #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
nrow(data)
data$Stressor<-factor(data$Stressor)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
#               estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA     -0.0114  0.0346  -0.3306  0.7410  -0.0792   0.0563      
# StressorA-W   -0.2554  0.1119  -2.2822  0.0225  -0.4748  -0.0361    * 
# StressorH     -0.3243  0.1107  -2.9300  0.0034  -0.5412  -0.1074   ** 
# StressorW     -0.1527  0.0435  -3.5107  0.0004  -0.2380  -0.0674  ***  

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("A"=1,"A-W"=1,"H"=1,"W"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#               Estimate Std. Error z value Pr(>|z|)   
# A-W - A == 0 -0.24401    0.11296  -2.160  0.03077 * 
# H - A == 0   -0.31287    0.11596  -2.698  0.00697 **
# W - A == 0   -0.14127    0.05351  -2.640  0.00829 **
# H - A-W == 0 -0.06886    0.15741  -0.437  0.66178   
# W - A-W == 0  0.10274    0.11525   0.891  0.37269   
# W - H == 0    0.17160    0.11892   1.443  0.14903   






############## Taxa ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
#                      estimate      se     zval    pval    ci.lb    ci.ub      
# TaxaAlgae+Protist     0.1587  0.0589   2.6947  0.0070   0.0433   0.2742   ** 
# TaxaDeuterostomia    -0.2117  0.0807  -2.6225  0.0087  -0.3700  -0.0535   ** 
# TaxaEcdysozoa        -0.2020  0.0614  -3.2910  0.0010  -0.3223  -0.0817  *** 
# TaxaLophotrochozoa   -0.0645  0.0387  -1.6680  0.0953  -0.1404   0.0113    . 
# TaxaRadiata          -0.2696  0.0832  -3.2414  0.0012  -0.4326  -0.1066   ** 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Algae+Protist"=1,"Deuterostomia"=1,"Ecdysozoa"=1,"Lophotrochozoa"=1,"Radiata"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                      Estimate Std. Error z value Pr(>|z|)    
# Deuterostomia - Algae+Protist == 0  -0.370437   0.099926  -3.707  0.00021 ***
# Ecdysozoa - Algae+Protist == 0      -0.360702   0.085015  -4.243 2.21e-05 ***
# Lophotrochozoa - Algae+Protist == 0 -0.223253   0.068680  -3.251  0.00115 ** 
# Radiata - Algae+Protist == 0        -0.428309   0.101916  -4.203 2.64e-05 ***
# Ecdysozoa - Deuterostomia == 0       0.009735   0.101139   0.096  0.92332    
# Lophotrochozoa - Deuterostomia == 0  0.147184   0.089407   1.646  0.09972 .  
# Radiata - Deuterostomia == 0        -0.057873   0.115913  -0.499  0.61759    
# Lophotrochozoa - Ecdysozoa == 0      0.137448   0.071855   1.913  0.05577 .  
# Radiata - Ecdysozoa == 0            -0.067608   0.103367  -0.654  0.51308    
# Radiata - Lophotrochozoa == 0       -0.205056   0.091733  -2.235  0.02539 *  

rm(list = ls())


############## Habitat ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
#                     estimate      se     zval    pval    ci.lb    ci.ub     
# HabitatFrigid       0.3240  0.1940   1.6698  0.0950  -0.0563   0.7043   . 
# HabitatTemperate   -0.0688  0.0489  -1.4076  0.1592  -0.1647   0.0270     
# HabitatTropical    -0.1062  0.0347  -3.0606  0.0022  -0.1742  -0.0382  ** 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Frigid"=1,"Temperate"=1,"Tropical"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                           Estimate Std. Error z value Pr(>|z|)  
# Temperate - Frigid == 0   -0.39283    0.20011  -1.963   0.0496 *
# Tropical - Frigid == 0    -0.43019    0.19712  -2.182   0.0291 *
# Tropical - Temperate == 0 -0.03736    0.05991  -0.624   0.5329  

rm(list = ls())


############## Trait ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
# estimate      se     zval    pval    ci.lb    ci.ub    
# TraitAccumulation   -0.1035  0.0408  -2.5340  0.0113  -0.1835  -0.0234  * 
# TraitBehaviour      -0.1344  0.1275  -1.0540  0.2919  -0.3844   0.1155    
# TraitGrowth         -0.0522  0.0529  -0.9873  0.3235  -0.1558   0.0514    
# TraitMetabolism     -0.0647  0.0415  -1.5602  0.1187  -0.1460   0.0166    
# TraitReproduction   -0.1526  0.0947  -1.6115  0.1071  -0.3382   0.0330    
# TraitSurvival       -0.1182  0.0838  -1.4112  0.1582  -0.2824   0.0460    

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Accumulation"=1,"Behaviour"=1,
                                     "Growth"=1,"Metabolism"=1,"Reproduction"=1,"Survival"=1), 
                                   type="Tukey"))), test=adjusted("none"))
# Estimate Std. Error z value Pr(>|z|)
# Behaviour - Accumulation == 0    -0.03096    0.13124  -0.236    0.813
# Growth - Accumulation == 0        0.05125    0.06169   0.831    0.406
# Metabolism - Accumulation == 0    0.03875    0.05235   0.740    0.459
# Reproduction - Accumulation == 0 -0.04913    0.09937  -0.494    0.621
# Survival - Accumulation == 0     -0.01477    0.08950  -0.165    0.869
# Growth - Behaviour == 0           0.08221    0.13358   0.615    0.538
# Metabolism - Behaviour == 0       0.06971    0.13243   0.526    0.599
# Reproduction - Behaviour == 0    -0.01817    0.15604  -0.116    0.907
# Survival - Behaviour == 0         0.01619    0.14789   0.109    0.913
# Metabolism - Growth == 0         -0.01250    0.06116  -0.204    0.838
# Reproduction - Growth == 0       -0.10038    0.09956  -1.008    0.313
# Survival - Growth == 0           -0.06601    0.09400  -0.702    0.483
# Reproduction - Metabolism == 0   -0.08788    0.10055  -0.874    0.382
# Survival - Metabolism == 0       -0.05352    0.09085  -0.589    0.556
# Survival - Reproduction == 0      0.03437    0.11878   0.289    0.772



############## Metal####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
# MetalArsenic      0.0601  0.1576   0.3815  0.7028  -0.2488   0.3691     
# MetalCadmium     -0.0488  0.0548  -0.8901  0.3734  -0.1563   0.0587     
# MetalCobalt      -0.0620  0.1175  -0.5275  0.5978  -0.2922   0.1683     
# MetalCopper      -0.1462  0.0458  -3.1919  0.0014  -0.2360  -0.0564  ** 
# MetalIron         0.1008  0.1273   0.7919  0.4284  -0.1486   0.3502     
# MetalLead        -0.3754  0.1141  -3.2888  0.0010  -0.5991  -0.1517  ** 
# MetalManganese   -0.1149  0.1345  -0.8543  0.3930  -0.3784   0.1487     
# MetalMercury     -0.0457  0.0761  -0.5999  0.5486  -0.1948   0.1035     
# MetalNickel       0.0000  0.1806   0.0001  0.9999  -0.3539   0.3540     
# MetalSilver      -0.1310  0.1128  -1.1617  0.2454  -0.3521   0.0900     
# MetalZinc         0.0256  0.0851   0.3005  0.7638  -0.1413   0.1924      

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Arsenic "=1,"Cadmium"=1,"Cobalt"=1,
                                     "Copper"=1,"Iron"=1,"Lead "=1,"Manganese"=1,
                                     "Mercury"=1,"Nickel"=1,"Silver"=1,"Zinc"=1), 
                                   type="Tukey"))), test=adjusted("none"))
# Estimate Std. Error z value Pr(>|z|)   
# Cadmium - Arsenic  == 0   -0.108947   0.164858  -0.661  0.50871   
# Cobalt - Arsenic  == 0    -0.122106   0.196376  -0.622  0.53408   
# Copper - Arsenic  == 0    -0.206376   0.163782  -1.260  0.20765   
# Iron - Arsenic  == 0       0.040644   0.202579   0.201  0.84099   
# Lead  - Arsenic  == 0     -0.435516   0.193344  -2.253  0.02429 * 
#   Manganese - Arsenic  == 0 -0.175008   0.206953  -0.846  0.39775   
# Mercury - Arsenic  == 0   -0.105791   0.174460  -0.606  0.54425   
# Nickel - Arsenic  == 0    -0.060121   0.239681  -0.251  0.80194   
# Silver - Arsenic  == 0    -0.191169   0.193197  -0.990  0.32242   
# Zinc - Arsenic  == 0      -0.034551   0.178329  -0.194  0.84637   
# Cobalt - Cadmium == 0     -0.013159   0.127699  -0.103  0.91793   
# Copper - Cadmium == 0     -0.097429   0.069387  -1.404  0.16028   
# Iron - Cadmium == 0        0.149591   0.138548   1.080  0.28028   
# Lead  - Cadmium == 0      -0.326569   0.125380  -2.605  0.00920 **
#   Manganese - Cadmium == 0  -0.066061   0.143279  -0.461  0.64475   
# Mercury - Cadmium == 0     0.003156   0.093314   0.034  0.97302   
# Nickel - Cadmium == 0      0.048826   0.188693   0.259  0.79582   
# Silver - Cadmium == 0     -0.082222   0.117577  -0.699  0.48436   
# Zinc - Cadmium == 0        0.074396   0.099525   0.748  0.45476   
# Copper - Cobalt == 0      -0.084270   0.125665  -0.671  0.50248   
# Iron - Cobalt == 0         0.162749   0.172646   0.943  0.34585   
# Lead  - Cobalt == 0       -0.313411   0.163592  -1.916  0.05539 . 
# Manganese - Cobalt == 0   -0.052902   0.140822  -0.376  0.70716   
# Mercury - Cobalt == 0      0.016314   0.138939   0.117  0.90653   
# Nickel - Cobalt == 0       0.061984   0.215309   0.288  0.77343   
# Silver - Cobalt == 0      -0.069064   0.157314  -0.439  0.66065   
# Zinc - Cobalt == 0         0.087555   0.140478   0.623  0.53311   
# Iron - Copper == 0         0.247019   0.135230   1.827  0.06775 . 
# Lead  - Copper == 0       -0.229140   0.121095  -1.892  0.05846 . 
# Manganese - Copper == 0    0.031368   0.141588   0.222  0.82467   
# Mercury - Copper == 0      0.100584   0.088613   1.135  0.25634   
# Nickel - Copper == 0       0.146254   0.186240   0.785  0.43228   
# Silver - Copper == 0       0.015207   0.120056   0.127  0.89921   
# Zinc - Copper == 0         0.171825   0.093389   1.840  0.06578 . 
# Lead  - Iron == 0         -0.476160   0.170934  -2.786  0.00534 **
#   Manganese - Iron == 0     -0.215651   0.184890  -1.166  0.24346   
# Mercury - Iron == 0       -0.146435   0.148262  -0.988  0.32331   
# Nickel - Iron == 0        -0.100765   0.220909  -0.456  0.64829   
# Silver - Iron == 0        -0.231813   0.169983  -1.364  0.17265   
# Zinc - Iron == 0          -0.075194   0.152643  -0.493  0.62228   
# Manganese - Lead  == 0     0.260509   0.176151   1.479  0.13917   
# Mercury - Lead  == 0       0.329725   0.136734   2.411  0.01589 * 
#   Nickel - Lead  == 0        0.375395   0.213597   1.757  0.07883 . 
# Silver - Lead  == 0        0.244347   0.159968   1.527  0.12664   
# Zinc - Lead  == 0          0.400965   0.140982   2.844  0.00445 **
#   Mercury - Manganese == 0   0.069216   0.153086   0.452  0.65117   
# Nickel - Manganese == 0    0.114886   0.224999   0.511  0.60962   
# Silver - Manganese == 0   -0.016161   0.169279  -0.095  0.92394   
# Zinc - Manganese == 0      0.140457   0.153733   0.914  0.36090   
# Nickel - Mercury == 0      0.045670   0.195940   0.233  0.81570   
# Silver - Mercury == 0     -0.085378   0.134736  -0.634  0.52630   
# Zinc - Mercury == 0        0.071240   0.113038   0.630  0.52854   
# Silver - Nickel == 0      -0.131048   0.212798  -0.616  0.53800   
# Zinc - Nickel == 0         0.025571   0.198149   0.129  0.89732   
# Zinc - Silver == 0         0.156618   0.136772   1.145  0.25216      





############## Trophic.level #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# Trophic.levelHerbivore          -0.0940  0.0359  -2.6168  0.0089  -0.1643  -0.0236   ** 
# Trophic.levelPredator           -0.1957  0.0518  -3.7749  0.0002  -0.2972  -0.0941  *** 
# Trophic.levelPrimary producer    0.0870  0.0573   1.5181  0.1290  -0.0253   0.1994      
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                   Estimate Std. Error z value Pr(>|z|)    
# Predator - Herbivore == 0         -0.10170    0.06235  -1.631 0.102891    
# Primary producer - Herbivore == 0  0.18099    0.06609   2.739 0.006169 ** 
# Primary producer - Predator == 0   0.28269    0.07563   3.738 0.000186 ***





############## Stressor:Taxa #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:TaxaAlgae+Protist       0.1513  0.0601   2.5159  0.0119   0.0334   0.2691    * 
# StressorA-W:TaxaAlgae+Protist    -0.1399  0.2929  -0.4778  0.6328  -0.7140   0.4342      
# StressorW:TaxaAlgae+Protist       0.2540  0.1409   1.8027  0.0714  -0.0222   0.5302    . 
# StressorA:TaxaDeuterostomia      -0.1211  0.0986  -1.2286  0.2192  -0.3143   0.0721      
# StressorA-W:TaxaDeuterostomia    -0.2207  0.1923  -1.1475  0.2512  -0.5975   0.1562      
# StressorW:TaxaDeuterostomia      -0.3445  0.1139  -3.0247  0.0025  -0.5678  -0.1213   ** 
# StressorA:TaxaEcdysozoa           0.0494  0.0871   0.5671  0.5706  -0.1214   0.2202      
# StressorH:TaxaEcdysozoa          -0.4562  0.1797  -2.5386  0.0111  -0.8084  -0.1040    * 
# StressorW:TaxaEcdysozoa          -0.3462  0.0842  -4.1097  <.0001  -0.5113  -0.1811  *** 
# StressorA:TaxaLophotrochozoa     -0.0545  0.0460  -1.1839  0.2365  -0.1446   0.0357      
# StressorA-W:TaxaLophotrochozoa   -0.2135  0.1670  -1.2785  0.2011  -0.5408   0.1138      
# StressorH:TaxaLophotrochozoa     -0.2878  0.1356  -2.1221  0.0338  -0.5536  -0.0220    * 
# StressorW:TaxaLophotrochozoa     -0.0335  0.0600  -0.5576  0.5771  -0.1511   0.0842      
# StressorA:TaxaRadiata            -0.1929  0.1220  -1.5811  0.1139  -0.4320   0.0462      
# StressorW:TaxaRadiata            -0.3282  0.1048  -3.1319  0.0017  -0.5335  -0.1228   ** 



############## Stressor:Habitat #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
# StressorW:HabitatFrigid         0.4314  0.2051   2.1031  0.0355   0.0294   0.8334    * 
# StressorA:HabitatTemperate     -0.0627  0.0577  -1.0866  0.2772  -0.1758   0.0504      
# StressorA-W:HabitatTemperate   -0.0600  0.1976  -0.3039  0.7612  -0.4473   0.3272      
# StressorH:HabitatTemperate     -0.4057  0.1110  -3.6551  0.0003  -0.6232  -0.1882  *** 
# StressorW:HabitatTemperate      0.0951  0.0718   1.3252  0.1851  -0.0456   0.2359      
# StressorA:HabitatTropical       0.0192  0.0380   0.5040  0.6142  -0.0554   0.0938      
# StressorA-W:HabitatTropical    -0.3432  0.1312  -2.6156  0.0089  -0.6003  -0.0860   ** 
# StressorH:HabitatTropical       0.2318  0.2869   0.8079  0.4191  -0.3305   0.7940      
# StressorW:HabitatTropical      -0.3153  0.0502  -6.2844  <.0001  -0.4136  -0.2170  *** 


############## Stressor:Habitat:Taxa #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
#                                                   estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:HabitatFrigid:TaxaAlgae+Protist          0.5860  0.2891   2.0270  0.0427   0.0194   1.1526    * 
# StressorA:HabitatTemperate:TaxaAlgae+Protist       0.1103  0.1291   0.8547  0.3927  -0.1427   0.3633      
# StressorW:HabitatTemperate:TaxaAlgae+Protist       0.2936  0.2036   1.4422  0.1493  -0.1054   0.6927      
# StressorA:HabitatTropical:TaxaAlgae+Protist        0.1631  0.0670   2.4339  0.0149   0.0318   0.2944    * 
# StressorA-W:HabitatTropical:TaxaAlgae+Protist     -0.2126  0.2933  -0.7248  0.4686  -0.7874   0.3623      
# StressorW:HabitatTropical:TaxaAlgae+Protist       -0.0395  0.2413  -0.1636  0.8700  -0.5124   0.4335      
# StressorA:HabitatTemperate:TaxaDeuterostomia      -0.3190  0.1364  -2.3384  0.0194  -0.5863  -0.0516    * 
# StressorA:HabitatTropical:TaxaDeuterostomia        0.0696  0.1329   0.5237  0.6005  -0.1909   0.3301      
# StressorA-W:HabitatTropical:TaxaDeuterostomia     -0.1408  0.1926  -0.7313  0.4646  -0.5183   0.2366      
# StressorW:HabitatTropical:TaxaDeuterostomia       -0.3157  0.1121  -2.8168  0.0048  -0.5353  -0.0960   ** 
# StressorA:HabitatTemperate:TaxaEcdysozoa          -0.1837  0.2273  -0.8081  0.4191  -0.6292   0.2618      
# StressorH:HabitatTemperate:TaxaEcdysozoa          -0.4488  0.1751  -2.5630  0.0104  -0.7920  -0.1056    * 
# StressorW:HabitatTemperate:TaxaEcdysozoa          -0.1727  0.2041  -0.8460  0.3976  -0.5728   0.2274      
# StressorA:HabitatTropical:TaxaEcdysozoa            0.0900  0.0901   0.9988  0.3179  -0.0866   0.2666      
# StressorW:HabitatTropical:TaxaEcdysozoa           -0.3716  0.0885  -4.1979  <.0001  -0.5451  -0.1981  *** 
# StressorW:HabitatFrigid:TaxaLophotrochozoa         0.2925  0.2721   1.0749  0.2824  -0.2408   0.8259      
# StressorA:HabitatTemperate:TaxaLophotrochozoa     -0.0314  0.0723  -0.4352  0.6634  -0.1731   0.1102      
# StressorA-W:HabitatTemperate:TaxaLophotrochozoa   -0.0467  0.1967  -0.2372  0.8125  -0.4321   0.3388      
# StressorH:HabitatTemperate:TaxaLophotrochozoa     -0.4401  0.1501  -2.9321  0.0034  -0.7343  -0.1459   ** 
# StressorW:HabitatTemperate:TaxaLophotrochozoa      0.0934  0.0765   1.2214  0.2219  -0.0565   0.2433      
# StressorA:HabitatTropical:TaxaLophotrochozoa      -0.0606  0.0564  -1.0756  0.2821  -0.1711   0.0498      
# StressorA-W:HabitatTropical:TaxaLophotrochozoa    -0.5812  0.3133  -1.8551  0.0636  -1.1952   0.0329    . 
# StressorH:HabitatTropical:TaxaLophotrochozoa       0.2329  0.2808   0.8294  0.4069  -0.3175   0.7833      
# StressorW:HabitatTropical:TaxaLophotrochozoa      -0.2687  0.0978  -2.7486  0.0060  -0.4603  -0.0771   ** 
# StressorA:HabitatTropical:TaxaRadiata             -0.1916  0.1190  -1.6101  0.1074  -0.4249   0.0416      
# StressorW:HabitatTropical:TaxaRadiata             -0.3305  0.1024  -3.2275  0.0012  -0.5312  -0.1298   ** 






############## Stressor:Trophic.level #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
#                                           estimate      se     zval    pval    ci.lb    ci.ub     
# StressorA:Trophic.levelHerbivore          -0.0150  0.0478  -0.3133  0.7540  -0.1086   0.0787     
# StressorA-W:Trophic.levelHerbivore        -0.2631  0.1864  -1.4114  0.1581  -0.6284   0.1022     
# StressorH:Trophic.levelHerbivore          -0.2485  0.1370  -1.8137  0.0697  -0.5169   0.0200   . 
# StressorW:Trophic.levelHerbivore          -0.1683  0.0554  -3.0390  0.0024  -0.2769  -0.0598  ** 
# StressorA:Trophic.levelPredator           -0.1183  0.0664  -1.7819  0.0748  -0.2485   0.0118   . 
# StressorA-W:Trophic.levelPredator         -0.3176  0.1419  -2.2379  0.0252  -0.5958  -0.0394   * 
# StressorH:Trophic.levelPredator           -0.4586  0.1827  -2.5095  0.0121  -0.8167  -0.1004   * 
# StressorW:Trophic.levelPredator           -0.2285  0.0742  -3.0776  0.0021  -0.3740  -0.0830  ** 
# StressorA:Trophic.levelPrimary producer    0.0924  0.0632   1.4614  0.1439  -0.0315   0.2162     
# StressorW:Trophic.levelPrimary producer    0.1053  0.1285   0.8196  0.4124  -0.1465   0.3570     

############## Habitat:Trophic.level #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
#                                                  estimate      se     zval    pval 
# HabitatFrigid:Trophic.levelHerbivore              0.2950  0.3084   0.9566  0.3388 
# HabitatTemperate:Trophic.levelHerbivore          -0.0903  0.0615  -1.4681  0.1421 
# HabitatTropical:Trophic.levelHerbivore           -0.1034  0.0446  -2.3160  0.0206 
# HabitatTemperate:Trophic.levelPredator           -0.1050  0.0843  -1.2449  0.2132 
# HabitatTropical:Trophic.levelPredator            -0.2521  0.0657  -3.8336  0.0001 
# HabitatFrigid:Trophic.levelPrimary producer       0.3427  0.2385   1.4366  0.1508 
# HabitatTemperate:Trophic.levelPrimary producer    0.1043  0.1275   0.8184  0.4131 
# HabitatTropical:Trophic.levelPrimary producer     0.0584  0.0668   0.8746  0.3818 
#                                                  ci.lb    ci.ub      
# HabitatFrigid:Trophic.levelHerbivore            -0.3094   0.8994      
# HabitatTemperate:Trophic.levelHerbivore         -0.2108   0.0302      
# HabitatTropical:Trophic.levelHerbivore          -0.1909  -0.0159    * 
# HabitatTemperate:Trophic.levelPredator          -0.2703   0.0603      
# HabitatTropical:Trophic.levelPredator           -0.3809  -0.1232  *** 
# HabitatFrigid:Trophic.levelPrimary producer     -0.1249   0.8102      
# HabitatTemperate:Trophic.levelPrimary producer  -0.1455   0.3542      
# HabitatTropical:Trophic.levelPrimary producer   -0.0725   0.1893      



############## Trophic.level:Metal #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
#                                               estimate      se     zval    pval 
# Trophic.levelHerbivore:MetalArsenic           0.1746  0.1782   0.9799  0.3271 
# Trophic.levelPrimary producer:MetalArsenic    0.1738  0.1909   0.9105  0.3626 
# Trophic.levelHerbivore:MetalCadmium          -0.1405  0.0619  -2.2676  0.0234 
# Trophic.levelPredator:MetalCadmium           -0.0408  0.1243  -0.3283  0.7427 
# Trophic.levelPrimary producer:MetalCadmium    0.4028  0.1418   2.8416  0.0045 
# Trophic.levelHerbivore:MetalCobalt           -0.0425  0.1822  -0.2334  0.8155 
# Trophic.levelPredator:MetalCobalt            -0.1379  0.1538  -0.8971  0.3697 
# Trophic.levelPrimary producer:MetalCobalt     0.1238  0.2902   0.4268  0.6696 
# Trophic.levelHerbivore:MetalCopper           -0.1742  0.0579  -3.0097  0.0026 
# Trophic.levelPredator:MetalCopper            -0.2332  0.0824  -2.8316  0.0046 
# Trophic.levelPrimary producer:MetalCopper    -0.0260  0.0864  -0.3011  0.7633 
# Trophic.levelPrimary producer:MetalIron       0.0991  0.1184   0.8374  0.4024 
# Trophic.levelHerbivore:MetalLead             -0.1336  0.1152  -1.1597  0.2462 
# Trophic.levelPredator:MetalLead              -2.8951  0.3893  -7.4370  <.0001 
# Trophic.levelHerbivore:MetalManganese        -0.0003  0.2620  -0.0012  0.9991 
# Trophic.levelPredator:MetalManganese         -0.1946  0.1553  -1.2528  0.2103 
# Trophic.levelHerbivore:MetalMercury           0.0341  0.0887   0.3844  0.7007 
# Trophic.levelPredator:MetalMercury           -0.1701  0.1146  -1.4852  0.1375 
# Trophic.levelHerbivore:MetalNickel            0.0105  0.1782   0.0591  0.9528 
# Trophic.levelHerbivore:MetalSilver           -0.1666  0.1342  -1.2416  0.2144 
# Trophic.levelPredator:MetalSilver            -0.1631  0.1878  -0.8688  0.3850 
# Trophic.levelHerbivore:MetalZinc              0.0584  0.1318   0.4433  0.6575 
# Trophic.levelPredator:MetalZinc              -0.0883  0.1315  -0.6711  0.5022 
# Trophic.levelPrimary producer:MetalZinc       0.1681  0.1663   1.0111  0.3119 
#                                              ci.lb    ci.ub      
# Trophic.levelHerbivore:MetalArsenic         -0.1746   0.5238      
# Trophic.levelPrimary producer:MetalArsenic  -0.2004   0.5480      
# Trophic.levelHerbivore:MetalCadmium         -0.2619  -0.0191    * 
# Trophic.levelPredator:MetalCadmium          -0.2845   0.2028      
# Trophic.levelPrimary producer:MetalCadmium   0.1250   0.6807   ** 
# Trophic.levelHerbivore:MetalCobalt          -0.3997   0.3147      
# Trophic.levelPredator:MetalCobalt           -0.4393   0.1634      
# Trophic.levelPrimary producer:MetalCobalt   -0.4449   0.6926      
# Trophic.levelHerbivore:MetalCopper          -0.2876  -0.0608   ** 
# Trophic.levelPredator:MetalCopper           -0.3947  -0.0718   ** 
# Trophic.levelPrimary producer:MetalCopper   -0.1953   0.1433      
# Trophic.levelPrimary producer:MetalIron     -0.1329   0.3311      
# Trophic.levelHerbivore:MetalLead            -0.3595   0.0922      
# Trophic.levelPredator:MetalLead             -3.6580  -2.1321  *** 
# Trophic.levelHerbivore:MetalManganese       -0.5138   0.5132      
# Trophic.levelPredator:MetalManganese        -0.4989   0.1098      
# Trophic.levelHerbivore:MetalMercury         -0.1398   0.2080      
# Trophic.levelPredator:MetalMercury          -0.3947   0.0544      
# Trophic.levelHerbivore:MetalNickel          -0.3387   0.3597      
# Trophic.levelHerbivore:MetalSilver          -0.4297   0.0964      
# Trophic.levelPredator:MetalSilver           -0.5312   0.2049      
# Trophic.levelHerbivore:MetalZinc            -0.1999   0.3167      
# Trophic.levelPredator:MetalZinc             -0.3461   0.1696      
# Trophic.levelPrimary producer:MetalZinc     -0.1578   0.4940      




############## Stressor:Metal ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
#                            estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:MetalArsenic      0.3543  0.2340   1.5142  0.1300  -0.1043   0.8128      
# StressorW:MetalArsenic     -0.1482  0.2007  -0.7386  0.4602  -0.5417   0.2452      
# StressorA:MetalCadmium      0.0356  0.0663   0.5367  0.5914  -0.0944   0.1655      
# StressorA-W:MetalCadmium   -0.1954  0.2370  -0.8242  0.4098  -0.6600   0.2692      
# StressorH:MetalCadmium     -0.3566  0.1302  -2.7382  0.0062  -0.6118  -0.1013   ** 
# StressorW:MetalCadmium     -0.0294  0.1147  -0.2564  0.7976  -0.2542   0.1953      
# StressorA:MetalCobalt      -0.0841  0.1331  -0.6320  0.5274  -0.3451   0.1768      
# StressorW:MetalCobalt       0.2791  0.1876   1.4876  0.1368  -0.0886   0.6468      
# StressorA:MetalCopper      -0.1051  0.0566  -1.8557  0.0635  -0.2160   0.0059    . 
# StressorA-W:MetalCopper    -0.3251  0.2963  -1.0975  0.2724  -0.9058   0.2555      
# StressorW:MetalCopper      -0.2030  0.0705  -2.8774  0.0040  -0.3412  -0.0647   ** 
# StressorA:MetalIron         0.0129  0.1522   0.0849  0.9324  -0.2853   0.3111      
# StressorW:MetalIron         0.2629  0.2094   1.2559  0.2092  -0.1474   0.6733      
# StressorA:MetalLead        -0.1281  0.1720  -0.7449  0.4563  -0.4653   0.2090      
# StressorW:MetalLead        -0.5410  0.1458  -3.7112  0.0002  -0.8267  -0.2553  *** 
# StressorA:MetalManganese    0.1309  0.1702   0.7692  0.4418  -0.2026   0.4644      
# StressorW:MetalManganese    0.0920  0.2154   0.4273  0.6692  -0.3301   0.5142      
# StressorA:MetalMercury      0.1589  0.1007   1.5782  0.1145  -0.0384   0.3562      
# StressorA-W:MetalMercury   -0.1926  0.1703  -1.1307  0.2582  -0.5264   0.1412      
# StressorW:MetalMercury     -0.2047  0.0961  -2.1312  0.0331  -0.3930  -0.0165    * 
# StressorA:MetalNickel       0.1839  0.2439   0.7541  0.4508  -0.2941   0.6619      
# StressorW:MetalNickel      -0.1905  0.2467  -0.7722  0.4400  -0.6741   0.2930      
# StressorA:MetalSilver      -0.1062  0.1258  -0.8442  0.3985  -0.3528   0.1404      
# StressorW:MetalSilver      -0.0627  0.2383  -0.2633  0.7924  -0.5297   0.4043      
# StressorA:MetalZinc         0.0165  0.0922   0.1786  0.8583  -0.1642   0.1971      
# StressorW:MetalZinc        -0.0063  0.2425  -0.0258  0.9794  -0.4815   0.4690      



############## Stressor:Trait ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
# StressorA:TraitAccumulation     -0.0097  0.0468  -0.2077  0.8355  -0.1014   0.0820      
# StressorA-W:TraitAccumulation   -0.2251  0.1537  -1.4642  0.1431  -0.5264   0.0762      
# StressorH:TraitAccumulation     -0.2325  0.1370  -1.6970  0.0897  -0.5009   0.0360    . 
# StressorW:TraitAccumulation     -0.2612  0.0724  -3.6052  0.0003  -0.4031  -0.1192  *** 
# StressorA:TraitBehaviour         0.0092  0.1871   0.0489  0.9610  -0.3576   0.3759      
# StressorW:TraitBehaviour        -0.1621  0.1899  -0.8536  0.3933  -0.5342   0.2101      
# StressorA:TraitGrowth            0.0030  0.0606   0.0487  0.9611  -0.1158   0.1217      
# StressorW:TraitGrowth           -0.1098  0.0957  -1.1465  0.2516  -0.2974   0.0779      
# StressorA:TraitMetabolism       -0.0489  0.0526  -0.9291  0.3528  -0.1520   0.0542      
# StressorA-W:TraitMetabolism     -0.5463  0.2078  -2.6286  0.0086  -0.9537  -0.1390   ** 
# StressorH:TraitMetabolism       -0.4258  0.2007  -2.1217  0.0339  -0.8191  -0.0325    * 
# StressorW:TraitMetabolism       -0.0153  0.0643  -0.2384  0.8116  -0.1413   0.1107      
# StressorA:TraitReproduction      0.1366  0.1275   1.0707  0.2843  -0.1134   0.3865      
# StressorW:TraitReproduction     -0.4627  0.1360  -3.4019  0.0007  -0.7293  -0.1961  *** 
# StressorA:TraitSurvival          0.0479  0.1344   0.3566  0.7214  -0.2155   0.3114      
# StressorW:TraitSurvival         -0.2564  0.1063  -2.4127  0.0158  -0.4647  -0.0481    * 


############## Stressor:Trait:Habitat #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
#                                                   estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:TraitMetabolism:HabitatFrigid           0.6434  0.2317   2.7772  0.0055   0.1893   1.0975   ** 
# StressorA:TraitAccumulation:HabitatTemperate     -0.0645  0.0773  -0.8344  0.4041  -0.2161   0.0870      
# StressorA-W:TraitAccumulation:HabitatTemperate   -0.1529  0.2440  -0.6265  0.5310  -0.6312   0.3254      
# StressorH:TraitAccumulation:HabitatTemperate     -0.3856  0.1374  -2.8069  0.0050  -0.6549  -0.1164   ** 
# StressorW:TraitAccumulation:HabitatTemperate      0.1911  0.1123   1.7009  0.0890  -0.0291   0.4112    . 
# StressorA:TraitBehaviour:HabitatTemperate        -0.1473  0.2979  -0.4945  0.6210  -0.7312   0.4366      
# StressorA:TraitGrowth:HabitatTemperate           -0.0587  0.1345  -0.4364  0.6625  -0.3224   0.2050      
# StressorW:TraitGrowth:HabitatTemperate           -0.2978  0.1812  -1.6433  0.1003  -0.6529   0.0574      
# StressorA:TraitMetabolism:HabitatTemperate        0.0033  0.1172   0.0285  0.9772  -0.2264   0.2331      
# StressorH:TraitMetabolism:HabitatTemperate       -0.3832  0.2344  -1.6348  0.1021  -0.8426   0.0762      
# StressorW:TraitMetabolism:HabitatTemperate        0.2230  0.1077   2.0710  0.0384   0.0120   0.4340    * 
# StressorW:TraitSurvival:HabitatTemperate         -0.1999  0.2336  -0.8555  0.3923  -0.6577   0.2580      
# StressorA:TraitAccumulation:HabitatTropical       0.0437  0.0556   0.7863  0.4317  -0.0652   0.1526      
# StressorA-W:TraitAccumulation:HabitatTropical    -0.2316  0.1838  -1.2599  0.2077  -0.5918   0.1287      
# StressorW:TraitAccumulation:HabitatTropical      -0.5760  0.0917  -6.2786  <.0001  -0.7559  -0.3962  *** 
# StressorA:TraitBehaviour:HabitatTropical          0.0774  0.2269   0.3411  0.7330  -0.3673   0.5220      
# StressorW:TraitBehaviour:HabitatTropical         -0.2135  0.1952  -1.0937  0.2741  -0.5961   0.1691      
# StressorA:TraitGrowth:HabitatTropical             0.0198  0.0654   0.3025  0.7623  -0.1084   0.1479      
# StressorW:TraitGrowth:HabitatTropical            -0.0672  0.1061  -0.6332  0.5266  -0.2752   0.1408      
# StressorA:TraitMetabolism:HabitatTropical        -0.0567  0.0567  -1.0006  0.3170  -0.1679   0.0544      
# StressorA-W:TraitMetabolism:HabitatTropical      -0.5862  0.1994  -2.9397  0.0033  -0.9770  -0.1954   ** 
# StressorW:TraitMetabolism:HabitatTropical        -0.2161  0.0794  -2.7214  0.0065  -0.3717  -0.0605   ** 
# StressorA:TraitReproduction:HabitatTropical       0.1963  0.1274   1.5405  0.1234  -0.0535   0.4461      
# StressorW:TraitReproduction:HabitatTropical      -0.5035  0.1307  -3.8537  0.0001  -0.7596  -0.2474  *** 
# StressorA:TraitSurvival:HabitatTropical           0.1012  0.1299   0.7786  0.4362  -0.1535   0.3559      
# StressorW:TraitSurvival:HabitatTropical          -0.3068  0.1132  -2.7112  0.0067  -0.5286  -0.0850   ** 




############## Trophic.level-Accumulation #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data <- subset(data, Trait == "Accumulation")
nrow(data)
data<-subset(data,data$Trait!="Abundance")
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
#                                 estimate      se     zval    pval    ci.lb    ci.ub 
# Trophic.levelHerbivore          -0.0711  0.0613  -1.1598  0.2461  -0.1913   0.0491 
# Trophic.levelPredator           -0.2149  0.0854  -2.5175  0.0118  -0.3822  -0.0476 
# Trophic.levelPrimary producer    0.3402  0.1515   2.2458  0.0247   0.0433   0.6372 
# 
# Trophic.levelHerbivore           
# Trophic.levelPredator          * 
#   Trophic.levelPrimary producer  * 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                   Estimate Std. Error z value Pr(>|z|)   
# Predator - Herbivore == 0          -0.1438     0.1051  -1.368  0.17119   
# Primary producer - Herbivore == 0   0.4113     0.1553   2.648  0.00809 **
# Primary producer - Predator == 0    0.5552     0.1739   3.192  0.00141 **








### C. Removing extreme value ###

############## overall #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
# estimate      se     zval    pval    ci.lb    ci.ub     
# -0.0783  0.0276  -2.8358  0.0046  -0.1325  -0.0242  ** 


############## Stressor #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
data$Stressor<-factor(data$Stressor)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
#               estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA     -0.0061  0.0340  -0.1798  0.8573  -0.0728   0.0605      
# StressorA-W   -0.2456  0.1156  -2.1252  0.0336  -0.4722  -0.0191    * 
# StressorH     -0.3229  0.1121  -2.8797  0.0040  -0.5426  -0.1031   ** 
# StressorW     -0.1428  0.0426  -3.3508  0.0008  -0.2263  -0.0593  *** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("A"=1,"A-W"=1,"H"=1,"W"=1), 
                                   type="Tukey"))), test=adjusted("none"))

#             Estimate Std. Error z value Pr(>|z|)   
# A-W - A == 0 -0.23951    0.11704  -2.046  0.04072 * 
# H - A == 0   -0.31677    0.11717  -2.704  0.00686 **
# W - A == 0   -0.13667    0.05294  -2.581  0.00984 **
# H - A-W == 0 -0.07725    0.16103  -0.480  0.63141   
# W - A-W == 0  0.10285    0.11914   0.863  0.38801   
# W - H == 0    0.18010    0.11995   1.501  0.13323  


############## Taxa ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
#                      estimate      se     zval    pval    ci.lb    ci.ub     
# TaxaAlgae+Protist     0.1649  0.0575   2.8707  0.0041   0.0523   0.2775  ** 
# TaxaDeuterostomia    -0.2051  0.0807  -2.5421  0.0110  -0.3632  -0.0470   * 
# TaxaEcdysozoa        -0.1938  0.0597  -3.2441  0.0012  -0.3108  -0.0767  ** 
# TaxaLophotrochozoa   -0.0657  0.0385  -1.7048  0.0882  -0.1411   0.0098   . 
# TaxaRadiata          -0.2307  0.0802  -2.8770  0.0040  -0.3879  -0.0735  ** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Algae+Protist"=1,"Deuterostomia"=1,"Ecdysozoa"=1,"Lophotrochozoa"=1,"Radiata"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                     Estimate Std. Error z value Pr(>|z|)    
# Deuterostomia - Algae+Protist == 0  -0.37003    0.09904  -3.736 0.000187 ***
# Ecdysozoa - Algae+Protist == 0      -0.35869    0.08285  -4.330 1.49e-05 ***
# Lophotrochozoa - Algae+Protist == 0 -0.23059    0.06775  -3.404 0.000665 ***
# Radiata - Algae+Protist == 0        -0.39567    0.09866  -4.011 6.06e-05 ***
# Ecdysozoa - Deuterostomia == 0       0.01134    0.10020   0.113 0.909899    
# Lophotrochozoa - Deuterostomia == 0  0.13944    0.08932   1.561 0.118498    
# Radiata - Deuterostomia == 0        -0.02564    0.11376  -0.225 0.821670    
# Lophotrochozoa - Ecdysozoa == 0      0.12810    0.07054   1.816 0.069375 .  
# Radiata - Ecdysozoa == 0            -0.03698    0.10000  -0.370 0.711519    
# Radiata - Lophotrochozoa == 0       -0.16508    0.08897  -1.855 0.063530 . 


############## Habitat ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
# HabitatPolar        0.3917  0.1861   2.1052  0.0353   0.0270   0.7563   * 
# HabitatTemperate   -0.0496  0.0476  -1.0415  0.2976  -0.1428   0.0437     
# HabitatTropical    -0.1066  0.0333  -3.1993  0.0014  -0.1719  -0.0413  ** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Polar"=1,"Temperate"=1,"Tropical"=1), 
                                   type="Tukey"))), test=adjusted("none"))
# Estimate Std. Error z value Pr(>|z|)   
# Temperate - Polar == 0    -0.44124    0.19204  -2.298  0.02158 * 
# Tropical - Polar == 0     -0.49826    0.18901  -2.636  0.00839 **
# Tropical - Temperate == 0 -0.05702    0.05805  -0.982  0.32597   



############## Trait ####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
# TraitAbundance       0.0876  0.1130   0.7752  0.4382  -0.1338   0.3090    
# TraitAccumulation   -0.0940  0.0412  -2.2808  0.0226  -0.1747  -0.0132  * 
# TraitBehaviour      -0.1394  0.1330  -1.0482  0.2945  -0.4000   0.1212    
# TraitGrowth         -0.0478  0.0543  -0.8803  0.3787  -0.1542   0.0586    
# TraitMetabolism     -0.0650  0.0423  -1.5382  0.1240  -0.1478   0.0178    
# TraitReproduction   -0.1537  0.0995  -1.5457  0.1222  -0.3487   0.0412    
# TraitSurvival       -0.1693  0.0872  -1.9408  0.0523  -0.3403   0.0017  . 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Abundance"=1,"Accumulation"=1,"Behaviour"=1,
                                     "Growth"=1,"Metabolism"=1,"Reproduction"=1,"Survival"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                  Estimate Std. Error z value Pr(>|z|)  
# Accumulation - Abundance == 0    -0.18155    0.11914  -1.524   0.1276  
# Behaviour - Abundance == 0       -0.22695    0.17374  -1.306   0.1915  
# Growth - Abundance == 0          -0.13536    0.12145  -1.115   0.2650  
# Metabolism - Abundance == 0      -0.15258    0.11782  -1.295   0.1953  
# Reproduction - Abundance == 0    -0.24131    0.14928  -1.617   0.1060  
# Survival - Abundance == 0        -0.25686    0.14047  -1.829   0.0675 .
# Behaviour - Accumulation == 0    -0.04540    0.13711  -0.331   0.7406  
# Growth - Accumulation == 0        0.04619    0.06409   0.721   0.4711  
# Metabolism - Accumulation == 0    0.02897    0.05437   0.533   0.5941  
# Reproduction - Accumulation == 0 -0.05976    0.10471  -0.571   0.5682  
# Survival - Accumulation == 0     -0.07531    0.09359  -0.805   0.4210  
# Growth - Behaviour == 0           0.09159    0.13984   0.655   0.5125  
# Metabolism - Behaviour == 0       0.07437    0.13827   0.538   0.5907  
# Reproduction - Behaviour == 0    -0.01436    0.16392  -0.088   0.9302  
# Survival - Behaviour == 0        -0.02991    0.15523  -0.193   0.8472  
# Metabolism - Growth == 0         -0.01722    0.06393  -0.269   0.7877  
# Reproduction - Growth == 0       -0.10595    0.10569  -1.002   0.3161  
# Survival - Growth == 0           -0.12150    0.09872  -1.231   0.2184  
# Reproduction - Metabolism == 0   -0.08873    0.10596  -0.837   0.4023  
# Survival - Metabolism == 0       -0.10428    0.09499  -1.098   0.2723  
# Survival - Reproduction == 0     -0.01555    0.12575  -0.124   0.9016  


############## Metal####################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
# MetalArsenic      0.0993  0.1521   0.6527  0.5139  -0.1989   0.3975      
# MetalCadmium     -0.0543  0.0538  -1.0094  0.3128  -0.1597   0.0511      
# MetalCobalt      -0.0453  0.1186  -0.3818  0.7026  -0.2777   0.1872      
# MetalCopper      -0.1579  0.0447  -3.5366  0.0004  -0.2455  -0.0704  *** 
# MetalIron         0.2114  0.1174   1.8008  0.0717  -0.0187   0.4416    . 
# MetalLead        -0.2957  0.1069  -2.7673  0.0057  -0.5052  -0.0863   ** 
# MetalManganese   -0.1028  0.1382  -0.7441  0.4568  -0.3737   0.1680      
# MetalMercury     -0.0389  0.0736  -0.5281  0.5974  -0.1831   0.1054      
# MetalNickel      -0.0092  0.1721  -0.0536  0.9573  -0.3464   0.3280      
# MetalSilver      -0.1266  0.1167  -1.0854  0.2777  -0.3553   0.1020      
# MetalZinc         0.0324  0.0869   0.3734  0.7088  -0.1378   0.2027      

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Arsenic "=1,"Cadmium"=1,"Cobalt"=1,
                                     "Copper"=1,"Iron"=1,"Lead "=1,"Manganese"=1,
                                     "Mercury"=1,"Nickel"=1,"Silver"=1,"Zinc"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                           Estimate Std. Error z value Pr(>|z|)   
# Cadmium - Arsenic  == 0   -0.153576   0.160021  -0.960  0.33719   
# Cobalt - Arsenic  == 0    -0.144582   0.192799  -0.750  0.45331   
# Copper - Arsenic  == 0    -0.257245   0.158334  -1.625  0.10423   
# Iron - Arsenic  == 0       0.112139   0.192176   0.584  0.55954   
# Lead  - Arsenic  == 0     -0.395040   0.185140  -2.134  0.03286 * 
# Manganese - Arsenic  == 0 -0.202141   0.205422  -0.984  0.32510   
# Mercury - Arsenic  == 0   -0.138169   0.168621  -0.819  0.41256   
# Nickel - Arsenic  == 0    -0.108518   0.229660  -0.473  0.63656   
# Silver - Arsenic  == 0    -0.225937   0.191373  -1.181  0.23776   
# Zinc - Arsenic  == 0      -0.066865   0.174627  -0.383  0.70179   
# Cobalt - Cadmium == 0      0.008994   0.128931   0.070  0.94438   
# Copper - Cadmium == 0     -0.103669   0.068496  -1.513  0.13015   
# Iron - Cadmium == 0        0.265715   0.129133   2.058  0.03962 * 
# Lead  - Cadmium == 0      -0.241463   0.118831  -2.032  0.04216 * 
# Manganese - Cadmium == 0  -0.048565   0.146987  -0.330  0.74110   
# Mercury - Cadmium == 0     0.015408   0.090827   0.170  0.86529   
# Nickel - Cadmium == 0      0.045058   0.180244   0.250  0.80260   
# Silver - Cadmium == 0     -0.072360   0.122149  -0.592  0.55359   
# Zinc - Cadmium == 0        0.086712   0.101014   0.858  0.39067   
# Copper - Cobalt == 0      -0.112663   0.126491  -0.891  0.37310   
# Iron - Cobalt == 0         0.256721   0.166539   1.542  0.12319   
# Lead  - Cobalt == 0       -0.250458   0.159551  -1.570  0.11647   
# Manganese - Cobalt == 0   -0.057559   0.149397  -0.385  0.70003   
# Mercury - Cobalt == 0      0.006414   0.138869   0.046  0.96316   
# Nickel - Cobalt == 0       0.036064   0.208909   0.173  0.86294   
# Silver - Cobalt == 0      -0.081354   0.162211  -0.502  0.61600   
# Zinc - Cobalt == 0         0.077718   0.143684   0.541  0.58858   
# Iron - Copper == 0         0.369384   0.125613   2.941  0.00328 **
# Lead  - Copper == 0       -0.137794   0.114342  -1.205  0.22816   
# Manganese - Copper == 0    0.055104   0.144975   0.380  0.70388   
# Mercury - Copper == 0      0.119077   0.085956   1.385  0.16595   
# Nickel - Copper == 0       0.148727   0.177723   0.837  0.40268   
# Silver - Copper == 0       0.031309   0.123821   0.253  0.80038   
# Zinc - Copper == 0         0.190381   0.095247   1.999  0.04563 * 
# Lead  - Iron == 0         -0.507179   0.158766  -3.194  0.00140 **
# Manganese - Iron == 0     -0.314280   0.181210  -1.734  0.08286 . 
# Mercury - Iron == 0       -0.250308   0.138567  -1.806  0.07086 . 
# Nickel - Iron == 0        -0.220657   0.208297  -1.059  0.28944   
# Silver - Iron == 0        -0.338076   0.165492  -2.043  0.04107 * 
# Zinc - Iron == 0          -0.179004   0.145758  -1.228  0.21941   
# Manganese - Lead  == 0     0.192899   0.174596   1.105  0.26923   
# Mercury - Lead  == 0       0.256871   0.129480   1.984  0.04727 * 
# Nickel - Lead  == 0        0.286521   0.202528   1.415  0.15715   
# Silver - Lead  == 0        0.169103   0.157950   1.071  0.28434   
# Zinc - Lead  == 0          0.328175   0.136801   2.399  0.01644 * 
# Mercury - Manganese == 0   0.063973   0.155550   0.411  0.68088   
# Nickel - Manganese == 0    0.093623   0.220614   0.424  0.67129   
# Silver - Manganese == 0   -0.023795   0.176026  -0.135  0.89247   
# Zinc - Manganese == 0      0.135277   0.159199   0.850  0.39547   
# Nickel - Mercury == 0      0.029650   0.187120   0.158  0.87410   
# Silver - Mercury == 0     -0.087768   0.136981  -0.641  0.52170   
# Zinc - Mercury == 0        0.071304   0.113032   0.631  0.52815   
# Silver - Nickel == 0      -0.117418   0.207822  -0.565  0.57208   
# Zinc - Nickel == 0         0.041654   0.191766   0.217  0.82804   
# Zinc - Silver == 0         0.159072   0.142036   1.120  0.26274   



############## Trophic.level #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# Trophic.levelHerbivore          -0.1030  0.0349  -2.9501  0.0032  -0.1714  -0.0346   ** 
# Trophic.levelPredator           -0.1699  0.0496  -3.4224  0.0006  -0.2672  -0.0726  *** 
# Trophic.levelPrimary producer    0.1150  0.0553   2.0785  0.0377   0.0066   0.2235    * 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                   Estimate Std. Error z value Pr(>|z|)    
# Predator - Herbivore == 0         -0.06690    0.06015  -1.112 0.266009    
# Primary producer - Herbivore == 0  0.21803    0.06428   3.392 0.000694 ***
# Primary producer - Predator == 0   0.28493    0.07301   3.902 9.53e-05 ***


############## Stressor:Taxa #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
#                                  estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:TaxaAlgae+Protist       0.1570  0.0595   2.6406  0.0083   0.0405   0.2736   ** 
# StressorA-W:TaxaAlgae+Protist    -0.1098  0.2719  -0.4036  0.6865  -0.6427   0.4232      
# StressorW:TaxaAlgae+Protist       0.2450  0.1277   1.9185  0.0550  -0.0053   0.4953    . 
# StressorA:TaxaDeuterostomia      -0.1153  0.0996  -1.1574  0.2471  -0.3105   0.0799      
# StressorA-W:TaxaDeuterostomia    -0.1928  0.2009  -0.9600  0.3370  -0.5865   0.2008      
# StressorW:TaxaDeuterostomia      -0.3381  0.1167  -2.8976  0.0038  -0.5669  -0.1094   ** 
# StressorA:TaxaEcdysozoa           0.0556  0.0862   0.6445  0.5193  -0.1134   0.2245      
# StressorH:TaxaEcdysozoa          -0.4393  0.1832  -2.3985  0.0165  -0.7983  -0.0803    * 
# StressorW:TaxaEcdysozoa          -0.3305  0.0805  -4.1066  <.0001  -0.4882  -0.1727  *** 
# StressorA:TaxaLophotrochozoa     -0.0538  0.0458  -1.1757  0.2397  -0.1435   0.0359      
# StressorA-W:TaxaLophotrochozoa   -0.2163  0.1761  -1.2286  0.2192  -0.5615   0.1288      
# StressorH:TaxaLophotrochozoa     -0.2945  0.1383  -2.1303  0.0332  -0.5655  -0.0235    * 
# StressorW:TaxaLophotrochozoa     -0.0342  0.0610  -0.5604  0.5752  -0.1538   0.0854      
# StressorA:TaxaRadiata            -0.1732  0.1169  -1.4820  0.1383  -0.4023   0.0559      
# StressorW:TaxaRadiata            -0.2739  0.1018  -2.6915  0.0071  -0.4734  -0.0744   ** 


############## Stressor:Habitat #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:HabitatPolar          0.5124  0.1950   2.6276  0.0086   0.1302   0.8946   ** 
# StressorA:HabitatTemperate     -0.0327  0.0566  -0.5776  0.5635  -0.1436   0.0782      
# StressorA-W:HabitatTemperate   -0.0449  0.2092  -0.2148  0.8299  -0.4549   0.3651      
# StressorH:HabitatTemperate     -0.4021  0.1128  -3.5641  0.0004  -0.6233  -0.1810  *** 
# StressorW:HabitatTemperate      0.0915  0.0723   1.2656  0.2056  -0.0502   0.2333      
# StressorA:HabitatTropical       0.0135  0.0375   0.3597  0.7191  -0.0600   0.0870      
# StressorA-W:HabitatTropical    -0.3273  0.1337  -2.4486  0.0143  -0.5893  -0.0653    * 
# StressorH:HabitatTropical       0.2289  0.2945   0.7772  0.4370  -0.3483   0.8061      
# StressorW:HabitatTropical      -0.2931  0.0487  -6.0233  <.0001  -0.3885  -0.1977  *** 


############## Stressor:Habitat:Taxa #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
#                                                    estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:HabitatPolar:TaxaAlgae+Protist           0.7029  0.2612   2.6911  0.0071   0.1910   1.2149   ** 
# StressorA:HabitatTemperate:TaxaAlgae+Protist       0.1998  0.1189   1.6814  0.0927  -0.0331   0.4328    . 
# StressorW:HabitatTemperate:TaxaAlgae+Protist       0.2740  0.2030   1.3499  0.1771  -0.1238   0.6717      
# StressorA:HabitatTropical:TaxaAlgae+Protist        0.1431  0.0673   2.1277  0.0334   0.0113   0.2750    * 
# StressorA-W:HabitatTropical:TaxaAlgae+Protist     -0.1857  0.2701  -0.6874  0.4918  -0.7150   0.3437      
# StressorW:HabitatTropical:TaxaAlgae+Protist       -0.0433  0.1950  -0.2223  0.8241  -0.4255   0.3388      
# StressorA:HabitatTemperate:TaxaDeuterostomia      -0.3061  0.1365  -2.2430  0.0249  -0.5737  -0.0386    * 
# StressorA:HabitatTropical:TaxaDeuterostomia        0.0772  0.1359   0.5686  0.5696  -0.1890   0.3435      
# StressorA-W:HabitatTropical:TaxaDeuterostomia     -0.1232  0.1997  -0.6169  0.5373  -0.5146   0.2682      
# StressorW:HabitatTropical:TaxaDeuterostomia       -0.3148  0.1145  -2.7489  0.0060  -0.5393  -0.0903   ** 
# StressorA:HabitatTemperate:TaxaEcdysozoa          -0.1825  0.2319  -0.7869  0.4313  -0.6370   0.2720      
# StressorH:HabitatTemperate:TaxaEcdysozoa          -0.4324  0.1785  -2.4225  0.0154  -0.7823  -0.0826    * 
# StressorW:HabitatTemperate:TaxaEcdysozoa          -0.1854  0.2074  -0.8937  0.3715  -0.5920   0.2212      
# StressorA:HabitatTropical:TaxaEcdysozoa            0.0945  0.0887   1.0663  0.2863  -0.0792   0.2683      
# StressorW:HabitatTropical:TaxaEcdysozoa           -0.3494  0.0836  -4.1807  <.0001  -0.5132  -0.1856  *** 
# StressorW:HabitatPolar:TaxaLophotrochozoa          0.3046  0.2726   1.1173  0.2638  -0.2297   0.8389      
# StressorA:HabitatTemperate:TaxaLophotrochozoa     -0.0268  0.0719  -0.3731  0.7091  -0.1678   0.1141      
# StressorA-W:HabitatTemperate:TaxaLophotrochozoa   -0.0417  0.2074  -0.2009  0.8408  -0.4481   0.3648      
# StressorH:HabitatTemperate:TaxaLophotrochozoa     -0.4444  0.1527  -2.9101  0.0036  -0.7437  -0.1451   ** 
# StressorW:HabitatTemperate:TaxaLophotrochozoa      0.0912  0.0778   1.1730  0.2408  -0.0612   0.2437      
# StressorA:HabitatTropical:TaxaLophotrochozoa      -0.0645  0.0560  -1.1530  0.2489  -0.1742   0.0452      
# StressorA-W:HabitatTropical:TaxaLophotrochozoa    -0.5949  0.3253  -1.8285  0.0675  -1.2325   0.0428    . 
# StressorH:HabitatTropical:TaxaLophotrochozoa       0.2300  0.2889   0.7961  0.4260  -0.3362   0.7961      
# StressorW:HabitatTropical:TaxaLophotrochozoa      -0.2693  0.0988  -2.7244  0.0064  -0.4630  -0.0756   ** 
# StressorA:HabitatTropical:TaxaRadiata             -0.1722  0.1141  -1.5089  0.1313  -0.3959   0.0515      
# StressorW:HabitatTropical:TaxaRadiata             -0.2744  0.0992  -2.7652  0.0057  -0.4689  -0.0799   ** 


############## Stressor:Trophic.level #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
#                                           estimate      se     zval    pval    ci.lb    ci.ub     
# StressorA:Trophic.levelHerbivore          -0.0254  0.0465  -0.5462  0.5849  -0.1165   0.0657     
# StressorA-W:Trophic.levelHerbivore        -0.2588  0.1856  -1.3943  0.1632  -0.6225   0.1050     
# StressorH:Trophic.levelHerbivore          -0.2582  0.1377  -1.8751  0.0608  -0.5281   0.0117   . 
# StressorW:Trophic.levelHerbivore          -0.1761  0.0540  -3.2632  0.0011  -0.2819  -0.0703  ** 
# StressorA:Trophic.levelPredator           -0.1056  0.0653  -1.6173  0.1058  -0.2335   0.0224     
# StressorA-W:Trophic.levelPredator         -0.2889  0.1480  -1.9517  0.0510  -0.5790   0.0012   . 
# StressorH:Trophic.levelPredator           -0.4346  0.1823  -2.3837  0.0171  -0.7919  -0.0772   * 
# StressorW:Trophic.levelPredator           -0.1876  0.0717  -2.6171  0.0089  -0.3282  -0.0471  ** 
# StressorA:Trophic.levelPrimary producer    0.1227  0.0612   2.0033  0.0451   0.0027   0.2427   * 
# StressorW:Trophic.levelPrimary producer    0.1163  0.1207   0.9635  0.3353  -0.1203   0.3530     

############## Habitat:Trophic.level #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
#                                                estimate      se     zval    pval    ci.lb    ci.ub 
# HabitatPolar:Trophic.levelHerbivore              0.3072  0.3012   1.0199  0.3078  -0.2831   0.8975 
# HabitatTemperate:Trophic.levelHerbivore          -0.0881  0.0607  -1.4522  0.1464  -0.2070   0.0308 
# HabitatTropical:Trophic.levelHerbivore           -0.1197  0.0430  -2.7818  0.0054  -0.2041  -0.0354 
# HabitatTemperate:Trophic.levelPredator           -0.0952  0.0823  -1.1566  0.2474  -0.2566   0.0661 
# HabitatTropical:Trophic.levelPredator            -0.2146  0.0622  -3.4516  0.0006  -0.3365  -0.0927 
# HabitatPolar:Trophic.levelPrimary producer       0.4422  0.2249   1.9661  0.0493   0.0014   0.8831 
# HabitatTemperate:Trophic.levelPrimary producer    0.1878  0.1183   1.5876  0.1124  -0.0441   0.4197 
# HabitatTropical:Trophic.levelPrimary producer     0.0621  0.0652   0.9529  0.3406  -0.0656   0.1898 
# 
# HabitatPolar:Trophic.levelHerbivore                
# HabitatTemperate:Trophic.levelHerbivore             
# HabitatTropical:Trophic.levelHerbivore           ** 
# HabitatTemperate:Trophic.levelPredator              
# HabitatTropical:Trophic.levelPredator           *** 
# HabitatPolar:Trophic.levelPrimary producer       * 
# HabitatTemperate:Trophic.levelPrimary producer      
# HabitatTropical:Trophic.levelPrimary producer         



############## Trophic.level:Metal #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
#                                               estimate      se     zval    pval    ci.lb 
# Trophic.levelHerbivore:MetalArsenic           0.1834  0.1792   1.0234  0.3061  -0.1678 
# Trophic.levelPrimary producer:MetalArsenic    0.2078  0.1859   1.1180  0.2636  -0.1565 
# Trophic.levelHerbivore:MetalCadmium          -0.1485  0.0627  -2.3701  0.0178  -0.2714 
# Trophic.levelPredator:MetalCadmium           -0.0431  0.1224  -0.3522  0.7247  -0.2831 
# Trophic.levelPrimary producer:MetalCadmium    0.3751  0.1389   2.6994  0.0069   0.1027 
# Trophic.levelHerbivore:MetalCobalt           -0.0408  0.1923  -0.2121  0.8320  -0.4176 
# Trophic.levelPredator:MetalCobalt            -0.1004  0.1576  -0.6372  0.5240  -0.4094 
# Trophic.levelPrimary producer:MetalCobalt     0.1602  0.3059   0.5237  0.6005  -0.4394 
# Trophic.levelHerbivore:MetalCopper           -0.1920  0.0586  -3.2793  0.0010  -0.3068 
# Trophic.levelPredator:MetalCopper            -0.2190  0.0823  -2.6611  0.0078  -0.3803 
# Trophic.levelPrimary producer:MetalCopper    -0.0242  0.0869  -0.2789  0.7803  -0.1946 
# Trophic.levelPrimary producer:MetalIron       0.2046  0.1113   1.8380  0.0661  -0.0136 
# Trophic.levelHerbivore:MetalLead             -0.1622  0.1140  -1.4222  0.1550  -0.3856 
# Trophic.levelPredator:MetalLead              -0.8727  0.2518  -3.4656  0.0005  -1.3662 
# Trophic.levelHerbivore:MetalManganese        -0.0076  0.2816  -0.0268  0.9786  -0.5595 
# Trophic.levelPredator:MetalManganese         -0.1569  0.1602  -0.9800  0.3271  -0.4709 
# Trophic.levelHerbivore:MetalMercury           0.0424  0.0886   0.4786  0.6322  -0.1312 
# Trophic.levelPredator:MetalMercury           -0.1604  0.1124  -1.4269  0.1536  -0.3806 
# Trophic.levelHerbivore:MetalNickel           -0.0069  0.1714  -0.0403  0.9679  -0.3428 
# Trophic.levelHerbivore:MetalSilver           -0.1628  0.1417  -1.1490  0.2506  -0.4405 
# Trophic.levelPredator:MetalSilver            -0.1348  0.1989  -0.6776  0.4980  -0.5246 
# Trophic.levelHerbivore:MetalZinc              0.0690  0.1408   0.4905  0.6238  -0.2069 
# Trophic.levelPredator:MetalZinc              -0.0915  0.1366  -0.6696  0.5031  -0.3593 
# Trophic.levelPrimary producer:MetalZinc       0.1840  0.1705   1.0791  0.2805  -0.1502 
#                                               ci.ub      
# Trophic.levelHerbivore:MetalArsenic          0.5345      
# Trophic.levelPrimary producer:MetalArsenic   0.5721      
# Trophic.levelHerbivore:MetalCadmium         -0.0257    * 
# Trophic.levelPredator:MetalCadmium           0.1968      
# Trophic.levelPrimary producer:MetalCadmium   0.6474   ** 
# Trophic.levelHerbivore:MetalCobalt           0.3360      
# Trophic.levelPredator:MetalCobalt            0.2085      
# Trophic.levelPrimary producer:MetalCobalt    0.7599      
# Trophic.levelHerbivore:MetalCopper          -0.0772   ** 
# Trophic.levelPredator:MetalCopper           -0.0577   ** 
# Trophic.levelPrimary producer:MetalCopper    0.1462      
# Trophic.levelPrimary producer:MetalIron      0.4228    . 
# Trophic.levelHerbivore:MetalLead             0.0613      
# Trophic.levelPredator:MetalLead             -0.3791  *** 
# Trophic.levelHerbivore:MetalManganese        0.5443      
# Trophic.levelPredator:MetalManganese         0.1570      
# Trophic.levelHerbivore:MetalMercury          0.2160      
# Trophic.levelPredator:MetalMercury           0.0599      
# Trophic.levelHerbivore:MetalNickel           0.3290      
# Trophic.levelHerbivore:MetalSilver           0.1149      
# Trophic.levelPredator:MetalSilver            0.2550      
# Trophic.levelHerbivore:MetalZinc             0.3449      
# Trophic.levelPredator:MetalZinc              0.1763      
# Trophic.levelPrimary producer:MetalZinc      0.5183      




############## Stressor:Metal #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
# StressorA:MetalArsenic      0.3554  0.2069   1.7179  0.0858  -0.0501   0.7609   . 
# StressorW:MetalArsenic     -0.1271  0.2042  -0.6223  0.5338  -0.5273   0.2732     
# StressorA:MetalCadmium      0.0233  0.0647   0.3607  0.7183  -0.1035   0.1502     
# StressorA-W:MetalCadmium   -0.2207  0.2485  -0.8883  0.3744  -0.7078   0.2663     
# StressorH:MetalCadmium     -0.3527  0.1284  -2.7466  0.0060  -0.6045  -0.1010  ** 
# StressorW:MetalCadmium     -0.0299  0.1157  -0.2586  0.7960  -0.2567   0.1969     
# StressorA:MetalCobalt      -0.0762  0.1359  -0.5609  0.5749  -0.3425   0.1901     
# StressorW:MetalCobalt       0.2688  0.1929   1.3935  0.1635  -0.1093   0.6468     
# StressorA:MetalCopper      -0.1184  0.0559  -2.1188  0.0341  -0.2279  -0.0089   * 
# StressorA-W:MetalCopper    -0.2883  0.2759  -1.0450  0.2960  -0.8290   0.2524     
# StressorW:MetalCopper      -0.2115  0.0671  -3.1499  0.0016  -0.3430  -0.0799  ** 
# StressorA:MetalIron         0.1405  0.1391   1.0100  0.3125  -0.1321   0.4130     
# StressorW:MetalIron         0.3341  0.1925   1.7357  0.0826  -0.0432   0.7113   . 
# StressorA:MetalLead        -0.1279  0.1744  -0.7336  0.4632  -0.4697   0.2139     
# StressorW:MetalLead        -0.3797  0.1296  -2.9292  0.0034  -0.6337  -0.1256  ** 
# StressorA:MetalManganese    0.1361  0.1777   0.7660  0.4437  -0.2121   0.4843     
# StressorW:MetalManganese    0.0854  0.2247   0.3802  0.7038  -0.3550   0.5259     
# StressorA:MetalMercury      0.1597  0.0976   1.6361  0.1018  -0.0316   0.3511     
# StressorA-W:MetalMercury   -0.1961  0.1763  -1.1120  0.2662  -0.5416   0.1495     
# StressorW:MetalMercury     -0.1994  0.0960  -2.0776  0.0377  -0.3875  -0.0113   * 
# StressorA:MetalNickel       0.1993  0.2406   0.8283  0.4075  -0.2723   0.6709     
# StressorW:MetalNickel      -0.1954  0.2237  -0.8736  0.3823  -0.6337   0.2430     
# StressorA:MetalSilver      -0.0989  0.1298  -0.7616  0.4463  -0.3533   0.1556     
# StressorW:MetalSilver      -0.1017  0.2481  -0.4100  0.6818  -0.5879   0.3845     
# StressorA:MetalZinc         0.0208  0.0938   0.2218  0.8244  -0.1630   0.2046     
# StressorW:MetalZinc         0.0093  0.2508   0.0371  0.9704  -0.4822   0.5009     



############## Stressor:Trait #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
data$Stressor<-factor(data$Stressor)
data$Trait<-factor(data$Trait)
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
# StressorA:TraitAbundance         0.2317  0.1681   1.3784  0.1681  -0.0978   0.5612      
# StressorW:TraitAbundance        -0.0146  0.1503  -0.0974  0.9224  -0.3092   0.2800      
# StressorA:TraitAccumulation     -0.0087  0.0475  -0.1841  0.8539  -0.1018   0.0843      
# StressorA-W:TraitAccumulation   -0.2099  0.1633  -1.2849  0.1988  -0.5300   0.1103      
# StressorH:TraitAccumulation     -0.2373  0.1386  -1.7124  0.0868  -0.5089   0.0343    . 
# StressorW:TraitAccumulation     -0.2307  0.0751  -3.0734  0.0021  -0.3779  -0.0836   ** 
# StressorA:TraitBehaviour         0.0083  0.1971   0.0423  0.9663  -0.3780   0.3947      
# StressorW:TraitBehaviour        -0.1627  0.1991  -0.8174  0.4137  -0.5530   0.2275      
# StressorA:TraitGrowth            0.0172  0.0626   0.2752  0.7831  -0.1055   0.1400      
# StressorW:TraitGrowth           -0.1278  0.0994  -1.2860  0.1984  -0.3227   0.0670      
# StressorA:TraitMetabolism       -0.0480  0.0541  -0.8885  0.3743  -0.1540   0.0579      
# StressorA-W:TraitMetabolism     -0.5511  0.2196  -2.5098  0.0121  -0.9814  -0.1207    * 
# StressorH:TraitMetabolism       -0.4198  0.2078  -2.0203  0.0433  -0.8271  -0.0125    * 
# StressorW:TraitMetabolism       -0.0183  0.0661  -0.2774  0.7815  -0.1479   0.1112      
# StressorA:TraitReproduction      0.1452  0.1345   1.0794  0.2804  -0.1184   0.4088      
# StressorW:TraitReproduction     -0.4766  0.1432  -3.3280  0.0009  -0.7572  -0.1959  *** 
# StressorA:TraitSurvival         -0.0756  0.1447  -0.5220  0.6017  -0.3593   0.2081      
# StressorW:TraitSurvival         -0.2607  0.1093  -2.3841  0.0171  -0.4750  -0.0464    * 
  
############## Stressor:Trait:Habitat #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
#                                                   estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:TraitMetabolism:HabitatPolar            0.6562  0.2345   2.7982  0.0051   0.1966   1.1159   ** 
# StressorA:TraitAbundance:HabitatTemperate         0.3697  0.2027   1.8243  0.0681  -0.0275   0.7669    . 
# StressorA:TraitAccumulation:HabitatTemperate     -0.0600  0.0763  -0.7860  0.4318  -0.2096   0.0896      
# StressorA-W:TraitAccumulation:HabitatTemperate   -0.1557  0.2620  -0.5941  0.5525  -0.6692   0.3579      
# StressorH:TraitAccumulation:HabitatTemperate     -0.3771  0.1384  -2.7243  0.0064  -0.6485  -0.1058   ** 
# StressorW:TraitAccumulation:HabitatTemperate      0.1825  0.1164   1.5688  0.1167  -0.0455   0.4106      
# StressorA:TraitBehaviour:HabitatTemperate        -0.1214  0.3095  -0.3924  0.6948  -0.7280   0.4851      
# StressorA:TraitGrowth:HabitatTemperate           -0.0234  0.1380  -0.1699  0.8651  -0.2938   0.2470      
# StressorW:TraitGrowth:HabitatTemperate           -0.2985  0.1900  -1.5712  0.1161  -0.6709   0.0739      
# StressorA:TraitMetabolism:HabitatTemperate        0.0140  0.1203   0.1164  0.9073  -0.2219   0.2499      
# StressorH:TraitMetabolism:HabitatTemperate       -0.3957  0.2430  -1.6283  0.1035  -0.8721   0.0806      
# StressorW:TraitMetabolism:HabitatTemperate        0.2224  0.1106   2.0115  0.0443   0.0057   0.4391    * 
# StressorW:TraitSurvival:HabitatTemperate         -0.2054  0.2213  -0.9282  0.3533  -0.6393   0.2284      
# StressorA:TraitAbundance:HabitatTropical         -0.0883  0.2817  -0.3133  0.7540  -0.6403   0.4638      
# StressorW:TraitAbundance:HabitatTropical         -0.1572  0.1531  -1.0268  0.3045  -0.4571   0.1428      
# StressorA:TraitAccumulation:HabitatTropical       0.0401  0.0570   0.7026  0.4823  -0.0717   0.1519      
# StressorA-W:TraitAccumulation:HabitatTropical    -0.2084  0.1964  -1.0612  0.2886  -0.5933   0.1765      
# StressorW:TraitAccumulation:HabitatTropical      -0.5180  0.0955  -5.4262  <.0001  -0.7051  -0.3309  *** 
# StressorA:TraitBehaviour:HabitatTropical          0.0733  0.2427   0.3019  0.7627  -0.4024   0.5489      
# StressorW:TraitBehaviour:HabitatTropical         -0.2077  0.2058  -1.0091  0.3129  -0.6111   0.1957      
# StressorA:TraitGrowth:HabitatTropical             0.0366  0.0680   0.5378  0.5907  -0.0967   0.1698      
# StressorW:TraitGrowth:HabitatTropical            -0.0904  0.1103  -0.8194  0.4125  -0.3065   0.1258      
# StressorA:TraitMetabolism:HabitatTropical        -0.0590  0.0583  -1.0111  0.3120  -0.1733   0.0553      
# StressorA-W:TraitMetabolism:HabitatTropical      -0.5888  0.2120  -2.7777  0.0055  -1.0043  -0.1733   ** 
# StressorW:TraitMetabolism:HabitatTropical        -0.2239  0.0822  -2.7217  0.0065  -0.3851  -0.0627   ** 
# StressorA:TraitReproduction:HabitatTropical       0.2056  0.1352   1.5210  0.1283  -0.0593   0.4705      
# StressorW:TraitReproduction:HabitatTropical      -0.5108  0.1381  -3.6994  0.0002  -0.7815  -0.2402  *** 
# StressorA:TraitSurvival:HabitatTropical          -0.0165  0.1408  -0.1173  0.9066  -0.2925   0.2595      
# StressorW:TraitSurvival:HabitatTropical          -0.3051  0.1192  -2.5607  0.0104  -0.5387  -0.0716    * 


############## Trophic.level-Accumulation #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
#                                 estimate      se     zval    pval    ci.lb    ci.ub 
# Trophic.levelHerbivore          -0.0724  0.0589  -1.2285  0.2193  -0.1878   0.0431 
# Trophic.levelPredator           -0.1908  0.0817  -2.3353  0.0195  -0.3509  -0.0307 
# Trophic.levelPrimary producer    0.3459  0.1472   2.3504  0.0188   0.0575   0.6344 
# 
# Trophic.levelHerbivore           
# Trophic.levelPredator          * 
# Trophic.levelPrimary producer  * 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                    Estimate Std. Error z value Pr(>|z|)   
# Predator - Herbivore == 0          -0.1184     0.1007  -1.176  0.23969   
# Primary producer - Herbivore == 0   0.4183     0.1512   2.766  0.00568 **
# Primary producer - Predator == 0    0.5367     0.1683   3.188  0.00143 **

############## Trophic.level-Abundance #################### 
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data <- subset(data, Trait == "Abundance")
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
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
#                                  estimate      se     zval    pval    ci.lb   ci.ub 
# Trophic.levelHerbivore          -0.1820  0.4154  -0.4381  0.6613  -0.9961  0.6321 
# Trophic.levelPredator           -0.0873  0.3618  -0.2412  0.8094  -0.7964  0.6219 
# Trophic.levelPrimary producer    0.2546  0.2869   0.8875  0.3748  -0.3077  0.8169 
# 
# Trophic.levelHerbivore           
# Trophic.levelPredator            
# Trophic.levelPrimary producer    
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                    Estimate Std. Error z value Pr(>|z|)
# Predator - Herbivore == 0          0.09472    0.55085   0.172    0.863
# Primary producer - Herbivore == 0  0.43659    0.50479   0.865    0.387
# Primary producer - Predator == 0   0.34187    0.46175   0.740    0.459












