library(metafor)

data <- read.csv("Supplementary Data 1.csv", header = TRUE)
data$Stressor <- as.factor(data$Stressor)
data$Taxa <- as.factor(data$Taxa)
data$Trait <- as.factor(data$Trait)
data$Metal <- as.factor(data$Metal)
data$Habitat <- as.factor(data$Habitat)

data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
m1<-rma.mv(yi, vi, mods=~Stressor:Taxa,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#          df      AIC       BIC     AICc    logLik     LRT   pval         QE 
# Full    19 957.6598 1041.6087 958.9414 -459.8299                38426.5385 
# Reduced  3 978.7068  991.9619 978.7462 -486.3534 53.0470 <.0001 41706.3867 

# ΔAIC  -21.047       ΔBIC. 49.6468


data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
m1<-rma.mv(yi, vi, mods=~Stressor:Habitat,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#         df      AIC      BIC     AICc    logLik     LRT   pval         QE 
# Full    12 946.3622 999.3826 946.8822 -461.1811                39558.5057 
# Reduced  3 978.7068 991.9619 978.7462 -486.3534 50.3446 <.0001 41706.3867 

# ΔAIC  -32.3446       ΔBIC. 7.0285


data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
m1<-rma.mv(yi, vi, mods=~Stressor:Metal,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
# df      AIC       BIC     AICc    logLik     LRT   pval         QE 
# Full    33 980.5116 1126.3177 984.3873 -457.2558                35477.0049 
# Reduced  3 978.7068  991.9619 978.7462 -486.3534 58.1951 0.0015 41706.3867 

# ΔAIC  1.8048      ΔBIC. 134.3558


data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
m1<-rma.mv(yi, vi, mods=~Stressor:Trait,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#         df      AIC       BIC     AICc    logLik     LRT   pval         QE 
# Full    26 984.0408 1098.9183 986.4367 -466.0204                25910.8040 
# Reduced  3 978.7068  991.9619 978.7462 -486.3534 40.6660 0.0129 41706.3867 
# ΔAIC 5.334     ΔBIC. 106.9564







data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
m1<-rma.mv(yi, vi, mods=~Stressor:Taxa:Trait,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#          df      AIC       BIC      AICc    logLik      LRT   pval         QE 
# Full    76 986.2879 1322.0837 1008.1238 -417.1440                 15291.2390 
# Reduced  3 978.7068  991.9619  978.7462 -486.3534 138.4188 <.0001 41706.3867 
# ΔAIC  7.5811    ΔBIC. 330.1218



data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
m1<-rma.mv(yi, vi, mods=~Stressor:Taxa:Metal,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#          df      AIC       BIC     AICc    logLik      LRT   pval         QE 
# Full    72 972.8266 1290.9489 992.2933 -414.4133                 31120.5846 
# Reduced  3 978.7068  991.9619 978.7462 -486.3534 143.8801 <.0001 41706.3867 

# ΔAIC  -5.8802    ΔBIC. 298.987


data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
m1<-rma.mv(yi, vi, mods=~Stressor:Taxa:Habitat,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#          df      AIC       BIC     AICc    logLik     LRT   pval         QE 
# Full    31 953.4849 1090.4542 956.8997 -445.7425                36786.9917 
# Reduced  3 978.7068  991.9619 978.7462 -486.3534 81.2219 <.0001 41706.3867 

# ΔAIC  -25.2219   ΔBIC. 104.6553





data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
m1<-rma.mv(yi, vi, mods=~Stressor:Trait:Metal,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#         df      AIC       BIC      AICc    logLik      LRT   pval         QE 
# Full    103 989.9522 1445.0438 1032.0426 -391.9761                 14941.0458 
# Reduced   3 978.7068  991.9619  978.7462 -486.3534 188.7546 <.0001 41706.3867 

# ΔAIC  11.2454    ΔBIC.  453.0819

data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
m1<-rma.mv(yi, vi, mods=~Stressor:Trait:Habitat,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#          df      AIC       BIC     AICc    logLik      LRT   pval         QE 
# Full    48 955.2385 1167.3200 963.5789 -429.6192                 18923.7626 
# Reduced  3 978.7068  991.9619 978.7462 -486.3534 113.4683 <.0001 41706.3867 

# ΔAIC  -23.4683    ΔBIC.  175.2696


data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
m1<-rma.mv(yi, vi, mods=~Stressor:Metal:Habitat,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#          df      AIC       BIC     AICc    logLik      LRT   pval         QE 
# Full    51 971.7284 1197.0651 981.1830 -434.8642                 33625.9479 
# Reduced  3 978.7068  991.9619 978.7462 -486.3534 102.9783 <.0001 41706.3867 

# ΔAIC   -6.9784   ΔBIC. 204.9782