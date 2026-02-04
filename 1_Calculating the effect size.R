library(readxl)     
library(readr) 
################################################
#                                              #
#    Calculating the effect size Fig. 2.a      #
#                                              #
################################################
data<-read_excel ("Supplementary Data 1.xlsx")
nrow(data)
lnRR <- log(data$x_gc_m/data$x_m)
var_lnRR <- (data$sd_gc_m)^2/((data$n_gc_m)*(data$x_gc_m)^2) +
  (data$sd_m)^2/((data$n_m)*(data$x_m)^2)
data$yi <- lnRR + (1/2) * (((data$sd_gc_m)^2/((data$n_gc_m)*(data$x_gc_m)^2) - 
                              (data$sd_m)^2/((data$n_m)*(data$x_m)^2)))
data$vi <- var_lnRR + (1/2) * (((data$sd_gc_m)^4/((data$n_gc_m)^2*(data$x_gc_m)^4) + 
                                  (data$sd_m)^4/((data$n_m)^2*(data$x_m)^4))) 
data$yi[data$Fitness==1] <- -1*data$yi[data$Fitness==1]

