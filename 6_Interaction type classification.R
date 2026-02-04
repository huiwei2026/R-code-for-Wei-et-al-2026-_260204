library(readxl)      
library(readr) 
################################################
#                                              #
#       Calculation of interaction types       #
#                                              #
################################################
data<-read_excel ("Supplementary Data 2.xlsx")
nrow(data)
numeric_cols <- c("sd1", "sd2", "avg1", "avg2", "n1", "n2")
data[numeric_cols] <- lapply(data[numeric_cols], as.numeric)
# Calculate effect sizes (yi) and sampling variances (vi)
lnRR <- log(data$avg2/data$avg1)
var_lnRR <- (data$sd2)^2/((data$n2)*(data$avg2)^2) +
  (data$sd1)^2/((data$n1)*(data$avg1)^2)
data$yi <- lnRR + (1/2) * (((data$sd2)^2/((data$n2)*(data$avg2)^2) -
                              (data$sd1)^2/((data$n1)*(data$avg1)^2)))
data$vi <- var_lnRR + (1/2) * (((data$sd2)^4/((data$n2)^2*(data$avg2)^4) +
                                  (data$sd1)^4/((data$n1)^2*(data$avg1)^4)))
data$yi[data$Fitness==1] <- -1*data$yi[data$Fitness==1]


# Calculate individual_effect
data$individual_effect<- log(data$avg2/data$avg1)
n_groups <- nrow(data) / 3
results <- vector("list", n_groups)
data$interaction_effect <- NA
data$variance <- NA

for (i in 1:n_groups) {
  start_row <- (i-1)*3 + 1
  end_row <- i*3
  
  group_data <- data[start_row:end_row, ]
  
  # Calculate interaction_effect and variance
  a <- group_data$avg2[3]*group_data$avg1[3] 
  b <- group_data$avg2[1]*group_data$avg2[2]
  interaction_effect_3 <- log(a) - log(b)
  data$interaction_effect[start_row + 2] <- interaction_effect_3  
  
  
  # interaction_effect-variance
  a <- (group_data$sd2[1])^2/(((group_data$avg2[1])^2)*group_data$n2[1])
  b <- (group_data$sd2[2])^2/(((group_data$avg2[2])^2)*group_data$n2[2])
  c <- (group_data$sd2[3])^2/(((group_data$avg2[3])^2)*group_data$n2[3])
  d <- (group_data$sd1[3])^2/(((group_data$avg1[3])^2)*group_data$n1[3])
  variance_3 <- a+b+c+d
  data$variance[start_row + 2] <- variance_3
  
}
data$ci <- sqrt(data$variance)/sqrt(data$n2)
data$lci <- data$interaction_effect- (data$ci*1.96)
data$uci <- data$interaction_effect+ (data$ci*1.96)
write_excel_csv (data,"Supplementary Data 2.csv")


# Identify the interaction type
data <- read.csv("Supplementary Data 2.csv", header = TRUE)
nrow(data)

data$int_type <- NA
data$int_direction <- NA

n_groups <- nrow(data) / 3

for (i in 1:n_groups) {
  start_row <- (i-1)*3 + 1
  
  lci_row3 <- data$lci[start_row + 2]
  uci_row3 <- data$uci[start_row + 2]
  
  if (!is.na(lci_row3) && !is.na(uci_row3) && lci_row3 <= 0 && uci_row3 >= 0) {
    data$int_type[start_row + 2] <- "add" 
  }
  
  ind_eff_row1 <- data$individual_effect[start_row]
  ind_eff_row2 <- data$individual_effect[start_row + 1]
  
  sum_ind_effect <- ind_eff_row1 + ind_eff_row2
  if (!is.na(sum_ind_effect)) {
    if (sum_ind_effect < 0) {
      # If sum is negative, mark as "neg" in the third row
      data$int_direction[start_row + 2] <- "neg"
    } else {
      # If sum is positive or zero, mark as "pos" in the third row
      data$int_direction[start_row + 2] <- "pos"
    }
  }
}

write.csv(data, "Supplementary Data 2.csv", row.names = FALSE)








