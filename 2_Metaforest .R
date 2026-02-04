rm(list = ls())

library(metaforest)   
library(caret)        
library(tidyverse)    
library(metafor)      
library(readxl)       
library(metadat)       
library(data.table)    
library(ggplot2)     
library(Matrix)       
library(parallel)   
library(doParallel)    

################################################
#                                              #
#              Metaforest Fig. 2.c             #
#                                              #
################################################
data<-read.csv ("Supplementary Data 1.csv", header =T)
nrow(data)
data$Stressor <- as.factor(data$Stressor)
data$Stressor.time  <- as.factor(data$Stressor.time )
data$Taxa <- as.factor(data$Taxa)
data$Trophic.level <- as.factor(data$Trophic.level)
data$Calcific.feature <- as.factor(data$Calcific.feature)
data$Habitat <- as.factor(data$Habitat)
data$Metal <- as.factor(data$Metal)
data$Trait <- as.factor(data$Trait)
data$Mobility.and.mode.of.life <- as.factor(data$Mobility.and.mode.of.life)
data$Research.setting <- as.factor(data$Research.setting)


moderators <- c("Stressor","Stressor.time  ","Taxa","Trophic.level",
                "Calcific.feature ","Habitat","Metal","Trait ", "Mobility.and.mode.of.life", "Research.setting") 
newnames <- c("Stressor","Stressor.time  ","Taxa","Trophic.level",
              "Calcific.feature ","Habitat","Metal","Trait ", "Mobility.and.mode.of.life", "Research.setting")
names(newnames) <- moderators
set.seed(123)
# Check how many iterations metaforest needs to converge
check_conv<- MetaForest(as.formula(paste0("yi~",
                                          paste(moderators, collapse = "+"))),
                        data = data,
                        whichweights = "random",
                        num.trees = 5000)

# Plot convergence trajectory
plot(check_conv)
# Perform recursive preselection
set.seed(123)
pre.rel<- preselect(check_conv, replications =100,
                    algorithm = "recursive")
pre.rel.p <- plot(pre.rel, label_elements = newnames)
plot(pre.rel.p)
ggsave("preselected.pdf",pre.rel.p  , width = 6.875, units = "in")  




set.seed(123)
cv_folds <- trainControl(method = "cv", 10,allowParallel = TRUE)
tuning_grid <- expand.grid(whichweights = c("random", "fixed", "unif"),
                           mtry = 2:6,
                           min.node.size = 2:6)
X <- dplyr::select(data, "vi", preselect_vars(pre.rel, cutoff = .95))
set.seed(123)
mf_cv<- train(y = data$yi,
              x = X,
              method = ModelInfo_mf(),
              trControl = cv_folds,
              tuneGrid = tuning_grid,
              keep.inbag = TRUE,
              num.threads = 6,
              verbose=TRUE,
              num.trees = 5000)
saveRDS(mf_cv , "mf_cv.RData")
mf_cv <- readRDS("mf_cv.RData")
final<- mf_cv$finalModel
importance.export <- varImp(mf_cv)$importance
write.csv(importance.export,"VI.csv")
imp<- VarImpPlot(final,label_elements =newnames )
imp
ggsave("VI_metaforest.pdf", imp, width =  6.875, units = "in")
ggsave("VI_metaforest.png", imp, width =  6.875, units = "in")
