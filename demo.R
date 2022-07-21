library(data.table)
source('~/Documents/iknl/code/function/function.R')
source('~/Documents/iknl/code/function/policy_evaluation_function.R')
source('~/Documents/iknl/code/policy_evaluation/main_function.R')


data <- DataImportGeneral(year=20156)
data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                         ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))

target.index <- ((data$BRAF==1)|(data$BRAF==0)) & (data$pN==0) & (data$MSI!=9)
data <- data[target.index,]
data <- droplevels(data)

load('~/Documents/iknl/code/policy_evaluation/post_y1_y0_2fp_20156.RData')
length_test <- (dim(post_y1_y0)[2])/2
p1_ <- post_y1_y0[,c(1:length_test)]
p_1 <- post_y1_y0[,c((length_test+1):(2*length_test))]
p1_ <- p1_[,target.index]
p_1 <- p_1[,target.index]

library(IntpTrtRule)
# use additive loss function.
#' initialize an object of class TrtRule_addcost.
#' @param data a patient-level data.frame N \times p,
#' N is the number of patients, and p is the dimension
#' of variables. The data.frame contains covariates,
#' treatment variable, and outcome variable.
#' @param var_names a list with four named elements.
#' Input elements from users are xvars, trt, y.obs
#' xvars is a character containing variables for modeling,
#' trt is a character indicating treatment variable,
#' and y.obs is a character indicating outcome name.
#' @param y1.hat,y0.hat a d \times N data.frame,
#' d is the number of draws from bayesian method.
#' y1.hat is the predicted potential outcome under treatment,
#' y0.hat is the predicted potential outcome under control.
obj <- TrtRule_addcost$new(data = data,
                           var_names = list(xvars=c("BRAF","MSI","pT_num",
                                                    "age_at_diagnosis"),
                                            trt="combined_chemo",
                                            y.obs="vitstat"),
                           y1.hat = p1_,
                           y0.hat = p_1)

obj$fit(min.ndsz=10, pt=0.05, max.depth=6, cost= 0.03)
#' print subgroups in tree or tree_pruned.
#' @param pruned a logical, indicating print
#' self$tree or self$tree_pruned. FALSE by default.
obj$Print_leaf_nodes()
#' plot the tree
#' @param pruned a logical, indicating print
#' self$tree or self$tree_pruned. TRUE by default.
obj$plot_sankey(TRUE)
#' compute feature importance via permute feature,
#' grow a tree model,
#' and compute difference in score between
#' unpermute data and permuted data.
obj$feature_importance(bootloop=10)


rho <- 0
L1_00 <- 1.02; L1_01 <- 1.02; L1_10 <- 0.02; L1_11 <- 0.02;
L0_00 <- 1; L0_01 <- 0; L0_10 <- 1; L0_11 <- 0;
cost_t <- c(L1_00,L1_01,L1_10,L1_11)
cost_c <- c(L0_00,L0_01,L0_10,L0_11)
obj <- TrtRule_GenCost$new(data,
                           var_names = list(xvars=c("MSI","BRAF","pT"),
                                            trt = "combined_chemo",
                                            y.obs = "vitstat"),
                           p1_,p_1)
obj$fit(10, 0.1, 5, rho=0, cost_t, cost_c)
obj$print_leaf_nodes()
obj$plot_sankey()
obj$feature_importance(20)

# why the first one is more refined than



