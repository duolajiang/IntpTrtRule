TrtRule_GenCost <- R6::R6Class(
  "TrtRule_GenCost",
  inherit = TrtRule,

  public = list(

    cost = NA,

    initialize = function(data,var_names,y1.hat,y0.hat){
      super$initialize(data,var_names,y1.hat,y0.hat)
      self$data$y1.hat.mean <- apply(y1.hat,2,mean)
      self$data$y0.hat.mean <- apply(y0.hat,2,mean)
      self$data$y1.hat.var <- apply(y1.hat,2,var)
      self$data$y0.hat.var <- apply(y0.hat,2,var)
    },

    est_OptimalTreatment = function(rho,cost_t,cost_c){
      self$cost <- list(cost_t=cost_t,
                        cost_c=cost_c)
      p11 <- rho*sqrt(self$y1.hat*(1-self$y1.hat)*self$y0.hat*(1-self$y0.hat)) + self$y1.hat*self$y0.hat
      p10 <- self$y1.hat - p11
      p01 <- self$y0.hat - p11
      p00 <- 1- p11 - p10 - p01
      deltaCost <- cost_t-cost_c
      DeltaPosteriorCost <- deltaCost[1]*p00 + deltaCost[2]*p01 + deltaCost[3]*p10 + deltaCost[4]*p11
      OptimalDecisionProb <- apply(DeltaPosteriorCost<0, 2, mean)
      self$data$y <- OptimalDecisionProb
    },

    fit = function(min.ndsz, pt, max.depth){
      self$constraints$min.ndsz <- min.ndsz
      self$constraints$pt <- pt
      self$constraints$max.depth <- max.depth

      self$tree <- private$fit_implementation()
      tree_pruned <- self$tree
      tree.pruned$trt <- 1*(tree_pruned$mu > 0.5)
      self$tree_pruned <- private$prune_search(1,tree.pruned)
    }
  ),

  private = list(

    Score_node = function(person_id){
      y <- self$data[self$data$person_id %in% person_id, "y"]
      mu <- mean(y)
      score <- -sum(y*log(mu)+(1-y)*log(1-mu))
      return(score)
    }



  )
)
