TrtRule_addcost <- R6::R6Class(
  "TrtRule_addcost",
  inherit = TrtRule,

  public = list(

    cost = NA,

    initialize = function(data,var_names,y1.hat,y0.hat){
      super$initialize(data,var_names,y1.hat,y0.hat)
      self$data$y1.hat.mean <- apply(y1.hat,2,mean)
      self$data$y0.hat.mean <- apply(y0.hat,2,mean)
      self$data$y1.hat.var <- apply(y1.hat,2,var)
      self$data$y0.hat.var <- apply(y0.hat,2,var)

      self$data$y <- self$data$y1.hat.mean - self$data$y0.hat.mean
      self$data$y.var <- self$data$y1.hat.var + self$data$y0.hat.var
    },

    fit = function(min.ndsz, pt, max.depth, cost){
      if (is.na(tree)) {
        self$constraints$min.ndsz <- min.ndsz
        self$constraints$pt <- pt
        self$constraints$max.depth <- max.depth
        self$cost <- cost

        self$tree <- private$fit_implementation()
        tree_pruned <- self$tree
        tree.pruned$trt <- 1*(tree_pruned$trt.effect > self$cost)
        self$tree_pruned <- private$prune_search(1,tree.pruned)
      } else {
        self$cost <- cost
        tree_pruned <- self$tree
        tree.pruned$trt <- 1*(tree_pruned$trt.effect > self$cost)
        self$tree_pruned <- private$prune_search(1,tree.pruned)
      }
    }

  ),

  private = list(

    Score_node = function(person_id){
      n <- length(person_id)
      y <- self$data[self$data$person_id %in% person_id,"y"]
      y.mean <- mean(y)
      #compute variance of this group treatment effect
      wgss <- sum(self$data[self$data$person_id %in% person_id,"y.var"])
      bgss <- sum((y-y.mean)^2)
      node_ss <- wgss + bgss

      return(node_ss)

    }


  )
)
