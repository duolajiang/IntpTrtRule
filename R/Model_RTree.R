#' @export
Model_RTree <- R6::R6Class(
  "Model_RTree",
  inherit = Model_Tree,
  public = list(
    initialize = function(data,var_names){
      super$initialize(data,var_names)
    },

    goodness_of_fit = function(){
      n <- dim(self$data)[1]
      y <- self$data$y
      y.pred <- self$datay_pred
      wgss <- sum(self$data$y.var)
      bgss <- sum((y-y.pred)^2)
      tree_score <- (wgss + bgss)/n

      return(tree_score)
    },

    goodness_of_fit_test = function(y.pred,y,y.var){
      n <- length(y.pred)
      wgss <- sum(y.var)
      bgss <- sum((y-y.pred)^2)
      tree_score <- (wgss + bgss)/n
      return(tree_score)
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
