#' @export
Model_CTree <- R6::R6Class(
  "Model_CTree",
  inherit = Model_Tree,
  public = list(
    initialize = function(data,var_names){
      super$initialize(data,var_names)
    },

    goodness_of_fit = function(){
      y <- self$data$y
      y.pred <- self$data$y_pred
      tree_score <- -sum(y*log(y.pred)+(1-y)*log(1-y.pred))
      return(tree_score)
    },

    goodness_of_fit_test = function(y.pred,y){
      n <- length(y.pred)
      tree_score <- -sum(y*log(y.pred)+(1-y)*log(1-y.pred))/n
      return(tree_score)
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
