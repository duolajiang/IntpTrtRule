#' @export
TrtRule_GenCost <- R6::R6Class(
  "TrtRule_GenCost",
  inherit = TrtRule,

  public = list(

    cost = NA,

    initialize = function(data,var_names,y1.hat,y0.hat,y1.var,y0.var){
      super$initialize(data,var_names,y1.hat,y0.hat)
      if(is.null(dim(y1.hat))){
        self$data$y1.hat.mean <- y1.hat
        self$data$y0.hat.mean <- y0.hat
        self$data$y1.hat.var <- y1.var
        self$data$y0.hat.var <- y0.var
      } else {
        self$data$y1.hat.mean <- apply(y1.hat,2,mean)
        self$data$y0.hat.mean <- apply(y0.hat,2,mean)
        self$data$y1.hat.var <- apply(y1.hat,2,var)
        self$data$y0.hat.var <- apply(y0.hat,2,var)
      }
    },

    fit = function(min.ndsz, pt, max.depth, rho, cost_t, cost_c){
      self$constraints$min.ndsz <- min.ndsz
      self$constraints$pt <- pt
      self$constraints$max.depth <- max.depth
      self$cost <- list(cost_t=cost_t,
                        cost_c=cost_c)
      self$data$y <- private$est_OptimalTreatment(rho,cost_t,cost_c)

      treeinfo_topology <- private$fit_implementation()
      self$tree <- treeinfo_topology$out
      self$tree_topology <- treeinfo_topology$out_topology
      tree_pruned <- self$tree
      tree_pruned$trt <- 1*(tree_pruned$mu > 0.5)
      self$tree_pruned <- private$prune_search(1,tree_pruned)
    },

    feature_importance = function(bootloop){
      #browser()
      n_sam <- nrow(self$data)
      Importance.permute <- NULL

      pb <- txtProgressBar(min = 0, max = bootloop, style = 3)
      #set.seed(123)
      for (i in seq(bootloop)) {
        train_id <- data.frame(person_id=sample(self$data$person_id,n_sam,replace = TRUE))
        train <- merge(self$data, train_id, by="person_id")
        obj <- Model_CTree$new(data = train,
                               var_names = self$var_names)
        obj$fit(self$constraints$min.ndsz,
                self$constraints$pt,
                self$constraints$max.depth)
        #2. compute J(T_b)
        score.bench <- obj$goodness_of_fit()
        #3. permute x_{j} then compute J(T_b)_j, and compute J(T_b)-J(T_b)_j
        importance.permute <- NULL
        for (xvar in self$var_names$xvars) {
          test.permute <- train
          permute.var <- test.permute[,xvar]
          permute.var <- permute.var[sample(1:n_sam,n_sam, replace=T)]
          test.permute[,xvar] <- permute.var
          test.pred <- obj$predict(test.permute)
          score.permute <- obj$goodness_of_fit_test(y.pred=test.pred,
                                                    y=test.permute$y)
          importance.var <- (score.bench - score.permute)/abs(score.bench)
          importance.permute <- c(importance.permute,importance.var)
        }
        Importance.permute <- rbind(Importance.permute,importance.permute)
        setTxtProgressBar(pb, i)
      }
      close(pb)

      Importance.permute <- abs(apply(Importance.permute, 2, mean))
      Importance.permute <- data.frame(xvars=self$var_names$xvars,score=Importance.permute)
      Importance.permute$xvars <- with(Importance.permute,
                                       factor(xvars,levels=xvars[order(score)]))

      #browser()
      p <- ggplot2::ggplot(data=Importance.permute, aes(x=xvars, y=score)) +
        geom_bar(stat="identity") +
        coord_flip()

      print(Importance.permute)

      return(p)
    },

    goodness_of_fit = function(){
      n <- dim(self$data)[1]
      y <- self$data[, "y"]
      y.pred <- self$data[,"y_pred"]
      score <- -sum(y*log(y.pred)+(1-y)*log(1-y.pred))/n

      print(tree_score)
    }


  ),

  private = list(

    Score_node = function(person_id){
      y <- self$data[self$data$person_id %in% person_id, "y"]
      mu <- mean(y)
      score <- -mean(y*log(mu)+(1-y)*log(1-mu))
      return(score)
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
      return(OptimalDecisionProb)
    },

    Optimal_trt_pred = function(mu){
      trt <- (1*(mu > 0.5))
      return(trt)
    }


  )
)
