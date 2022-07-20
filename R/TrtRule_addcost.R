#' @export
#' @import ggplot2
TrtRule_addcost <- R6::R6Class(
  "TrtRule_addcost",
  inherit = TrtRule,

  public = list(

    cost = NA,
    # when is.null(dim(y1.hat)), then user must provide y1.var and y0.var
    initialize = function(data,var_names,y1.hat,y0.hat,y1.hat.var=NULL,y0.hat.var=NULL){
      super$initialize(data,var_names,y1.hat,y0.hat)
      if(is.null(dim(y1.hat))){
        if(is.null(y1.hat.var)) {
          stop("since you only provide the mean of y1.hat and y0.hat, you should provide variance of y1.hat and y0.hat!")
        }

        self$data$y1.hat.mean <- y1.hat
        self$data$y0.hat.mean <- y0.hat
        self$data$y1.hat.var <- y1.hat.var
        self$data$y0.hat.var <- y0.hat.var
      } else {
        self$data$y1.hat.mean <- apply(y1.hat,2,mean)
        self$data$y0.hat.mean <- apply(y0.hat,2,mean)
        self$data$y1.hat.var <- apply(y1.hat,2,var)
        self$data$y0.hat.var <- apply(y0.hat,2,var)
      }
      self$data$y <- self$data$y1.hat.mean - self$data$y0.hat.mean
      self$data$y.var <- self$data$y1.hat.var + self$data$y0.hat.var
    },

    fit = function(min.ndsz, pt, max.depth, cost){
      #browser()
      if ((is.na(self$cost))|!((self$constraints$min.ndsz==min.ndsz)&
                               (self$constraints$pt==pt)&
                               (self$constraints$max.depth==max.depth))){
        self$constraints$min.ndsz <- min.ndsz
        self$constraints$pt <- pt
        self$constraints$max.depth <- max.depth
        self$cost <- cost

        self$tree <- private$fit_implementation()
        tree_pruned <- self$tree
        tree_pruned$trt <- 1*(tree_pruned$trt.effect > self$cost)
        self$tree_pruned <- private$prune_search(1,tree_pruned)
      } else {
        self$cost <- cost
        tree_pruned <- self$tree
        tree_pruned$trt <- 1*(tree_pruned$trt.effect > self$cost)
        self$tree_pruned <- private$prune_search(1,tree_pruned)
      }
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
        obj <- Model_RTree$new(data = train,
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
                                                    y=test.permute$y,
                                                    y.var=test.permute$y.var)
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
      y <- self$data[,"y"]
      y.pred <- self$data[,"y_pred"]
      wgss <- sum(self$data$y.var)
      bgss <- sum((y-y.pred)^2)
      tree_score <- (wgss + bgss)/n

      print(tree_score)
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

    },

    Optimal_trt_pred = function(mu){
      trt <- 1*(mu > self$cost)
      return(trt)
    }


  )
)
