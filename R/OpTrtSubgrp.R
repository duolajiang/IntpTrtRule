OpTrtSubgrp <- R6::R6Class(
  "OpTrtSubgrp",
  public = list(
    data = data.frame(),
    xvars = character(),
    xvars.ctg = character(),

    y1.hat = data.frame(),
    y0.hat = data.frame(),

    initialize = function(x, trt, y.obs, y1.hat, y0.hat, xvars){
      self$data <- cbind(x, trt, y.obs)
      colnames(self$data) <- c(colnames(x), "trt", "y.obs")
      self$data$y1.hat.mean <- apply(y1.hat,2,mean)
      self$data$y0.hat.mean <- apply(y0.hat,2,mean)
      self$data$person_id <- seq(dim(self$data)[1])

      self$xvars <- xvars
      self$xvars_ctg <- self$data %>% select(xvars) %>% select_if(Negate(is.numeric)) %>% colnames()
    },

    est_optimalTrt = function(rho, cost){
      p11 <- rho*sqrt(y1.hat*(1-y1.hat)*y0.hat*(1-y0.hat)) + y1.hat*y0.hat
      p10 <- y1.hat - p11
      p01 <- y0.hat - p11
      p00 <- 1- p11 - p10 - p01
      costControl <- cost[1,]
      costTreatment <- cost[2,]
      deltaCost <- costTreatment-costControl
      DeltaPosteriorCost <- deltaCost[1]*p00 + deltaCost[2]*p01 + deltaCost[3]*p10 + deltaCost[4]*p11
      OptimalDecisionProb <- apply(DeltaPosteriorCost<0, 2, mean)
      self$data$ProTrtOp <- OptimalDecisionProb
    },

    fit = function(min.ndsz, pt, max.depth){

      split.var <- 1:length(xvars);ctg=6:length(xvars); mtry=length(split.var)
      # FLEXIBILTIY OF TREE STRUCTURE
      min.ndsz=min.ndsz # mininum number of observations in a node
      pt=pt # minimun number of treatment/control in a node
      max.depth=max.depth # max depth of tree

      tree <- grow.INT.ILL.OptimalDecision(dat,split.var,ctg,min.ndsz, pt,
                                           max.depth,mtry=length(split.var),alpha,loc)


      out <- list.nd <- temp.list <- temp.name <- NULL
      list.nd <- list(dat);
      name <- 1
      #i <- 1
      while (length(list.nd)!=0) {
        for (i in 1:length(list.nd)){
          print(i)
          if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){
            split <- private$partition(list.nd[[i]], name[i], min.ndsz=min.ndsz,
                                       pt=pt, split.var=split.var, ctg=ctg,
                                       max.depth=max.depth,mtry=mtry,alpha,loc)
            out <- rbind(out, split$info)
            if (!is.null(split$left)) {
              temp.list <- c(temp.list, list(split$left, split$right))
              temp.name <- c(temp.name, split$name.l, split$name.r)
            }
          }
        }
        list.nd <- temp.list; name <- temp.name
        temp.list <-  temp.name <- NULL
      }
      out$node <- as.character(out$node)
      out <- out[order(out$node),]

      self$tree <- out

    }



  ),
  private = list()
)
