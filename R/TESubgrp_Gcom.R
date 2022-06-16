#' @title regression tree
#' @param x a dataframe of patient-level covariates
#' @export
#' @import R6
TESubgrp_Gcom <- R6::R6Class(
  "TEsubgrp_Gcom",
  inherit = TESubgrp,
  #-------------------------public fields-----------------------------#
  public = list(

    initialize = function(x,trt,y1,y0,y,xvars){
      self$data <- cbind(x,trt,y1,y0,y)
      colnames(self$data) <- c(colnames(x),"trt","y1.hat","y0.hat","y")
      super$initialize(xvars)
    }
  ),

  private = list(
    IncrementLL = function(){
      trt <- dat$trt
      n <- nrow(dat)
      n.l <- sum(z==1)
      n.r <- sum(z==0)
      # n.l.t <- sum((z==1)&(trt==1))/sum((z==1))
      # n.l.c <- sum((z==1)&(trt==0))/sum((z==1))
      # n.r.t <- sum((z==0)&(trt==1))/sum((z==0))
      # n.r.c <- sum((z==0)&(trt==0))/sum((z==0))
      #if ((n.l<min.ndsz)|(n.r<min.ndsz)|(min(n.l.t,n.l.c,n.r.t,n.r.c) < pt)) {
      if ((n.l<min.ndsz)|(n.r<min.ndsz)) {
        score <- NA
        return(score)
      }

      y.l <- dat[z==1,'y']
      y.r <- dat[z==0,'y']
      y <- dat$y

      mu.l <- mean(y.l)
      mu.r <- mean(y.r)
      mu   <- mean(y)

      predvar.l <- dat[z==1,'predictive_var']
      predvar.r <- dat[z==0,'predictive_var']
      predvar   <- dat$predictive_var

      sigma.l <- (sum(predvar.l) + sum((y.l-mu.l)^2))/n.l
      sigma.r <- (sum(predvar.r) + sum((y.r-mu.r)^2))/n.r
      sigma   <- (sum(predvar) + sum((y-mu)^2))/n

      ll.l <- -1/2*n.l*log(sigma.l)
      ll.r <- -1/2*n.r*log(sigma.r)
      ll   <- -1/2*n*log(sigma)

      score <- ll.l+ll.r-ll

      return(score)
    },

    ATE_estimator = function(index){
      y1.hat <- self$data$y1.hat[index]
      y0.hat <- self$data$y0.hat[index]
      ate <- mean(y1.hat) - mean(y0.hat)

      return(ate)
    }
  )
)
