#' @export
Model_Tree <- R6::R6Class(
  "Model_Tree",
  public = list(
    data = data.frame(),

    tree = NA,

    var_names = list(xvars = character(),
                     xvars.ctg = character(),
                     trt = character(),
                     y.obs = character()),

    constraints = list(min.ndsz=NA,
                       pt = NA,
                       max.depth = NA),

    initialize = function(data,var_names){
      self$data <- data
      self$var_names <- var_names
    },

    fit = function(min.ndsz, pt, max.depth){
      self$constraints$min.ndsz <- min.ndsz
      self$constraints$pt <- pt
      self$constraints$max.depth <- max.depth
      self$tree <- private$fit_implementation()
    },

    predict = function(test){
      #browser()
      y <- NULL
      for (i in seq(dim(test)[1])){
        obs <- test[i,]
        node.id <- 1
        node <- self$tree[self$tree$node==node.id,]
        while (!is.na(node$vname)){
          varname <- as.character(node$vname)
          if(is.element(node$var,self$var_names$xvars.ctg)){
            cut.value <- as.character(node$cut)
            if(is.element(obs[,varname],strsplit(cut.value,' '))){
              node.id <- paste(node.id,'1',sep = '')
            } else{
              node.id <- paste(node.id,'2',sep = '')
            }
          } else{
            # cant use as.numeric(node$cut); below is efficient thant as.numeric(as.character(node$cut))
            cut.value <- as.numeric(levels(node$cut)[node$cut])
            if(obs[,varname]<=cut.value){
              node.id <- paste(node.id,'1',sep = '')
            }else{
              node.id <- paste(node.id,'2',sep = '')
            }
          }
          node <- self$tree[self$tree$node==node.id,]
        }
        y <- c(y,node$mu)
      }
      return(y)
    }


  ),
  private = list(

    fit_implementation = function(){
      #browser()
      out <- list.nd <- temp.list <- temp.name <- NULL
      name <- 1
      list.nd <- list(self$data$person_id)

      #i <- 1
      while (length(list.nd)!=0) {
        for (i in 1:length(list.nd)){
          #print(i)
          if (!is.null(list.nd[[i]]) && length(list.nd[[i]]) > 1){
            split <- private$partition(data_id=list.nd[[i]], name[i],
                                       min.ndsz=self$constraints$min.ndsz,
                                       pt=self$constraints$pt,
                                       max.depth=self$constraints$max.depth,
                                       split.var=self$var_names$xvars, ctg=self$var_names$xvars.ctg,
                                       mtry=length(self$var_names$xvars))

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
      return(out)
    },


    partition = function(data_id, name,
                         min.ndsz, pt, max.depth,
                         split.var, ctg,
                         mtry){
      # NOTE THAT CTG INDICATES THE COLUMNS FOR CATEGORICAL VARIABLES.
      #browser()
      out <- match.call(expand = F)
      out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
      name.l <- paste(name, 1, sep=""); name.r <- paste(name, 2, sep="")
      n <- length(data_id)
      var <- vname <- NA; cut <- NA; max.score <- 0;
      dat <- self$data[self$data$person_id %in% data_id, ]
      trt <- dat[,self$var_names$trt]
      trt.effect <- mean(dat$y1.hat.mean-dat$y0.hat.mean)
      y <- dat$y
      mu <- mean(y)

      # ============================================================
      # COMPUTE THE SCORE IN CURRENT NODE
      score_node <- private$Score_node(data_id)
      # ============================================================
      n.1 <- sum(trt==1); n.0 <- n-n.1 ; n0 <- n*pt
      # CONTROL THE MAX TREE DEPTH
      depth <- nchar(name)
      if (depth <= max.depth && n >= min.ndsz && min(n.1, n.0) >= n0) {
        m.try <- ifelse(is.null(mtry), length(split.var), mtry)
        for(var_name in split.var) {
          x <- dat[,var_name]; v.name <- var_name; temp <- sort(unique(x));
          if(length(temp) > 1) {
            if (is.element(var_name,ctg)) zcut <- power.set(temp) ############################ CLASS VARIABLE
            else zcut <- temp[-length(temp)]
            # print(i); print(temp); print(zcut)
            for(j in zcut) {
              score <- NA
              #### CLASS VARIABLE
              if (is.element(var_name,ctg)) {
                grp <- sign(is.element(x, j))
                cut1 <- paste(j, collapse=" ")
              } else {
                grp <- sign(x <= j)
                cut1 <- as.character(j)
              }
              #browser()
              if((sum(grp)>self$constraints$min.ndsz) & (n-sum(grp)>self$constraints$min.ndsz)){
                dat.l <- dat[grp==1,]
                dat.r <- dat[grp==0,]
                n.l <- sum(grp)
                n.r <- n-n.l
                score_cnode_l <- private$Score_node(dat.l$person_id)
                score_cnode_r <- private$Score_node(dat.r$person_id)
                score <- score_node - (score_cnode_l+score_cnode_r)
              }
              #print(cbind(var=i, cut=j, score=score,var=variance))
              if (!is.na(score) && score > max.score) {
                max.score <- score; var <- var_name; vname <- v.name; cut <- cut1; best.cut<-j
              }
            }
          }
        }
      }
      if (!is.na(var)) {
        out$name.l <- name.l; out$name.r <- name.r
        if (is.element(var,ctg)) { ############################
          out$left  <- dat[is.element(dat[,var], best.cut),"person_id"]
          out$right <- dat[!is.element(dat[,var], best.cut),"person_id"]}
        else {
          out$left  <- dat[dat[,var]<= best.cut,"person_id"]
          out$right <- dat[dat[,var]> best.cut, "person_id"]
        }
        out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=length(out$left),n.r=length(out$right),
                               trt.effect=trt.effect, mu=mu,
                               lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]],
                               var = var, vname=vname, cut= cut,
                               score=max.score, se=sqrt(score_node/n))
        #print(variance)
      }
      else{
        out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=NA, n.r=NA,
                               trt.effect=trt.effect, mu=mu,
                               lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]],
                               var=NA, vname=NA, cut=NA,
                               score=NA, se=sqrt(score_node/n))
        #print(variance)
      }

      self$data[self$data$person_id %in% data_id, "y_pred"] <- mu

      return(out)

    }

  )
)
