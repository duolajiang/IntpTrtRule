#' @title regression tree
#' @param x a dataframe of patient-level covariates
#' @param z a vector of binary treatment, 1 treatment and 0 observation
#' @param xvars a character vector indicating the regressor variables names.
#' @export
#' @import dplyr
#' @import R6
TESubgrp <- R6::R6Class(
  "TEsubgrp",
  #-------------------------public fields-----------------------------#
  public = list(
    data = data.frame(),
    xvars = character(),
    xvars.ctg = character(),

    tree = NA,
    constraints = list(min.ndsz=NA,
                        pt = NA,
                        max.depth = NA),

    initialize = function(x,trt,y.obs,
                          y1.hat.mean, y0.hat.mean,
                          y1.hat.var=0,y0.hat.var=0,
                          xvars){
      self$data <- cbind(x,trt,y.obs,
                         y1.hat.mean,y0.hat.mean,
                         y1.hat.var,y0.hat.var)
      colnames(self$data) <- c(colnames(x),"trt","y.obs",
                               "y1.hat.mean","y0.hat.mean",
                               "y1.hat.var","y0.hat.var")
      self$xvars <- xvars
      self$xvars.ctg <- self$data %>% select(xvars) %>% select_if(Negate(is.numeric)) %>% colnames()
      self$data$person_id <- seq(dim(self$data)[1])
    },


    # nd is short for "node"
    fit = function(min.ndsz, pt, max.depth){
      #browser()
      self$constraints$min.ndsz <- min.ndsz
      self$constraints$pt <- pt
      self$constraints$max.depth <- max.depth

      out <- list.nd <- temp.list <- temp.name <- NULL
      list.nd <- list(self$data$person_id);
      name <- 1

      #i <- 1
      while (length(list.nd)!=0) {
        for (i in 1:length(list.nd)){
          #print(i)
          if (!is.null(list.nd[[i]]) && length(list.nd[[i]]) > 1){
            split <- private$partition(data_id=list.nd[[i]], name[i],
                                       min.ndsz=self$constraints$min.ndsz,
                                       pt=self$constraints$pt,
                                       max.depth=self$constraints$max.depth,
                                       split.var=self$xvars, ctg=self$xvars.ctg,
                                       mtry=length(xvars))

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

    },


    print_leaf_nodes = function(i=1,path=NULL,data=self$tree,dat=self$data,cate=self$xvars.ctg){
      #browser()
      if(!is.element(i,data$node)){
        return(NULL)
      }

      if (i!=1){
        parent.node <- substr(i,1,nchar(i)-1)
        vname <- as.character(data[data$node==parent.node,'vname'])
        if(substr(i,nchar(i),nchar(i))==1){
          if(is.element(vname,cate)){
            cut <- as.character(unlist(strsplit(as.character(data[data$node==parent.node,'cut']),split=" ")))
            cut.index <- is.element(levels(dat[,vname]),cut)
            cut.val <- levels(dat[,vname])[cut.index]
            cut <- NULL
            for(j in 1:length(cut.val)){
              cut <- paste(cut,cut.val[j],sep = ' ')
            }
            vname.con <- paste(vname,'in', cut, sep = ' ')
            if (parent.node==1){
              path <- paste(path,'{',vname.con,sep = ' ')
            } else{
              path <- paste(path,vname.con,sep = ' AND ')
            }
          } else{
            cut <- data[data$node==parent.node,'cut']
            vname.con <- paste(vname,'<=', cut, sep = ' ')
            path <- ifelse(parent.node==1,paste(path,'{',vname.con,sep = ' '),paste(path,vname.con,sep = ' AND '))
          }
        } else{
          if(is.element(vname,cate)){
            cut <- as.character(unlist(strsplit(as.character(data[data$node==parent.node,'cut']),split=" ")))
            cut.index <- !is.element(levels(dat[,vname]),cut)
            cut.val <- levels(dat[,vname])[cut.index]
            cut <- NULL
            for(j in 1:length(cut.val)){
              cut <- paste(cut,cut.val[j],sep = ' ')
            }
            vname.con <- paste(vname,'in', cut, sep = ' ')
            if (parent.node==1){
              path <- paste(path,'{',vname.con,sep = ' ')
            } else{
              path <- paste(path,vname.con,sep = ' AND ')
            }
          } else{
            cut <- data[data$node==parent.node,'cut']
            vname.con <- paste(vname,'>', cut, sep = ' ')
            path <- ifelse(parent.node==1,paste(path,'{',vname.con,sep = ' '),paste(path,vname.con,sep = ' AND '))
          }
        }
      }

      if((!is.element(paste(i,1,sep = ''),data$node))&(!is.element(paste(i,2,sep = ''),data$node))){

        path <- paste(path,'}',
                      'ate=',round(data[data$node==i,'trt.effect'],3),
                      # #'(',round(data[data$node==i,'lower'],3),',',round(data[data$node==i,'upper'],3),')',
                      'size=',data[data$node==i,'size'],
                      't%=', round(data[data$node==i,'n.t']/data[data$node==i,'size'],3)*100,
                      'se=', format(data[data$node==i,'se'], digits=3),
                      sep = ' ')
        print(path)
      }

      self$print_leaf_nodes(paste(i,1,sep = ''),path,data,dat,cate)
      self$print_leaf_nodes(paste(i,2,sep = ''),path,data,dat,cate)

    },


    SankeyNetworkPlot = function(tree=self$tree,dat=self$data,cate=self$xvars.ctg,link_group=FALSE){
      rownames(tree) <- seq(0,dim(tree)[1]-1,1)
      Link <- data.frame(source=character(),target=character(),IDsource=numeric(),IDtarget=numeric(),stringsAsFactors = FALSE)
      i <- 1
      for (j in seq(2,dim(tree)[1],1)) {
        nodeid <- tree[j,'node']
        isleft <- (substr(tree[j,'node'],start = nchar(tree[j,'node']),stop = nchar(tree[j,'node']))==1)
        node.parent.id <- substr(tree[j,'node'],start = 1,stop = nchar(tree[j,'node'])-1)
        node.parent <- tree[tree$node==node.parent.id,]
        var.name <- as.character(node.parent[,'vname'])
        if(is.element(var.name,cate)){
          cut <- as.character(unlist(strsplit(as.character(node.parent[,'cut']),split=" ")))
          cut.index <- is.element(levels(dat[,var.name]),cut)
          cut.val.left <- levels(dat[,var.name])[cut.index];   cut.val.left <- CombineLevelsString(cut.val.left)
          cut.val.right <- levels(dat[,var.name])[!cut.index]; cut.val.right <- CombineLevelsString(cut.val.right)
          target <- ifelse(isleft==TRUE,paste(var.name,'=',cut.val.left,sep = ''),paste(var.name,'=',cut.val.right,sep = ''))
          source <- ifelse(node.parent.id==1,'all',GetNodeName(tree,dat,cate,node.parent.id))
        } else{
          cut.val <- as.character(node.parent[,'cut'])
          target <- ifelse(isleft==TRUE,paste(var.name,'<=',cut.val,sep = ''),paste(var.name,'>',cut.val,sep = ''))
          source <- ifelse(node.parent.id==1,'all',GetNodeName(tree,dat,cate,node.parent.id))
        }
        Link[i,'source'] <- source; Link[i,'target'] <- target; Link[i,'value'] <- tree[j,'size']
        Link[i,'IDsource'] <- rownames(node.parent); Link[i,'IDtarget'] <- rownames(tree[j,])
        if(link_group==TRUE){
          Link[i,'group'] <- ifelse(tree[j,'trt']==1,'1','0')
        }
        i <- i+1
      }


      n_occur <- data.frame(table(Link$target))
      n_occur <- n_occur[n_occur$Freq > 1,]
      if(dim(n_occur)[1]!=0){
        for (i in seq(dim(n_occur)[1])){
          target <- n_occur[i,'Var1']
          target.data <- Link[Link$target==target,]
          k <- NULL
          for (j in seq(dim(target.data)[1])) {
            k <- paste(k,' ',sep = '')
            id.row <- rownames(target.data[j,])
            Link[id.row,'target'] <- paste(Link[id.row,'target'],k,sep = '')
            id.target <- Link[id.row,'IDtarget']
            if(is.element(id.target,Link$IDsource)){
              Link[Link$IDsource==id.target,'source'] <- Link[id.row,'target']
            }
          }
        }
      }


      # From these flows we need to create a node data frame: it lists every entities involved in the flow
      nodes <- data.frame(
        name=c(as.character(Link$source),
               as.character(Link$target)) %>% unique()
      )


      # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
      Link$IDsource <- match(Link$source, nodes$name)-1
      Link$IDtarget <- match(Link$target, nodes$name)-1

      if(link_group==TRUE){
        nodes$group <- c('nodes.group')
        my_color <- 'd3.scaleOrdinal() .domain(["1", "0","nodes.group"]) .range(["red", "silver","silver"])'

        p <- sankeyNetwork(Links = Link, Nodes = nodes,
                           Source = "IDsource", Target = "IDtarget",
                           Value = "value", NodeID = "name",LinkGroup = 'group',
                           NodeGroup = 'group',colourScale = my_color,fontSize = 15,
                           sinksRight=FALSE)
        p
      } else{
        p <- sankeyNetwork(Links = Link, Nodes = nodes,
                           Source = "IDsource", Target = "IDtarget",
                           Value = "value", NodeID = "name",fontSize = 15,
                           sinksRight=FALSE)
        p
      }
    },

    PlotCostSensTree = function(tree.pruned=self$tree,dat=self$data, cate=self$xvars.ctg,cost=0){
      #browser()
      tree.pruned$trt <- 1*(tree.pruned$trt.effect > cost)
      tree.pruned <- private$prune_search(1,tree.pruned)
      path <- NULL
      self$print_leaf_nodes(1,path,tree.pruned,dat,cate)
      self$SankeyNetworkPlot(tree.pruned,dat,cate,link_group = TRUE)
    }
  ),


  private = list(
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
        dat <- self$data[data_id,]

        trt <- dat$trt
        trt.effect <- private$ATE_estimator(data_id)

        # ============================================================
        # COMPUTE THE SCORE IN CURRENT NODE
        score_node <- private$SS_score(dat$y1.hat.mean, dat$y0.hat.mean, dat$y1.hat.var, dat$y0.hat.var)
        # ============================================================
        n.1 <- sum(trt==1); n.0 <- n-n.1 ; n0 <- n*pt
        # CONTROL THE MAX TREE DEPTH
        depth <- nchar(name)
        if (depth <= max.depth && n >= min.ndsz && min(n.1, n.0) >= n0) {
          print(c(n,min.ndsz))
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
                  score_cnode_l <- private$SS_score(dat.l$y1.hat.mean, dat.l$y0.hat.mean, dat.l$y1.hat.var, dat.l$y0.hat.var)
                  score_cnode_r <- private$SS_score(dat.r$y1.hat.mean, dat.r$y0.hat.mean, dat.r$y1.hat.var, dat.r$y0.hat.var)
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
                                 trt.effect=trt.effect,
                                 #lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]],
                                 var = var, vname=vname, cut= cut,
                                 score=max.score, se=sqrt(score_node/n))
          #print(variance)
        }
        else{
          out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=NA, n.r=NA, trt.effect=trt.effect,
                                 #lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]],
                                 var=NA, vname=NA, cut=NA, score=NA, se=sqrt(score_node/n))
          #print(variance)
        }
        return(out)

    },


    # sum of square of the node
    SS_score = function(y1.hat.mean, y0.hat.mean, y1.hat.var, y0.hat.var){
      n <- length(y1.hat.mean)
      y <- y1.hat.mean - y0.hat.mean
      y.mean <- mean(y)
      #compute variance of this group treatment effect
      wgss <- sum(y1.hat.var + y0.hat.var)
      bgss <- sum((y-y.mean)^2)
      node_ss <- wgss + bgss

      return(node_ss)
    },


    # maximum likelihood score
    MLL_score = function(y1.hat.mean, y0.hat.mean, y1.hat.var, y0.hat.var){
      n <- length(y1.hat.mean)
      y <- y1.hat.mean - y0.hat.mean
      #compute variance of this group treatment effect
      wgvar <- mean(y1.hat.var + y0.hat.var)
      bgvar <- var(y)
      node_mll <- -1/2*n*log(wgvar + bgvar)
      return(node_mll)
    },


    ATE_estimator = function(index) {
      #browser()
      y1.hat.mean <- self$data$y1.hat.mean[index]
      y0.hat.mean <- self$data$y0.hat.mean[index]
      est <- mean(y1.hat.mean - y0.hat.mean)
      return(est)
    },


    prune_search = function(i,data){
      #browser()
      if(!is.element(i,data$node)){
        #print(i)
        return(NULL)
      }
      #print(i)

      left.subtree <- prune_search(paste(i,1,sep = ''),data)
      right.subtree <- prune_search(paste(i,2,sep = ''),data)
      tree <- rbind(data[data$node==i,],left.subtree,right.subtree)

      if((!is.null(left.subtree))&(!is.null(right.subtree))){
        n.l <- dim(left.subtree)[1]; n.r <- dim(right.subtree)[1]
        if(n.l==n.r){
          n.t <- sum(left.subtree$trt == right.subtree$trt)
          if (n.t==n.l){
            tree <- data[data$node==i,]
            tree[,c('n.l','n.r','var','vname','cut','score')] <- NA
            return(tree)
          } else{
            return(tree)
          }
        } else{
          return(tree)
        }
      } else{
        return(tree)
      }
    }
  )
)




