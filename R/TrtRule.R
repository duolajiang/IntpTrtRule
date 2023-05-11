#' @export
#' @import dplyr
#'
TrtRule <- R6::R6Class(
  "TrtRule",
  public = list(
    data = data.frame(),

    var_names = list(xvars = character(),
                     xvars.ctg = character(),
                     trt = character(),
                     y.obs = character()),

    y1.hat = 0,
    y0.hat = 0,

    tree = NA,
    tree_pruned = NA,
    tree_topology = NA,

    constraints = list(min.ndsz=NA,
                       pt = NA,
                       max.depth = NA,
                       score.threhold = NA),

    initialize = function(data,var_names,
                          y1.hat,y0.hat){
      self$data <- data
      self$var_names <- var_names
      self$var_names$xvars.ctg <- self$data %>%
        select(var_names$xvars) %>%
        select_if(Negate(is.numeric)) %>%
        colnames()
      self$data$person_id <- seq(dim(self$data)[1])
      self$y1.hat <- y1.hat
      self$y0.hat <- y0.hat
      self$data$y_pred <- 0
      self$data$y_pred_trt <- -1
      self$data$node_id <- NA
    },

    fit = function(){},

    # ==============================================
    # PREDICT TREATMENT EFFECT
    # test: data for prediction
    # ctg: categorical covariate index
    # ==============================================
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
    },

    est_loss_decision = function(){},

    est_outcome_decision = function(){},

    est_treatment_decision = function(){},


    Print_leaf_nodes = function(pruned=FALSE){
      if(isTRUE(pruned)){
        private$print_leaf_nodes(i=1,
                                 path=NULL,
                                 data=self$tree_pruned,
                                 dat=self$data,
                                 cate=self$var_names$xvars.ctg)
      } else {
        private$print_leaf_nodes(i=1,
                                 path=NULL,
                                 data=self$tree,
                                 dat=self$data,
                                 cate=self$var_names$xvars.ctg)
      }
    },

    plot_sankey = function(pruned=TRUE){
      if(isTRUE(pruned)){
        if(dim(self$tree_pruned)[1]==1) {

          stop("no subgroup found given the cost")
          } else{
            private$SankeyNetworkPlot(tree=self$tree_pruned,
                                    dat=self$data,
                                    cate=self$var_names$xvars.ctg,
                                    link_group=TRUE)
            }
      } else{
        private$SankeyNetworkPlot(tree=self$tree,
                                  dat=self$data,
                                  cate=self$var_names$xvars.ctg,
                                  link_group=FALSE)

      }
    },

    feature_importance = function(){}


  ),

  private = list(

    fit_implementation = function(data_id=NA){
      #browser()

      out <- out_topology <- list.nd <- temp.list <- temp.name <- NULL
      name <- 1

      if(is.na(data_id)){
        list.nd <- list(self$data$person_id)
      } else{
        list.nd <- list(data_id)
      }

      #i <- 1
      while (length(list.nd)!=0) {
        for (i in 1:length(list.nd)){
          #print(i)
          if (!is.null(list.nd[[i]]) && length(list.nd[[i]]) > 1){
            split <- private$partition(data_id=list.nd[[i]], name[i],
                                       min.ndsz=self$constraints$min.ndsz,
                                       pt=self$constraints$pt,
                                       max.depth=self$constraints$max.depth,
                                       score.threhold = self$constraints$score.threhold,
                                       split.var=self$var_names$xvars, ctg=self$var_names$xvars.ctg,
                                       mtry=length(self$var_names$xvars))

            out <- rbind(out, split$info)
            out_topology <- rbind(out_topology, split$topology)
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
      out_topology <- out_topology %>% mutate_if(is.factor, as.character)
      return(list(out=out,
                  out_topology=out_topology))
    },

    SankeyNetworkPlot = function(tree,dat,cate,link_group){
      #browser()
      rownames(tree) <- seq(0,dim(tree)[1]-1,1)
      Link <- data.frame(source=character(),target=character(),IDsource=numeric(),IDtarget=numeric(),IDrowsource = character(),IDrowtarget = character(),stringsAsFactors = FALSE)
      node <- data.frame(id= character(),name = character(), group = numeric(), stringsAsFactors = FALSE)
      node[1,"id"] <- "1"
      node[1,"name"] <- "all"
      node[1,"group"] <- tree[1,"mu"]
      for (j in seq(2,dim(tree)[1],1)) {
        nodeid <- tree[j,'node']
        isleft <- (substr(tree[j,'node'],start = nchar(tree[j,'node']),stop = nchar(tree[j,'node']))==1)
        node.parent.id <- substr(tree[j,'node'],start = 1,stop = nchar(tree[j,'node'])-1)
        node.parent <- tree[tree$node==node.parent.id,]
        var.name <- as.character(node.parent[,'vname'])
        if(is.element(var.name,cate)){
          #cut <- as.character(unlist(strsplit(as.character(node.parent[,'cut']),split=" ")))
          #cut.index <- is.element(levels(dat[,var.name]),cut)
          #cut.val.left <- levels(dat[,var.name])[cut.index];   cut.val.left <- private$CombineLevelsString(cut.val.left)
          #cut.val.right <- levels(dat[,var.name])[!cut.index]; cut.val.right <- private$CombineLevelsString(cut.val.right)
          cut.val.left <- self$tree_topology[self$tree_topology$node == node.parent.id, "left.node.var"]
          cut.val.right <- self$tree_topology[self$tree_topology$node == node.parent.id, "right.node.var"]
          target <- ifelse(isleft==TRUE,paste(var.name,'=',cut.val.left,sep = ''),paste(var.name,'=',cut.val.right,sep = ''))
          source <- ifelse(node.parent.id==1,'all',private$GetNodeName(tree,cate,node.parent.id))
        } else{
          cut.val <- as.character(node.parent[,'cut'])
          target <- ifelse(isleft==TRUE,paste(var.name,'<=',cut.val,sep = ''),paste(var.name,'>',cut.val,sep = ''))
          source <- ifelse(node.parent.id==1,'all',private$GetNodeName(tree,cate,node.parent.id))
        }
        Link[j-1,'source'] <- source; Link[j-1,'target'] <- target; Link[j-1,'value'] <- tree[j,'size']
        Link[j-1,'IDsource'] <- node.parent$node; Link[j-1,'IDtarget'] <- tree[j,'node']
        Link[j-1,'group'] <- tree[j,'mu']

        node[j,"id"] <- tree[j,'node']
        node[j,"name"] <- target

        if(link_group==TRUE){
          node_child <- paste(tree[j,"node"],"1",sep = "")
          if(!(node_child %in% tree$node)) {
            node[j,'group'] <- ifelse(tree[j,'trt']==1,'T','C')
          } else {
            node[j,'group'] <- tree[j,'mu']
          }
        } else {
          node[j,"group"] <- tree[j,"mu"]
        }
      }

      node$row_id <- seq(dim(node)[1])

      for (i in seq(dim(Link)[1])) {
        Link_i_IDrowsource <- node %>% filter(id == Link[i,"IDsource"]) %>% select(row_id)
        Link[i,"IDrowsource"] <- Link_i_IDrowsource$row_id -1
        Link_i_IDrowtarget <- node %>% filter(id == Link[i,"IDtarget"]) %>% select(row_id)
        Link[i,"IDrowtarget"] <- Link_i_IDrowtarget$row_id -1
      }

      # n_occur <- data.frame(table(Link$target))
      # n_occur <- n_occur[n_occur$Freq > 1,]
      # if(dim(n_occur)[1]!=0){
      #   for (i in seq(dim(n_occur)[1])){
      #     target <- n_occur[i,'Var1']
      #     target.data <- Link[Link$target==target,]
      #     k <- NULL
      #     for (j in seq(dim(target.data)[1])) {
      #       k <- paste(k,' ',sep = '')
      #       id.row <- rownames(target.data[j,])
      #       Link[id.row,'target'] <- paste(Link[id.row,'target'],k,sep = '')
      #       id.target <- Link[id.row,'IDtarget']
      #       if(is.element(id.target,Link$IDsource)){
      #         Link[Link$IDsource==id.target,'source'] <- Link[id.row,'target']
      #       }
      #     }
      #   }
      # }


      # From these flows we need to create a node data frame: it lists every entities involved in the flow
      # nodes <- data.frame(
      #   name=c(as.character(Link$source),
      #          as.character(Link$target)) %>% unique()
      # )

      # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
      # Link$IDsource <- match(Link$source, nodes$name)-1
      # Link$IDtarget <- match(Link$target, nodes$name)-1

      Link$group <- cut(Link$group, breaks = 10, labels = seq(10))
      for(i in seq(dim(node)[1])){
        if(i==1) {
          node[i,'group'] <- 1
        } else {
          if(link_group==TRUE){
            if(!node[i,'group'] %in% c("C","T")){
              node[i,'group'] <- Link[Link$IDtarget==node[i,'id'],'group']
            }
          } else{
            node[i,'group'] <- Link[Link$IDtarget==node[i,'id'],'group']
          }
        }
      }

      #NOTE, must convert ID to numeric, otherwise, plots error, the same for group.
      cols <- colorRampPalette(c("grey", "red"))(10)
      Link$IDrowsource <- as.numeric(Link$IDrowsource)
      Link$IDrowtarget <- as.numeric(Link$IDrowtarget)
      Link$group <- as.character(Link$group)
      node$group <- as.character(node$group)

      if(link_group==TRUE){
        my_color <- paste("d3.scaleOrdinal() .domain([1,2,3,4,5,6,7,8,9,10,'C','T']).range([",
                          paste('"',cols,'"', collapse = ","),",'white','red'])")
        p <- networkD3::sankeyNetwork(Links = Link, Nodes = node,
                                      Source = "IDrowsource", Target = "IDrowtarget",
                                      Value = "value", NodeID = "name",LinkGroup = 'group',
                                      NodeGroup = 'group',
                                      colourScale = my_color,fontSize = 15,
                                      sinksRight=FALSE)
        p
      } else{
        my_color <- paste("d3.scaleOrdinal() .domain([1,2,3,4,5,6,7,8,9,10]).range([",
                          paste('"',cols,'"', collapse = ","),"])")
        p <- networkD3::sankeyNetwork(Links = Link, Nodes = node,
                                      Source = "IDrowsource", Target = "IDrowtarget",
                                      Value = "value", NodeID = "name",LinkGroup = 'group',
                                      NodeGroup = 'group',
                                      colourScale = my_color,fontSize = 15,
                                      sinksRight=FALSE)
        p
      }

      # if(link_group==TRUE){
      #   p <- networkD3::sankeyNetwork(Links = Link, Nodes = node,
      #                                 Source = "IDrowsource", Target = "IDrowtarget",
      #                                 Value = "value", NodeID = "name",LinkGroup = 'group',
      #                                 NodeGroup = 'group',
      #                                 colourScale = my_color,fontSize = 15,
      #                                 sinksRight=FALSE)
      #   p
      # } else{
      #   p <- networkD3::sankeyNetwork(Links = Link, Nodes = node,
      #                                 Source = "IDrowsource", Target = "IDrowtarget",
      #                                 LinkGroup = 'group',NodeGroup = 'group',
      #                                 Value = "value", NodeID = "name",fontSize = 15,
      #                                 sinksRight=FALSE)
      #   p
      # }
    },

    prune_search = function(i,data){
      #browser()
      if(!is.element(i,data$node)){
        #print(i)
        return(NULL)
      }
      #print(i)

      left.subtree <- private$prune_search(paste(i,1,sep = ''),data)
      right.subtree <- private$prune_search(paste(i,2,sep = ''),data)
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
    },

    partition = function(data_id, name,
                         min.ndsz, pt, max.depth, score.threhold,
                         split.var, ctg,
                         mtry){
      # NOTE THAT CTG INDICATES THE COLUMNS FOR CATEGORICAL VARIABLES.
      #browser()
      out <- match.call(expand = F)
      out$info <- out$name.l <- out$name.r <- out$topology <- out$left <- out$right <- out$... <- NULL
      name.l <- paste(name, 1, sep=""); name.r <- paste(name, 2, sep="")
      n <- length(data_id)
      var <- vname <- NA; cut <- NA; max.score <- score.threhold;
      dat <- self$data[self$data$person_id %in% data_id, ]
      trt <- dat[,self$var_names$trt]
      trt.effect <- mean(dat$y1.hat.mean-dat$y0.hat.mean)
      y <- dat$y
      mu <- mean(y)

      # y.obs.bar <- mean(dat[,self$var_names$y.obs])
      # y1.hat.bar <- mean(dat[,'y1.hat.mean'])
      # y0.hat.bar <- mean(dat[,'y1.hat.mean'])
      # ============================================================
      # COMPUTE THE SCORE IN CURRENT NODE
      score_node <- private$Score_node(data_id)
      # ============================================================
      n.1 <- sum(trt==1); n.0 <- n-n.1 ; n0 <- n*pt
      # CONTROL THE MAX TREE DEPTH
      depth <- nchar(name)
      if (depth <= max.depth && n >= min.ndsz && min(n.1, n.0) >= n0) {
        #print(c(n,min.ndsz))
        m.try <- ifelse(is.null(mtry), length(split.var), mtry)
        for(var_name in split.var) {
          x <- dat[,var_name]; v.name <- var_name; temp <- sort(unique(x));
          if(length(temp) > 1) {
            if (is.element(var_name,ctg)) zcut <- private$power.set(temp) ############################ CLASS VARIABLE
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
          out$right <- dat[!is.element(dat[,var], best.cut),"person_id"]
          cut.var.l <- private$CombineLevelsString(unique(dat[is.element(dat[,var], best.cut),var]))
          cut.var.r <- private$CombineLevelsString(unique(dat[!is.element(dat[,var], best.cut),var]))
          }
        else {
          out$left  <- dat[dat[,var]<= best.cut,"person_id"]
          out$right <- dat[dat[,var]> best.cut, "person_id"]
          cut.var.l <- paste(var, "<=", best.cut, sep="")
          cut.var.r <- paste(var, ">", best.cut, sep="")
        }
        out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=length(out$left),n.r=length(out$right),
                               trt.effect=trt.effect, mu=mu,
                               lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]],
                               var = var, vname=vname, cut= cut,
                               score=max.score, se=sqrt(score_node/n))

        out$topology <- data.frame(node=name, left.node=name.l, right.node=name.r,
                                   left.node.var=cut.var.l, right.node.var=cut.var.r)
        #print(variance)
      }
      else{
        out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=NA, n.r=NA,
                               trt.effect=trt.effect, mu=mu,
                               lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]],
                               var=NA, vname=NA, cut=NA,
                               score=NA, se=sqrt(score_node/n))

        out$topology <- data.frame(node=name, left.node=NA, right.node=NA,
                                   left.node.var=NA, right.node.var=NA)
        #print(variance)
      }

      self$data[self$data$person_id %in% data_id, "y_pred"] <- mu
      self$data[self$data$person_id %in% data_id, "y_pred_trt"] <- private$Optimal_trt_pred(mu)
      self$data[self$data$person_id %in% data_id, "node_id"] <- as.character(name)

      return(out)

    },


    print_leaf_nodes = function(i,path,data,dat,cate){
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
                      'mu=',round(data[data$node==i,'mu'],3),
                      # #'(',round(data[data$node==i,'lower'],3),',',round(data[data$node==i,'upper'],3),')',
                      'size=',data[data$node==i,'size'],
                      't%=', round(data[data$node==i,'n.t']/data[data$node==i,'size'],3)*100,
                      'se=', format(data[data$node==i,'se'], digits=3),
                      sep = ' ')
        print(path)
      }

      private$print_leaf_nodes(paste(i,1,sep = ''),path,data,dat,cate)
      private$print_leaf_nodes(paste(i,2,sep = ''),path,data,dat,cate)

    },



    Score_node = function(person_id){},

    Optimal_trt_pred = function(mu){},


    # ===========================================================================
    # THE power.set() FUNCTION PROVIDES THE POWER SET FOR A CATEGORICAL VARIABLE
    # ===========================================================================
    power.set = function(x) {
      if(length(x) == 0) return(vector(mode(x), 0))
      x <- sort(unique(x)); n <- length(x); K <- NULL
      for(m in x) K <- rbind(cbind(K, FALSE), cbind(K, TRUE))
      out <- apply(K, 1, function(x, s) s[x], s = x)
      out <- out[-c(1, length(out))]
      l <- length(out); i <- 1
      out[!sapply(out, length)>=ceiling(n/2+.5)]
    },


    # =============================================================
    # for printing
    # =============================================================
    GetNodeName = function(tree,cate,nodeid){
      #browser()
      isleft <- (substr(nodeid,start = nchar(nodeid),stop = nchar(nodeid))==1)
      node.parent.id <- substr(nodeid,start = 1,stop = nchar(nodeid)-1)
      node.parent <- tree[tree$node==node.parent.id,]
      var.name <- as.character(node.parent[,'vname'])
      if(is.element(var.name,cate)){
        #cut <- as.character(unlist(strsplit(as.character(node.parent[,'cut']),split=" ")))
        #cut.index <- is.element(levels(dat[,var.name]),cut)
        #cut.val.left <- levels(dat[,var.name])[cut.index];   cut.val.left <- private$CombineLevelsString(cut.val.left)
        #cut.val.right <- levels(dat[,var.name])[!cut.index]; cut.val.right <- private$CombineLevelsString(cut.val.right)
        cut.val.left <- self$tree_topology[self$tree_topology$node == node.parent.id, "left.node.var"]
        cut.val.right <- self$tree_topology[self$tree_topology$node == node.parent.id, "right.node.var"]
        target <- ifelse(isleft==TRUE,paste(var.name,'=',cut.val.left,sep = ''),paste(var.name,'=',cut.val.right,sep = ''))
      } else{
        cut.val <- as.character(node.parent[,'cut'])
        target <- ifelse(isleft==TRUE,paste(var.name,'<=',cut.val,sep = ''),paste(var.name,'>',cut.val,sep = ''))
      }
      return(target)
    },


    CombineLevelsString = function(stringarray){
      single <- NULL
      for(i in seq(length(stringarray))){
        single <- paste(single,as.character(stringarray[i]),'')
      }
      return(single)
    }

  )


)
