#' PlotEvent
#'
#' PLot data
#'
#' @param compare
#' @param dataset The dataset that you want to depict. It can be a single character or a numeric value refers to the element of the compare object.
#' @param stat A character or a numeric value that refers to the statistics (elemnt of the dataset) you want ot depict.
#' @param data A character or numeric vector defining the clusters or group of clusters you want to depict.
#' @param cladogram
#' @param sizemeth
#' @param comp.size
#' @param size
#' @param node_pch_opt
#' @param tip_pch_opt
#' @param node_text_opt
#' @param tip_text_opt
#' @param ...
#'
#' @return
#' @export
#' @import ape readr plyr stringr tidyr
#'
#' @examples
#'
#'
#'
###########
### meg át kell írni mindent hogy ne peldaul  if(tip_pch_opt != FALSE) legyen hanem if(is.list(tip_pch_opt))
###########
PlotEvent <- function(compare, dataset = NULL, stat = NULL, data = NULL, cladogram = FALSE,
                      sizemeth = c("standard","comparable", "log"), comp.size = c(0,1000), base,
                      size = 1, node_pch_opt, tip_pch_opt, node_text_opt, tip_text_opt, ...){
  if(is.character(dataset)){
    dataset <- which(names(compare) %in% dataset)
  }else{
    if(!is.numeric(dataset)){
      stop(paste0("You need to give a character or a numeric value defining the dataset you want to depict!\n"))
    }
  }
  if(is.character(stat)){
    stat <- which(names(compare[[dataset]]) %in% stat)
  }else{
    if(!is.numeric(stat)){
      stop(paste0("You need to give a character or a numeric value defining the statistics you want to depict!\n"))
    }
  }
  switch(sizemeth,
         standard = {compare[[dataset]]$transformed <- (compare[[dataset]][[stat]] - min(compare[[dataset]][[stat]])) / (max(compare[[dataset]][[stat]]) - min(compare[[dataset]][[stat]]))}, # ez jo igy?
         comparable = {compare[[dataset]]$transformed <- (abs(compare[[dataset]][[stat]]) - min(comp.size)) / (max(max(comp.size) - min(comp.size)))},
         log = {compare[[dataset]]$transformed <- log(na.omit(compare[[dataset]][[stat]]), base = base)})
  if(nrow(compare[[dataset]][[stat]]) != compare$tree$Nnode + length(compare$tree$tip.label)){
    stop("Compare object was corrupted. The number of nodes in tree and in the event data is different!")
  }
  if(length(data) == 1 & !is.null(data)){
    if(is.character(data)){
      data <- which(colnames(compare[[dataset]][[stat]]) %in% data)
    }else{
      if(!is.numeric(data)){
        stop(paste0("You need to give a character or a numeric vector defining the clusters or group of clusters you want to depict!\n"))
      }
    }
    tree_nodes <- length(compare$tree$tip.label) + 1 : compare$tree$Nnode
    tree_tips <- 1:length(compare$tree$tip.label)
    if(cladogram){
      compare$tree$edge.length <- NULL
      params <- list(...)
      if(length(params) == 0){
      plot.phylo(compare$tree, label.offset = 0.55 + abs(-0.7) + 0.2,
                   node.depth = 2, show.tip.label = TRUE,
                   cex = 1, edge.width = 2)
      }else{
        plot.phylo(compare$tree, ...)
      }
      if(missing(node_pch_opt)){
        nodelabels(pch = 16, col = "red", node = tree_nodes,
                   cex = compare[[dataset]][["transformed"]][tree_nodes, data] * size)
      }else{
        if(is.list(node_pch_opt)){
          do.call(nodelabels, c(list(node = tree_nodes, cex = compare[[dataset]][["transformed"]][tree_nodes,data] * size),
                              node_pch_opt))
        }
      }
      if(missing(tip_pch_opt)){
        tiplabels(pch = 16, col = "red", tip = tree_tips,
                  cex = compare[[dataset]][["transformed"]][tree_tips, data] * size, adj = c(0.55, 0.5))
      }else{
        if(is.list(tip_pch_opt)){
          do.call(tiplabels, c(list(tip = tree_tips, cex = compare[[dataset]][["transformed"]][tree_tips, data] * size),
                             tip_pch_opt))
        }
      }
      if(missing(node_text_opt)){
        nodelabels(text = compare[[dataset]][[stat]][tree_nodes, data], node = tree_nodes,
                   adj = c(0.5, 0.5), bg = NULL, frame = "none", cex = 1)
      }else{
        if(is.list(node_text_opt)){
          do.call(nodelabels, c(list(text = compare[[dataset]][[stat]][tree_nodes, data], node = tree_nodes),
                              node_text_opt))
        }
      }
      if(missing(tip_text_opt)){
        tiplabels(text = compare[[dataset]][[stat]][tree_tips, data], tip = tree_tips, cex = 1,
                  bg = NULL, frame = "none", adj = -0.7)
      }else{
        if(is.list(tip_text_opt)){
          do.call(tiplabels, c(list(text = compare[[dataset]][[stat]][tree_tips, data], tip = tree_tips),
                             tip_text_opt))
        }
      }
    }else{
      params <- list(...)
      if(length(params) == 0){
        plot.phylo(compare$tree, label.offset = 0.55 + abs(-0.7) + 0.2,
                   show.tip.label = TRUE, cex = 1, edge.width = 2)
      }else{
        plot.phylo(compare$tree, ...)
      }
      if(missing(node_pch_opt)){
        nodelabels(pch = 16, col = "red", node = tree_nodes,
                   cex = compare[[dataset]][["transformed"]][tree_nodes, data] * size)
      }else{
        if(is.list(node_pch_opt)){
          do.call(nodelabels, c(list(node = tree_nodes, cex = compare[[dataset]][["transformed"]][tree_nodes,data] * size),
                              node_pch_opt))
        }
      }
      if(missing(tip_pch_opt)){
        tiplabels(pch = 16, col = "red", tip = tree_tips,
                  cex = compare[[dataset]][["transformed"]][tree_tips, data] * size, adj = c(0.55, 0.5))
      }else{
        if(is.list(tip_pch_opt)){
          do.call(tiplabels, c(list(node = tree_tips, cex = compare[[dataset]][["transformed"]][tree_tips, data] * size),
                             tip_pch_opt))
        }
      }
      if(missing(node_text_opt)){
        nodelabels(text = compare[[dataset]][[stat]][tree_nodes, data], node = tree_nodes,
                   adj = c(0.5, 0.5), bg = NULL, frame = "none", cex = 1)
      }else{
        if(is.list(node_text_opt)){
          do.call(nodelabels, c(list(text = compare[[dataset]][[stat]][tree_nodes, data], node = tree_nodes),
                              node_text_opt))
        }
      }
      if(missing(tip_text_opt)){
        tiplabels(text = compare[[dataset]][[stat]][tree_tips, data], tip = tree_tips, cex = 1,
                  bg = NULL, frame = "none", adj = -0.7)
      }else{
        if(is.list(tip_text_opt)){
          do.call(tiplabels, c(list(text = compare[[dataset]][[stat]][tree_tips, data], tip = tree_tips),
                             tip_text_opt))
        }
      }
    }
  }else{
    if(is.character(data)){
      data <- which(colnames(compare[[dataset]][[stat]]) %in% data)
    }else{
      if(!is.numeric(data)){
        stop(paste0("You need to give a character or a numeric vector defining the clusters or group of clusters you want to depict!\n"))
      }
    }

    tree_nodes <- length(compare$tree$tip.label) + 1 : compare$tree$Nnode
    tree_tips <- 1:length(compare$tree$tip.label)
    for(i in 1:length(data)){
      if(cladogram){
        compare$tree$edge.length <- NULL
        params <- list(...)
        if(length(params) == 0){
          plot.phylo(compare$tree, label.offset = 0.55 + abs(-0.7) + 0.2,
                     node.depth = 2, show.tip.label = TRUE,
                     cex = 1, edge.width = 2)
        }else{
          plot.phylo(compare$tree, main = colnames(compare[[dataset]][["transformed"]])[i],  ...)
        }
        if(missing(node_pch_opt)){
          nodelabels(pch = 16, col = "red", node = tree_nodes,
                     cex = compare[[dataset]][["transformed"]][tree_nodes, data[i]] * size)
        }else{
          if(is.list(node_pch_opt)){
            do.call(nodelabels, c(list(node = tree_nodes, cex = compare[[dataset]][["transformed"]][tree_nodes,data[i]] * size),
                                node_pch_opt))
          }
        }
        if(missing(tip_pch_opt)){
          tiplabels(pch = 16, col = "red", tip = tree_tips,
                    cex = compare[[dataset]][["transformed"]][tree_tips, data[i]] * size, adj = c(0.55, 0.5))
        }else{
          if(is.list(tip_pch_opt)){
            do.call(tiplabels, c(list(tip = tree_tips, cex = compare[[dataset]][["transformed"]][tree_tips, data[i]] * size),
                               tip_pch_opt))
          }
        }
        if(missing(node_text_opt)){
          nodelabels(text = compare[[dataset]][[stat]][tree_nodes, data[i]], node = tree_nodes,
                     adj = c(0.5, 0.5), bg = NULL, frame = "none", cex = 1)
        }else{
          if(is.list(node_text_opt)){
            do.call(nodelabels, c(list(text = compare[[dataset]][[stat]][tree_nodes, data[i]], node = tree_nodes),
                                node_text_opt))
          }
        }
        if(missing(tip_text_opt)){
          tiplabels(text = compare[[dataset]][[stat]][tree_tips, data[i]], tip = tree_tips, cex = 1,
                    bg = NULL, frame = "none", adj = -0.7)
        }else{
          if(is.list(tip_text_opt)){
            do.call(tiplabels, c(list(text = compare[[dataset]][[stat]][tree_tips, data[i]], tip = tree_tips),
                               tip_text_opt))
          }
        }
      }else{
        params <- list(...)
        if(length(params) == 0){
          plot.phylo(compare$tree, label.offset = 0.55 + abs(-0.7) + 0.2,
                     show.tip.label = TRUE, cex = 1, edge.width = 2, main = colnames(compare[[dataset]][["transformed"]])[i])
        }else{
          plot.phylo(compare$tree, ...)
        }
        if(missing(node_pch_opt)){
          nodelabels(pch = 16, col = "red", node = tree_nodes,
                     cex = compare[[dataset]][["transformed"]][tree_nodes, data[i]] * size)
        }else{
          if(is.list(node_pch_opt)){
            do.call(nodelabels, c(list(node = tree_nodes, cex = compare[[dataset]][["transformed"]][tree_nodes,data[i]] * size),
                                node_pch_opt))
          }
        }
        if(missing(tip_pch_opt)){
          tiplabels(pch = 16, col = "red", tip = tree_tips,
                    cex = compare[[dataset]][["transformed"]][tree_tips, data[i]] * size, adj = c(0.55, 0.5))
        }else{
          if(is.list(tip_pch_opt)){
            do.call(tiplabels, c(list(tip = tree_tips, cex = compare[[dataset]][["transformed"]][tree_tips, data[i]] * size),
                               tip_pch_opt))
          }
        }
        if(missing(node_text_opt)){
          nodelabels(text = compare[[dataset]][[stat]][tree_nodes, data[i]], node = tree_nodes,
                     adj = c(0.5, 0.5), bg = NULL, frame = "none", cex = 1)
        }else{
          if(is.list(node_text_opt)){
            do.call(nodelabels, c(list(text = compare[[dataset]][[stat]][tree_nodes, data[i]], node = tree_nodes),
                                node_text_opt))
          }
        }
        if(missing(tip_text_opt)){
          tiplabels(text = compare[[dataset]][[stat]][tree_tips, data[i]], tip = tree_tips, cex = 1,
                    bg = NULL, frame = "none", adj = -0.7)
        }else{
          if(is.list(tip_text_opt)){
            do.call(tiplabels, c(list(text = compare[[dataset]][[stat]][tree_tips, data[i]], tip = tree_tips),
                               tip_text_opt))
          }
        }
      }
    }
  }
}
