#' Perform hierachical clustering of clusters base on their events...
#'
#' @param compare a compare object
#' @param dataset The dataset that you want to depict. It can be a single character or a numeric value refers to the element of the compare object.
#' @param stat A character or a numeric value that refers to the statistics (elemnt of the dataset) you want ot depict.
#' @param cl_method A character string defining the agglomeration method to be used. See hclust {stats} function for details.
#' @param ... arguments passed to the hclust.vector function
#'
#' @return this function returns a compare object retain only just part of the original clusters
#' @export
#' @import fastcluster
#'
#' @examples
############
## be kene epiteni meg a pvclust fuggvenyeket, hogy bootstrap ertekeket lehessen szamolni,
# de elsokent ha korrelaciot lehetne csinalni az is jo lenne...
############
# source("~/R/R/DolloImport3.R")
# source("~/R/R/retrieve_event3.R")
# library(ape)
# library(readr)
# library(stringr)
# library(plyr)
# library(fastcluster)
# tree_path <- system.file("extdata", "71g_species.tree", package = "CompareTools")
# tree <- read.tree(tree_path)
# file_path <- system.file("extdata", "Dolloout", package = "CompareTools")
# files <- file_path
# probe1 <- DolloImport3(tree, files = files)
# files <- "~/R/Dolloout_all_cluster"
# probe2 <- DolloImport3(tree, files = files)
# is.list(probe1$event_data$gains)
#
# str(probe1$event_data$gains)
# is.list(probe1)
# load("~/R/Big_dendro.RData")
HclustComp <- function(compare, dataset = NULL, stat = NULL,
                    cl_method = "ward", ...){
  if(is.matrix(compare)){
    return(hclust.vector(t(compare), method = cl_method, ...))
  }else{
    if(!is.list(compare)){
      stop(paste0("Give a compare object(list) or a matrix!"))
    }
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
    return(hclust.vector(t(compare[[dataset]][[stat]]), method = cl_method, ...))
  }
}
# dendro <- HclustComp(probe1$event_data$gains)
# dendro <- CorComp(probe2$event_data$gains)
# save(list = ls(), file = "~/R/Big_dendro.RData")
# dendro$height <- dendro$height^2
# plot(dendro, labels = FALSE)
#         if(plot){
#         if(is.null(group_n)){
#             plot(hclust.vector(t(compare[[dataset]][[stat]]), ...),
#                  main = paste0("Cluster Dendrogram of the ", names(compare)[dataset], " dataset and the ", names(compare[[dataset]])[stat], " statistics"),
#                  xlab = "Cluster names")
#           }
#         }else{
#           h <- hclust(as.dist(1-compare[[dataset]][[stat]]), method = cl_method, ...)
#           groups <- cutree(h, k = group_n)
#           plot(h,
#                main = paste0("Cluster Dendrogram of the ", names(compare)[dataset], " dataset and the ", names(compare[[dataset]])[stat], " statistics"),
#                xlab = "Cluster names")
#           rect.hclust(h, k = group_n, border = "red")
#         }
#         }else{
#           if(!is.null(group_n)){
#             h <- hclust(as.dist(1-compare[[dataset]][[stat]]), method = cl_method, ...)
#             groups <- cutree(h, k = group_n)
#             df <- data.frame(clusters = names(groups), groups = groups)
#             return(df)
#           }else{
#             stop(paste0("You need to define the group_n argumentum!"))
#           }
#         }
#
#   }
# }
#
