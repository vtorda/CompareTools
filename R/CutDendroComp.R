#' Cut a dendrogram into a predefined number of groups
#'
#' @param dendro an object of class hclust
#' @param n_groups The dataset that you want to depict. It can be a single character or a numeric value refers to the element of the compare object.
#' @param plot Logical. If FALSE you have to define the group_n argumentum and a
#' dataframe is returned with cluster names and groups a according hierachical clustering.
#' This dataframe can be inputed to the GoupComp function. Otherwise you can plot a dendrogram with our withour a heatmap
#' @param min_num Integer. Minimum size of a group.
#' @param rect_opt argument passed to the rect.hclust function as a list
#'@param ... arguments passed to the plot function
#'
#' @return this function returns a compare object retain only just part of the original clusters
#' @export
#' @import fastcluster
#'
#' @examples
############
## be kene epiteni meg a pvclust fuggvenyeket, hogy bootstrap ertekeket lehessen szamolni
############
# source("~/R/R/DolloImport3.R")
# source("~/R/R/retrieve_event3.R")
# library(ape)
# library(readr)
# library(stringr)
# library(plyr)
# library(fastcluster)
# majd meg meg kene vizsgalni hogy helyes dendro objektumot adunk be
# plot ba beepiteni a ... argumentumot mert alapbol nem megy ezt kapom: Error: '...' used in an incorrect context
# be kene epiteni hogy addig vegjon es vonogasson ossze mig eleg nagy csoportok nem lesznek
# n_groups <- 25
# min_num <- 30
CutDendroComp <- function(dendro, n_groups = NULL, height = NULL, plot = FALSE, min_num = NULL,
                          rect_opt, ...){
  if(is.null(n_groups) && is.null(height)){
    stop("Either n_groups or height should be specified")
  }
  if(plot){
    groups <- cutree(dendro, k = n_groups, h = height)
    if(is.null(min_num)){
    n_groups <- length(unique(groups))
    plot(dendro, ...)
    if(missing(rect_opt)){
      rect.hclust(dendro, k = n_groups, cluster = groups, border = "red")
    }else{
      do.call(rect.hclust, c(list(tree = dendro, k = n_groups,
                                  cluster = groups, rect_opt)))
    }
    return(groups)
    }else{
      if(!is.numeric(min_num)){
        stop("Give a numeric number to the min_num argumentum...")
      }
      group_tb <- table(groups)
      merge <- names(group_tb)[group_tb < min_num]
      groups[groups %in% merge] <- "merged"
      n_groups <- length(unique(groups))
      plot(dendro, ...)
      if(missing(rect_opt)){
        rect.hclust(dendro, k = n_groups-1, cluster = groups[!groups %in% "merged"], border = "red")
      }else{
        do.call(rect.hclust, c(list(tree = dendro, k = n_groups-1,
                                    cluster = groups[!groups %in% "merged"], rect_opt)))
      }
      return(data.frame(clusters = names(groups), groups = groups, stringsAsFactors = FALSE))
    }
  }else{
    groups <- cutree(dendro, k = n_groups, h = height)
    if(is.null(min_num)){
    return(groups)
    }else{
    group_tb <- table(groups)
    merge <- names(group_tb)[group_tb < min_num]
    groups[groups %in% merge] <- "merged"
    return(data.frame(clusters = names(groups), groups = groups, stringsAsFactors = FALSE))
    }
  }
}
# gr <- CutDendroComp(dendro = dendro, n_groups = 25, plot = TRUE, min_num = 30, labels = FALSE)
#
# df <- data.frame(clusters = names(groups), groups = groups)
# table(groups)
# g <- CutDendroComp(dendro, n_groups = 10, plot = TRUE)
# png(filename = "~/R/benchmark/dendro_probe.png")
# plot(dendro)
# dev.off()
# png(filename = "~/R/benchmark/dendro_probe.png")
# CutDendroComp(dendro, n_groups = 10, plot = TRUE)
# dev.off()
# png(filename = "~/R/benchmark/dendro_probe2.png")
# CutDendroComp(dendro, n_groups = 30, plot = TRUE)
# dev.off()
