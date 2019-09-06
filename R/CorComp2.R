#' Perform pair-wise correlation of clusters based on chosen statistics
#'
#' @param compare a compare object
#' @param dataset The dataset that you want to depict. It can be a single character or a numeric value refers to the element of the compare object.
#' @param stat A character or a numeric value that refers to the statistics (elemnt of the dataset) you want ot depict.
#' @param plot Logical. If FALSE you have to define the group_n argumentum and a
#' dataframe is returned with cluster names and groups a according hierachical clustering.
#' This dataframe can be inputed to the GoupComp function. Otherwise you can plot a dendrogram with our withour a heatmap
#' @param heatmap Logical. Whether a heatmap should be ploted. Only when heatmap is FALSE you can draw rectangles around hierarchical cluster.
#' @param group_n Numeric. Defining the number of groups which into the clusters should be sorted
#' @param cor_method A character string defining the  correlation coefficient you want to use: "pearson", "kendall" or "spearman". The deafault is "pearson".
#' @param cl_method A character string defining the agglomeration method to be used. See hclust {stats} function for details.
#' @param ... argument passed to the hclust function
#'
#'
#' @return this function returns a compare object retain only just part of the original clusters
#' @export
#' @import ape readr plyr stringr tidyr
#'
#' @examples
############
## be kene epiteni meg a pvclust fuggvenyeket, hogy bootstrap ertekeket lehessen szamolni
############
CorComp <- function(compare, dataset = NULL, stat = NULL, plot = TRUE, heatmap = TRUE,
                    group_n = NULL, cor_method = "pearson", cl_method = "ward.D2", ...){
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
  compare[[dataset]][[stat]] <- cor(compare[[dataset]][[stat]], method = cor_method)
  if(plot){
  if(is.null(group_n)){
    if(heatmap){
      heatmap(compare[[dataset]][[stat]], keep.dendro = FALSE, hclustfun = function(x) hclust(as.dist(1-compare[[dataset]][[stat]]), method = cl_method, ...))
    }else{
      plot(hclust(as.dist(1-compare[[dataset]][[stat]]), method = cl_method, ...),
           main = paste0("Cluster Dendrogram of the ", names(compare)[dataset], " dataset and the ", names(compare[[dataset]])[stat], " statistics"),
           xlab = "Cluster names")
    }
  }else{
    h <- hclust(as.dist(1-compare[[dataset]][[stat]]), method = cl_method, ...)
    groups <- cutree(h, k = group_n)
    plot(h,
         main = paste0("Cluster Dendrogram of the ", names(compare)[dataset], " dataset and the ", names(compare[[dataset]])[stat], " statistics"),
         xlab = "Cluster names")
    rect.hclust(h, k = group_n, border = "red")
  }
  }else{
    if(!is.null(group_n)){
      h <- hclust(as.dist(1-compare[[dataset]][[stat]]), method = cl_method, ...)
      groups <- cutree(h, k = group_n)
      df <- data.frame(clusters = names(groups), groups = groups)
      return(df)
    }else{
      stop(paste0("You need to define the group_n argumentum!"))
    }
  }
}

