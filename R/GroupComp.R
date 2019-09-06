#' Creare groups among clusters
#'
#' @param compare a compare object
#' @param groups A data.frame object with cluster names in the first column and group names in the second column
#' @param event_data A character or a numeric vector defining datasets to modify
#'
#' @return this function returns a compare object with grouped clusters
#' @export
#' @import ape readr plyr stringr tidyr
#'
#' @examples
#######
## fontos hogy kell egy parameter amivel suffixokat lehet adni
##meg kell oldani hogy ha node-okat fajokat kizarunk akkor is tudjon mukodni
#######
# compare <- probe1_sliced2
# groups <- group_df
# event_data <- c(2,3)
# k <- 1
GroupComp <- function(compare, groups = NULL, event_data = NULL){
  if(is.null(groups)){
    stop(paste0("Please give a data.frame which contains cluster names and group names!\n"))
  }
  if(is.null(event_data)){
    if(!"event_data" %in% names(compare)){
      stop(paste0("There is no data called event_data!\n
                  Please define names or numbers refere to compare event data!\n"))
    }
    if(any(!as.character(unlist(groups[,1])) %in% colnames(compare$event_data$gains))){
      warning(paste0("some of user defined cluster names are not in the data: ",
                     as.character(unlist(groups[,1]))[!as.character(unlist(groups[,1])) %in% colnames(compare$event_data$gains)], "\n"))
    }
      group_names <- unique(groups[[2]])
      group_cl_sep <- vector("list", length = length(group_names))
      for(i in seq_along(group_names)){
        group_cl_sep[[i]] <- unlist(groups[unlist(groups[,2]) %in% group_names[i], 1])
      }
      grouped_data <- vector("list", length = length(compare$event_data))
      for(j in 1:length(compare$event_data)){
        grouped_event <- sapply(group_cl_sep, function(x) rowSums(compare$event_data[[j]][,colnames(compare$event_data[[j]]) %in% as.character(x)]))
        colnames(grouped_event) <- group_names
        grouped_data[[j]] <- grouped_event
      }
      names(grouped_data) <- names(compare$event_data)
      old_names <- names(compare)
      compare <- c(compare[-length(compare)], list(grouped_data),
                   compare["raw_data"])
      names(compare) <- c(old_names[-length(old_names)],
                          "event_data_GROUPED", "raw_data")
      return(compare)
  }else{
    if(is.character(event_data)){
      event_data <- which(names(compare) %in% event_data)
    }else{
      if(!is.numeric(event_data)){
        stop(paste0("You need to give a character or a numeric vector defininf the data you want to modifiy!\n"))
      }
    }
    for(k in 1:length(event_data)){
    if(any(!as.character(unlist(groups[,1])) %in% colnames(compare[[event_data[k]]]$gains))){
        warning(paste0("some of user defined cluster names are not in the data: ",
                       as.character(unlist(groups[,1]))[!as.character(unlist(groups[,1])) %in% colnames(compare[[event_data[k]]]$gains)], "\n"))
      }
      group_names <- unique(groups[[2]])
      group_cl_sep <- vector("list", length = length(group_names))
      for(i in seq_along(group_names)){
        group_cl_sep[[i]] <- unlist(groups[unlist(groups[,2]) %in% group_names[i], 1])
      }
      grouped_data <- vector("list", length = length(compare[[event_data[k]]]))
      for(j in 1:length(compare[[event_data[k]]])){
        grouped_event <- sapply(group_cl_sep, function(x) rowSums(compare[[event_data[k]]][[j]][,colnames(compare[[event_data[k]]][[j]]) %in% as.character(x)]))
        colnames(grouped_event) <- group_names
        grouped_data[[j]] <- grouped_event
      }
      names(grouped_data) <- names(compare[[event_data[k]]])
      old_names <- names(compare)
      compare <- c(compare[-length(compare)], list(grouped_data),
                   compare["raw_data"])
      names(compare) <- c(old_names[-length(old_names)],
                          paste0(old_names[event_data[k]], "_GROUPED"), "raw_data")
    }
    return(compare)
    }
}
