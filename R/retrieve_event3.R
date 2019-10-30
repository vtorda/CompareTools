#' Retrieve event3
#'
#' @param compare a compare object
#' @param partitioning logical, does the clusters handled separately (TRUE) or all together (FALSE).
#'     If TRUE and no groups is given, all clusters will be handeled separately.
#'     If a group file is given clusteres will be grouped according to that file.
#' @param groups a Tab Separated file containing the cluster names and their position in predefined groups
#' @param part_name
#'
#' @return this function returns a compare object with additional event matrices
#' @export
#' @import ape readr plyr stringr tidyr
#'
#' @examples
#######
## majd meg beepiteni hogy a part_name-eket csekkolja, szerepel-e mar ez a nev a listaba
## aztan hogy itt is ellenorizze hogy minden node meg van a faban es a raw databan..
#######
retrieve_event3 <- function(compare,  new_object = FALSE){
  if(any(names(compare) %in% "event_data")){
    stop(paste0("The compare object have already had an element names event_data!"))
  }
  raw_data <- compare[["raw_data"]]
  tree <- compare$tree
  node_n <- nrow(raw_data)
  #first retrieve all events (pattern = " ") then split it to the clusters and orthogroups (pattern "/")
  cluster_names <- na.omit(unique(str_split(unique(unlist(str_split(raw_data,
                                                                    pattern = " "))),
                                            pattern = "/", simplify = TRUE)[,1]))
  all_nodeg_sep <- vector("list", length = node_n+1)
  all_nodel_sep <- vector("list", length = node_n+1)
  for(l in 1:node_n){
    clfreq <- table(sapply(str_split(unlist(str_split(raw_data[l, 1],
                                                      " ")), "/"), "[", 1))
    all_nodeg_sep[[l]] <- clfreq
    clfreq2 <- table(sapply(str_split(unlist(str_split(raw_data[l, 2],
                                                       " ")), "/"), "[", 1))
    all_nodel_sep[[l]] <- clfreq2
  }

  all_nodeg_sep[[node_n+1]] <- table(cluster_names)
  all_nodel_sep[[node_n+1]] <- table(cluster_names)

  Gxx <- ldply(all_nodeg_sep, rbind) %>% replace(is.na(.), 0)
  Gxx <- as.matrix(Gxx[-(node_n+1),])
  Gxx <- Gxx[,order(colnames(Gxx))]
  Lxx <-ldply(all_nodel_sep,rbind) %>% replace(is.na(.), 0)
  Lxx <- as.matrix(Lxx[-(node_n+1),])
  Lxx <- Lxx[,order(colnames(Lxx))]
  all_cl_sep <- list(Gxx, Lxx)
  names(all_cl_sep) <- c("gains", "losses")
  if(new_object){
    compare_ls2 <- c(compare[1], list(all_cl_sep),
                     compare["raw_data"])
    names(compare_ls2) <- c(names(compare)[1],
                            "event_data", "raw_data")
    return(compare_ls2)
  }else{
    compare_ls2 <- c(compare[1:(length(compare) - 1)], list(all_cl_sep),
                   compare["raw_data"])
    names(compare_ls2) <- c(names(compare)[1: (length(compare) - 1)],
                          "event_data", "raw_data")
    return(compare_ls2)
  }
}
