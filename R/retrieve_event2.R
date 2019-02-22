#' Retrieve event2
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
#######

retrieve_event2 <- function(compare, partitioning = FALSE, groups = NULL,
                           part_name = NULL){
  raw_data <- compare[["raw_data"]]
  tree <- compare$tree
  node_n <- nrow(raw_data)
  #first retrieve all events (pattern = " ") then split it to the clusters and orthogroups (pattern "/")
  cluster_names <- na.omit(unique(str_split(unique(unlist(str_split(raw_data,
                                                                    pattern = " "))),
                                            pattern = "/", simplify = TRUE)[,1]))
# create node paths
root_no <- getMRCA(tree, tree$tip.label)
path <- vector(mode = "character", length = node_n)
  for(j in 1:length(path)){
    node <- j
    node_init <- node
    node_path <- vector()
    count <- 1
    while(node != root_no){
      node_path[count] <- tree$edge[tree$edge[,2] == node,1]
      node <- node_path[count]
      count <- count + 1
    }
    path[j] <- str_c(c(node_init, node_path), collapse = ";")
  }
node_ls <- vector("list", node_n)
  for(j in 1:node_n){
    path_temp <- path[j]
    node_ls[[j]] <- as.numeric(str_split(path_temp, ";", simplify = TRUE)[1,])
  }
if(partitioning){
  if(is.null(groups)){
    all_cl_sep <- vector("list", length = length(cluster_names))
    for(i in seq_along(cluster_names)){
      gains <- str_count(raw_data[,1], paste0(cluster_names[i], "/\\d"))
      losses <- str_count(raw_data[,2], paste0(cluster_names[i], "/\\d"))
      gains[is.na(gains)] <- 0
      losses[is.na(losses)] <- 0
      all_cl_sep[[i]] <- matrix(c(gains, losses), nrow = node_n, ncol = 2)
     }
    names(all_cl_sep) <- cluster_names
    # calculate the commulative values
    for(k in 1:length(all_cl_sep)){
      copy_num1 <- vector(mode = "numeric", length = node_n)
      copy_num2 <- vector(mode = "numeric", length = node_n)
      for(j in 1:node_n){
        copy_num1[j] <- sum(all_cl_sep[[k]][node_ls[[j]], 1])
        copy_num2[j] <- sum(all_cl_sep[[k]][node_ls[[j]], 2])
      }
      all_cl_sep[[k]] <- cbind(all_cl_sep[[k]], copy_num1)
      all_cl_sep[[k]] <- cbind(all_cl_sep[[k]], copy_num2)
      all_cl_sep[[k]] <- unname(all_cl_sep[[k]])
    }
    compare_ls2 <- c(compare[1:(length(compare) - 1)], list(all_cl_sep),
                       compare["raw_data"])
    names(compare_ls2) <- c(names(compare)[1: (length(compare) - 1)],
                              part_name, "raw_data")
    return(compare_ls2)
  }else{
  ################################
  ######write a warning if the group file is not in the correct form, dim check########
  ##### does it contain every cluster?
  ################################
    if(is.character(groups) & length(groups) > 1){ #### here the states should be defined more precisely!!!
      stop("wrongly defined group path\nPlease give a single character with the absolute path of the group file!")
    }else{
      if(is.character(groups) & length(groups) == 1){
        groups_df <- read_tsv(groups)
      }
      if(!is.data.frame(groups_df)){
        stop("Provide a Tab Separated Value (TSV) file in the right form. See man page")
      }else{
        ###
        ### hogy lehet hogy a groups file-ban nem az osszes cluster nev szerepel? mig a cun van van 350! ugy hogy nem feltetel hogy minend clustert belevonjunk egy-egy groupolasba
        ###
      }
    }
  group_names <- unique(groups_df[[2]])
  group_cl_sep <- vector("list", length = length(group_names))
  for(i in seq_along(group_names)){
    group1 <- unlist(groups_df[unlist(groups_df[,2]) %in% group_names[i], 1])
    clust_seq <- str_c(group1, "/\\d", collapse = "|")
    gains <- str_count(raw_data[,1], clust_seq)
    losses <- str_count(raw_data[,2], clust_seq)
    gains[is.na(gains)] <- 0
    losses[is.na(losses)] <- 0
    group_cl_sep[[i]] <- matrix(c(gains, losses), nrow = node_n, ncol = 2)
  }
  for(k in 1:length(group_cl_sep)){
    copy_num1 <- vector(mode = "numeric", length = node_n)
    copy_num2 <- vector(mode = "numeric", length = node_n)
    for(j in 1:node_n){
      copy_num1[j] <- sum(group_cl_sep[[k]][node_ls[[j]], 1])
      copy_num2[j] <- sum(group_cl_sep[[k]][node_ls[[j]], 2])
    }
    group_cl_sep[[k]] <- cbind(group_cl_sep[[k]], copy_num1)
    group_cl_sep[[k]] <- cbind(group_cl_sep[[k]], copy_num2)
    group_cl_sep[[k]] <- unname(group_cl_sep[[k]])
  }
  names(group_cl_sep) <- group_names
  compare_ls2 <- c(compare[1:(length(compare) - 1)], list(group_cl_sep),
                       compare["raw_data"])
  names(compare_ls2) <- c(names(compare)[1: (length(compare) - 1)],
                              part_name, "raw_data")
  return(compare_ls2)
  }
}else{
  all_cl <- list()
  all_cl <- vector("list", length = 1)
  gains <- str_count(raw_data[,1], "\\d/\\d")
  losses <- str_count(raw_data[,2], "\\d/\\d")
  gains[is.na(gains)] <- 0
  losses[is.na(losses)] <- 0
  all_cl[[1]] <- matrix(c(gains, losses), nrow = node_n, ncol = 2)
  copy_num1 <- vector(mode = "numeric", length = node_n)
  copy_num2 <- vector(mode = "numeric", length = node_n)
  for(j in 1:node_n){
    copy_num1[j] <- sum(all_cl[[1]][node_ls[[j]], 1])
    copy_num2[j] <- sum(all_cl[[1]][node_ls[[j]], 2])
  }
  all_cl[[1]] <- cbind(all_cl[[1]], copy_num1)
  all_cl[[1]] <- cbind(all_cl[[1]], copy_num2)
  all_cl[[1]] <- unname(all_cl[[1]])
  names(all_cl) <- "all_cluster"
  compare_ls2 <- c(compare[1:(length(compare) - 1)], list(all_cl),
                   compare["raw_data"])
  names(compare_ls2) <- c(names(compare)[1: (length(compare) - 1)],
                          part_name, "raw_data")
  return(compare_ls2)
  }
}
