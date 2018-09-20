#' Retrieve event
#'
#' @param compare
#' @param partitioning
#' @param groups
#' @param part_name
#'
#' @return this function returns a compare object with additional event matrices
#' @export
#'
#' @examples
#######
## majd meg beepiteni hogy a part_name-eket csekkolja, szerepel-e mar ez a nev a listaba
#######
# compare <- compare_ls
# group_path <- system.file("extdata", "groups.txt", package = "CompareTools")
# groups <- c(group_path, "bla")
# groups <- read_tsv(group_path)
# part_title <- "test"
# # tehat akkor eloszor
# part_name <- "probe"
# compare <- probe1
# partitioning <- TRUE
# groups <- groups
# part_name <- "grouping1"
retrieve_event <- function(compare, partitioning = FALSE, groups = NULL, part_name = NULL){
  raw_data <- compare[["raw_data"]]
  tree <- compare$tree
  node_n <- nrow(raw_data)
  #first retrieve all events (pattern = " ") then split it to the clusters and orthogroups (pattern "/")
  cluster_names <- na.omit(unique(str_split(unique(unlist(str_split(raw_data, pattern = " "))), pattern = "/", simplify = TRUE)[,1]))
  all_cl_sep <- vector("list", length = length(cluster_names))
  for(i in seq_along(cluster_names)){
    gains <- str_count(raw_data[,1], paste0(cluster_names[i], "/\\d"))
    losses <- str_count(raw_data[,2], paste0(cluster_names[i], "/\\d"))
    gains[is.na(gains)] <- 0
    losses[is.na(losses)] <- 0
    all_cl_sep[[i]] <- matrix(c(gains, losses), nrow = node_n, ncol = 2)
  }
  names(all_cl_sep) <- cluster_names
  # create node paths
  j <- NULL
  i <- NULL
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
  j <- NULL
  i <- NULL
# calculate the commulative gains and losses
#most meg van a path akkor mindengyiket ki kell bontogatni es osszeadni a dolgokat
  j <- 1
  k <- 1
  node_ls <- vector("list", node_n)
  for(j in 1:node_n){
    path_temp <- path[j]
    node_ls[[j]] <- as.numeric(str_split(path_temp, ";", simplify = TRUE)[1,])
  }


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
  #ezek utan mar csak ossze kell addni az adott matrixokat es meg lesz a particionalt matrix
  if(partitioning){
    if(is.null(groups)){
      compare_ls2 <- c(compare[1: (length(compare) - 1)], list(all_cl_sep), compare["raw_data"])
      names(compare_ls2) <- c(names(compare)[1: (length(compare) - 1)], part_name, "raw_data")
      return(compare_ls2)
    }else{
      ################################
      ######write a warning if the group file is not in the correct form########
      ##### does it contain every cluster?
      ################################
      if(is.character(groups) & length(groups) > 1){
        stop("wrongly defined group path\nPlease give a single character with the absolute path of the group file!")
      }else{
        if(is.character(groups) & length(groups) == 1){
          groups <- read_tsv(groups)
        }
        if(!is.data.frame(groups)){
          stop("provide a data.frame object")
        }else{
          ###
          ### hogy lehet hogy a groups file-ban nem az osszes cluster nev szerepel? mig a cun van van 350!
          ###
        }
      }
      # mivel mar minden ki van szamolva, ezert csak ossze kell adogatni
      group_names <- unique(groups[[2]])
      #i <- 1
      group_cl_sep <- list()
      for(i in 1:length(group_names)){
        clusters_names <- groups[groups[[2]] %in% group_names[i],][[1]]
        pr <- all_cl_sep[which(names(all_cl_sep) %in% as.character(clusters_names))]
        group_cl_sep[[i]] <- Reduce("+", pr)
      }
      names(group_cl_sep) <- group_names
      #length(group_cl_sep)
      compare_ls2 <- c(compare[1: (length(compare) - 1)], list(group_cl_sep), compare["raw_data"])
      names(compare_ls2) <- c(names(compare)[1: (length(compare) - 1)], part_name, "raw_data")
      return(compare_ls2)
    }

  }else{
    all_cl <- list()
    all_cl[[1]] <- Reduce("+", all_cl_sep)
    compare_ls2 <- c(compare[1: (length(compare) - 1)], list(all_cl), compare["raw_data"])
    names(compare_ls2) <- c(names(compare)[1: (length(compare) - 1)], part_name, "raw_data")
    return(compare_ls2)
  }
}
# as.numeric(str_split(path_temp, ";", simplify = TRUE)[1,])
# probe1 <- retrieve_event(compare_ls)
# probe2 <- retrieve_event(compare_ls, partitioning = TRUE)
# probe3 <- retrieve_event(compare_ls, partitioning = TRUE, groups = groups)
# str_count(raw_data[,1], paste0(cluster_names[1], "/\\d"))
# sapply(raw_data[,1], function(x) str_count(x, paste0(cluster_names[1], "/\\d")), USE.NAMES = FALSE)
# gains <- sapply(raw_data[,1], function(x) str_count(x, paste0(cluster_names[i], "/\\d")), USE.NAMES = FALSE)
# losses <- sapply(raw_data[,2], function(x) str_count(x, paste0(cluster_names[i], "/\\d")), USE.NAMES = FALSE)
# gains[is.na(gains)] <- 0
# losses[is.na(losses)] <- 0
# all_cl_sep[[i]] <- matrix(c(gains, losses), nrow = node_n, ncol = 2)
