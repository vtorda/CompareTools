#' Transform compare object
#'
#' @param compare a compare object
#' @param builtin A character of One of the four built in transformation: cummulative gain (CumGain)
#' cummulative losses (CumLoss), netto gain (NetGain) and netto cummulative gain (NetCumGain).
#' @param custom A user defined function
#'
#' @return this function returns a compare object retain only just part of the original clusters
#' @export
#' @import ape readr plyr stringr tidyr
#'
#' @examples
#######
## meg kell oldani hogy ha node-okat fajokat kizarunk akkor is tudjon mukodni
#######
TransformComp <- function(compare, builtin = c("CumGain", "CumLoss", "NetGain", "NetCumGain"), custom = NULL){
  if(is.null(custom)){
    switch(builtin,
           NetGain = {
             for(i in 1:(length(compare)-2)){
               element_n <- length(compare[[i+1]])
               compare[[i+1]][[element_n + 1]] <- compare[[i+1]][["gains"]] - compare[[i+1]][["losses"]]
               names(compare[[i+1]])[element_n + 1] <- "NetGain"
             }
            return(compare)
           },
          CumGain = {
            tree <- compare$tree
            node_n <- nrow(compare$raw_data)
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
            for(k in 1:(length(compare)-2)){
              gains <- compare[[k+1]]$gains
              iter_v <- c(1:nrow(gains))[-root_no]
              for(j in iter_v){
                gains[j, ] <- colSums(compare[[k+1]]$gains[node_ls[[j]], ])
              }
              element_n <- length(compare[[k+1]])
              compare[[k+1]][[element_n+1]] <- gains
              names(compare[[k+1]])[element_n+1] <- "CumGain"
            }
            return(compare)
          },
          CumLoss = {
            tree <- compare$tree
            node_n <- nrow(compare$raw_data)
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
            for(k in 1:(length(compare)-2)){
              losses <- compare[[k+1]]$losses
              iter_v <- c(1:nrow(losses))[-root_no]
              for(j in iter_v){
                losses[j, ] <- colSums(compare[[k+1]]$losses[node_ls[[j]], ])
              }
              element_n <- length(compare[[k+1]])
              compare[[k+1]][[element_n+1]] <- losses
              names(compare[[k+1]])[element_n+1] <- "CumLoss"
            }
            return(compare)
          },
          NetCumGain = {
            tree <- compare$tree
            node_n <- nrow(compare$raw_data)
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
            for(k in 1:(length(compare)-2)){
              losses <- compare[[k+1]]$losses
              gains <- compare[[k+1]]$gains
              iter_v <- c(1:nrow(losses))[-root_no]
              for(j in iter_v){
                losses[j, ] <- colSums(compare[[k+1]]$losses[node_ls[[j]], ])
                gains[j, ] <- colSums(compare[[k+1]]$gains[node_ls[[j]], ])
              }
              element_n <- length(compare[[k+1]])
              compare[[k+1]][[element_n+1]] <- gains - losses
              names(compare[[k+1]])[element_n+1] <- "NetCumGain"
            }
            return(compare)
          })
  }else{

  }
}
