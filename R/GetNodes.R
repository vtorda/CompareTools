#' Get node numbers of subclades
#'
#' @param compare a compare object
#' @param species Character vector. At least two species which specifiy a monophyletic clade
#' @param tips whether include species or not
#'
#' @return this function returns a compare object with grouped clusters
#' @export
#' @import ape readr plyr stringr tidyr
#'
#' @examples
GetNodes <- function(compare, species, tips = TRUE){
  mrca <- getMRCA(compare$tree, species)
  all_tips <- extract.clade(compare$tree, node = mrca)$tip.label
  w_count <- 1
  node_v <- vector()
  for(i in 1:length(all_tips)){
    tipn <- which(compare$tree$tip.label %in% all_tips[i])
    node_v[w_count] <- tipn
    w_count <- w_count + 1
    while(tipn != mrca){
    tipn <- compare$tree$edge[compare$tree$edge[,2] %in% tipn, 1]
    node_v[w_count] <- tipn
    w_count <- w_count + 1
    }
  }
  node_v <- unique(node_v)
  if(tips){
    return(node_v)
  }else{
    excl <- which(compare$tree$tip.label %in% all_tips)
    node_v <- node_v[!node_v %in% excl]
    return(node_v)
  }
}
