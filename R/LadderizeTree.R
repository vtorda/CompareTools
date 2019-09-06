#' ladderize tree
#'
#' ladderizing a phylo object by not just indexing but reordering the tip labels
#'
#' @param tree a \code{phylo} object
#' @param temp_file a character specifiying a a temporary directory see details
#' @param orientation
#'
#' @import ape
#'
#' @return this function returns \code{phylo} object which is ladderized with new tip label order
#' @export
#'
#'
#' @examples
#'

# tree <- rtree(15)
# plot.phylo(tree)
# tree_lad_ape <- ladderize(tree, right = FALSE)
# plot.phylo(tree_lad_ape)
# #comparing the object structures: edge matris has changed but not the tip.label order
# identical(tree, tree_lad_ape)
# identical(tree$tip.label, tree_lad_ape$tip.label)
# tree_lad_compare <- LadderizeTree(tree)
# plot.phylo(tree_lad_ape, main = "ape function")
# plot.phylo(tree_lad_compare, main = "modified ape function")
# identical(tree_lad_ape, tree_lad_compare)
# identical(tree_lad_ape$tip.label, tree_lad_compare$tip.label)
# # now the tip label order will represent the order you see on the plot
# @author Torda Varga
# @seealso \link[ape]{ladderize}
LadderizeTree <- function(tree, temp_file = "temp", orientation = "left"){
  if(file.exists(paste0("./", temp_file))){
    stop("The chosen temporary file exists! Please choose an other temp_file name")
  }
  if(orientation == "left"){
    right <- FALSE
  }else{
    right <- TRUE
  }
  tree_temp <- ladderize(tree, right = right)
  write.tree(tree_temp, file = paste0("./", temp_file, ".tre"))
  tree_lad <- read.tree(paste0("./", temp_file, ".tre"))
  file.remove(paste0("./", temp_file, ".tre"))
  return(tree_lad)
}
