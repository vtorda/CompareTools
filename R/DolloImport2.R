#' DolloImport2
#'
#' Import data from compare analysis
#'
#' @param tree a phylo object describing the species tree used in the COMPARE analysis
#' @param files a character containing the path with the name of the the output file from the step Dollo parsimony of the COMPARE analysis.
#' @param partitioning
#' @param groups
#' @param part_name
#'
#' @return this function returns a \code{list} object with the following elements: tree, Compare matrices and raw_data. raw data is a matrix with two columns (1 = gains, 2 = losses) and as many rows as many nodes + tips the species tree has. These are aligned.
#' @export
#' @import ape readr plyr stringr tidyr
#'
#' @examples
#'
#'
DolloImport2 <- function(tree, files = NULL, partitioning = FALSE, groups = NULL,
                         part_name = NULL){
    if(is.null(files)){
      stop("You need to give a path with the file name!")
      }
    #import raw data
    raw_file <- read_lines(files)
    rows_needed <- c(which(str_detect(raw_file, "GAIN") == TRUE), which(str_detect(raw_file, "LOSSES") == TRUE))
    raw_df <- raw_file[rows_needed]
    raw_df2 <- str_split(raw_df, "\t")
    raw_df2 <- ldply(raw_df2, rbind)
    raw_df2 <- raw_df2[,-c(2,4)]
    colnames(raw_df2) <- c("event", "node", "chars")
    # itt eloszor le kene gyartani az objektumot a tree-vel es a raw file-al
    #match the phylo object nodes and the compare output nodes:
    node_def <- which(str_detect(raw_file, "node\\d") == TRUE)
    node_df <- raw_file[node_def]
    node_df2 <- str_split(node_df, "\t")
    node_df2 <- ldply(node_df2, rbind)
    colnames(node_df2) <- c("nodes", "species")
    node_df2$nodes <- str_replace(node_df2$nodes, "node", "")
    node_df2$species <- str_sub(node_df2$species, start = 2) # spacet kitakaritom
    dollo_species <- unique(unlist(sapply(node_df2$species, function(x) str_split(x, " "))))
    if(any(dollo_species %in% tree$tip.label == FALSE)){
      wrong_names_d <- which(dollo_species %in% tree$tip.label == FALSE)
      warning(paste0("The following species names don't match!\t", "Species in Dollo file:\t", dollo_species[wrong_names]))
    }
    node_df2$tree_nodes <- rep(NA, nrow(node_df2))
    for(i in 1:nrow(node_df2)){
      species <- unlist(str_split(node_df2[i,"species"], " "))
      if(length(species) > 1){
        node_df2[i, "tree_nodes"] <- getMRCA(tree, species)
      }else{
        node_df2[i, "tree_nodes"] <- which(tree$tip.label == species)
      }
    }
    compare_ls <- list()
    compare_ls[[1]] <- list()
    compare_ls[[1]] <- tree
    #create the raw matrix
    node_df2$gain_event <- rep(NA, nrow(node_df2))
    node_df2$loss_event <- rep(NA, nrow(node_df2))
    gain_split <- raw_df2[raw_df2$event == "GAIN", ]
    loss_split <- raw_df2[raw_df2$event == "LOSSES", ]
    node_df2$gain_event <- gain_split$chars[match(node_df2$nodes, gain_split$node)]
    node_df2$loss_event <- loss_split$chars[match(node_df2$nodes, loss_split$node)]
    node_df2 <- node_df2[,-c(1,2)]
    node_df2 <- node_df2[order(node_df2$tree_nodes),]
    raw_matrix <- as.matrix(node_df2[,2:3])
    raw_matrix <- unname(raw_matrix)
    compare_ls[[2]] <- raw_matrix
    names(compare_ls) <- c("tree", "raw_data")
    # use the retrieve_event function
    events_ls <- retrieve_event2(compare = compare_ls, partitioning = partitioning, groups = groups, part_name = part_name)
    return(events_ls)
}
# file_path <- system.file("extdata", "Dolloout", package = "CompareTools")
# files <- file_path
# tree <- fungi_tree_lad
# group_path <- system.file("extdata", "groups.txt", package = "CompareTools")
# groups <- group_path
# groups <- read_tsv(group_path)
# groups <- NULL
# partitioning <- TRUE
# part_name <- "probe"
# probe1 <- DolloImport(tree = tree, files = files, partitioning = FALSE, part_name = "no_partitioning")
# probe2 <- DolloImport(tree = tree, files = files, partitioning = TRUE, part_name = "every_cl_part")
# probe3 <- DolloImport(tree = tree, files = files, partitioning = TRUE, part_name = "grouping1", groups = groups)
# probe4 <- CompareTools::retrieve_event(probe1, partitioning = TRUE, groups = groups, part_name = "grouping1")
