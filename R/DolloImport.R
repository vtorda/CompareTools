#' DolloImport
#'
#' Import data from compare analysis
#'
#' @param tree
#' @param files
#' @param path
#' @param pattern
#' @param partitioning
#' @param groups
#'
#' @return this function returns a \code{list} object with...
#' @export
#' @import ape readr plyr stringr
#'
#' @examples
#'
#'
files <- file_path
tree <- fungi_tree_lad
group_path <- system.file("extdata", "groups.txt", package = "CompareTools")
groups <- group_path
groups <- read_tsv(group_path)
partitioning <- TRUE
DolloImport <- function(tree, files = NULL, partitioning = FALSE, groups = NULL){
    if(is.null(files)){
      warning(paste0("You need to give a file name!")) # TODO stopra atirni
    }
    raw_file <- read_lines(files)
    rows_needed <- c(which(str_detect(raw_file, "GAIN") == TRUE), which(str_detect(raw_file, "LOSSES") == TRUE))
    raw_df <- raw_file[rows_needed]
    raw_df2 <- str_split(raw_df, "\t")
    raw_df2 <- data.frame(matrix(unlist(raw_df2), nrow = length(raw_df2), byrow = TRUE), stringsAsFactors = FALSE)
    ######################Zsolti
    if(partitioning){
           cat("Partitioning: ")
           raw_df2z <- raw_df2[,-c(ncol(raw_df2) - 1, ncol(raw_df2) - 3)]
           colnames(raw_df2z) <- c("changes", "event", "node", "chars")
           cutter <- separate_rows(raw_df2z, chars, sep = " ")
           cutter <- separate(cutter, chars, sep = "/", into = c("clusters","orthogroups"))
           cun <- sort(unique(cutter$clusters))
             if(is.null(groups)){
               allcl <- lapply(cun, function(x) cutter[which(cutter$clusters == x),])
             }else{
               ################################
               ######write a warning if the group file is not in the correct form########
               ##### does it contain every cluster?
               ################################
               partition <- unique(groups[[2]])
               allcl <- lapply(partition, function(x) cutter[which(cutter$clusters %in% groups[groups[,2] == x, 1][[1]]),])
               cun <- partition
             }
           all1 <- allcl[[1]]
           df_listC <- list()
           akt_listC <- list()
           for (jkl in 1:length(allcl)){
             #jkl=1
             cat(jkl," ")
             akt <- allcl[[jkl]]
             glcounts <- table(akt$node,akt$event)

               if(any(colnames(glcounts) %in% "GAIN")){
                 akt[akt$event=="GAIN","changes"] <- glcounts[match(akt[akt$event=="GAIN","node"], rownames(glcounts)),1] # we modify the number of changes (gain) in the separated cluster table
               }
               if(any(colnames(glcounts) %in% "LOSSES")){
                 akt[akt$event=="LOSSES","changes"]<-glcounts[match(akt[akt$event=="LOSSES","node"], rownames(glcounts)),2] # we modify the number of changes (loss) in the separated cluster table
               }
             raw_df2 <- akt[,c(1:3)]
             nodes <- unique(raw_df2$node)
             df_split <- data.frame(nodes = nodes, gains = rep(NA, length(nodes)),
                                    losses = rep(NA, length(nodes)), stringsAsFactors = FALSE)
             df_split_gain <- raw_df2[raw_df2$event %in% "GAIN",]
             df_split_losses <-  raw_df2[raw_df2$event %in% "LOSSES",]
             gains_match <- match(df_split$nodes, df_split_gain$node)
             losses_match <- match(df_split$nodes, df_split_losses$node)
             df_split$gains[!is.na(gains_match)] <- df_split_gain[na.omit(gains_match), "changes"]
             df_split$losses[!is.na(losses_match)] <- df_split_losses[na.omit(losses_match), "changes"]
             node_def <- which(str_detect(raw_file, "0\tnode") == TRUE)
             node_df <- raw_file[node_def]
             node_df2 <- str_split(node_df, "\t")
             node_df2 <- ldply(node_df2, rbind)
             #node_df2y <- data.frame(matrix(unlist(node_df2), nrow=length(node_df2), byrow=T),stringsAsFactors=FALSE)
             node_df2 <- node_df2[,-1]
             colnames(node_df2) <- c("nodes", "species")
             node_df2$nodes <- str_replace(node_df2$nodes, "node", "")
             node_df2$species <- str_sub(node_df2$species, start = 2) # spacet kitakaritom

             dollo_species <- unique(unlist(sapply(node_df2$species, function(x) str_split(x, " "))))
             if(any(dollo_species %in% tree$tip.label == FALSE)){
               wrong_names_d <- which(dollo_species %in% tree$tip.label == FALSE)
               warning(paste0("Species names don't match!\t", "Species in Dollo file:\t", dollo_species[wrong_names]))
             }else{
               warning(paste0("Species names match!"))
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
             node_df2$gains <- rep(NA, nrow(node_df2))
             node_df2$losses <- rep(NA, nrow(node_df2))
             node_match <- match(node_df2$nodes, df_split$nodes)
             node_df2$gains[!is.na(node_match)] <- df_split[na.omit(node_match), "gains"]
             node_df2$losses[!is.na(node_match)] <- df_split[na.omit(node_match), "losses"]
             node_df2[is.na(node_df2$losses), "losses"] <- 0
             node_df2[is.na(node_df2$gains), "gains"] <- 0
             node_df2$tree_nodes <- as.numeric(node_df2$tree_nodes)
             node_df2$gains <- as.numeric(node_df2$gains)
             node_df2$losses <- as.numeric(node_df2$losses)
             node_df2$net_gains <- node_df2$gains - node_df2$losses
             j <- NULL
             i <- NULL
             root_no <- getMRCA(tree, tree$tip.label)
             node_df2$path <- rep(NA, nrow(node_df2))
             for(j in 1:nrow(node_df2)){
               node <- node_df2[j,"tree_nodes"]
               node_init <- node
               node_path <- vector()
               count <- 1
               while(node != root_no){
                 node_path[count] <- tree$edge[tree$edge[,2] == node,1]
                 node <- node_path[count]
                 count <- count + 1
               }
               node_df2[j, "path"] <- str_c(c(node_init, node_path), collapse = ";")
             }
             j <- NULL
             i <- NULL
             #most meg van a path akkor mindengyiket ki kell bontogatni es osszeadni a dolgokat
             node_df2$copy_num <- rep(NA, nrow(node_df2))
             for(j in 1:nrow(node_df2)){
               path <- node_df2[j, "path"]
               v <- as.numeric(unlist(str_split(path, ";")))
               v_l <- length(v)
               nums <- 0
               for(i in 0:(v_l-1)){
                 nums <- nums + node_df2[node_df2$tree_nodes %in% v[v_l - i], "net_gains"]
               }
               node_df2[j, "copy_num"] <- nums
             }
             node_df2 <- node_df2[,!colnames(node_df2) %in% "path"]
             attr(node_df2, "file") <- files

             df_listC[[jkl]] <- node_df2
             akt_listC[[jkl]] <- akt
           }
           names(df_listC) <- cun
           names(akt_listC) <- cun
           return(list(df_listC, akt_listC))
         }else{
           cat("No partitioning")
           raw_df2 <- raw_df2[,-c(ncol(raw_df2), ncol(raw_df2) - 1, ncol(raw_df2) - 3)]
           colnames(raw_df2) <- c("changes", "event", "node")
           #nodes <- unique(as.vector(tree$edge))
           nodes <- unique(raw_df2$node)
           df_split <- data.frame(nodes = nodes, gains = rep(NA, length(nodes)),
                                  losses = rep(NA, length(nodes)), stringsAsFactors = FALSE)
           df_split_gain <- raw_df2[raw_df2$event %in% "GAIN",]
           df_split_losses <-  raw_df2[raw_df2$event %in% "LOSSES",]
           gains_match <- match(df_split$nodes, df_split_gain$node)
           losses_match <- match(df_split$nodes, df_split_losses$node)
           df_split$gains[!is.na(gains_match)] <- df_split_gain[na.omit(gains_match), "changes"]
           df_split$losses[!is.na(losses_match)] <- df_split_losses[na.omit(losses_match), "changes"]
           node_def <- which(str_detect(raw_file, "0\tnode") == TRUE)
           node_df <- raw_file[node_def]
           node_df2 <- str_split(node_df, "\t")
           node_df2 <- data.frame(matrix(unlist(node_df2), nrow=length(node_df2), byrow=T),stringsAsFactors=FALSE)
           node_df2 <- node_df2[,-1]
           colnames(node_df2) <- c("nodes", "species")
           node_df2$nodes <- str_replace(node_df2$nodes, "node", "")
           node_df2$species <- str_sub(node_df2$species, start = 2) # spacet kitakaritom

           dollo_species <- unique(unlist(sapply(node_df2$species, function(x) str_split(x, " "))))
             if(any(dollo_species %in% tree$tip.label == FALSE)){
               wrong_names_d <- which(dollo_species %in% tree$tip.label == FALSE)
               warning(paste0("Species names don't match!\t", "Species in Dollo file:\t", dollo_species[wrong_names]))
             }else{
               warning(paste0("Species names match!"))
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
           node_df2$gains <- rep(NA, nrow(node_df2))
           node_df2$losses <- rep(NA, nrow(node_df2))
           node_match <- match(node_df2$nodes, df_split$nodes)
           node_df2$gains[!is.na(node_match)] <- df_split[na.omit(node_match), "gains"]
           node_df2$losses[!is.na(node_match)] <- df_split[na.omit(node_match), "losses"]
           node_df2[is.na(node_df2$losses), "losses"] <- 0
           node_df2[is.na(node_df2$gains), "gains"] <- 0
           node_df2$tree_nodes <- as.numeric(node_df2$tree_nodes)
           node_df2$gains <- as.numeric(node_df2$gains)
           node_df2$losses <- as.numeric(node_df2$losses)
           node_df2$net_gains <- node_df2$gains - node_df2$losses
           j <- NULL
           i <- NULL
           root_no <- getMRCA(tree, tree$tip.label)
           node_df2$path <- rep(NA, nrow(node_df2))
             for(j in 1:nrow(node_df2)){
               node <- node_df2[j,"tree_nodes"]
               node_init <- node
               node_path <- vector()
               count <- 1
               while(node != root_no){
                 node_path[count] <- tree$edge[tree$edge[,2] == node,1]
                 node <- node_path[count]
                 count <- count + 1
               }
               node_df2[j, "path"] <- str_c(c(node_init, node_path), collapse = ";")
             }
           j <- NULL
           i <- NULL
           #most meg van a path akkor mindengyiket ki kell bontogatni es osszeadni a dolgokat
           node_df2$copy_num <- rep(NA, nrow(node_df2))
             for(j in 1:nrow(node_df2)){
               path <- node_df2[j, "path"]
               v <- as.numeric(unlist(str_split(path, ";")))
               v_l <- length(v)
               nums <- 0
               for(i in 0:(v_l-1)){
                 nums <- nums + node_df2[node_df2$tree_nodes %in% v[v_l - i], "net_gains"]
               }
               node_df2[j, "copy_num"] <- nums
             }
           node_df2 <- node_df2[,!colnames(node_df2) %in% "path"]
           attr(node_df2, "file") <- files
           return(node_df2)
           }
} # v5
