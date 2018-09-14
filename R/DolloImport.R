DolloOutImport <- function(tree, files = NULL, path = NULL, pattern = NULL, partitioning = "NO", groups = NULL){
  if(!is.null(path) & is.null(pattern)){
    if(is.null(files)){
      warning(paste0("You need to give at least one file name or a pattern!"))
    }
    files <- paste0(path, files)
  }
  if(!is.null(pattern) & is.null(path)){
    files <- list.files("./", pattern = pattern)
  }
  if(!is.null(pattern) & !is.null(path)){
    files <- list.files(path = path, pattern = pattern)
  }
  # itt most kette valik hogy tobb filet hivok vagy nem...
  if(length(files) == 1){
    raw_file <- read_lines(files)
    rows_needed <- c(which(str_detect(raw_file, "GAIN") == TRUE), which(str_detect(raw_file, "LOSSES") == TRUE))
    raw_df <- raw_file[rows_needed]
    raw_df2 <- str_split(raw_df, "\t")
    raw_df2 <- data.frame(matrix(unlist(raw_df2), nrow = length(raw_df2), byrow = TRUE), stringsAsFactors = FALSE)
    ######################Zsolti
    switch(partitioning,
           YES={
             cat("Partitioning: ")
             raw_df2z <- raw_df2[,-c(ncol(raw_df2) - 1, ncol(raw_df2) - 3)]
             colnames(raw_df2z) <- c("changes", "event", "node", "chars")
             cutter<-separate_rows(raw_df2z,chars, sep = " ")
             cutter<-separate(cutter,chars, sep = "/", into=c("clusters","orthogroups"))
             cun<-sort(unique(cutter$clusters))
             if(is.null(groups)){
               allcl<-lapply(cun,function(x) cutter[which(cutter$clusters==x),])
             }else{
               csop<-unique(groups$csop)
               csopcl<-lapply(csop,function(x) cutter[which(cutter$clusters %in% groups[groups[,2]==x,1][[1]]),])
               cun<-csop
               allcl<-csopcl
             }
             df_listC<-list()
             akt_listC<-list()
             for (jkl in 1:length(allcl)){
               #jkl=1
               cat(jkl," ")
               akt<-allcl[[jkl]]
               glcounts<-table(akt$node,akt$event)
               if(isTRUE(grep("GAIN",colnames(glcounts))>0)){
                 akt[akt$event=="GAIN","changes"]<-glcounts[match(akt[akt$event=="GAIN","node"], rownames(glcounts)),1] # we modify the number of changes (gain) in the separated cluster table
               }
               if(isTRUE(grep("LOSSES",colnames(glcounts))>0)){
                 akt[akt$event=="LOSSES","changes"]<-glcounts[match(akt[akt$event=="LOSSES","node"], rownames(glcounts)),2] # we modify the number of changes (loss) in the separated cluster table
               }
               raw_df2<-akt[,c(1:3)]
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
               node_dfX<-plyr::ldply(node_df2,rbind)
               #node_df2 <- data.frame(matrix(unlist(node_df2), nrow=length(node_df2), byrow=T),stringsAsFactors=FALSE)
               node_df2 <- node_dfX[,-1]
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
             names(df_listC)<-cun
             names(akt_listC)<-cun
             return(list(df_listC, akt_listC))
           }, #------------------------------------------------------------------------------------regiscript
           NO={cat("No partitioning")
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
    )
    ###############################################################################
    if(length(files) > 1){
      df_list <- list()
      for(j in 1:length(files)){
        warning(paste0("process file:\t", files[j], "\n"))
        raw_file <- read_lines(files[j])
        rows_needed <- c(which(str_detect(raw_file, "GAIN") == TRUE), which(str_detect(raw_file, "LOSSES") == TRUE))
        raw_df <- raw_file[rows_needed]
        raw_df2 <- str_split(raw_df, "\t")
        raw_df2 <- data.frame(matrix(unlist(raw_df2), nrow = length(raw_df2), byrow = TRUE), stringsAsFactors = FALSE)
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
          warning(paste0("Species names match!\n"))
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
        i <- NULL
        root_no <- getMRCA(tree, tree$tip.label)
        node_df2$path <- rep(NA, nrow(node_df2))
        for(k in 1:nrow(node_df2)){
          node <- node_df2[k,"tree_nodes"]
          node_init <- node
          node_path <- vector()
          count <- 1
          while(node != root_no){
            node_path[count] <- tree$edge[tree$edge[,2] == node,1]
            node <- node_path[count]
            count <- count + 1
          }
          node_df2[k, "path"] <- str_c(c(node_init, node_path), collapse = ";")
        }
        i <- NULL
        k <- NULL
        #most meg van a path akkor mindengyiket ki kell bontogatni es osszeadni a dolgokat
        node_df2$copy_num <- rep(NA, nrow(node_df2))
        for(k in 1:nrow(node_df2)){
          path <- node_df2[k, "path"]
          v <- as.numeric(unlist(str_split(path, ";")))
          v_l <- length(v)
          nums <- 0
          for(i in 0:(v_l-1)){
            nums <- nums + node_df2[node_df2$tree_nodes %in% v[v_l - i], "net_gains"]
          }
          node_df2[k, "copy_num"] <- nums
        }
        node_df2 <- node_df2[,!colnames(node_df2) %in% "path"]

        df_list[[j]] <- node_df2
      }
      names(df_list) <- files
      return(df_list)
    }
  }
} # v5
