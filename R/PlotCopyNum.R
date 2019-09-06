PlotCopyNum <- function(tree, df, stat = c("gains", "losses", "copy_num", "net_gains"), sizemeth = c("standard","comparable","log10sqrt","log10"),
                        size = 5, pch = 16, tip_pch_adj = 0.55, label.offset = 0.2, edge.width = 2, comp.size=c(0,1000), gain_col= "steelblue1", loss_col="red",
                        cladogram = FALSE, pch_col = "blue", node_lab_adj = c(0.5, 0.5), term_adj = -0.7,
                        plot_copy_numbers = FALSE, node_lab_cex = 1, show_species = TRUE, species_size = 1){
  if(length(stat) > 1){
    stop("please choose only one stat value: gains, losses, net_gains or copy_num")
  }
  if(inherits(df, "data.frame")){
    #cat("data_frame_if\n")
    if(length(unique(as.vector(tree$edge))) != nrow(df)){
      stop("not an identical number of nodes and tips in the tree and in the data frame!")
    }
    coln <- which(colnames(df) == stat) #<------------------------------------------------
    if(length(coln) == 0){
      warning(paste0("please provide a column called copy_num\n"))
    }
    if(length(coln) > 1){
      warning(paste0("more than one columns called copy_num\nplease provied only one column with a name copy_num\n"))
    }
    switch(sizemeth,
           standard={ df$nums_trans <- (df[,coln] - min(df[,coln])) / (max(df[,coln]) - min(df[,coln])) },
           comparable={ df$nums_trans <- (abs(df[,coln]) - min(comp.size)) / (max(max(comp.size) - min(comp.size))) },
           log10sqrt={ df$nums_trans <- log10(sqrt(na.omit(df[,coln]))) },
           log10={ df$nums_trans <- log10(na.omit(df[,coln])) },
           sqrtd10={ df$nums_trans <- sqrt(na.omit(df[,coln]))/100 }
    )
    if(plot_copy_numbers == FALSE){
      df[,stat] <- rep(NA, nrow(df))
    }
    df$losses<-df$losses*-1
    df_tips <- df[df$tree_nodes %in% 1:length(tree$tip.label),]
    pch_col_tip<-ifelse(df_tips[,coln]>=0,gain_col,loss_col)
    df_innodes <- df[!df$tree_nodes %in% 1:length(tree$tip.label),]
    pch_col_node<-ifelse(df_innodes[,coln]>=0,gain_col,loss_col)
    coln2 <- which(colnames(df_tips) == stat)
    df_tipsM<-str_pad(df_tips[,coln2], max(nchar(df_tips[,coln2])), side="left", pad=" ")
    if(cladogram){
      tree$edge.length <- NULL
      plot.phylo(tree, label.offset = tip_pch_adj+abs(term_adj)+label.offset, node.depth = 2, show.tip.label = show_species, cex = species_size, edge.width=edge.width)
      nodelabels(pch = pch, col = pch_col_node, node = df_innodes$tree_nodes, cex = df_innodes$nums_trans * size)
      tiplabels(pch = pch, col = pch_col_tip, tip = df_tips$tree_nodes, cex = df_tips$nums_trans * size, adj = c(tip_pch_adj, 0.5))
      nodelabels(text = df_innodes[,stat], node = df_innodes$tree_nodes, adj = node_lab_adj, bg = NULL, frame = "none", cex = node_lab_cex)
      tiplabels(text = df_tipsM, tip = df_tips$tree_nodes, cex = node_lab_cex,
                bg = NULL, frame = "none", adj = term_adj)
    }else{
      plot.phylo(tree, label.offset = tip_pch_adj+abs(term_adj)+label.offset, show.tip.label = show_species, cex = species_size, edge.width=edge.width)
      nodelabels(pch = pch, col = pch_col_node, node = df_innodes$tree_nodes, cex = df_innodes$nums_trans * size)
      tiplabels(pch = pch, col = pch_col_tip, tip = df_tips$tree_nodes, cex = df_tips$nums_trans * size, adj = c(tip_pch_adj, 0.5))
      nodelabels(text = df_innodes[,stat], node = df_innodes$tree_nodes, adj = node_lab_adj, bg = NULL, frame = "none", cex = node_lab_cex)
      tiplabels(text = df_tipsM, tip = df_tips$tree_nodes, cex = node_lab_cex,
                bg = NULL, frame = "none", adj =  term_adj)
    }
  }else{
    if(inherits(df, "list")){
      #cat("list_if\n")
      df_list <- df
      for(i in 1:length(df_list)){
        df <- df_list[[i]]
        if(length(unique(as.vector(tree$edge))) != nrow(df)){
          stop("not an identical number of nodes and tips in the tree and in the data frame!")
        }
        coln <- which(colnames(df) == stat)
        if(length(coln) == 0){
          warning(paste0("please provide a column called copy_num\n"))
        }
        if(length(coln) > 1){
          warning(paste0("more than one columns called copy_num\nplease provied only one column with a name copy_num\n"))
        }
        switch(sizemeth,
               standard={ df$nums_trans <- (df[,coln] - min(df[,coln])) / (max(df[,coln] - min(df[,coln]))) },
               comparable={ df$nums_trans <- (abs(df[,coln]) - min(comp.size)) / (max(max(comp.size) - min(comp.size))) },
               log10sqrt={ df$nums_trans <- log10(sqrt(na.omit(df[,coln]))) },
               log10={ df$nums_trans <- log10(na.omit(df[,coln])) },
               sqrtd10={ df$nums_trans <- sqrt(na.omit(df[,coln]))/100 }
        )
        df$losses<-df$losses*-1
        df_tips <- df[df$tree_nodes %in% 1:length(tree$tip.label),]
        pch_col_tip<-ifelse(df_tips[,coln]>=0,gain_col,loss_col)
        df_innodes <- df[!df$tree_nodes %in% 1:length(tree$tip.label),]
        pch_col_node<-ifelse(df_innodes[,coln]>=0,gain_col,loss_col)
        coln2 <- which(colnames(df_tips) == stat)
        df_tipsM<-str_pad(df_tips[,coln2], max(nchar(df_tips[,coln2])), side="left", pad=" ")
        if(cladogram){
          tree$edge.length <- NULL
          plot.phylo(tree, label.offset = tip_pch_adj+abs(term_adj)+label.offset, node.depth = 2, show.tip.label = show_species, cex = species_size, edge.width=edge.width)
          nodelabels(pch = pch, col = pch_col_node, node = df_innodes$tree_nodes, cex = df_innodes$nums_trans * size)
          tiplabels(pch = pch, col = pch_col_tip, tip = df_tips$tree_nodes, cex = df_tips$nums_trans * size, adj = c(tip_pch_adj, 0.5))
          nodelabels(text = df_innodes[,stat], node = df_innodes$tree_nodes, adj = node_lab_adj, bg = NULL, frame = "none", cex = node_lab_cex)
          tiplabels(text = df_tipsM, tip = df_tips$tree_nodes, cex = node_lab_cex,
                    bg = NULL, frame = "none", adj =  term_adj)
        }else{
          plot.phylo(tree, label.offset = tip_pch_adj+abs(term_adj)+label.offset, show.tip.label = show_species, cex = species_size, edge.width=edge.width)
          nodelabels(pch = pch, col = pch_col_node, node = df_innodes$tree_nodes, cex = df_innodes$nums_trans * size)
          tiplabels(pch = pch, col = pch_col_tip, tip = df_tips$tree_nodes, cex = df_tips$nums_trans * size, adj = c(tip_pch_adj, 0.5))
          nodelabels(text = df_innodes[,stat], node = df_innodes$tree_nodes, adj = node_lab_adj, bg = NULL,
                     frame = "none", cex = node_lab_cex)
          tiplabels(text = df_tips$copy_num, tip = df_tips$tree_nodes, cex = node_lab_cex,
                    bg = NULL, frame = "none", adj =  term_adj)
        }
      }
    }else{
      warning("the input is neither a data.frame nor a list!\n please proved a list or a data.frame object!\n")
    }
  }
}

