#' Slice compare object
#'
#' @param compare a compare object
#' @param clusters a character vetor. With this vector you input cluster names you want to keep in the data.
#' @param data a character or a numeric value or a vector. With this parameter the required data can be sliced.
#' Default is the event_data, which is the original dataset imported in R.
#' @param new_name_suf suffix to add to the original data name(s)
#' @param overwrite should the original data be replaced with the sliced version?
#'
#' @return this function returns a compare object retain only just part of the original clusters
#' @export
#' @import ape readr plyr stringr tidyr
#'
#' @examples

SliceComp <- function(compare, clusters = NULL, data = "event_data",
                      new_name_suf = "SLICED", overwrite = TRUE){
  if(is.character(data)){
    data <- which(names(compare) %in% data)
  }
  if(!is.numeric(data) & !is.character(data)){
    stop(paste0("Data parameter should have a character or numeric value!\n"))
  }
  if(!is.character(clusters)){
    stop(paste0("cluster parameter should be a character vector!\n"))
  }
  data_obj <- compare[data]
  old_names <- names(compare)
  for(i in 1:length(data_obj)){
    data_obj[[i]] <- lapply(data_obj[[i]], function(x) x[,colnames(x) %in% clusters])
  }
  if(overwrite){
    compare[data] <- data_obj
    new_names <- old_names
    new_names[data] <- str_c(old_names[data], new_name_suf, sep = "_")
    names(compare) <- new_names
    return(compare)
  }else{
    compare <- c(compare[-length(compare)], data_obj,
                     compare["raw_data"])
    data_names <- str_c(old_names[data], new_name_suf, sep = "_")
    names(compare) <- c(old_names[-length(old_names)],
                            data_names, "raw_data")
    return(compare)
  }
}

