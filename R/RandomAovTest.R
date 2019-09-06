#' Perform random oneway anova test
#'
#' @param compare a compare object
#' @param dataset
#' @param stat
#' @param aov_groups
#' @param s
#'
#' @return this function returns a compare object with grouped clusters
#' @export
#' @import ape readr plyr stringr tidyr
#'
#' @examples
######
## Nem kene p adjustolni az aov eredmenyeket?
######
# compare <- comp_obj
# dataset <- "event_data_GROUPED"
# stat <- "gains"
RandomAovTest <- function(compare, dataset = NULL, stat = NULL, aov_groups, s) {
  aov_groups <- aov_groups[order(aov_groups)]
  data <- as.matrix(compare[[dataset]][[stat]][aov_groups,])
  fact <- names(aov_groups)
  output <- NULL
  for (i in 1:ncol(data)) {
    # Single Factor Randomization ANOVA.
    # Output is a vector with DF treatment, DF error, MS error, F statistic, and p value.
    null_sum <- summary(aov(data[,i] ~ fact))[[1]]
    Fstat <- null_sum$F[1]
    MSe <- null_sum$Mean[2]
    DFt <- null_sum$Df[1]
    DFe <- null_sum$Df[2]
    perm_resp <- rep(0,length(data[,i])*s)
    dim(perm_resp) <- c(length(data[,i]),s)
    for (j in 1:s) {
      perm_resp[,j] <- sample(data[,i],length(data[,i]), replace=FALSE)}

    Fnull <- rep(0,s)
    for (j in 1:s) {
      temp <- summary(aov(perm_resp[,j] ~ fact))
      Fnull[j] <- temp[[1]]$F[1]}

    pval <- sum(Fnull > Fstat) / length(Fnull)

    ranova <- c(DFt,DFe,MSe,Fstat,pval)
    names(ranova) <- c("DFt","DFe","MSe","F","p")

    output <- rbind(output,ranova)
  }
  rownames(output) <- colnames(data)
  return(output)
}
