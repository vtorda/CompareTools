library(CompareTools)
library(ape)
library(microbenchmark)
file_path <- "~/Documents/MUNKA/PhD/CompareTools/inst/extdata/Dollogeneratum/"
treefiles <- list.files(file_path, pattern = "^Tree")
treefiles <- treefiles[-4]
dollofiles <- list.files("~/Documents/MUNKA/PhD/CompareTools/inst/extdata/Dollogeneratum/", pattern = "^dollo")
reorder_index <- c(12,13,14,15,8,9,10,11,4,5,6,7,1,2,3)
treefiles <- treefiles[reorder_index]
dollofiles <- dollofiles[reorder_index]
for(i in 1:15){
  cat(i, "\t", dollofiles[i], "\n")
  tree <- read.tree(paste0(file_path, treefiles[i]))
  file <- paste0(file_path, dollofiles[i])
  start_memory <- gc(reset = TRUE)
  micro_log <- microbenchmark(DolloImport3(tree = tree, files = file), times = 1)
  end_memory <- gc()
  save(list = c("start_memory", "end_memory", "micro_log"), file = paste0(file_path, "benchmark_", dollofiles[i], ".RData"))
}
