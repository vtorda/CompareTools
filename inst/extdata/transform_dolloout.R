
files <- "~/Documents/MUNKA/PhD/CompareTools/inst/extdata/Dolloout_whole_genome"
raw_file <- read_lines(files)
raw_file2 <- str_replace(raw_file, "^\\d+\t", "")
sink(file = "Dolloout_all_cluster")
cat(raw_file2, sep = "\n")
sink()
