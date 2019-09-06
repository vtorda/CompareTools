load("~/Documents/MUNKA/PhD/CompareTools/benchmark/comp_obj.RData")
library(CompareTools)
A <- GetNodes(comp_obj, c("Rozella allomycis", "Laccaria bicolor"), tips = FALSE)
all_n <- c(1:141)
all_n2 <- all_n[!all_n %in% A]
all_n3 <- all_n2[!all_n2 %in% 1:71]
all_n3 <- all_n3[-1]
aov_groups <- c(A, all_n3)
names(aov_groups) <- c(rep("Fungi", length(A)), rep("Animals_Protists", length(all_n3)))
all_cl <- RandomAovTest(comp_obj, dataset = "event_data", "gains", aov_groups = aov_groups, 100)
save(all_cl, file = "all_aov_test.RData")
