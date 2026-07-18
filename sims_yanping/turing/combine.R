library(dplyr)
library(tidyr)

setwd("/home/ypei/obsReLOOP")

# read and combine task results
files <- list.files(
  "results",
  pattern = "^sim_results_task_.*\\.rds$",
  full.names = TRUE
)

print(files)

sim_results <- do.call(rbind, lapply(files, readRDS))

sim_results <- sim_results[order(
  sim_results$sim,
  sim_results$scenario,
  sim_results$functional_form,
  sim_results$matching_quality
), ]

saveRDS(sim_results, file = "results/sim_results_all500y0hat.rds")
write.csv(sim_results, file = "results/sim_results_all500y0hat.csv", row.names = FALSE)
