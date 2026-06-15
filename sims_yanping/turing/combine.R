setwd("/home/ypei/obsReLOOP/results")
files <- list.files("yhatRF", pattern = "^sim_results_task_.*\\.rds$", full.names = TRUE)

sim_results <- do.call(rbind, lapply(files, readRDS))
sim_results <- sim_results[order(sim_results$sim, sim_results$scenario), ]

saveRDS(sim_results, file = "yhatRF/yhatRF_sim_results_all.rds")
write.csv(sim_results, file = "yhatRF/yhatRF_sim_results_all.csv", row.names = FALSE)
