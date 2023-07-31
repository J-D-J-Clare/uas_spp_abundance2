library(MCMCvis)

models <- list.files("detection/models/power_analysis/", pattern = "power_0", full.names = TRUE)

lout <- list()

for(i in 1:length(models)){
	f <- models[[i]]
	r <- readRDS(f)[[1]]
	ll <- list()
	
	ll[[1]] <- MCMCsummary(r, params = c("beta0FP", "carsd", "sigma", "p", "psi", "beta", "theta1", "theta2"))
	ll[[2]] <- MCMCsummary(r, params = "phi")
	ll[[3]] <- MCMCsummary(r, params = "NTP")
	lout[[i]] <- ll

}

write.csv(as.data.frame(models), "detection/power_model_names.csv")
saveRDS(lout, "detection/power_results.rds")

