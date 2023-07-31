library(tidyverse)

r <- readRDS("~/../../Volumes/az_drive/detection/power_results.rds")
names <- read.csv("~/../../Volumes/az_drive/detection/power_model_names.csv") |> 
  mutate(X = rep( c(seq(10, 80, by = 10), 100, 90), 2), 
         model = sort(rep(c("Adult", "Juv"), 10)) , 
         idx = paste0(model, "_", X))

lpar <- lapply(1:length(r), function(x) { r[[x]][[1]] |> mutate(prop = names$X[x], model = names$model[x], par = rownames(r[[x]][[1]])) }) |> bind_rows()
lntp <- lapply(1:length(r), function(x) { r[[x]][[3]] |> mutate(prop = names$X[x], model = names$model[x]) }) |> bind_rows()

# === total count 
lntp |> 
  filter(prop > 10) |>
  filter(prop > 10) |>
  group_by(model, prop) |> summarize(total = sum(mean), 
                               lower = sum(`2.5%`), 
                               upper = sum(`97.5%`)) |> ungroup() |>
  ggplot() + geom_point(aes(x = prop, y = total)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = prop), width = .25) + 
  theme_bw() + facet_wrap(.~model, scales = "free")

# === Probability of detection: site-level
lpar |> 
  filter(prop > 10) |>
  filter(grepl("p", par) & !grepl("psi", par) ) |> 
  ggplot(aes(x = prop, y = mean, group = par, fill = par)) + geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = par, group = par), alpha = .5, linetype = "dotted", linewidth = .25) + 
  labs(x = "Propotion of validated cells", y = "P(detection)", fill = "Site") +
  facet_wrap(.~ model)  + 
  theme_bw() 

# === Probability of detection, covariates
lpar |> 
  filter(prop > 10) |>
  filter(grepl("psi", par) ) |> 
  ggplot(aes(x = prop, y = mean, group = par, fill = par)) + geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = par, group = par), alpha = .5, linetype = "dotted", linewidth = .25) + 
  labs(x = "Propotion of validated cells", y = "Coefficient", fill = "Site") +
  scale_fill_manual(values = 1:4, labels = c('Intercept','Wind','Height', "AreaSeg")) +
  facet_wrap(.~ model)  + 
  theme_bw() 

# === spatial term coefficients
lpar |> 
  filter(prop > 10) |>
  filter(grepl("theta", par) ) |> 
  ggplot(aes(x = prop, y = mean, group = par, fill = par)) + geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = par, group = par), alpha = .5, linetype = "dotted", linewidth = .25) + 
  labs(x = "Propotion of validated cells", y = "Coefficient", fill = "Site") +
  facet_wrap(.~ model)  + 
  theme_bw() 

# === Beta Coefficients
lpar |> 
  filter(prop > 10) |>
  filter(grepl("beta", par) ) |> 
  ggplot(aes(x = prop, y = mean, group = par, fill = par, group = par)) + geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = par, group = par), alpha = .5, linetype = "dotted", linewidth = .25) + 
  scale_fill_manual(values = 1:7, 
                    labels = c('FP_TRI','TP_TRI','FP_TPI', 'TP_TPI', 'FP_HLI', 'TP_HLI', "FP_Intercept")) +
  facet_wrap(.~ model)  + 
  theme_bw() 

# === Cell-level abundances (mae or rmse)
# 1. Comparable are only the cells that are included in the full model [rest need to be subset]
dats <- lapply(list.files("~/../../Volumes/az_drive/detection/power_analysis/", pattern = "power_data", full.names = TRUE), readRDS)
refnames <- filter(dats[[10]][[2]], val == 0) |> pull(ucid)

ll <- list()
for(i in 1:20){
  if (i<=10) {mod = "Adult"} else {mod = "Juv"}
  j <- rep(1:10, 2)[i]
  dats[[i]][[2]] |> filter(val == 0) |> pull(ucid) -> id
  
  lntp |> filter(model == !!mod, prop == j*10) |> 
    mutate(ucid = id) -> tmp
  ll[[i]] <- tmp
}
ll |> bind_rows() |> filter(ucid %in% refnames) -> dfntp

dfntp|>
  mutate(base = c(rep(filter(dfntp, model == 'Adult', prop == 100) |> pull(mean), 10),
                  rep(filter(dfntp, model == 'Juv', prop == 100) |> pull(mean), 10)), 
         mae = ((mean - base)) ) |> #-> tmp
  filter(prop > 10) |>
  ggplot() + geom_jitter(aes(x = prop, y = mae), alpha = .2) +
  # geom_errorbar(aes(ymin = lower, ymax = upper, x = prop), width = .25) + 
  theme_bw() + facet_wrap(.~model, scales = "free")

