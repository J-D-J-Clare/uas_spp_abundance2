library(tidyverse)

r <- readRDS("~/../../Volumes/az_drive/detection/power_results.rds")
names <- read.csv("~/../../Volumes/az_drive/detection/power_model_names.csv") |> 
  mutate(X = rep( c(seq(10, 80, by = 10), 100, 90), 2), 
         model = sort(rep(c("Adult", "Juv"), 10)) , 
         idx = paste0(model, "_", X))

lpar <- lapply(1:length(r), function(x) { r[[x]][[1]] |> mutate(prop = names$X[x], model = names$model[x], par = rownames(r[[x]][[1]])) }) |> bind_rows()
lntp <- lapply(1:length(r), function(x) { r[[x]][[3]] |> mutate(prop = names$X[x], model = names$model[x]) }) |> bind_rows()
lntp <- lapply(1:length(r), function(x) { r[[x]][[3]] |> mutate(prop = names$X[x], model = names$model[x]) }) |> bind_rows()


lntp |> 
  filter(prop > 10) |>
  group_by(model, prop) |> summarize(total = sum(mean), 
                               lower = sum(`2.5%`), 
                               upper = sum(`97.5%`)) |> ungroup() |>
  ggplot() + geom_point(aes(x = prop, y = total)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = prop), width = .25) + 
  theme_bw() + facet_wrap(.~model, scales = "free")

lpar |> 
  filter(prop > 10) |>
  filter(grepl("p", par) & !grepl("psi", par) ) |> 
  ggplot(aes(x = prop, y = mean, group = par, fill = par)) + geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = par, group = par), alpha = .5, linetype = "dotted", linewidth = .25) + 
  #geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, x = par), width = .25, position = position_dodge(width = 3)) + 
  facet_wrap(.~ model)  + 
  theme_bw() 

lpar |> 
  filter(prop > 10) |>
  filter(grepl("psi", par) ) |> 
  ggplot(aes(x = prop, y = mean, group = par, fill = par)) + geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = par, group = par), alpha = .5, linetype = "dotted", linewidth = .25) + 
  #geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, x = par), width = .25, position = position_dodge(width = 3)) + 
  facet_wrap(.~ model)  + 
  theme_bw() 

lpar |> 
  filter(grepl("theta", par) ) |> 
  ggplot(aes(x = prop, y = mean, group = par, fill = par)) + geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = par, group = par), alpha = .5, linetype = "dotted", linewidth = .25) + 
  #geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, x = par), width = .25, position = position_dodge(width = 3)) + 
  facet_wrap(.~ model)  + 
  theme_bw() 

lpar |> 
  filter(grepl("beta", par) ) |> 
  ggplot(aes(x = prop, y = mean, group = par, fill = par, group = par)) + geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = par, group = par), alpha = .5, linetype = "dotted", linewidth = .25) + 
  #geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, x = par), width = .25, position = position_dodge(width = 3)) + 
  facet_wrap(.~ model)  + 
  theme_bw() 


lntp |> 
  filter(prop > 10) |>
  mutate(base = rep(filter(lntp, prop == 100) |> pull(mean), 10))
  group_by(model, prop) |> summarize(total = sum(mean), 
                                     lower = sum(`2.5%`), 
                                     upper = sum(`97.5%`)) |> ungroup() |>
  ggplot() + geom_point(aes(x = prop, y = total)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = prop), width = .25) + 
  theme_bw() + facet_wrap(.~model, scales = "free")

