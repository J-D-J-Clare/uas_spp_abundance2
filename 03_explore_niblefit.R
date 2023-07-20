pkgs <- c("tidyverse", "sf", "ggplot2", "MCMCvis")
sapply(pkgs, require, character.only = TRUE)
source("00_helper_fn.R")


# ======================================================================
# === Adult or Juvenile plants 
# ======================================================================
post <- readRDS("~/../../Volumes/az_drive/detection/models/nf2.2.1nb_posterior.rds")
samp <- map(post, as.data.frame) |> bind_rows() 
dat <- st_read("data/celldat_covar.geojson") |> 
  filter(!is.na(tpi) | !is.na(slope))
# ===
DC <- readRDS("~/../../Volumes/az_drive/detection/models/nf2.2nb_DataConstants.rds")
Data <- DC[['Data']]
Constants <- DC[['Constants']]
# ===

MCMCtrace(post, pdf = FALSE, params = "NTP[1]", ISB = FALSE)
MCMCsummary(post, params = c("beta0FP", "carsd", "sigma", "p", "psi", "beta", "theta1", "theta2"))

phisumm <- MCMCsummary(post, params = "phi")
hist(phisumm$Rhat)
ntpsumm <- MCMCsummary(post, params = "NTP")
hist(ntpsumm$Rhat)
# musumm <- MCMCsummary(post, params = "mu")
# musumm2 <- MCMCsummary(post, params = "mu2")
# 
# 
# idx <- which(is.na(ntpsumm$Rhat))
# ntpsumm[idx,]
# MCMCtrace(post, pdf = FALSE, params = "NTP[1077]", ISB = FALSE)
# idxr <- which(ntpsumm$Rhat > 1.1)
# ntpsumm[idxr,]
# MCMCtrace(post, pdf = FALSE, params = "NTP[530]", ISB = FALSE)
# 
# 
# 
# badx <- Data$ucellidxNV[idx]
# musumm[idx,]
# musumm[nrow(ntpsumm)+idx,]
# Data$n[idx]
# 
# 
# dCompoundNMix(x = 0, N = rnbinom(1, mu = 0.006, size = 0.1), p=.1, lamFP=0.007,
#                  log = 1)
# 
# ll <- sapply(1:1e4, function(x){dCompoundNMix(x = 0, N = rnbinom(1, mu = 0.006, size = 0.1), 
#                                               p=.1, lamFP=0.007, log = 1)})
# table(ll)


# === Visualize parameters
betas <- samp |>
  select(contains('beta')) |> select(1:6)
names(betas) <- c("tri.FP", "tpi.FP", "hli.FP", 
                  "tri.TP", "tpi.TP", "hli.TP")
betas |> 
  pivot_longer(everything()) |> mutate(par = "beta") -> betas
# ---
ps <- samp |>
  select(starts_with('p['))
names(ps) <- unique(dat$site)
ps |> 
  pivot_longer(everything()) |> mutate(par = "p") -> ps
# ---
phis <- samp |>
  select(contains("phi")) |> 
  pivot_longer(everything()) |> group_by(name) |> 
  summarize(value = mean(value)) |> ungroup() |>
  mutate(name = as.numeric(str_extract_all(name, "(?<=\\[).+?(?=\\])") |> unlist())) |> 
  arrange(name) |> mutate(par = "phi")
# ---
psis <- samp |>
  select(contains('psi')) 
names(psis) <- c('Int', 'Wind_max', 'Height_avg', "Area_segmented")
psis |> 
  pivot_longer(everything()) |> mutate(par = "psi") -> psis
# ---
sigma <- data.frame(name = "sigma", value = pull(samp, sigma), par = "sigma")
# ---

betas |> 
  bind_rows(ps, psis, phis |> mutate(name = as.character(name)), sigma) |> 
  group_by(par, name) |> 
  summarize(mean = mean(value), 
            lower = quantile(value, 0.025), 
            upper = quantile(value, 0.975)) |>
  ungroup() |> # write.csv("results/beta_p_psi_phi_sigma_juv.csv", row.names = FALSE)
  filter(par == "beta") |>
  ggplot(aes(y = name, x = mean)) + geom_point() +
  geom_errorbar(aes(x = mean, y = name, xmin = lower, xmax = upper), width = .1) + 
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray") +
  labs(x = "Coefficient", y = "") + 
  facet_grid(.~par) + 
  theme_bw()

# --- theta
thetas <- samp |>
  select(contains('theta')) 
D <- seq(0, 1.5, l = 100)
eff <- sapply(D, function(x) { 
  apply(thetas, 1, function(y) { y[1] * exp(-x * y[2]) }) 
})
eff |> as.data.frame() |> 
  pivot_longer(everything()) |> 
  mutate(name = rep(D, 8000)) |> 
  group_by(name) |> summarize(mean = mean(value), 
                              lower = quantile(value, 0.025), 
                              upper = quantile(value, 0.975)) |>
  ungroup() |> # write.csv("results/distL_juv.csv", row.names = FALSE)
  ggplot(aes(x = name, y = mean)) + geom_line(linewidth = .5, linetype = "dotted") + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) + 
  scale_x_continuous(labels = function(x){x*100}) + 
  labs(x = "Distance from fire edge, m", y = "Effect size") +
  theme_bw()

# --- posterior predictive checks
pp <- predTP(post = samp, data = Data, counts = TRUE, val = TRUE)
ppmu <- apply(pp, 2, mean)

cat("TP: average R2:", cor(ppmu, Data$nTP)^2)
cat("TP: average MAE:", mae(ppmu, Data$nTP)/25)

maeseq <- apply(pp, 1, function(x) {mae(x, Data$nTP)/25 })

pp |> 
  as.data.frame() |> 
  pivot_longer(everything()) |> 
  group_by(name) |>
  summarize(Estimate = mean(value), Q2.5 = quantile(value, 0.025), 
            Q97.5 = quantile(value, 0.975)) |> ungroup() |> 
  arrange(as.numeric(name)) |>
  mutate(obs = Data$nTP) |> # write.csv("results/ppc_juv.csv", row.names = FALSE)
  ggplot(aes(x = obs, y = Estimate)) + geom_point(alpha = 0.25) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0, alpha = .25) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw()


  # --- posterior predictive checks FP
  pp <- predFP(post = samp, data = Data, counts = TRUE, val = TRUE) |> t()
  ppmu <- apply(pp, 2, mean)
  
  cat("FP: average R2:", cor(ppmu, Data$FP)^2)
  cat("FP: average MAE:", mae(ppmu, Data$FP)/25)
  
  maeseq <- apply(pp, 1, function(x) {mae(x, Data$FP)/25 })
  
  pp |> 
    as.data.frame() |> 
    pivot_longer(everything()) |> mutate(name = rep(1:729, 8000)) |>
    group_by(name) |>
    summarize(Estimate = mean(value), Q2.5 = quantile(value, 0.025), 
              Q97.5 = quantile(value, 0.975)) |> ungroup() |> 
    # arrange(as.numeric(name)) |>
    mutate(obs = Data$FP) |> # write.csv("results/ppcFP_juv.csv", row.names = FALSE)
    ggplot(aes(x = obs, y = Estimate)) + geom_point(alpha = 0.25) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0, alpha = .25) + 
    geom_abline(intercept = 0, slope = 1) + 
    theme_bw()
  
  # --- posterior predictive checks TP+FP
  pp <- predUnk(post = samp, data = Data) 
  ppmu <- apply(pp, 2, mean)
  
  cat("Unk: average R2:", cor(ppmu, Data$n)^2)
  cat("Unk: average MAE:", mae(ppmu, Data$n)/25)
  
  maeseq <- apply(pp, 1, function(x) {mae(x, Data$n)/25 })
  
  pp |> 
    as.data.frame() |> 
    pivot_longer(everything()) |> mutate(name = rep(1:length(Data$n), 8000)) |>
    group_by(name) |>
    summarize(Estimate = mean(value), Q2.5 = quantile(value, 0.025), 
              Q97.5 = quantile(value, 0.975)) |> ungroup() |> 
    # arrange(as.numeric(name)) |>
    mutate(obs = Data$n) |> # write.csv("results/ppcUnk_juv.csv", row.names = FALSE)
    ggplot(aes(x = obs, y = Estimate)) + geom_point(alpha = 0.25) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0, alpha = .25) + 
    geom_abline(intercept = 0, slope = 1) + 
    theme_bw()
