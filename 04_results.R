pkgs <- c("tidyverse", "sf", "terra", 'stars', 'ggspatial')
sapply(pkgs, require, character.only = TRUE)
source("00_helper_fn.R")
# ===========================


geo <- st_read("data/celldat_covar.geojson")
# ===========================
dat <- geo |> select(ucid) |> 
  right_join(read.csv("data/celldat_input.csv"))
# ======================
dat025 <- geo |> select(ucid) |> 
  right_join(read.csv("data/celldat_input025.csv"))
# ======================
datsite <- read.csv("data/site_level_data.csv")
# ======================
# ---

# === p: probability of detection
post <- readRDS("~/../../Volumes/az_drive/detection/models/nf1.2.1nb_pars.rds")
post |> 
  select(starts_with('p[')) |> rename_with(.fn = function(x){unique(dat$site)}) |> 
  pivot_longer(cols = everything()) |> 
  group_by(name) |> 
  summarize(mean = mean(value), l95 = quantile(value, 0.025), u95 = quantile(value, 0.975)) |> 
  ungroup() # -> ps

# === psi: probability of detection - covariates
# --- summary reports dP(d) when covariate is increased by 2 sd
post <- readRDS("~/../../Volumes/az_drive/detection/models/nf2.2.1nb_pars.rds")
post |> 
  select(contains('psi')) |> rename_with(.fn = function(x){c('Intercept', 'Wind_max', 'Height_avg', 'Area_segment')}) |> 
  mutate(dWind_max = plogis(Intercept + Wind_max) - plogis(Intercept),
         dHeight_avg = plogis(Intercept + Height_avg) - plogis(Intercept),
         dArea_segment = plogis(Intercept + Area_segment) - plogis(Intercept)) |>
  pivot_longer(cols = everything()) |> 
  group_by(name) |> 
  summarize(mean = mean(value), l95 = quantile(value, 0.025), u95 = quantile(value, 0.975)) |> 
  ungroup()

aInt <- 0.561
jInt <- -2.24
cat("Wind: Relative reduction in p(adult): ", 1 - (-0.186  + aInt)/aInt )
cat("Wind: Relative reduction in p(juv): ", 1 - (0.0522 + jInt )/jInt )

cat("Height: Relative reduction in p(adult): ", plogis(0.689 + aInt)/plogis(aInt))
cat("Height: Relative reduction in p(juv): ", (1.89 + jInt )/jInt)

cat("PArea: Relative reduction in p(adult): ", 1 - (0.599 + aInt)/aInt)
cat("PArea: Relative reduction in p(juv): ", 1 - (0.165 + jInt )/jInt)

# === FP ppc
post <- readRDS("~/../../Volumes/az_drive/detection/models/nf2.2.1nb_pars.rds")
DC <- readRDS("~/../../Volumes/az_drive/detection/models/nf2.2nb_DataConstants.rds")
Data <- DC[['Data']]
Constants <- DC[['Constants']]

ppcfp <- t(predFP(post = post, data = Data, counts = TRUE, val = FALSE))
colnames(ppcfp) <- 1:ncol(ppcfp)

ppcfpmu <- apply(ppcfp, 2, mean)

cat("FP: mean count: ", mean(ppcfpmu/25) )
cat("FP: count SD: ", sd(ppcfpmu/25) )
cat("FP: count 95%CI: ", quantile(ppcfpmu, probs = c(0.025, 0.975))/25)

cat("FP: average R2", cor(ppcfpmu[Data$ucellidxV], Data$FP)^2)
cat("FP: average MAE", mae(ppcfpmu[Data$ucellidxV], Data$FP)/25)

maeseq <- apply(ppcfp[,Data$ucellidxV], 1, function(x) { mae(x, Data$FP) })/25
cat("FP: mean absolute error", mean(maeseq))

ppcfp <- t(predFP(post = post, data = Data, counts = TRUE, val = TRUE))
colnames(ppcfp) <- 1:ncol(ppcfp)

ppcfp |> as.data.frame() |> 
  pivot_longer(cols = everything()) |> 
  group_by(name) |> summarize(mean = mean(value), l95 = quantile(value, .025), u95 = quantile(value, 0.975)) |>
  ungroup() |> mutate(name = as.numeric(name)) |> arrange(name) |> 
  mutate(site = Data$siteV) |> 
  group_by(site) |> summarize(smean = mean(mean/25), sd = sd(mean/25)) |> 
  ungroup() # -> fp

# === TP~FP comparison
mu <- readRDS("~/../../Volumes/az_drive/detection/models/nf1.2.1nb_mu.rds")

mufp <- apply(mu[, 1:nrow(dat)], 2, mean)
mutp <- apply(mu[, (nrow(dat)+1):ncol(mu)], 2, mean)
dat |> select(site) |> 
  mutate(mufp = mufp, mutp = mutp, p = ps$mean[Data$sitefull]) |> 
  group_by(site) |> 
  summarize(ntpfp = mean(mufp)/mean(mutp), 
                              tpfp = mean(mufp)/mean(mutp*p)) |> 
  ungroup() -> tpfp
  
# === FP betas
post <- readRDS("~/../../Volumes/az_drive/detection/models/nf2.2nb_pars.rds")
post |> 
  select(contains('beta') & contains("1]"), "beta0FP") |> rename_with(.fn = function(x){c('TRI', 'TPI', 'HLI', 'Intercept')}) |> 
  mutate(dTRI = exp(Intercept + TRI) - exp(Intercept),
         dTPI = exp(Intercept + TPI) - exp(Intercept),
         dHLI = exp(Intercept + HLI) - exp(Intercept)) |>
  pivot_longer(cols = everything()) |> 
  group_by(name) |> 
  summarize(mean = mean(value), l95 = quantile(value, 0.025), u95 = quantile(value, 0.975)) |> 
  ungroup()

aInt <- exp(-1.09)
jInt <- exp(0.0383)

cat("TRI: Change in expected abundance, adult: ", (1.71 + aInt)/aInt - 1 )
cat("TRI: Change in expected abundance, juv: ", (1.25 + jInt)/jInt -1 )

cat("TPI: Change in expected abundance, adult: ", (0.118 + aInt)/aInt - 1)
cat("TPI: Change in expected abundance, juv: ", (0.172 + jInt )/jInt - 1)

cat("HLI: Change in expected abundance, adult: ", (-0.833 + aInt)/aInt - 1)
cat("HLI: Change in expected abundance, juv: ", (-0.481 + jInt)/jInt - 1)


# === nTP ppc
post <- readRDS("~/../../Volumes/az_drive/detection/models/nf2.2.1nb_pars.rds")
DC <- readRDS("~/../../Volumes/az_drive/detection/models/nf2.2nb_DataConstants.rds")
Data <- DC[['Data']]
Constants <- DC[['Constants']]

ppctp <- predTP(post = post, data = Data, counts = TRUE, val = FALSE)

ppctpmu <- apply(ppctp, 2, mean)

cat("NTP: mean count: ", mean(ppctpmu/25) )
cat("NTP: count SD: ", sd(ppctpmu/25) )
cat("NTP: count 95%CI: ", quantile(ppctpmu, probs = c(0.025, 0.975))/25)

cat("NTP: average R2", cor(ppctpmu[Data$ucellidxV], Data$nTP)^2)
cat("NTP: average MAE", mae(ppctpmu[Data$ucellidxV], Data$nTP)/25)

maeseq <- apply(ppctp[,Data$ucellidxV], 1, function(x) { mae(x, Data$nTP) })/25
cat("NTP: mean absolute error", mean(maeseq))

ppctp <- predTP(post = post, data = Data, counts = TRUE, val = TRUE)
colnames(ppctp) <- 1:ncol(ppctp)

ppctp |> as.data.frame() |> 
  pivot_longer(cols = everything()) |> 
  group_by(name) |> summarize(mean = mean(value), l95 = quantile(value, .025), u95 = quantile(value, 0.975)) |>
  ungroup() |> mutate(name = as.numeric(name)) |> arrange(name) |> 
  mutate(site = Data$sitefull) |> 
  group_by(site) |> summarize(smean = mean(mean/25), sd = sd(mean/25)) |> 
  ungroup()

# === TP betas
post <- readRDS("~/../../Volumes/az_drive/detection/models/nf2.2.1nb_pars.rds")
remu <- apply(post |> select(contains("phi")) |> as.matrix(), 1, mean)

post |> 
  select(contains('beta') & contains("2]"),) |> rename_with(.fn = function(x){c('TRI', 'TPI', 'HLI')}) |> 
  mutate(Intercept = remu) |>
  mutate(dTRI = exp(Intercept + TRI) - exp(Intercept),
         dTPI = exp(Intercept + TPI) - exp(Intercept),
         dHLI = exp(Intercept + HLI) - exp(Intercept)) |>
  pivot_longer(cols = everything()) |> 
  group_by(name) |> 
  summarize(mean = mean(value), l95 = quantile(value, 0.025), u95 = quantile(value, 0.975)) |> 
  ungroup()

aInt <- exp(-1.06 )
jInt <- exp(-0.0552)

cat("TRI: Change in expected abundance, adult: ", (0.473 + aInt)/aInt - 1 )
cat("TRI: Change in expected abundance, juv: (overlaps zero)", (0.520  + jInt)/jInt -1 )

cat("TPI: Change in expected abundance, adult: ", (-0.0311 + aInt)/aInt - 1)
cat("TPI: Change in expected abundance, juv: ", (0.155 + jInt )/jInt - 1)

cat("HLI: Change in expected abundance, adult: ", (0.178 + aInt)/aInt - 1)
cat("HLI: Change in expected abundance, juv: (ovarlaps zero)", (0.0225 + jInt)/jInt - 1)


# === TP edge effect
ll <- list()
for(i in 1:2) {
  post <- readRDS(paste0("~/../../Volumes/az_drive/detection/models/nf", i,".2.1nb_pars.rds"))
  thetas <- post |> select(contains('theta')) 
  sigma <- post |> select('sigma') |> pull(sigma)
  Int <- c(exp(-0.0886), exp(2.01))[i]
  
  D <- seq(0, 1.5, by = .05)
  eff <- sapply(D, function(x) { apply(thetas, 1, function(y) {  y[1] * exp(-x * y[2]) }) })
  
  ll[[i]] <- sapply(1:nrow(thetas), function(i) { rnbinom(length(D), mu =  0.25*(exp(log(Int) + eff[i,]) - Int), size = sigma[i]) }) |> t()
}

map(ll, as.data.frame) |> bind_rows() |>
  pivot_longer(everything()) |>  
  mutate(name = rep(D, 2*8000), 
         Model = sort(rep(c("Adult", "Juvenile"), 31*8000))) |> 
  group_by(Model, name) |> summarize(mean = mean(value), 
                              lower = quantile(value, 0.025), 
                              upper = quantile(value, 0.975)) |>
  ungroup() |> filter(name %in% c(.1, .5, 1)) |>
  mutate(Model = factor(Model, levels = c("Juvenile", "Adult"))) #|>
  ggplot(aes(x = name, y = mean, group = Model, colour = Model, fill = Model)) + 
  geom_line(linewidth = 1.5, linetype = "solid") + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25, linetype = "dotted", linewidth = .25) + 
  scale_x_continuous(labels = function(x){x*100}) + 
  geom_vline(xintercept = c(.1, .5, 1), alpha = .33, linetype = "dashed") +
  scale_fill_viridis_d(begin = .2, end = 0.8, direction = -1) +
  scale_colour_viridis_d(begin = .2, end = 0.8, direction = -1) +
  labs(x = "Distance from fire edge, m", 
       y = expression(paste("Effect size, [", plants~m^{-2}, "]")) ) +
  theme_bw()

ggsave("figures/dist_ppc_marginal.png", width = 12, height = 8, units = "cm")

# === Predicted counts
ll <- list()
for(i in 1:2){
  class <- c("adult", "juvenile")[i]
  
  NTP <- readRDS(paste0("~/../../Volumes/az_drive/detection/models/nf", i,".2.1nb_NTP.rds"))
  DC <- readRDS(paste0("~/../../Volumes/az_drive/detection/models/nf", i, ".2nb_DataConstants.rds"))
  Data <- DC[['Data']]
  Constants <- DC[['Constants']]
  
  NTPmu <- apply(NTP, 2, mean)
  NTPlower <- apply(NTP, 2, quantile, probs = 0.025)
  NTPupper <- apply(NTP, 2, quantile, probs = 0.975)
  
  df <- if(i == 1) { dat } else { dat025 } 
  df$Abundance <- df$nTP
  df$lower <- df$nTP
  df$upper <- df$nTP
  df[Data$ucellidxNV,"Abundance"] <- NTPmu
  df[Data$ucellidxNV,"lower"] <- NTPlower
  df[Data$ucellidxNV,"upper"] <- NTPupper
  df$Model <- if(i == 1) { "Adult" } else { "Juvenile" } 
  
  ll[[i]] <- df
}

df <- ll |> bind_rows()

# --- Area per site in m^2
sareas <- table(df$site)*25

# --- Area by burnt status: 
table(df$burnt)*25/20000

# --- Predicted N: 
df |> group_by(Model, burnt) |> 
  summarize(N = sum(Abundance), Nlower = sum(lower), Nupper = sum(upper)) |> ungroup() |> 
  arrange(Model, N) |> print(n = Inf)

# --- Predicted Density:
df |> 
  group_by(Model, site) |> 
  summarize(N = sum(Abundance), Nlower = sum(lower), Nupper = sum(upper)) |> ungroup() |> 
  arrange(Model, N) |> 
  mutate(N = N/rep(sareas, 2), Nlower = Nlower/rep(sareas, 2), Nupper = Nupper/rep(sareas, 2)) |>print(n = Inf)

# --- Predicted N by site: 
# Table S3: preds
df |> group_by(Model, site) |> 
  summarize(N = sum(Abundance), Nlower = sum(lower), Nupper = sum(upper)) |> 
  ungroup() |> as.data.frame() |> select(-geometry) |> 
  mutate(site_area = rep(sareas, 2), 
         density = N/site_area, densityL = Nlower/site_area, densityU = Nupper/site_area) -> ctsbysite

# --- N from OBIA
ll <- list.files("data/plantdat_artr_site/", full.names = TRUE)

artr <- lapply(ll, function(i) {
  tmp <- read.csv(i) |> filter(is.na(FN) | FN == 0) |> 
    mutate(Model = c("Adult", "Juvenile")[grepl("025", i) + 1 ]) 
  }) |> bind_rows()

ctsbysite$prop <- as.numeric(table(artr$site, artr$Model))/ctsbysite$N
ctsbysite$propL <- as.numeric(table(artr$site, artr$Model))/ctsbysite$Nlower
ctsbysite$propU <- as.numeric(table(artr$site, artr$Model))/ctsbysite$Nupper

ctsbysite[,c("prop", "propL", "propU")]-1
props <- ctsbysite[,c("prop", "propL", "propU")]
props$overlap <- sign(props$propL-1) + sign(props$propU-1)
sum(props$overlap == 0)/nrow(props)

# Table S3:OBIA
table(artr$site, artr$Model)
table(artr$site, artr$Model)/rep(sareas, 2)
