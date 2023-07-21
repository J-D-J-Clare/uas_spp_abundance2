pkgs <- c("tidyverse", "sf", "terra", 'ggspatial')
sapply(pkgs, require, character.only = TRUE)
source("00_helper_fn.R")

geo <- st_read("data/celldat_covar.geojson")
# ===========================
dat <- geo |> select(ucid) |> 
  right_join(read.csv("data/celldat_input.csv"))
# ======================
dat025 <- geo |> select(ucid) |> 
  right_join(read.csv("data/celldat_input025.csv"))
# ======================

# === betas
ll <- list.files("results/", pattern = "^beta_", full.names = TRUE)

lapply(ll, function(x) {read.csv(x) |> 
    mutate(mod = str_split(x, pattern = "_", simplify = TRUE)[1,6])}) |> 
  bind_rows() |> 
  filter(par == "beta", !grepl("Intercept", name), !grepl("Int", name)
         #!grepl("FP", name)
         ) |> 
  mutate(name = sub("TP", "(TP)", name),
         name = sub("FP", "(FP)", name), 
         name = sub("tpi", "TPI", name), 
         name = sub("tri", "TRI", name),
         name = sub("hli", "HLI", name),
         mod = sub(".csv", "", mod),
         mod = ifelse(mod == "adult", "Adult", "Juvenile")) |> 
  mutate(fp = ifelse(grepl("FP", name) == 1, "FP", "TP"),
         name = as.factor(gsub('\\.(.*)$', "", name))) |> 
  mutate(
    # name = fct_relevel(name, rev(c('CorralsTrail', 'SouthTrail', 'LowerDryCreek', 'InitialPointBurn', 'NorthHam', 'Cold', 'PonyComplex1', 'PonyComplex2', 'SodaNaturalArea2', 'SodaNaturalArea1'))), 
    fpf = factor(fp, levels=c('TP','FP'))) |> 
  ggplot(aes(y = name, x = mean, colour = mod)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x = mean, y = name, xmin = lower, xmax = upper), width = .1, position = position_dodge(width = .5)) + 
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray") +
  scale_colour_discrete(breaks=c('juv', 'adult')) +
  scale_colour_viridis_d(begin = .2, end = 0.8) +
  # scale_y_discrete(labels = c("Height_avg", "Area segmented", "Wind_max") ) +
  labs(x = "Effect size", y = "", colour = "Size class") + 
  theme_bw() +  
  facet_wrap(.~fpf, nrow = 2)

# ggsave("figures/p.png", width = 12, height = 8, units = "cm")
# ggsave("figures/psi.png", width = 12, height = 8, units = "cm")
# ggsave("figures/beta.png", width = 12, height = 12, units = "cm")
  

# === dist
ll <- list.files("results/", pattern = "^distL_", full.names = TRUE)
lapply(ll, function(x) {read.csv(x) |> 
    mutate(mod = str_split(x, pattern = "_", simplify = TRUE)[1,2])}) |> 
  bind_rows() |> 
  mutate(mod = sub(".csv", "", mod),
         mod = ifelse(mod == "adult", "Adult", "Juvenile")) -> out
out |>
  ggplot(aes(x = name, y = mean, group = mod, colour = mod)) + 
  geom_line(linewidth = 1.5, linetype = "solid") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = mod, group = mod, colour = mod), alpha = .15, linetype = "dotted", linewidth = .25) + 
  scale_x_continuous(labels = function(x){x*100}) +
  scale_y_continuous(labels = function(x){x/25}) +
  scale_colour_viridis_d(begin = .2, end = 0.8) +
  scale_fill_viridis_d(begin = .2, end = 0.8) +
  labs(x = "Distance from fire edge, m", 
       y = expression(paste("Effect size,", ~Delta,"E[",log(lambda),"]"~m^{-2})), 
       colour = "Size class", fill = "Size class") +
  theme_bw() 

ggsave("figures/dist.png", width = 12, height = 8, units = "cm")


# === PPC TP 
ll <- list.files("results/", pattern = "^ppc_", full.names = TRUE)

lapply(ll, function(x) {read.csv(x) |> 
    mutate(mod = str_split(x, pattern = "_", simplify = TRUE)[1,2]) |> 
    mutate(mod = sub(".csv", "", mod))}) |> 
  bind_rows() |> 
  mutate(mod = ifelse(mod == "adult", "Adult size class", "Juvenile size class")) |>
  # filter(mod == "juv") |>
  ggplot(aes(x = obs, y = Estimate)) + 
  geom_point(alpha = 0.25) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0, alpha = .25) + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed counts", y = "Estimated counts") +
  theme_bw() +
  facet_wrap(.~mod, scales = "free")

ggsave("figures/ppcnTP.png", width = 18, height = 10, units = "cm")


# === PPC FP
ll <- list.files("results/", pattern = "^ppcFP_", full.names = TRUE)

lapply(ll, function(x) {read.csv(x) |> 
    mutate(mod = str_split(x, pattern = "_", simplify = TRUE)[1,2]) |> 
    mutate(mod = sub(".csv", "", mod))}) |> 
  bind_rows() |> 
  mutate(mod = ifelse(mod == "adult", "Adult size class", "Juvenile size class")) |>
  # filter(mod == "juv") |>
  ggplot(aes(x = obs, y = Estimate)) + 
  geom_point(alpha = 0.25) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0, alpha = .25) + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed counts", y = "Estimated counts") +
  theme_bw() +
  facet_wrap(.~mod, scales = "free")

ggsave("figures/ppcFP.png", width = 18, height = 10, units = "cm")


# === PPC Unknowns
ll <- list.files("results/", pattern = "^ppcUnk_", full.names = TRUE)

lapply(ll, function(x) {read.csv(x) |> 
    mutate(mod = str_split(x, pattern = "_", simplify = TRUE)[1,2]) |> 
    mutate(mod = sub(".csv", "", mod))}) |> 
  bind_rows() |> 
  mutate(mod = ifelse(mod == "adult", "Adult size class", "Juvenile size class")) |>
  # filter(mod == "juv") |>
  ggplot(aes(x = obs, y = Estimate)) + 
  geom_point(alpha = 0.25) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0, alpha = .25) + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed counts", y = "Estimated counts") +
  theme_bw() +
  facet_wrap(.~mod, scales = "free")

ggsave("figures/ppcUnk.png", width = 18, height = 10, units = "cm")

# === Predictive maps
blines <- st_read("~/../../Volumes/az_drive/FieldData/FieldData_2022/burn_lines.geojson") |> st_transform(32611)
ll <- list()
for(k in 1:2){
  class <- c("adult", "juvenile")[k]
  
  NTP <- readRDS(paste0("~/../../Volumes/az_drive/detection/models/nf", k,".2nb_NTP.rds"))
  DC <- readRDS(paste0("~/../../Volumes/az_drive/detection/models/nf", k, ".2nb_DataConstants.rds"))
  Data <- DC[['Data']]
  Constants <- DC[['Constants']]

  NTPhat <- apply(NTP, 2, sd)

  df <- if(k == 1) { dat } else { dat025 } 
  df$Abundance <- df$nTP
  df[Data$ucellidxNV, "Abundance"] <- NTPhat
  ll[[k]] <- df
}
for(k in 1:2){
  class <- c("adult", "juvenile")[k]
  
  sites <- unique(ll[[k]]$site)
  for(i in sites){
    maxct <- bind_rows(ll) |> filter(site == !!i, val == 0) |> pull(Abundance) |> max()
    df <- ll[[k]] |> mutate(Validated = as.factor(val)) |>
      filter(site == !!i) |> mutate(Abundance = ifelse(val == 1, 0, Abundance))
    df |>
      ggplot() +
      # geom_sf(aes(fill = Abundance, colour = Validated), linewidth = .33) +
      geom_sf(aes(fill = Abundance)) +
      geom_sf(data = blines |> st_crop(df |> filter(site == !!i)), 
              mapping = aes(), fill = NA, colour = "white", linewidth = 1) +
      annotation_scale(location = "br", pad_x = unit(1, "cm"), pad_y = unit(0.01, "cm") ) +
      scale_fill_viridis_c(option = "viridis", trans = "sqrt", 
                           limits = c(0, maxct)) +
      #scale_color_manual(breaks = "1", values = c("1" = "black")) +
      labs(x = 'Easting', y = "Northing", title = paste0(i, ", ", class),
           fill = "Abundance, SD") +
      #theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.background = element_blank())

    ggsave(paste0("figures/pred_maps/", i, "_sd_", class,".png")) #, width = 12, units = "cm"
  }
}

# Figure 1c
dat |> 
  filter(site == "InitialPointBurn") |>
  select(site, ucid, val, TP, FN, FP, Unk) |> 
  pivot_longer(cols = c(TP, FN, FP, Unk)) |> 
  mutate(value = ifelse(name == "Unk" & val == 1, NA, value), 
         name = ifelse(name == "Unk", "TP + FP", name)) |>
  ggplot() + 
  geom_sf(aes(fill = value)) +
  scale_fill_viridis_c(option = "viridis", trans = "sqrt") +
  labs(x = "Easting", y = "Northing", fill = "Count") +
  facet_wrap(.~name) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank(), text = element_text(size = 14)) 
  
ggsave("figures/countsV.png")
