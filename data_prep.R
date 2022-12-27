# This code relies on an external data storage device. 

pkgs <- c("tidyverse", "terra", "sf", "stars")
sapply(pkgs, require, character.only = TRUE)


# --- load data
# - field
plants <- st_read("~/../../Volumes/az_drive/FieldData/FieldData_2022/plant_data.geojson") |> 
  filter(site != "DuncanSaddle") |> 
  st_transform(32611)
cells <- st_read("~/../../Volumes/az_drive/FieldData/FieldData_2022/cells.geojson") |> 
  filter(site != "DuncanSaddle") |> select(site, cid) |>
  st_transform(32611)
blines <- st_read("~/../../Volumes/az_drive/FieldData/FieldData_2022/burn_lines.geojson") |>
  filter(site != "DuncanSaddle") |> select(burn) |>
  st_transform(32611)

# - rasters: template, 5m
l0 <- list.files("~/../../Volumes/az_drive/uas_2022_5m/", full.names = TRUE)
rnull <- lapply(l0[-grep("DuncanSaddle", l0)], function(x) {
  rast(x) |> set.crs("epsg:32611") })
# - rasters: chms
lchm <- list.files("~/../../Volumes/az_drive/uas_2022/", pattern = "_chm_", full.names = TRUE, recursive = TRUE) 
chmMax <- lapply(1:10, function(x) {
  rast(lchm[-grep("DuncanSaddle", lchm)][x]) |> crop(rnull[[x]]) |> resample(rnull[[x]], method = "max")})
chmAvg <- lapply(1:10, function(x) {
  rast(lchm[-grep("DuncanSaddle", lchm)][x]) |> crop(rnull[[x]]) |> resample(rnull[[x]], method = "average")})
# - rasters: heterogeneity (scales 3, 6, 8)
l3 <- list.files("~/../../Volumes/az_drive/wave_chms/", pattern = "L3", full.names = TRUE) 
l6 <- list.files("~/../../Volumes/az_drive/wave_chms/", pattern = "L6", full.names = TRUE) 
l8 <- list.files("~/../../Volumes/az_drive/wave_chms/", pattern = "L8", full.names = TRUE) 
het3 <- lapply(1:10, function(x) {
  rast(l3[-grep("DuncanSaddle", l3)][x]) |> crop(rnull[[x]]) |> resample(rnull[[x]])})
het6 <- lapply(1:10, function(x) {
  rast(l6[-grep("DuncanSaddle", l6)][x]) |> crop(rnull[[x]]) |> resample(rnull[[x]])})
het8 <- lapply(1:10, function(x) {
  rast(l8[-grep("DuncanSaddle", l8)][x]) |> crop(rnull[[x]]) |> resample(rnull[[x]])})
# - rasters: dtm
ll <- list.files("~/../../Volumes/az_drive/uas_2022/", pattern = "_dtm_", full.names = TRUE, recursive = TRUE) 
elev <- lapply(1:10, function(x) {
  rast(ll[[x]]) |> crop(rnull[[x]]) |> resample(rnull[[x]], method = "average")}) 
tpis <- lapply(1:10, function(x) {
  tmp <- rast(ll[x]) |> set.crs("epsg:32611") |> 
    aggregate(fact = 166) |> terrain(v = "TPI", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })
slopes <- lapply(1:10, function(x) {
  tmp <- rast(ll[x]) |> set.crs("epsg:32611") |> 
    aggregate(fact = 166) |> terrain(v = "slope", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })

cdat <- data.frame(site = NULL, het3 = NULL, het6 = NULL, het8 = NULL, 
                   tpi = NULL, slope = NULL, elev = NULL,
                   max.ht = NULL, avg.ht = NULL)
for(i in 1:length(unique(cells$site))) {
  j <- unique(cells$site)[i]
  df <- data.frame(site = rep(j, table(cells$site)[i]), 
                   het3 = NA, het6 = NA, het8 = NA, 
                   tpi = NA, slope = NA, elev = NA,
                   max.ht = NA, avg.ht = NA)
  df[,2] <- extract(het3[[i]], cells[cells$site == j,])[,2]
  df[,3] <- extract(het6[[i]], cells[cells$site == j,])[,2]
  df[,4] <- extract(het8[[i]], cells[cells$site == j,])[,2]
  df[,5] <- extract(tpis[[i]], cells[cells$site == j,])[,2]
  df[,6] <- extract(slopes[[i]], cells[cells$site == j,])[,2]
  df[,7] <- extract(elev[[i]], cells[cells$site == j,])[,2]
  df[,8] <- extract(chmMax[[i]], cells[cells$site == j,])[,2]
  df[,9] <- extract(chmAvg[[i]], cells[cells$site == j,])[,2]
  cdat <- bind_rows(cdat, df)
}
# -----

# --- combine raster values with field
# - add: distance to fire edge, easting, northing
cells |> 
  mutate(burnt = st_intersection(cells, blines)$burn, 
         east = as.numeric(st_coordinates(geometry)[,1]), 
         north = as.numeric(st_coordinates(geometry)[,2]), 
         ucid = as.numeric(paste0(as.numeric(as.factor(site)), "00",cid)),
         dist = 0) |> 
  bind_cols(cdat |> select(-site)) -> celldf
# add distance
dmat <- st_distance(celldf, by_element = FALSE)
idx <- which(celldf$burnt == 1)
idxb <- dmat[-idx,]  
for(i in idx) {celldf$dist[i] = min(idxb[,i])}


# add counts
rcounts <- lapply(1:10, function(x) {
  tmp <- st_as_sf(st_as_stars(rnull[[x]]), as.points = FALSE, merge = TRUE) |> select()
  tmp |> 
    mutate(n = lengths(st_intersects(geometry, plants |>
                                            filter(Species %in% c("ARTR", "ARTRW", "ARTR4", "ARTRV")), 
                                     sparse = TRUE)),
           shrubs_n = lengths(st_intersects(geometry, plants |>
                                       filter(Class == "Shrub"), 
                                     sparse = TRUE)),
           gt25_n = lengths(st_intersects(geometry, plants |>
                                       filter(Species %in% c("ARTR", "ARTRW", "ARTR4", "ARTRV"), Ht_gt_25 == 1), 
                                     sparse = TRUE))) 
  }) |> bind_rows()

celldf |>
  mutate(n = st_intersection(celldf, rcounts)$n,
         shrubs_n = st_intersection(celldf, rcounts)$shrubs_n,
         gt25_n = st_intersection(celldf, rcounts)$gt25_n) |>
  as.data.frame() |> select(-geometry) -> celldf

write.csv(celldf, "data/celldat.csv", row.names = FALSE)
# --- end cell-level data
