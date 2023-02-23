# This code relies on an external data storage device. 

pkgs <- c("tidyverse", "terra", "sf", "stars")
sapply(pkgs, require, character.only = TRUE)

`%notin%` <- negate(`%in%`)

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

# - rasters: templates, 5m resolution
l0 <- list.files("~/../../Volumes/az_drive/uas_2022_5m/", full.names = TRUE)
rnull <- lapply(l0[-grep("DuncanSaddle", l0)], function(x) {
  tmp <- rast(x) 
  set.crs(tmp, "epsg:32611")
  return(tmp)
  })
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
# - rasters: dtm + topographic derivatives
ll <- list.files("~/../../Volumes/az_drive/uas_2022/", pattern = "_dtm_", full.names = TRUE, recursive = TRUE) 
elev <- lapply(1:10, function(x) {
  rast(ll[[x]]) |> crop(rnull[[x]]) |> resample(rnull[[x]], method = "average")}) 
tpis <- lapply(1:10, function(x) {
  tmp <- rast(ll[x]) 
  set.crs(tmp, "epsg:32611") 
  tmp |>
    aggregate(fact = 166) |> terrain(v = "TPI", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })
slopes <- lapply(1:10, function(x) {
  tmp <- rast(ll[x])
  set.crs(tmp, "epsg:32611") 
  tmp |>
    aggregate(fact = 166) |> terrain(v = "slope", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })

# - combine raster info into a data frame
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

# --- combine raster values with field data and site chars
# - add: distance to fire edge, easting, northing
cells |> 
  st_intersection(blines) |>
  mutate(east = as.numeric(st_coordinates(geometry)[,1]), 
         north = as.numeric(st_coordinates(geometry)[,2]), 
         ucid = paste0(site, "_", cid),
         dist = 0) |> 
  bind_cols(cdat |> select(-site)) -> celldf
# add distance
dmatb <- st_distance(celldf, filter(blines, burn == 0))
dmatnb <- st_distance(celldf, filter(blines, burn == 1))
idb <- which(celldf$burn == 1)
idnb <- which(celldf$burn == 0)
for(i in idb) {celldf$dist[i] = min(dmatb[i,])}
for(i in idnb) {celldf$dist[i] = -min(dmatnb[i,])}

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
  st_intersection(rcounts) |>
  as.data.frame() |> select(-geometry) -> celldfout

write.csv(celldfout, "data/celldat.csv", row.names = FALSE)

# --- full raster dataset for predictions
rcell <- rcounts |> select()
s = lapply(1:length(rnull), function(x) {
  j <- unique(cells$site)[x]
  rep(j, length(as.matrix(rnull[[x]])))
}) |> unlist()

rcell <- rcounts |> select() |> 
  mutate(site = s,
         cid = NA,
         het3 = NA, het6 = NA, het8 = NA, 
         tpi = NA, slope = NA, elev = NA, 
         max.ht = NA, avg.ht = NA,
         east = as.numeric(st_coordinates(st_centroid(geometry))[,1]), 
         north = as.numeric(st_coordinates(st_centroid(geometry))[,2])) 

for(i in 1:length(rnull)) {
  j <- unique(cells$site)[i]
  idx <- which(rcell$site == j)
  
  rcell[idx,"cid"] <- extract(rnull[[i]], rcell[idx,])[,2] 
  rcell[idx,"het3"] <- extract(het3[[i]], rcell[idx,])[,2]
  rcell[idx,"het6"] <- extract(het6[[i]], rcell[idx,])[,2]
  rcell[idx,"het8"] <- extract(het8[[i]], rcell[idx,])[,2]
  rcell[idx,"tpi"] <- extract(tpis[[i]], rcell[idx,],)[,2]
  rcell[idx,"slope"] <- extract(slopes[[i]], rcell[idx,])[,2]
  rcell[idx,"elev"] <- extract(elev[[i]], rcell[idx,])[,2]
  rcell[idx,"max.ht"] <- extract(chmMax[[i]], rcell[idx,])[,2]
  rcell[idx,"avg.ht"] <- extract(chmAvg[[i]], rcell[idx,])[,2]
}

rcell |> 
  mutate(ucid = paste0(site, "_", cid)) |> 
  bind_cols(rcounts |> as.data.frame() |> select(-geometry)) |> 
  mutate(burn = lengths(st_intersects(geometry, blines |> 
                                        filter(burn == 1))),
         val = ifelse(ucid %in% celldfout$ucid, 1, 0),
         dist = 0) -> rcellout
# --- add distance
dmatb <- st_distance(rcellout, filter(blines, burn == 0))
dmatnb <- st_distance(rcellout, filter(blines, burn == 1))
idb <- which(rcellout$burn == 1)
idnb <- which(rcellout$burn == 0)
for(i in idb) {rcellout$dist[i] = min(dmatb[i,])}
for(i in idnb) {rcellout$dist[i] = -min(dmatnb[i,])}
# ---

st_write(rcellout, "data/celldat_full.geojson")

# --- end cell-level data


# === Add Individual-level data
# inputs: 
# plant_data.geojson - field data; 
# {site}_5m_raster.tif - define ROI and CID attribute for each plant
# celldat_full.geojson - cell-level dataset with plant counts
# predicted_plants_rf_v1_{site}.geojson - predicted labels + P(ARTR) for segmented layer

# outputs: plantdat_artr_full.csv, plantdat_artr_full.geojson
ll <- list()
sites <- unique(cells$site)
for(i in 1:length(sites)) {
  site <- sites[i]
  # --------------------------------
  dfg <- st_read("data/celldat_full.geojson") |> 
    filter(site == !!site) 
  
  # specify paths
  path <- "~/../../Volumes/az_drive/FieldData/FieldData_2022/predicted_plants_rf_v1_sites/"
  # load segments with predictions
  dfi <- st_read(paste0(path, "predicted_plants_rf_v1_", site, ".geojson")) |>
    select(ID, Species, Class, site, plant_type_preds, prob_ARTR, uid) |> 
    # re-code all sagebrush species
    mutate(Species0 = ifelse(Species %in% c('ARTR', 'ARTR4', 'ARTRV', 'ARTRW'), 'ARTR', Species)) |>
    st_crop(ext(rnull[[i]])) 
  
  st_agr(dfi) <- 'constant'
  st_agr(dfg) <- 'constant'
  
  dfi |> 
    select(-ID, -uid) |>
    mutate(ucid = st_intersection(dfi |>
                                    st_centroid(), dfg[,"ucid"]) |> as.data.frame() |> pull(ucid), 
           # assign segments field ID, Species by a single overlap. Otherwise NA
           id = sapply(st_intersects(geometry, st_buffer(plants, .1)), 
                       function(x) { tmp = plants$ID[x]
                       out = ifelse(length(tmp) == 0 | length(tmp) > 1, NA, tmp)}),
           # assign whether the matching segments are below or above 25cm
           Ht_gt_25 = sapply(st_intersects(geometry, st_buffer(plants, .1)), 
                             function(x) { tmp = plants$Ht_gt_25[x]
                             out = ifelse(length(tmp) == 0 | length(tmp) > 1, NA, tmp)}),
           # assign 1 if a plant is detected but not counted as so because there are two plants overlapping the segment
           detMult = sapply(st_intersects(geometry, st_buffer(plants, .1)), 
                             function(x) { tmp = plants$Ht_gt_25[x]
                             out = ifelse(length(tmp) > 1, length(tmp), 0)})) |> 
    # remove marginal matched plants that are outside the central 5m cell, ie overlap by polygon margins
    filter( !(!is.na(Species0) & is.na(id)), !(is.na(Species0) & !is.na(id)) ) -> a
  
  st_agr(plants) <- 'constant'
  
  plants |> 
    filter(site == !!site) |>
    mutate(Species0 = ifelse(Species %in% c('ARTR', 'ARTR4', 'ARTRV', 'ARTRW'), 'ARTR', Species)) |>
    rename("id" = "ID") |> 
    select(id, CID, Class, Species, Species0, site, Ht_gt_25) |> 
    # missed ID
    mutate(plant_type_preds = NA, prob_ARTR = NA, 
           ucid = paste0(site, "_", CID)) |> 
    # filter only field plants that were not detected to avoid duplicates
    filter(id %notin% unique(a$id)) -> tmp
  
  st_agr(tmp) <- 'constant'
  
  out <- a |>
    st_centroid() |>
    mutate(z = NA, y = 1) |>
    mutate(z = ifelse(Species0 == 'ARTR', 1, z), # present 
           z = ifelse(Species0 != 'ARTR' & !is.na(id), 0, z), # present not ARTR
           # if z is not field data & not validated as true/false presence -> unknown 
           ) |>  
    # append field data that was not detected
    bind_rows(tmp |>
                mutate(z = ifelse(Species0 == 'ARTR', 1, 0), # field spp ARTR
                       y = ifelse(is.na(plant_type_preds), 0, 1)) ) |>
    # filter only plants that were classified as ARTR or field-recorded as ARTR
    # Species0 reduces all sage shrubs variants of interest to ARTR 
    filter(Species0 == 'ARTR' | plant_type_preds == 'ARTR')  
  
  ll[[i]] <- out
}

ll |> bind_rows() |>
  write_sf("data/plantdat_artr_full.geojson", delete_dsn = TRUE)
ll |> bind_rows() |> 
  as.data.frame() |> select(-geometry) |> 
  write.csv("data/plantdat_artr_full.csv", row.names = FALSE)

# === end Individual-level data

# === Add precieved counts to the cell-level dataset
plantdat_artr <- st_read("data/plantdat_artr_full.geojson") 

dfg <- st_read('data/celldat_full.geojson') |>
  select(-gt25_n, -shrubs_n) |>
  mutate(# n_obs: number of drone detected plants 
         n_obs = lengths(st_intersects(geometry, plantdat_artr |> filter(y == 1), sparse = TRUE)), 
         # n: number of field mapped plants
         n = ifelse(val == 1, n, NA),
         # ngt25: number of field mapped plants taller than 25 cm
         ngt25 = lengths(st_intersects(geometry, filter(plantdat_artr, Ht_gt_25 == 1 & z == 1), sparse = TRUE)),
         ngt25 = ifelse(val == 1, ngt25, NA), 
         # nTP: number of field mapped + detected with UAS
         nTP = lengths(st_intersects(geometry, filter(plantdat_artr, y == 1 & z == 1), sparse = TRUE)),
         nTP = ifelse(val == 1, nTP, NA), 
         # nTPM: number of field mapped + detected + not NTP b/c 1 segment = 2 field points
         nTPM = lengths(st_intersects(geometry, filter(plantdat_artr, detMult > 0), sparse = TRUE)),
         nTPM = ifelse(val == 1, nTPM, NA), 
         # nFP: number of field mapped + detected + classified as one but not in fact ARTR
         # can be assumed to be absent, since we have only 17 plants across all sites
         nFP = lengths(st_intersects(geometry, filter(plantdat_artr, y == 1 & z == 0), sparse = TRUE)),
         nFP = ifelse(val == 1, nFP, NA),
         # nM: number of field plants not detected with UAS
         nM = lengths(st_intersects(geometry, filter(plantdat_artr, y == 0 & z == 1), sparse = TRUE)),
         nM = ifelse(val == 1, nM, NA),
         # nU: number of plants detected w/ UAS, but unknown validation status
         nU = lengths(st_intersects(geometry, filter(plantdat_artr, detMult == 0 & is.na(z)), sparse = TRUE)), 
         nU = ifelse(val == 1, nU, NA))

write_sf(dfg, "data/celldat_full.geojson", delete_dsn = TRUE)
dfg |> 
  as.data.frame() |> select(-geometry) |> 
  write.csv("data/celldat_full.csv")
# ===

