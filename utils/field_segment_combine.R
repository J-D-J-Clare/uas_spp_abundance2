pkgs <- c("tidyverse", "terra", "sf", "exactextractr")
sapply(pkgs, require, character.only = TRUE)

# === Content
# === 0. Load data
# === 1. Match field/gps with segments and spectral info
# === 2. Match classified layer with inputs

# ================================
path <- "D:/Andrii/ML_data_combine/"
plant_df <- st_read(paste0(path, "plant_data.geojson"))
cell_df <- st_read(paste0(path, "cells.geojson")) |> 
  st_transform(32611)

sites <- sapply(list.files(paste0(path, "chms/")), function(x){
  str_split(x, "_")[[1]][[1]]
}) |> as.character()

# --- create output list
ll <- list()

for(i in 1:length(sites)) {

  print(paste("Working on:", i, "/ 10 sites"))
# =========== Match field and gps data
site <- sites[i]

# === 1. Match Field/gps data with segments and extract CHM attributes

ttops <- st_read(paste0(path, "segmentation_outputs/", site,"_ttops.geojson")) 
cro <- st_read(paste0(path, "segmentation_outputs/", site, "_crowns.geojson"))

pldf <- plant_df |> filter(site == !!site)

cro |>
  mutate(Species = sapply(st_intersects(geometry, st_buffer(pldf, .05)), 
                          function(x) { tmp = pldf$Species[x]
                          out = ifelse(length(tmp) == 0 | length(tmp) > 1, NA, tmp)}),
         Class = sapply(st_intersects(geometry, st_buffer(pldf, .05)), 
                        function(x) { tmp = pldf$Class[x]
                        out = ifelse(length(tmp) == 0 | length(tmp) > 1, NA, tmp)})) |>
  rename(ID = treeID) -> ml_df

# --- extract CHM attributes
chm <- rast(list.files(paste0(path, "chms/"), pattern = site, full.names = TRUE))
ortho.rgb <- rast(list.files(paste0(path, "orthos_rgb/"), pattern = site, full.names = TRUE)) 
ortho.multispec <- rast(list.files(paste0(path, "orthos_multispec/"), pattern = site, full.names = TRUE))

ml_df |> 
  mutate(site = !!site) |>
  mutate(area = as.numeric(st_area(geometry)), 
         perimeter = as.numeric(st_length(st_cast(geometry,"MULTILINESTRING"))),
         edge_to_area_ratio = perimeter/area) |> 
  bind_cols(exactextractr::exact_extract(chm, ml_df, fun = c("stdev", "mean", "max", "min"))) |> 
  rename_with(function(x){paste0("chm_", x)}, c("stdev", "mean", "max", "min")) -> b


# --- extract RGB values
b |> 
  bind_cols(exactextractr::exact_extract(ortho.rgb, ml_df, fun = c("stdev", "mean")) |> 
              rename_with(function(x){paste0("rgb_", x)})) |> 
  select(-contains(!!site)) -> b

# --- extract multispec values
b |> 
  bind_cols(exactextractr::exact_extract(ortho.multispec, ml_df, fun = c("stdev", "mean")) |> 
              rename_with(function(x){paste0("multispec_", x)})) |> 
  select(-contains(!!site)) -> b

# --- add spectral indices
nir <- ortho.multispec[["NIR"]]
red <- ortho.multispec[["Red"]]
blue <- ortho.multispec[["Blue"]]
RE <- ortho.multispec[["RedEdge"]]

ndvi <- ((nir - red)/(nir + red))
msavi <- (1/2) * (2 * (nir +1) - sqrt((2 * nir + 1)^2 - 8*(nir - red))) 

sidx <- c(ndvi, msavi)
names(sidx) <- c("ndvi", "msavi")

b |> 
  bind_cols(exactextractr::exact_extract(sidx, ml_df, fun = c("stdev", "mean")) |> 
              rename_with(function(x){paste0("spec_", x)})) |> 
  mutate(uid = paste0(site, "_", ID)) -> b

# --- store in a list
ll[[i]] <- b

rm(b)
}

tmp <- lapply(ll, function(x) {
  x  |> st_cast('MULTIPOLYGON') |> 
    st_intersection(cell_df |> filter(site == x$site[1]) |> st_buffer(10) |> st_bbox() |> st_as_sfc() ) |>
    mutate(uid = paste0(site, "_", ID) )
})
lapply(tmp, function(x) {x |> write_sf(paste0(path, "ml_dataset_sites/ml_dataset_", x$site[1], ".geojson"), delete_dsm = TRUE)})


dfc <- do.call("bind_rows", tmp)
dfc |> 
  # select(-Z) |> st_cast('MULTIPOLYGON') |> write_sf(paste0(path, "ml_dataset.geojson"), delete_dsn = TRUE)
  as.data.frame() |> 
  select(-Z, -geometry) |>
  write.csv(paste0(path, "ml_dataset.csv"), row.names = FALSE)
# ===============================================================================
  
# === Reconstruct 5m raster from a sample of cell-level data (centroids)
cells <- st_read("~/Desktop/UAV_demography/FieldData/SurveyPts_shapes_2022/SoutTrail/SourthTrail_pts_wgs84.shp") |> 
  st_transform(32611)

cells |>
  mutate(X = st_coordinates(geometry)[,1],
         Y = st_coordinates(geometry)[,2], 
         Z = cid) |> as.data.frame() |>
  dplyr::select(-cid, -geometry) -> xydat
  
nr <- diff(range(xydat$Y))/5 + 1
nc <- diff(range(xydat$X))/5 + 1
ncell <- nr*nc
rp = rast(nrows = nr, ncols = nc, vals = sample(runif(ncell)), 
         extent = ext(cells) + 2.5) 

plot(rp)
plot(cells$geometry, add = TRUE, pch = 19)
plot(st_buffer(cells, 5)$geometry, add = TRUE, pch = 19)

# export out the raster: z values are place holders
writeRaster(rp, 
            filename = "~/Desktop/UAV_demography/FieldData/2022_campaign_field_data/cells_data/SouthTrail_5m_raster.tif")
st_write(cells, "~/Desktop/UAV_demography/FieldData/2022_campaign_field_data/cells_data/SouthTrail_cells.geojson")
# === end



# === 2. Merge classified data with the geospatial layer
preds <- read.csv("D:/Andrii/ML_data_combine/predicted_plants_rf_v2.csv")
dfg <- st_read("D:/Andrii/ML_data_combine/ml_dataset.geojson")
# --- reproduce Trevor's workflow for predictions: NA's are excluded
# - use unique id to merge: uid = site + ID
mldf <- read.csv("D:/Andrii/ML_data_combine/ml_dataset.csv") |> mutate(uid = paste0(site, "_", ID))
# nas <- rowSums(mldf[,5:31])
# pdat <- mldf[-which(is.na(nas == TRUE)), ]
# ---
dfg |> 
  left_join(preds, by = "uid") -> out

fs <- unique(out$site)
lapply(1:length(fs), function(x) {
  tmp <- filter(out, site == fs[x])
  st_write(tmp, paste0("D:/Andrii/ML_data_combine/predicted_plants_rf_v2_sites/predicted_plants_rf_v2_", fs[x], ".geojson"))
})

st_write(out, "D:/Andrii/ML_data_combine/predicted_plants_rf_v2.geojson")
# === end of merge





