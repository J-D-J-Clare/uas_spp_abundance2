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
  filter(site != "DuncanSaddle") |> select(site, cid) |> mutate(ucid = paste0(site, "_", cid)) |>
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
tris <- lapply(1:10, function(x) {
  tmp <- rast(ll[x]) 
  set.crs(tmp, "epsg:32611") 
  tmp |>
    aggregate(fact = 166) |> terrain(v = "TRI", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })
tris.sm <- lapply(1:10, function(x) {
  tmp <- rast(ll[x]) 
  set.crs(tmp, "epsg:32611") 
  tmp |>
    terrain(v = "TRI", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })
tpis <- lapply(1:10, function(x) {
  tmp <- rast(ll[x]) 
  set.crs(tmp, "epsg:32611") 
  tmp |>
    aggregate(fact = 166) |> terrain(v = "TPI", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })
tpis.sm <- lapply(1:10, function(x) {
  tmp <- rast(ll[x]) 
  set.crs(tmp, "epsg:32611") 
  tmp |> 
    terrain(v = "TPI", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })
slopes <- lapply(1:10, function(x) {
  tmp <- rast(ll[x])
  set.crs(tmp, "epsg:32611") 
  tmp |>
    aggregate(fact = 166) |> terrain(v = "slope", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })
slopes.sm <- lapply(1:10, function(x) {
  tmp <- rast(ll[x])
  set.crs(tmp, "epsg:32611") 
  tmp |>
    terrain(v = "slope", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })
flowdirs <- lapply(1:10, function(x) {
  tmp <- rast(ll[x])
  set.crs(tmp, "epsg:32611") 
  tmp |>
    aggregate(fact = 166) |> terrain(v = "flowdir", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })
flowdirs.sm <- lapply(1:10, function(x) {
  tmp <- rast(ll[x])
  set.crs(tmp, "epsg:32611") 
  tmp |>
    terrain(v = "flowdir", neighbors = 8) |>
    crop(rnull[[x]]) |> resample(rnull[[x]]) })
# - rasters: view counts
ll <- list.files("~/../../Volumes/az_drive/uas_2022/", pattern = "_view_counts_", full.names = TRUE, recursive = TRUE) 
viewcts <- lapply(1:10, function(x) {
  rast(ll[[x]]) |> crop(rnull[[x]]) |> resample(rnull[[x]], method = "med")}) 
# end data import


# === extract true plant counts (aka field/ground counts)
rcounts <- lapply(1:10, function(x) {
  s <- unique(cells$site)[x]
  tmp <- st_as_sf(st_as_stars(rnull[[x]]), as.points = FALSE, merge = TRUE) |> select()
  tmp |> 
    mutate(cid = extract(rnull[[x]], tmp)[,2], 
           site = s) |> 
    mutate(nARTR = lengths(st_intersects(geometry, plants |>
                                           filter(Species %in% c("ARTR", "ARTRW", "ARTRV")), 
                                         sparse = TRUE)),
           nNonARTR = lengths(st_intersects(geometry, plants |>
                                              filter(Class == "Shrub", Species %notin% c("ARTR", "ARTRW", "ARTRV")), 
                                            sparse = TRUE)),
           nARTRgt25 = lengths(st_intersects(geometry, plants |>
                                               filter(Species %in% c("ARTR", "ARTRW", "ARTRV"), Ht_gt_25 == 1), 
                                             sparse = TRUE)),
           nNonARTRgt25 = lengths(st_intersects(geometry, plants |>
                                                  filter(Class == "Shrub", Species %notin% c("ARTR", "ARTRW", "ARTRV"), Ht_gt_25 == 0), 
                                                sparse = TRUE)), 
           nShrubs = lengths(st_intersects(geometry, plants |> 
                                             filter(Class == "Shrub"))), 
           val = lengths(st_intersects(geometry, cells|> 
                                         filter(site == s), sparse = TRUE)), 
           ucid = paste0(site, "_", cid)) 
}) |> bind_rows()  

# === export field counts
rcounts |> # write_sf("data/celldat_fcounts.geojson", delete_dsn = TRUE)
  as.data.frame() |> select(-geometry) |> write.csv("data/celldat_fcounts.csv", row.names = FALSE)
# ===

# === extract cell-level covariates 
rcell <- rcounts |> select(site, cid, ucid, val) |> 
  mutate(het3 = NA, het6 = NA, het8 = NA, 
         tri = NA, tri.sm = NA, 
         tpi = NA, tpi.sm = NA, 
         slope = NA, slope.sm = NA, 
         flowdir = NA, flowdir.sm = NA,
         elev = NA,
         max.ht = NA, avg.ht = NA,
         viewct = NA,
         east = as.numeric(st_coordinates(st_centroid(geometry))[,1]), 
         north = as.numeric(st_coordinates(st_centroid(geometry))[,2]), 
         burnt = abs(lengths(st_intersects(geometry, blines |> filter(burn == 0))) -1) ,
         dist = 0) 

for(i in 1:length(rnull)) {
  j <- unique(cells$site)[i]
  idx <- which(rcell$site == j)
  
  rcell[idx,"cid"] <- extract(rnull[[i]], rcell[idx,])[,2] 
  rcell[idx,"het3"] <- extract(het3[[i]], rcell[idx,])[,2]
  rcell[idx,"het6"] <- extract(het6[[i]], rcell[idx,])[,2]
  rcell[idx,"het8"] <- extract(het8[[i]], rcell[idx,])[,2]
  rcell[idx,"tri"] <- extract(tris[[i]], rcell[idx,],)[,2]
  rcell[idx,"tri.sm"] <- extract(tris.sm[[i]], rcell[idx,],)[,2]
  rcell[idx,"tpi"] <- extract(tpis[[i]], rcell[idx,],)[,2]
  rcell[idx,"tpi.sm"] <- extract(tpis.sm[[i]], rcell[idx,],)[,2]
  rcell[idx,"slope"] <- extract(slopes[[i]], rcell[idx,])[,2]
  rcell[idx,"slope.sm"] <- extract(slopes.sm[[i]], rcell[idx,])[,2]
  rcell[idx,"flowdir"] <- extract(flowdirs[[i]], rcell[idx,])[,2]
  rcell[idx,"flowdir.sm"] <- extract(flowdirs.sm[[i]], rcell[idx,])[,2]
  rcell[idx,"elev"] <- extract(elev[[i]], rcell[idx,])[,2]
  rcell[idx,"max.ht"] <- extract(chmMax[[i]], rcell[idx,])[,2]
  rcell[idx,"avg.ht"] <- extract(chmAvg[[i]], rcell[idx,])[,2]
  rcell[idx,"viewct"] <- extract(viewcts[[i]], rcell[idx,])[,2]
}

# --- add distance
dmatb <- st_distance(rcell, filter(blines, burn == 0))
dmatnb <- st_distance(rcell, filter(blines, burn == 1))
idb <- which(rcell$burn == 1)
idnb <- which(rcell$burn == 0)
for(i in idb) {rcell$dist[i] = min(dmatb[i,])}
for(i in idnb) {rcell$dist[i] = -min(dmatnb[i,])}
# ---

# === export cell-level covariates
rcell |> # write_sf("data/celldat_covar.geojson", delete_dsn = TRUE) 
  as.data.frame() |> select(-geometry) |> write.csv("data/celldat_covar.csv", row.names = FALSE)
# ===
rcell <- st_read("data/celldat_covar.geojson")


# === Add Individual-level data

for(k in 1:length(unique(rcell$site))) {

site <- unique(rcell$site)[k]  
# this is pretty convoluted stuff. takes time to follow what each does.
  
segpol <- st_read(paste0("~/../../Volumes/az_drive/FieldData/FieldData_2022/predicted_plants_rf_v1_sites/predicted_plants_rf_v1_", site, ".geojson") ) |> 
  select(uid, site, chm_max, area, plant_type_preds, prob_ARTR) 
st_agr(segpol) <- "constant"
# classified polygons
segpol_nonartr <- segpol |> 
  filter(chm_max <= .5) |> # for medium plants only
  filter(plant_type_preds != 'ARTR')
segpol_artr <- segpol |> 
  filter(chm_max <= .5) |> # for medium plants only
  filter(plant_type_preds == 'ARTR')
# field points
artr <- plants |> 
  filter(!grepl("same", notes),
    site == !!site, Ht_gt_25 == 1, 
    Species %in% c('ARTR', 'ARTRW', 'ARTRV'))
st_agr(artr) <- "constant"
nonartr <- plants |> 
  filter(Species %notin% c('ARTR', 'ARTRW', 'ARTRV'))
st_agr(nonartr) <- "constant"

# === Level 1: perfect detections (one to one match)
l1idx <- lengths(st_intersects(segpol_artr, artr, sparse = TRUE) )

l1empty <- segpol_artr[l1idx == 0,]
# summarize certain L1
l1TP <- segpol_artr[l1idx == 1,] #TP pols
l1artr <- artr |> 
  mutate(var = lengths(st_intersects(geometry, segpol_artr[l1idx < 2,]))) |>   
  filter(var == 1) |> select(-var) |>
  mutate(TP = 1, FP = 0, FN = 0)  |> # TP pts
  st_intersection(segpol_artr[l1idx < 2, c('chm_max', 'area','uid','geometry')])

# === Level 1a: 1+ pts over 1 polygon (result: some FN and some TP)
l1artrFNtmp <- artr |> 
  st_intersection(segpol_artr[l1idx > 1, c('chm_max', 'area','uid', 'geometry')]) |>   
  mutate(dist = apply(
    as.matrix(st_distance(geometry, st_centroid(segpol_artr[l1idx > 1, ]) )), 
    1, min) ) |> 
  group_by(uid) |>
  filter(dist != min(dist)) |> ungroup() |> select(-dist)

notFN <- l1empty |> st_intersection(l1artrFNtmp |> st_buffer(.15) |> pull(geometry))
notFNpts <- l1artrFNtmp |> st_buffer(.15) |> st_intersection(l1empty[,'geometry']) |> st_centroid()

# --- summarize L1a
l1empty |> filter(uid %notin%  notFN$uid) -> l1empty
l1TP |> bind_rows(segpol_artr |> 
                    filter(uid %in% l1artrFNtmp$uid)) -> l1TP
l1artr |> 
  bind_rows(artr |> 
              st_intersection(segpol_artr[l1idx > 1, c('chm_max', 'area','uid', 'geometry')]) |>   
              mutate(dist = apply(
                as.matrix(st_distance(geometry, st_centroid(segpol_artr[l1idx > 1, ]) )), 
                1, min) ) |> 
              group_by(uid) |>
              filter(dist == min(dist)) |> ungroup() |> select(-dist)) |> 
  bind_rows(notFNpts) |>
  mutate(TP = 1, FP = 0, FN = 0) |> 
  st_centroid() -> l1artr

l1artrFN <- l1artrFNtmp |> filter(uid %notin% notFNpts$uid) |>
  mutate(TP = 0, FP = 0, FN = 1, chm_max = NA, area = NA)
# --- end L1a

# = Level2: check l1empty overlap with non-detected plants (radius = 0.15 m)
l2artrtmp <- artr |> 
  # plants not overlapping w/ anything
  mutate(var = lengths(st_intersects(geometry, segpol_artr))) |>   
  filter(var == 0) |> select(-var) |> st_buffer(.15)
# plants + 15cm buffer that overlap w/ empty pols 
l2idx <- lengths(st_intersects(l1empty, l2artrtmp, sparse = TRUE))

# a pol in l1empty can overlap 1+ points with 15cm buffers
# pol is a valid detection; overlapping pts except closest are FN
l2artrtmp |> 
  st_intersection(l1empty[l2idx > 0, c('chm_max', 'area','uid','geometry')]) |>
  mutate(dist = apply(
    as.matrix(st_distance(geometry, st_centroid(l1empty[l2idx > 0, ]) )), 
    1, min) ) |> 
  group_by(uid) |>
  filter(dist != min(dist)) |> ungroup() |> select(-dist) |> 
  st_centroid() -> l2artrFNadd

# init summary
l2TP <- l1empty[l2idx > 0, ]
l2artr <- l2artrtmp |> filter(ID %notin% l2artrFNadd$ID) |>
  st_intersection(l2TP[,c('chm_max', 'area', 'uid','geometry')]) |> 
  mutate(TP = 1, FP = 0, FN = 0) |>
  st_centroid()

# summarize FN pts at L2
l2artrFN <- l2artrtmp |> 
  mutate(var = lengths(st_intersects(geometry, l1empty)), uid = NA) |> 
  filter(var == 0) |> select(-var) |> 
  bind_rows(l2artrFNadd) |> st_centroid() |>
  mutate(TP = 0, FP = 0, FN = 1, chm_max = NA, area = NA)

# summarize FP pols at L2
l2FP <- l1empty |> 
  filter(l2idx == 0) |> 
  mutate(var = lengths(st_intersects(geometry, plants |>
                               filter(!grepl("same", notes),
                                      site == !!site, Ht_gt_25 == 0,
                                      Species %in% c('ARTR', 'ARTRW', 'ARTRV')))) ) |>
  filter(var == 0) |> select(-var) |>
  mutate(ID = NA, Class = NA, Species = NA, Ht_gt_25 = NA, notes = NA, elev = NA,
         TP = 0, FP = 1, FN = 0) 
# add CID to each FP polygon; NA if polygon outside validated cells
l2FP |>
  mutate(CID = extract(rnull[[k]], l2FP |> st_centroid())[,2], 
         FP = ifelse(CID %in% c(rcell$cid[rcell$site == !!site & rcell$val == 1]), 
                     FP, NA),
         TP = ifelse(is.na(FP), NA, TP), FN = ifelse(is.na(FP), NA, FN)) |>
  st_intersection(rcell[,'geometry']) -> l2FP
# --- end of L2

# === combine sets

l2FP |>
  st_centroid() |>
  bind_rows(l1artr, l2artr) |>
  bind_rows(l1artrFN, l2artrFN)  -> ca

ca |> write_sf( paste0("../plantdat_artr_site/plantdat_artr_lt50_", site, ".geojson"), delete_dsn = TRUE )
ca |> 
  as.data.frame() |> select(-geometry) |> 
  write.csv( paste0("data/plantdat_artr_site/plantdat_artr_lt50_", site, ".csv"),row.names = FALSE)

}

# === Combine individual TP, FN, FP indices into cell-level counts
rcell <- st_read("data/celldat_covar.geojson")
ll <- list.files("data/plantdat_artr_site/", pattern = "_lt50_", full.names = TRUE)

map(ll, read.csv) |> 
  bind_rows() |> 
  filter(!is.na(CID)) |> 
  mutate(Unk = ifelse(is.na(TP) & is.na(FN) & is.na(FP), 1, 0)) |>
  select(site, TP, FN, FP, Unk, CID) |>
  group_by(site, CID) |> summarize_all(sum) |> ungroup() |>
  mutate(ucid = paste0(site, "_", CID)) -> cts

rcell |> 
  select(site, cid, ucid, val) |> 
  left_join(cts |> select(-site, -CID), by = 'ucid' ) |> 
  mutate(Unk = ifelse(val == 0 & is.na(Unk), 0, Unk),
         Unk = ifelse(val == 1 & is.na(Unk), 0, Unk), 
         TP = ifelse(val == 1 & is.na(TP), 0, TP), 
         FN = ifelse(val == 1 & is.na(FN), 0, FN), 
         FP = ifelse(val == 1 & is.na(FP), 0, FP)) -> out
  
out |> #write_sf("data/celldat_match_counts_lt50.geojson", delete_dsn = TRUE)  
  as.data.frame() |> select(-geometry) |> 
  write.csv("data/celldat_match_counts_lt50.csv", row.names = FALSE)


# === Visualize detection patterns
# === individual-level
ll <- list.files("data/plantdat_artr_site/", pattern = "lt50", full.names = TRUE)

map(ll, read.csv) |> 
  bind_rows() -> idf

idf |> 
  mutate(TP = as.factor(TP)) |>
  ggplot(aes(chm_max, colour = TP ) ) +
  stat_ecdf(geom = 'step', pad = FALSE, linewidth = 1.5) +
  theme_bw()

idf |> 
  mutate(TP = as.factor(TP)) |>
  ggplot(aes(x = TP, y = chm_max)) + 
  geom_violin() +
  theme_bw()

idf |> 
  filter(!is.na(TP), FP != 1) |> 
  mutate(ucid = paste0(site, "_", CID)) |>
  filter(ucid %notin% outs) -> fdf

1 - table(fdf$TP)[1]/table(fdf$TP)[2]


# === cell-level
cols <- viridis::viridis(10)
cts <- read.csv("data/celldat_match_counts.csv")

leg <- c("FP" = cols[3], "FN" = cols[8])
cts |> 
  ggplot() + 
  geom_jitter(aes(TP, FP, colour = "FP"), size = 1, alpha = .75) + 
  geom_jitter(aes(TP, FN, colour = "FN"), size = 1, alpha = .75) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(name = "Point", values = leg) +
  labs(y = "") +
  theme_bw()

cts |> 
  ggplot(aes(x = site, y = FN)) + geom_violin()

a <- cts |> filter(!is.na(TP))

table(a$site)

filter(cts, FN > 10) |> pull(ucid) -> outs



