pkgs <- c("tidyverse", "sf")
sapply(pkgs, require, character.only = TRUE)

# === Merge classified data with the geospatial layer
preds <- read.csv("~/../../Volumes/az_drive/FieldData/FieldData_2022/ML_classification_segment/predicted_plants_v1.csv")
dfg <- st_read("~/../../Volumes/az_drive/FieldData/FieldData_2022/ML_classification_segment/ml_dataset_v1.geojson")
# - use unique id to merge: uid = site + ID [ID as in treeID]
mldf <- read.csv("~/../../Volumes/az_drive/FieldData/FieldData_2022/ML_classification_segment/ml_dataset_v1.csv") #|> mutate(uid = paste0(site, "_", ID))

# ---
dfg |> 
  left_join(preds, by = "uid") -> out

fs <- unique(out$site)
lapply(1:length(fs), function(x) {
  tmp <- filter(out, site == fs[x])
  st_write(tmp, paste0("~/../../Volumes/az_drive/FieldData/FieldData_2022/ML_classification_segment/predicted_plants_v1_sites/predicted_plants_v1_", fs[x], ".geojson"))
})

st_write(out, "~/../../Volumes/az_drive/FieldData/FieldData_2022/ML_classification_segment/predicted_plants_v1.geojson")
# === end of merge
