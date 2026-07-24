
# =============================================================================
# USGS_uv_mean_temp_do_map.R
#
# Purpose:
#   1. Read the aggregated USGS unit-value dataset
#   2. Calculate mean water temperature and dissolved oxygen by site
#   3. Join site coordinates
#   4. Map:
#        - Point color = mean water temperature
#        - Point size  = mean dissolved oxygen
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load packages
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(maps)
library(readr)
library(viridis)


# -----------------------------------------------------------------------------
# 2. File paths
# -----------------------------------------------------------------------------

UV_RDS_FILE <- paste0(
  "/Users/kellyloria/Documents/DRI/Stream_WQ_predictions/",
  "post-fire-stream-water-quality-ml/USGS_uv_all_sites.rds"
)

# This CSV should contain:
# site_no, site_name, state_cd, lat, and lon
SITE_METADATA_FILE <- paste0(
  "/Users/kellyloria/Documents/DRI/Stream_WQ_predictions/",
  "post-fire-stream-water-quality-ml/USGS_qualifying_sites.csv"
)

OUTPUT_FILE <- paste0(
  "/Users/kellyloria/Documents/DRI/Stream_WQ_predictions/",
  "post-fire-stream-water-quality-ml/",
  "USGS_mean_temperature_DO_map.png"
)


# -----------------------------------------------------------------------------
# 3. Confirm that input files exist
# -----------------------------------------------------------------------------

if (!file.exists(UV_RDS_FILE)) {
  stop(
    "Aggregated UV data file was not found:\n",
    UV_RDS_FILE
  )
}

if (!file.exists(SITE_METADATA_FILE)) {
  stop(
    "Site metadata file was not found:\n",
    SITE_METADATA_FILE,
    "\n\nThe metadata file is needed to provide site latitude and longitude."
  )
}


# -----------------------------------------------------------------------------
# 4. Read aggregated unit-value data
# -----------------------------------------------------------------------------

uv_all <- readRDS(UV_RDS_FILE)

message(
  "Rows read from aggregated UV dataset: ",
  format(nrow(uv_all), big.mark = ",")
)

message(
  "Columns in aggregated UV dataset: ",
  paste(names(uv_all), collapse = ", ")
)


# -----------------------------------------------------------------------------
# 5. Check required water-quality columns
# -----------------------------------------------------------------------------

required_uv_columns <- c(
  "site_no",
  "Wtemp_Inst",
  "DO_Inst"
)

missing_uv_columns <- setdiff(
  required_uv_columns,
  names(uv_all)
)

if (length(missing_uv_columns) > 0) {
  stop(
    "The aggregated UV dataset is missing required columns: ",
    paste(missing_uv_columns, collapse = ", ")
  )
}


# -----------------------------------------------------------------------------
# 6. Calculate mean temperature and DO by site
# -----------------------------------------------------------------------------

site_water_quality <- uv_all %>%
  mutate(
    site_no = as.character(site_no),
    Wtemp_Inst = as.numeric(Wtemp_Inst),
    DO_Inst = as.numeric(DO_Inst)
  ) %>%
  group_by(site_no) %>%
  summarise(
    n_records = n(),
    
    n_temp = sum(
      !is.na(Wtemp_Inst)
    ),
    
    n_do = sum(
      !is.na(DO_Inst)
    ),
    
    mean_temp_c = if_else(
      n_temp > 0,
      mean(Wtemp_Inst, na.rm = TRUE),
      NA_real_
    ),
    
    mean_do_mg_l = if_else(
      n_do > 0,
      mean(DO_Inst, na.rm = TRUE),
      NA_real_
    ),
    
    .groups = "drop"
  )

message(
  "Sites summarized: ",
  nrow(site_water_quality)
)


# -----------------------------------------------------------------------------
# 7. Read site metadata
# -----------------------------------------------------------------------------

site_metadata <- read_csv(
  SITE_METADATA_FILE,
  col_types = cols(
    site_no = col_character(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

required_metadata_columns <- c(
  "site_no",
  "lat",
  "lon"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  names(site_metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "The site metadata file is missing required columns: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}


# Keep one metadata row per site
site_metadata <- site_metadata %>%
  mutate(
    site_no = as.character(site_no)
  ) %>%
  distinct(
    site_no,
    .keep_all = TRUE
  )


# -----------------------------------------------------------------------------
# 8. Join site means to coordinates
# -----------------------------------------------------------------------------

map_sites <- site_water_quality %>%
  left_join(
    site_metadata,
    by = "site_no"
  ) %>%
  filter(
    !is.na(lon),
    !is.na(lat),
    !is.na(mean_temp_c),
    !is.na(mean_do_mg_l)
  )

message(
  "Sites available for mapping: ",
  nrow(map_sites)
)

if (nrow(map_sites) == 0) {
  stop(
    "No sites have coordinates plus non-missing mean temperature and DO."
  )
}


# -----------------------------------------------------------------------------
# 9. Create western-state background map
# -----------------------------------------------------------------------------

western_states <- c(
  "california",
  "washington",
  "oregon",
  #"utah",
  #"colorado",
  "idaho",
  "nevada"
 # "arizona"
 # "new mexico",
  #"montana",
  #"wyoming"
)

western_map <- map_data("state") %>%
  filter(
    region %in% western_states
  )


# -----------------------------------------------------------------------------
# 10. Create map
# -----------------------------------------------------------------------------

mean_temp_do_map <- ggplot() +
  
  geom_polygon(
    data = western_map,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "grey95",
    color = "grey60",
    linewidth = 0.3
  ) +
  
  geom_point(
    data = map_sites,
    aes(
      x = lon,
      y = lat,
      color = mean_temp_c,
      size = mean_do_mg_l
    ),
    alpha = 0.8
  ) +
  
  scale_color_viridis_c(
    option = "plasma",
    direction = 1,
    name = expression(
      "Mean temperature (" * degree * "C)"
    )
  ) +
  
  scale_size_continuous(
    name = expression(
      "Mean DO (mg L"^{-1} * ")"
    ),
    range = c(2, 8)
  ) +
  
  coord_fixed(
    ratio = 1.3,
    xlim = c(-125, -102),
    ylim = c(31, 50),
    expand = FALSE
  ) +
  
  theme_minimal() +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  
  labs(
    title = "Mean water temperature and dissolved oxygen at USGS gages",
    subtitle = paste0(
      "Point color represents mean temperature; ",
      "point size represents mean dissolved oxygen; n = ",
      nrow(map_sites),
      " sites"
    ),
    x = "Longitude",
    y = "Latitude"
  )


print(mean_temp_do_map)


# -----------------------------------------------------------------------------
# 11. Save map
# -----------------------------------------------------------------------------

ggsave(
  filename = ,
  plot = mean_temp_do_map,
  width = 10,
  height = 7,
  dpi = 300
)

message(
  "Map saved to:\n",
  OUTPUT_FILE
)


# -----------------------------------------------------------------------------
# 12. Optional: save the site-level values used in the map
# -----------------------------------------------------------------------------

MAP_DATA_FILE <- paste0(
  "/Users/kellyloria/Documents/DRI/Stream_WQ_predictions/",
  "post-fire-stream-water-quality-ml/",
  "USGS_site_mean_temperature_DO.csv"
)

write_csv(
  map_sites,
  MAP_DATA_FILE
)

message(
  "Site-level map data saved to:\n",
  MAP_DATA_FILE
)
