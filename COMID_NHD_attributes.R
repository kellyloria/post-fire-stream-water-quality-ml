library(readr)
library(dplyr)
library(purrr)
library(nhdplusTools)

sites <- read_csv("/Users/kellyloria/Documents/DRI/Stream_WQ_predictions/data/USGS_sites_eligible_for_download.csv")

sites <- sites %>%
  mutate(
    COMID = map_int(site_no, \(x)
                    discover_nhdplus_id(
                      nldi_feature = list(
                        featureSource = "nwissite",
                        featureID = paste0("USGS-", x)
                      )
                    ))
  )


inputDir <- "/Users/kellyloria/Documents/UNR/Course work/Fall2020Projects/NHD_Tools/NHD_AttributeFiles"


# =============================================================================
# Retrieve COMIDs
# =============================================================================

# First try direct NLDI lookup using the USGS gage ID
dat <- sites %>%
  mutate(
    COMID_nwis = map_chr(
      site_no,
      get_comid_from_nwis
    )
  )


# Use coordinates when the direct USGS lookup does not return a COMID
dat <- dat %>%
  mutate(
    COMID_coordinates = pmap_chr(
      list(lon, lat),
      get_comid_from_coordinates
    ),
    COMID = coalesce(
      COMID_nwis,
      COMID_coordinates
    )
  )


# Examine COMID lookup success
dat %>%
  count(
    comid_found = !is.na(COMID)
  )


# Review sites that did not receive a COMID
sites_without_comid <- sites %>%
  filter(is.na(COMID)) 

print(sites_without_comid)



# Use coordinates when the direct USGS lookup does not return a COMID
dat_WCI <- sites_without_comid %>%
  mutate(
    COMID_nwis = pmap_chr(
      list(lon, lat),
      get_comid_from_coordinates),
    COMID_coordinates = pmap_chr(
      list(lon, lat),
      get_comid_from_coordinates
    ))


all_dat <- rbind(dat, dat_WCI)

all_dat <- all_dat %>%
  mutate(
    COMID = coalesce(COMID_nwis, COMID_coordinates),
    COMID = as.character(COMID)
  )

# =============================================================================
# Helper functions
# =============================================================================

# Read an NHDPlus attribute file while preserving original column names
read_nhd_attribute <- function(filename) {
  
  filepath <- file.path(inputDir, filename)
  
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  
  read.csv(
    filepath,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


# Find the first matching column name from a set of possible names
find_column <- function(data, possible_names) {
  
  matched <- possible_names[possible_names %in% names(data)]
  
  if (length(matched) == 0) {
    return(NA_character_)
  }
  
  matched[[1]]
}


# Convert COMID values to a consistent type before joining
standardize_comid <- function(data) {
  
  if (!"COMID" %in% names(data)) {
    stop("The supplied data frame does not contain a COMID column.")
  }
  
  data %>%
    mutate(COMID = as.character(COMID))
}


# Safely retrieve a COMID from an NWIS site
get_comid_from_nwis <- function(site_no) {
  
  site_no <- as.character(site_no)
  
  result <- tryCatch(
    {
      nhdplusTools::discover_nhdplus_id(
        nldi_feature = list(
          featureSource = "nwissite",
          featureID = paste0("USGS-", site_no)
        )
      )
    },
    error = function(e) {
      message(
        "Direct NWIS lookup failed for site ",
        site_no,
        ": ",
        conditionMessage(e)
      )
      
      NA
    }
  )
  
  # discover_nhdplus_id() may return a number, vector, or data frame,
  # depending on the package version and lookup method.
  if (length(result) == 0 || all(is.na(result))) {
    return(NA_character_)
  }
  
  if (is.data.frame(result)) {
    
    comid_col <- find_column(
      result,
      c("comid", "COMID", "nhdplus_comid", "identifier")
    )
    
    if (is.na(comid_col)) {
      return(NA_character_)
    }
    
    return(as.character(result[[comid_col]][1]))
  }
  
  as.character(result[[1]])
}


# Coordinate-based fallback for sites that fail the NWIS lookup
get_comid_from_coordinates <- function(lon, lat) {
  
  if (is.na(lon) || is.na(lat)) {
    return(NA_character_)
  }
  
  result <- tryCatch(
    {
      point <- sf::st_sfc(
        sf::st_point(c(as.numeric(lon), as.numeric(lat))),
        crs = 4326
      )
      
      nhdplusTools::discover_nhdplus_id(point)
    },
    error = function(e) {
      message(
        "Coordinate lookup failed at ",
        lon,
        ", ",
        lat,
        ": ",
        conditionMessage(e)
      )
      
      NA
    }
  )
  
  if (length(result) == 0 || all(is.na(result))) {
    return(NA_character_)
  }
  
  if (is.data.frame(result)) {
    
    comid_col <- find_column(
      result,
      c("comid", "COMID", "nhdplus_comid", "identifier")
    )
    
    if (is.na(comid_col)) {
      return(NA_character_)
    }
    
    return(as.character(result[[comid_col]][1]))
  }
  
  as.character(result[[1]])
}


library(tidyverse)

inputDir <- paste0(
  "/Users/kellyloria/Documents/UNR/Course work/",
  "Fall2020Projects/NHD_Tools/NHD_AttributeFiles"
)

dat <- all_dat %>%
  mutate(
    site_no = as.character(site_no),
    COMID = coalesce(
      as.character(COMID_nwis),
      as.character(COMID_coordinates)
    )
  )


dat <- dat %>%
  mutate(
    COMID_source = case_when(
      !is.na(COMID_nwis) ~ "NWIS site lookup",
      is.na(COMID_nwis) & !is.na(COMID_coordinates) ~ "Coordinate lookup",
      TRUE ~ "No COMID found"
    )
  )

dat %>%
  count(COMID_source)

# =============================================================================
# Olson geology
# =============================================================================

olson_attributes <- read.csv(
  file.path(inputDir, "OLSON_CAT_CONUS.TXT"),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    COMID = as.character(COMID),
    
    olson_compressive_strength_MPa =
      as.numeric(CAT_OLSON_UCS),
    
    olson_hydraulic_conductivity_um_s =
      as.numeric(CAT_OLSON_PERM)
  ) %>%
  distinct(COMID, .keep_all = TRUE)


# =============================================================================
# Generalized geology
# =============================================================================

geology_attributes <- read.csv(
  file.path(inputDir, "BUSHREED_ACC_CONUS.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    COMID = as.character(COMID),
    
    GGT_granitic_pct =
      as.numeric(ACC_BUSHREED2),
    
    GGT_sedimentary_pct =
      as.numeric(ACC_BUSHREED5),
    
    GGT_volcanic_pct =
      as.numeric(ACC_BUSHREED6),
    
    GGT_water_pct =
      as.numeric(ACC_BUSHREED8)
  ) %>%
  distinct(COMID, .keep_all = TRUE)


# =============================================================================
# Precipitation
# =============================================================================

precip_attributes <- read.csv(
  file.path(inputDir, "PPT2015_ANN_CONUS.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    COMID = as.character(COMID),
    
    local_mean_annual_precip_mm =
      as.numeric(CAT_PPT2015_ANN),
    
    upstream_mean_annual_precip_mm =
      as.numeric(ACC_PPT2015_ANN),
    
    total_upstream_annual_precip_mm =
      as.numeric(TOT_PPT2015_ANN)
  ) %>%
  distinct(COMID, .keep_all = TRUE)


# =============================================================================
# Basin characteristics
# =============================================================================

basin_attributes <- read.csv(
  file.path(inputDir, "BASIN_CHAR_CAT_CONUS.TXT"),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    COMID = as.character(COMID),
    
    local_catchment_area_km2 =
      as.numeric(CAT_BASIN_AREA),
    
    local_catchment_slope_pct =
      as.numeric(CAT_BASIN_SLOPE),
    
    local_mean_elevation_m =
      as.numeric(CAT_ELEV_MEAN),
    
    local_min_elevation_m =
      as.numeric(CAT_ELEV_MIN),
    
    local_max_elevation_m =
      as.numeric(CAT_ELEV_MAX),
    
    stream_slope_pct =
      as.numeric(CAT_STREAM_SLOPE),
    
    stream_length_km =
      as.numeric(CAT_STREAM_LENGTH)
  ) %>%
  distinct(COMID, .keep_all = TRUE)


# =============================================================================
# Bedrock permeability
# =============================================================================

bedperm_attributes <- read.csv(
  file.path(inputDir, "BEDPERM_ACC_CONUS.TXT"),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    COMID = as.character(COMID),
    
    upstream_not_principal_aquifer_pct =
      as.numeric(ACC_BEDPERM_1),
    
    upstream_sand_gravel_aquifer_pct =
      as.numeric(ACC_BEDPERM_6)
  ) %>%
  distinct(COMID, .keep_all = TRUE)


# =============================================================================
# Baseflow index
# =============================================================================

bfi_raw <- read.csv(
  file.path(inputDir, "BFI_CONUS.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Inspect the fields if needed
names(bfi_raw)

bfi_attributes <- bfi_raw %>%
  mutate(
    COMID = as.character(COMID)
  ) %>%
  select(
    COMID,
    matches("BFI", ignore.case = TRUE)
  ) %>%
  rename_with(
    ~ paste0("bfi_", .x),
    -COMID
  ) %>%
  distinct(COMID, .keep_all = TRUE)


# =============================================================================
# Canopy cover within 100-m stream buffer
# =============================================================================

canopy_raw <- read.csv(
  file.path(inputDir, "CNPY11_BUFF100_CONUS.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Inspect the fields if needed
names(canopy_raw)

canopy_attributes <- canopy_raw %>%
  mutate(
    COMID = as.character(COMID)
  ) %>%
  select(
    COMID,
    matches(
      "CNPY|CANOPY|BUFF",
      ignore.case = TRUE
    )
  ) %>%
  rename_with(
    ~ paste0("canopy_", .x),
    -COMID
  ) %>%
  distinct(COMID, .keep_all = TRUE)


# =============================================================================
# NLCD 2016 accumulated upstream land cover
# =============================================================================

nlcd_raw <- read.csv(
  file.path(inputDir, "NLCD16_ACC_CONUS.TXT"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Inspect the fields if needed
names(nlcd_raw)

nlcd_attributes <- nlcd_raw %>%
  mutate(
    COMID = as.character(COMID)
  ) %>%
  select(
    COMID,
    starts_with("ACC_")
  ) %>%
  rename_with(
    ~ paste0("nlcd16_", .x),
    -COMID
  ) %>%
  distinct(COMID, .keep_all = TRUE)


# =============================================================================
# Road density
# =============================================================================

road_raw <- read.csv(
  file.path(inputDir, "ROAD_DENS_CAT_CONUS.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Inspect the fields if needed
names(road_raw)

road_attributes <- road_raw %>%
  mutate(
    COMID = as.character(COMID)
  ) %>%
  select(
    COMID,
    matches(
      "ROAD|DENS",
      ignore.case = TRUE
    )
  ) %>%
  rename_with(
    ~ paste0("road_", .x),
    -COMID
  ) %>%
  distinct(COMID, .keep_all = TRUE)


# Road density ---------------------------------------------------------------

road_raw <- read.csv(
  file.path(inputDir, "ROAD_DENS_CAT_CONUS.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Remove columns with blank or NA names
road_raw <- road_raw[
  ,
  !is.na(names(road_raw)) & names(road_raw) != "",
  drop = FALSE
]

road_attributes <- road_raw %>%
  transmute(
    COMID = as.character(COMID),
    road_density_total =
      as.numeric(CAT_TOTAL_ROAD_DENS)
  ) %>%
  distinct(COMID, .keep_all = TRUE)

# =============================================================================
# Merge all attributes
# =============================================================================

all_dat_attributes <- dat %>%
  left_join(
    olson_attributes,
    by = "COMID",
    relationship = "many-to-one"
  ) %>%
  left_join(
    geology_attributes,
    by = "COMID",
    relationship = "many-to-one"
  ) %>%
  left_join(
    precip_attributes,
    by = "COMID",
    relationship = "many-to-one"
  ) %>%
  left_join(
    basin_attributes,
    by = "COMID",
    relationship = "many-to-one"
  ) %>%
  left_join(
    bedperm_attributes,
    by = "COMID",
    relationship = "many-to-one"
  ) %>%
  left_join(
    bfi_attributes,
    by = "COMID",
    relationship = "many-to-one"
  ) %>%
  left_join(
    canopy_attributes,
    by = "COMID",
    relationship = "many-to-one"
  ) %>%
  left_join(
    nlcd_attributes,
    by = "COMID",
    relationship = "many-to-one"
  ) %>%
  left_join(
    road_attributes,
    by = "COMID",
    relationship = "many-to-one"
  )


# =============================================================================
# Inspect result
# =============================================================================

glimpse(all_dat_attributes)

print(
  all_dat_attributes,
  n = Inf,
  width = Inf
)


  # =============================================================================
# Identify representative fields for match diagnostics
# =============================================================================

bfi_check_col <- names(all_dat_attributes)[
  str_detect(names(all_dat_attributes), "^bfi_")
][1]

canopy_check_col <- names(all_dat_attributes)[
  str_detect(names(all_dat_attributes), "^canopy_")
][1]

nlcd_check_col <- names(all_dat_attributes)[
  str_detect(names(all_dat_attributes), "^nlcd16_")
][1]

road_check_col <- names(all_dat_attributes)[
  str_detect(names(all_dat_attributes), "^road_")
][1]


# =============================================================================
# Match summary
# =============================================================================

attribute_match_summary <- all_dat_attributes %>%
  summarise(
    total_sites = n(),
    
    sites_with_COMID =
      sum(!is.na(COMID)),
    
    sites_without_COMID =
      sum(is.na(COMID)),
    
    olson_matches =
      sum(!is.na(olson_compressive_strength_MPa)),
    
    geology_matches =
      sum(!is.na(GGT_granitic_pct)),
    
    precipitation_matches =
      sum(!is.na(local_mean_annual_precip_mm)),
    
    basin_matches =
      sum(!is.na(local_catchment_area_km2)),
    
    bedrock_matches =
      sum(!is.na(upstream_not_principal_aquifer_pct)),
    
    BFI_matches =
      sum(!is.na(.data[[bfi_check_col]])),
    
    canopy_matches =
      sum(!is.na(.data[[canopy_check_col]])),
    
    NLCD_matches =
      sum(!is.na(.data[[nlcd_check_col]])),
    
    road_density_matches =
      sum(!is.na(.data[[road_check_col]]))
  )

attribute_match_summary


# =============================================================================
# Rename NHDPlus attributes and create summary variables
# =============================================================================

all_dat_attributes <- all_dat_attributes %>%
  rename(
    
    # ---------- Baseflow Index ----------
    local_BFI = bfi_CAT_BFI,
    upstream_BFI = bfi_ACC_BFI,
    upstream_BFI_total = bfi_TOT_BFI,
    
    # ---------- Riparian Canopy (% within 100 m stream buffer) ----------
    local_canopy_pct = canopy_CAT_CNPY11_BUFF100,
    upstream_canopy_pct = canopy_ACC_CNPY11_BUFF100,
    upstream_canopy_total_pct = canopy_TOT_CNPY11_BUFF100,
    
    # ---------- NLCD 2016 Land Cover ----------
    upstream_open_water_pct          = nlcd16_ACC_NLCD16_11,
    upstream_ice_snow_pct            = nlcd16_ACC_NLCD16_12,
    upstream_developed_open_pct      = nlcd16_ACC_NLCD16_21,
    upstream_developed_low_pct       = nlcd16_ACC_NLCD16_22,
    upstream_developed_medium_pct    = nlcd16_ACC_NLCD16_23,
    upstream_developed_high_pct      = nlcd16_ACC_NLCD16_24,
    upstream_barren_pct              = nlcd16_ACC_NLCD16_31,
    upstream_deciduous_forest_pct    = nlcd16_ACC_NLCD16_41,
    upstream_evergreen_forest_pct    = nlcd16_ACC_NLCD16_42,
    upstream_mixed_forest_pct        = nlcd16_ACC_NLCD16_43,
    upstream_shrub_pct               = nlcd16_ACC_NLCD16_52,
    upstream_grassland_pct           = nlcd16_ACC_NLCD16_71,
    upstream_pasture_pct             = nlcd16_ACC_NLCD16_81,
    upstream_crops_pct               = nlcd16_ACC_NLCD16_82,
    upstream_woody_wetland_pct       = nlcd16_ACC_NLCD16_90,
    upstream_emergent_wetland_pct    = nlcd16_ACC_NLCD16_95,
    upstream_nodata_pct              = nlcd16_ACC_NODATA) %>%
  mutate(
    
    # ---------- Developed land ----------
    upstream_developed_pct =
      upstream_developed_open_pct +
      upstream_developed_low_pct +
      upstream_developed_medium_pct +
      upstream_developed_high_pct,
    
    # ---------- Forest ----------
    upstream_forest_pct =
      upstream_deciduous_forest_pct +
      upstream_evergreen_forest_pct +
      upstream_mixed_forest_pct,
    
    # ---------- Agriculture ----------
    upstream_agriculture_pct =
      upstream_pasture_pct +
      upstream_crops_pct,
    
    # ---------- Wetlands ----------
    upstream_wetland_pct =
      upstream_woody_wetland_pct +
      upstream_emergent_wetland_pct,
    
    # ---------- Open natural landscapes ----------
    upstream_open_natural_pct =
      upstream_barren_pct +
      upstream_shrub_pct +
      upstream_grassland_pct,
    
    # ---------- Anthropogenic land ----------
    upstream_human_land_pct =
      upstream_developed_pct +
      upstream_agriculture_pct,
    
    # ---------- Vegetated land ----------
    upstream_vegetated_pct =
      upstream_forest_pct +
      upstream_shrub_pct +
      upstream_grassland_pct +
      upstream_pasture_pct,
    
    # ---------- Aquatic area ----------
    upstream_water_wetland_pct =
      upstream_open_water_pct +
      upstream_wetland_pct,
    
    # ---------- Check NLCD percentages ----------
    upstream_landcover_sum =
      upstream_open_water_pct +
      upstream_ice_snow_pct +
      upstream_developed_pct +
      upstream_barren_pct +
      upstream_forest_pct +
      upstream_shrub_pct +
      upstream_grassland_pct +
      upstream_agriculture_pct +
      upstream_wetland_pct +
      upstream_nodata_pct
  )

### quick plot 

library(tidyverse)

landcover_vars <- c(
  "upstream_open_water_pct",
  "upstream_developed_pct",
  "upstream_barren_pct",
  "upstream_forest_pct",
  "upstream_shrub_pct",
  "upstream_grassland_pct",
  "upstream_agriculture_pct",
  "upstream_wetland_pct"
)

landcover_pca_data <- all_dat_attributes %>%
  select(site_no, site_name, all_of(landcover_vars)) %>%
  drop_na()

landcover_pca <- prcomp(
  landcover_pca_data %>%
    select(-site_no, -site_name),
  center = TRUE,
  scale. = TRUE
)


library(tidyverse)

landcover_pca_data <- all_dat_attributes %>%
  select(site_no, site_name, all_of(landcover_vars)) %>%
  drop_na()

landcover_pca <- prcomp(
  landcover_pca_data %>%
    select(-site_no, -site_name),
  center = TRUE,
  scale. = TRUE
)

summary(landcover_pca)

variance_explained <- tibble(
  PC = paste0("PC", 1:length(landcover_pca$sdev)),
  Variance = landcover_pca$sdev^2 /
    sum(landcover_pca$sdev^2)
)

variance_explained


loadings <- landcover_pca$rotation %>%
  as.data.frame() %>%
  rownames_to_column("Variable")

loadings

site_scores <- landcover_pca$x %>%
  as.data.frame() %>%
  mutate(
    site_no = landcover_pca_data$site_no,
    site_name = landcover_pca_data$site_name
  )

all_dat_attributes <- left_join(
  all_dat_attributes,
  site_scores,
  by = c("site_no", "site_name")
)

library(ggplot2)

ggplot(site_scores,
       aes(PC1, PC2,
           label = site_no)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3) +
  theme_classic()