
# =============================================================================
# USGS_uv_read_and_bind.R
#
# Purpose:
#   Read all site-level Parquet files from a folder and combine them by rows
#   into a single sub-daily time-series dataset.
#
# Expected columns:
#   site_no
#   dateTime
#   Flow_Inst
#   DO_Inst
#   Wtemp_Inst
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load packages
# -----------------------------------------------------------------------------

library(arrow)
library(dplyr)
library(purrr)


# -----------------------------------------------------------------------------
# 2. Specify the folder containing the Parquet files
# -----------------------------------------------------------------------------

# Use the relative path below when uv_data_by_site is inside the current
# R working directory.

UV_DATA_DIR <- "/Users/kellyloria/Documents/DRI/Stream_WQ_predictions/post-fire-stream-water-quality-ml/uv_data_by_site/"


# Alternatively, use the complete folder path:
#
# UV_DATA_DIR <- paste0(
#   "/Users/kellyloria/Documents/DRI/",
#   "USGS_subdaily_analysis/uv_data_by_site"
# )


# Display the current working directory
getwd()


# Display the full path R will use
normalizePath(
  UV_DATA_DIR,
  mustWork = FALSE
)


# -----------------------------------------------------------------------------
# 3. Find all Parquet files
# -----------------------------------------------------------------------------

parquet_files <- list.files(
  path = UV_DATA_DIR,
  pattern = "\\.parquet$",
  full.names = TRUE
)


if (length(parquet_files) == 0) {
  stop(
    "No Parquet files were found in: ",
    normalizePath(UV_DATA_DIR, mustWork = FALSE)
  )
}


message(
  "Number of Parquet files found: ",
  length(parquet_files)
)


print(
  head(parquet_files)
)


# -----------------------------------------------------------------------------
# 4. Read each Parquet file and bind all rows
# -----------------------------------------------------------------------------

# map_dfr() reads each file and binds the resulting data frames by rows.

uv_all <- parquet_files %>%
  set_names(basename(.)) %>%
  map_dfr(
    function(file_path) {
      
      message(
        "Reading: ",
        basename(file_path)
      )
      
      read_parquet(
        file_path
      )
    },
    .id = "source_file"
  )


# -----------------------------------------------------------------------------
# 5. Standardize and arrange the combined time series
# -----------------------------------------------------------------------------

uv_all <- uv_all %>%
  mutate(
    site_no = as.character(site_no),
    
    dateTime = as.POSIXct(
      dateTime,
      tz = "UTC"
    )
  ) %>%
  arrange(
    site_no,
    dateTime
  )


# -----------------------------------------------------------------------------
# 6. Review the combined dataset
# -----------------------------------------------------------------------------

print(
  uv_all
)


glimpse(
  uv_all
)


message(
  "Total rows: ",
  format(
    nrow(uv_all),
    big.mark = ","
  )
)


message(
  "Total sites: ",
  n_distinct(uv_all$site_no)
)


message(
  "First observation: ",
  min(
    uv_all$dateTime,
    na.rm = TRUE
  )
)


message(
  "Last observation: ",
  max(
    uv_all$dateTime,
    na.rm = TRUE
  )
)


# -----------------------------------------------------------------------------
# 7. Optional: remove the source-file column
# -----------------------------------------------------------------------------

# The source_file column identifies the Parquet file from which each row came.
# Remove it if it is not needed.

# uv_all <- uv_all %>%
#   select(-source_file)


# -----------------------------------------------------------------------------
# 8. Optional: save the combined time-series dataset
# -----------------------------------------------------------------------------

# Saving as Parquet is much more efficient than saving a very large CSV.

write_rds(
  uv_all,
  "USGS_uv_all_sites.rds"
)


# A CSV can be created, but it may be very large and slow to read and write.

# readr::write_csv(
#   uv_all,
#   "USGS_uv_all_sites.csv"
# )

unique(uv_all$site_no)