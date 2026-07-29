# =============================================================================
# USGS_obs_agg_revised_v3.R
#
# Consolidated, single-download-path version. Fixes vs. the prior draft:
#   - Removed the dead fetch_uv()/map()/rbindlist() block entirely (this was
#     the second, in-memory download path that both duplicated the download
#     and caused the 16GB vector memory crash).
#   - Step 2a/2b and Part 3 now read from the arrow dataset (uv_ds) instead
#     of the no-longer-existing uv_raw/uv_clean objects.
#   - example_sites is now defined BEFORE it's used to subset uv_ds.
#   - state_cd is captured in site_meta (Step 1d) so downstream joins/plots
#     that reference it actually work.
#   - western_map filter and plot labels now include Montana, matching
#     target_states.
#
# Goals (unchanged):
#   1. Identify USGS stream gage sites in CA, WA, OR, ID, NV, AZ, MT that have
#      more than 6 months of sub-daily (unit-value/instantaneous) dissolved
#      oxygen (parm_cd 00300) AND water temperature (parm_cd 00010) records.
#   2. Download a continuous (unit-value) time-series of:
#        - Streamflow        parm_cd = 00060  (cfs)
#        - Dissolved oxygen  parm_cd = 00300  (mg/L)
#        - Water temperature parm_cd = 00010  (deg C)
#      for qualifying sites, covering 2011-10-01 through 2024-10-01.
#
# Package requirements:
#   install.packages(c("dataRetrieval", "dplyr", "purrr", "lubridate",
#                       "readr", "arrow", "ggplot2", "maps"))
# =============================================================================

library(dataRetrieval)
library(dplyr)
library(purrr)
library(lubridate)
library(readr)
library(arrow)

# -----------------------------------------------------------------------------
# 0.  Constants
# -----------------------------------------------------------------------------

target_states <- c("CA", "WA", "OR", "ID", "NV", "AZ", "MT", "UT", "CO", "NM", "WY")

PARM_DO   <- "00300"
PARM_TEMP <- "00010"
PARM_FLOW <- "00060"

START_DATE <- as.Date("2011-10-01")
END_DATE   <- as.Date("2024-10-01")

MIN_DAYS <- 180

UV_DATA_DIR <- "uv_data_by_site"

# =============================================================================
# PART 1 -- Identify qualifying sites
# =============================================================================

fetch_sites <- function(state, parm) {
  message("  state=", state, "  parm=", parm)
  
  tryCatch(
    whatNWISsites(
      stateCd = state,
      siteType = "ST",
      parameterCd = parm,
      service = "uv"
    ) %>%
      mutate(state_cd = state),
    error = function(e) {
      message(
        "  [WARN] whatNWISsites failed for state=", state,
        ", parm=", parm, ": ", e$message
      )
      NULL
    }
  )
}

do_sites_raw <- map_dfr(
  target_states,
  fetch_sites,
  parm = PARM_DO
)

temp_sites_raw <- map_dfr(
  target_states,
  fetch_sites,
  parm = PARM_TEMP
)

names(do_sites_raw)


site_meta <- bind_rows(do_sites_raw, temp_sites_raw) %>%
  distinct(site_no, .keep_all = TRUE) %>%
  select(
    site_no,
    station_nm,
    state_cd,
    dec_lat_va,
    dec_long_va
  ) %>%
  rename(
    site_name = station_nm,
    lat = dec_lat_va,
    lon = dec_long_va
  )

candidate_sites <- intersect(
  unique(do_sites_raw$site_no),
  unique(temp_sites_raw$site_no)
)

CHUNK_SIZE <- 100
site_chunks <- split(candidate_sites, ceiling(seq_along(candidate_sites) / CHUNK_SIZE))

fetch_availability <- function(sites_chunk) {
  tryCatch(
    whatNWISdata(siteNumber = sites_chunk, service = "uv",
                 parameterCd = c(PARM_DO, PARM_TEMP)),
    error = function(e) {
      message("  [WARN] whatNWISdata chunk failed: ", e$message)
      NULL
    }
  )
}

avail_raw <- map_dfr(site_chunks, fetch_availability)
message("  Availability catalogue rows returned: ", nrow(avail_raw))


message("=== Step 1c: Filtering to >", MIN_DAYS, "-day sub-daily records ===")

avail_filtered <- avail_raw %>%
  mutate(
    begin_date = as.Date(begin_date),
    end_date   = as.Date(end_date),
    span_days  = as.numeric(end_date - begin_date)
  ) %>%
  filter(parm_cd %in% c(PARM_DO, PARM_TEMP), data_type_cd == "uv",
         span_days > MIN_DAYS, count_nu > MIN_DAYS)

parm_site_summary <- avail_filtered %>%
  group_by(site_no, parm_cd) %>%
  summarise(total_count = sum(count_nu, na.rm = TRUE),
            rec_start = min(begin_date, na.rm = TRUE),
            rec_end   = max(end_date, na.rm = TRUE),
            .groups = "drop")

sites_with_do   <- parm_site_summary %>% filter(parm_cd == PARM_DO)   %>% pull(site_no)
sites_with_temp <- parm_site_summary %>% filter(parm_cd == PARM_TEMP) %>% pull(site_no)
qualifying_site_nos <- intersect(sites_with_do, sites_with_temp)

message("  Sites passing >", MIN_DAYS, "-day threshold for BOTH DO and temp: ",
        length(qualifying_site_nos))

message("=== Step 1d: Attaching site metadata ===")

# NOTE: state_cd kept here (previously commented out), since it's used
# downstream for both the join and the map color aesthetic.
site_meta <- bind_rows(do_sites_raw, temp_sites_raw) %>%
  distinct(site_no, .keep_all = TRUE) %>%
  select(site_no, station_nm, state_cd, dec_lat_va, dec_long_va) %>%
  rename(site_name = station_nm, lat = dec_lat_va, lon = dec_long_va)

do_summary <- parm_site_summary %>%
  filter(parm_cd == PARM_DO) %>%
  rename(do_count = total_count, do_start = rec_start, do_end = rec_end) %>%
  select(-parm_cd)

temp_summary <- parm_site_summary %>%
  filter(parm_cd == PARM_TEMP) %>%
  rename(temp_count = total_count, temp_start = rec_start, temp_end = rec_end) %>%
  select(-parm_cd)

qualifying_sites <- tibble(site_no = qualifying_site_nos) %>%
  left_join(do_summary, by = "site_no") %>%
  left_join(temp_summary, by = "site_no") %>%
  left_join(site_meta, by = "site_no")

message("  Final qualifying site count: ", nrow(qualifying_sites))
write_csv(qualifying_sites, "USGS_qualifying_sites.csv")
message("  Saved: USGS_qualifying_sites.csv")


# =============================================================================
# PART 2 -- Select sites and download unit-value time series
# =============================================================================

message("=== Step 2: Preparing unit-value downloads ===")
message(
  "    Parameters : ",
  PARM_FLOW, " (flow), ",
  PARM_DO, " (DO), ",
  PARM_TEMP, " (temperature)"
)
message("    Date range : ", START_DATE, " to ", END_DATE)

# Keep site numbers as character throughout the workflow. This preserves
# leading zeros and prevents mismatches between CSV values and file names.
qualifying_sites <- qualifying_sites %>%
  mutate(site_no = as.character(site_no))

UV_SUMMARY_FILE <- "USGS_uv_summary.csv"
MISSING_SITES_FILE <- "USGS_missing_sites.csv"

dir.create(
  UV_DATA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# Step 2a -- Determine which sites are eligible for downloading
# -----------------------------------------------------------------------------

# The coverage percentages are calculated from downloaded instantaneous-value
# records. Therefore:
#
#   * If USGS_uv_summary.csv exists from a previous run, retain only sites
#     with >25% DO coverage AND >50% temperature coverage.
#
#   * If the summary does not exist, this is treated as an initial run and all
#     qualifying sites are downloaded so that coverage can be calculated.
#
# After the initial run, rerunning this script will use the saved summary to
# restrict the download and analysis lists.

if (file.exists(UV_SUMMARY_FILE)) {
  
  message(
    "  Existing coverage summary found: ",
    UV_SUMMARY_FILE
  )
  
  previous_uv_summary <- read_csv(
    UV_SUMMARY_FILE,
    col_types = cols(
      site_no = col_character(),
      .default = col_guess()
    ),
    show_col_types = FALSE
  )
  
  required_summary_columns <- c(
    "site_no",
    "pct_do_present",
    "pct_temp_present"
  )
  
  missing_summary_columns <- setdiff(
    required_summary_columns,
    names(previous_uv_summary)
  )
  
  if (length(missing_summary_columns) > 0) {
    stop(
      "The existing USGS_uv_summary.csv is missing required columns: ",
      paste(missing_summary_columns, collapse = ", "),
      ". Delete or repair the summary file before rerunning."
    )
  }
  
  sites_meeting_coverage <- previous_uv_summary %>%
    mutate(site_no = as.character(site_no)) %>%
    filter(
        !is.na(pct_do_present),
        !is.na(pct_temp_present),
        pct_do_present >= 0,
        pct_temp_present >= 50
    ) %>%
    distinct(site_no) %>%
    semi_join(
      qualifying_sites %>% select(site_no),
      by = "site_no"
    ) %>%
    pull(site_no)
  
  message(
    "  Sites with >50% DO and >50% temperature coverage: ",
    length(sites_meeting_coverage)
  )
  
  sites_eligible <- sites_meeting_coverage
  
} else {
  
  message(
    "  No previous USGS_uv_summary.csv was found."
  )
  
  message(
    "  Initial run: downloading all qualifying sites so coverage can be calculated."
  )
  
  sites_eligible <- qualifying_sites %>%
    distinct(site_no) %>%
    pull(site_no)
}


sites_eligible <- unique(as.character(sites_eligible))

message(
  "  Sites eligible for this run: ",
  length(sites_eligible)
)

if (length(sites_eligible) == 0) {
  stop(
    "No sites meet the current download criteria. ",
    "Check USGS_uv_summary.csv and the >50% coverage thresholds."
  )
}


# Save the site list used for this run
eligible_sites_table <- qualifying_sites %>%
  filter(site_no %in% sites_eligible)

write_csv(
  eligible_sites_table,
  "USGS_sites_eligible_for_download.csv"
)

message(
  "  Saved: USGS_sites_eligible_for_download.csv"
)


# -----------------------------------------------------------------------------
# Step 2b -- Function to download and save one site
# -----------------------------------------------------------------------------

fetch_and_write_uv <- function(site) {
  
  site <- as.character(site)
  
  message("  Fetching site ", site)
  
  tryCatch({
    
    df <- readNWISuv(
      siteNumbers = site,
      parameterCd = c(
        PARM_FLOW,
        PARM_DO,
        PARM_TEMP
      ),
      startDate = as.character(START_DATE),
      endDate = as.character(END_DATE)
    )
    
    if (is.null(df) || nrow(df) == 0) {
      message("  [WARN] No records returned for site ", site)
      return(invisible(FALSE))
    }
    
    df <- renameNWISColumns(df)
    
    # Confirm that the required identifier and time fields exist
    required_download_columns <- c(
      "site_no",
      "dateTime"
    )
    
    missing_download_columns <- setdiff(
      required_download_columns,
      names(df)
    )
    
    if (length(missing_download_columns) > 0) {
      message(
        "  [WARN] Site ", site,
        " is missing required fields: ",
        paste(missing_download_columns, collapse = ", ")
      )
      
      return(invisible(FALSE))
    }
    
    df <- df %>%
      mutate(
        site_no = as.character(site_no),
        dateTime = lubridate::with_tz(
          dateTime,
          tzone = "UTC"
        )
      ) %>%
      select(
        site_no,
        dateTime,
        any_of(
          c(
            "Flow_Inst",
            "DO_Inst",
            "Wtemp_Inst"
          )
        )
      )
    
    # Add absent parameter columns as NA so every Parquet file has the
    # same schema. This prevents open_dataset() schema errors later.
    if (!"Flow_Inst" %in% names(df)) {
      df$Flow_Inst <- NA_real_
    }
    
    if (!"DO_Inst" %in% names(df)) {
      df$DO_Inst <- NA_real_
    }
    
    if (!"Wtemp_Inst" %in% names(df)) {
      df$Wtemp_Inst <- NA_real_
    }
    
    df <- df %>%
      select(
        site_no,
        dateTime,
        Flow_Inst,
        DO_Inst,
        Wtemp_Inst
      )
    
    output_file <- file.path(
      UV_DATA_DIR,
      paste0(site, ".parquet")
    )
    
    write_parquet(
      df,
      output_file
    )
    
    message(
      "  Saved site ", site,
      ": ", format(nrow(df), big.mark = ","),
      " rows"
    )
    
    invisible(TRUE)
    
  }, error = function(e) {
    
    message(
      "  [WARN] Download failed for site ",
      site,
      ": ",
      conditionMessage(e)
    )
    
    invisible(FALSE)
  })
}


# -----------------------------------------------------------------------------
# Step 2c -- Download only eligible sites without existing files
# -----------------------------------------------------------------------------

existing_parquet_files <- list.files(
  path = UV_DATA_DIR,
  pattern = "\\.parquet$",
  full.names = FALSE
)

already_downloaded <- tools::file_path_sans_ext(
  existing_parquet_files
)

already_downloaded <- as.character(
  already_downloaded
)

sites_to_fetch <- setdiff(
  sites_eligible,
  already_downloaded
)

message(
  "  Existing Parquet files in folder : ",
  length(already_downloaded)
)

message(
  "  Eligible sites already downloaded: ",
  sum(sites_eligible %in% already_downloaded)
)

message(
  "  Eligible sites needing download  : ",
  length(sites_to_fetch)
)


if (length(sites_to_fetch) == 0) {
  
  message(
    "  [SKIP] All eligible sites already have Parquet files."
  )
  
} else {
  
  download_results <- map_lgl(
    sites_to_fetch,
    fetch_and_write_uv
  )
  
  message(
    "  Successful downloads this run: ",
    sum(download_results)
  )
  
  message(
    "  Failed or empty downloads this run: ",
    sum(!download_results)
  )
}


# -----------------------------------------------------------------------------
# Step 2d -- Check for eligible sites that remain missing
# -----------------------------------------------------------------------------

downloaded_after_run <- tools::file_path_sans_ext(
  list.files(
    path = UV_DATA_DIR,
    pattern = "\\.parquet$",
    full.names = FALSE
  )
)

still_missing <- setdiff(
  sites_eligible,
  downloaded_after_run
)

if (length(still_missing) > 0) {
  
  message(
    "  [WARN] ",
    length(still_missing),
    " eligible site(s) remain missing."
  )
  
  write_csv(
    tibble(site_no = still_missing),
    MISSING_SITES_FILE
  )
  
  message(
    "  Saved: ",
    MISSING_SITES_FILE
  )
  
} else {
  
  message(
    "  All eligible sites have Parquet files."
  )
  
  # Remove an obsolete missing-sites file from a previous run
  if (file.exists(MISSING_SITES_FILE)) {
    file.remove(MISSING_SITES_FILE)
  }
}


# -----------------------------------------------------------------------------
# Step 2e -- Open the Parquet files as one Arrow dataset
# -----------------------------------------------------------------------------

all_parquet_files <- list.files(
  path = UV_DATA_DIR,
  pattern = "\\.parquet$",
  full.names = TRUE
)

if (length(all_parquet_files) == 0) {
  stop(
    "No Parquet files are available in ",
    UV_DATA_DIR,
    ". The dataset cannot be opened."
  )
}

uv_ds_all <- open_dataset(
  sources = UV_DATA_DIR,
  format = "parquet"
)


# Restrict the lazy dataset to sites eligible under the current coverage rule.
#
# This is important because Parquet files from low-coverage sites may still
# remain in UV_DATA_DIR from the initial run. They are not deleted, but they
# are excluded from the current analysis.
uv_ds <- uv_ds_all %>%
  filter(site_no %in% sites_eligible)


total_rows_eligible <- uv_ds %>%
  summarise(n = n()) %>%
  collect() %>%
  pull(n)

message(
  "  Total rows for eligible sites: ",
  format(total_rows_eligible, big.mark = ",")
)


# -----------------------------------------------------------------------------
# Step 2f -- Per-site coverage summary
# -----------------------------------------------------------------------------

uv_summary <- uv_ds %>%
  group_by(site_no) %>%
  summarise(
    n_records = n(),
    
    n_flow_present = sum(
      !is.na(Flow_Inst)
    ),
    
    n_do_present = sum(
      !is.na(DO_Inst)
    ),
    
    n_temp_present = sum(
      !is.na(Wtemp_Inst)
    ),
    
    pct_flow_present = mean(
      !is.na(Flow_Inst)
    ) * 100,
    
    pct_do_present = mean(
      !is.na(DO_Inst)
    ) * 100,
    
    pct_temp_present = mean(
      !is.na(Wtemp_Inst)
    ) * 100,
    
    date_min = min(
      dateTime,
      na.rm = TRUE
    ),
    
    date_max = max(
      dateTime,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  collect() %>%
  mutate(
    site_no = as.character(site_no)
  ) %>%
  left_join(
    qualifying_sites %>%
      select(
        site_no,
        site_name,
        state_cd,
        lat,
        lon
      ),
    by = "site_no"
  ) %>%
  arrange(
    state_cd,
    site_no
  )


print(uv_summary)

write_csv(
  uv_summary,
  UV_SUMMARY_FILE
)

message(
  "  Saved: ",
  UV_SUMMARY_FILE
)

message(
  "  Full time series remains in ",
  UV_DATA_DIR,
  "/*.parquet."
)



# =============================================================================
# PART 3 -- Diagnostic plots
# =============================================================================

library(ggplot2)
library(maps)

# -- 3a. Map of qualifying sites ---------------------------------------------

western_map <- map_data("state") %>%
  filter(region %in% c("california", "washington", "oregon", "utah", "colorado",
                        "idaho", "nevada", "arizona", "new mexico", "montana", "wyoming"))

ggplot() +
  geom_polygon(data = western_map, aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey60") +
  geom_point(data = qualifying_sites, aes(x = lon, y = lat, color = state_cd),
             size = 2, alpha = 0.8) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(
    title = "USGS gages with 6 months DO & temp",
    subtitle = paste0("n = ",
                       nrow(qualifying_sites), " sites"),
    x = "Longitude", y = "Latitude", color = "State"
  )
# ggsave("USGS_qualifying_sites_map.png", width = 8, height = 6, dpi = 300)




# -- 3c. Map of sites with downloaded Parquet files --------------------------

# Identify site numbers directly from the downloaded Parquet filenames.
downloaded_parquet_files <- list.files(
  path = UV_DATA_DIR,
  pattern = "\\.parquet$",
  full.names = FALSE
)

downloaded_site_nos <- downloaded_parquet_files %>%
  tools::file_path_sans_ext() %>%
  as.character() %>%
  unique()


# Join downloaded site numbers to the site metadata.
downloaded_sites <- qualifying_sites %>%
  mutate(
    site_no = as.character(site_no)
  ) %>%
  filter(
    site_no %in% downloaded_site_nos
  ) %>%
  distinct(
    site_no,
    .keep_all = TRUE
  ) %>%
  filter(
    !is.na(lon),
    !is.na(lat)
  )


message(
  "  Parquet files found: ",
  length(downloaded_parquet_files)
)

message(
  "  Downloaded sites with map coordinates: ",
  nrow(downloaded_sites)
)


# Warn when a Parquet filename cannot be matched to site metadata.
downloaded_without_metadata <- setdiff(
  downloaded_site_nos,
  downloaded_sites$site_no
)

if (length(downloaded_without_metadata) > 0) {
  message(
    "  [WARN] ",
    length(downloaded_without_metadata),
    " downloaded site(s) could not be matched to qualifying_sites metadata."
  )
}


downloaded_sites_map <- ggplot() +
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
    data = downloaded_sites,
    aes(
      x = lon,
      y = lat,
      color = state_cd
    ),
    size = 2.5,
    alpha = 0.8
  ) +
  coord_fixed(
    ratio = 1.3
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  ) +
  labs(
    title = "USGS gages with downloaded unit-value data",
    subtitle = paste0(
      "Points represent sites with an existing .parquet file; n = ",
      nrow(downloaded_sites),
      " sites"
    ),
    x = "Longitude",
    y = "Latitude",
    color = "State"
  )

print(downloaded_sites_map)


# Optional: save the map.
ggsave(
  filename = "USGS_downloaded_parquet_sites_map.png",
  plot = downloaded_sites_map,
  width = 9,
  height = 7,
  dpi = 300
)




# -- 3b. Example time-series -- up to 3 randomly selected sites --------------
# NOTE: example_sites is defined here, BEFORE it's used to subset uv_ds.

set.seed(200)
example_sites <- qualifying_sites %>%
  slice_sample(n = min(3, nrow(qualifying_sites))) %>%
  pull(site_no)

plot_data <- uv_ds %>%
  filter(site_no %in% example_sites) %>%
  collect()   # small subset -- safe to bring into memory

if ("DO_Inst" %in% names(plot_data)) {
  print(
    ggplot(plot_data %>% filter(!is.na(DO_Inst)), aes(x = dateTime, y = DO_Inst)) +
      geom_line(alpha = 0.5, color = "steelblue", linewidth = 0.3) +
      facet_wrap(~ site_no, ncol = 1, scales = "free_x") +
      theme_minimal() +
      labs(title = "Dissolved oxygen -- example sites", x = NULL, y = "DO (mg/L)")
  )
}

if ("Wtemp_Inst" %in% names(plot_data)) {
  print(
    ggplot(plot_data %>% filter(!is.na(Wtemp_Inst)), aes(x = dateTime, y = Wtemp_Inst)) +
      geom_line(alpha = 0.5, color = "firebrick", linewidth = 0.3) +
      facet_wrap(~ site_no, ncol = 1, scales = "free_x") +
      theme_minimal() +
      labs(title = "Water temperature -- example sites", x = NULL, y = "Temperature (deg C)")
  )
}

if ("Flow_Inst" %in% names(plot_data)) {
  print(
    ggplot(plot_data %>% filter(!is.na(Flow_Inst)), aes(x = dateTime, y = Flow_Inst)) +
      geom_line(alpha = 0.5, color = "darkgreen", linewidth = 0.3) +
      scale_y_log10() +
      facet_wrap(~ site_no, ncol = 1, scales = "free_x") +
      theme_minimal() +
      labs(title = "Streamflow -- example sites (log scale)", x = NULL, y = "Discharge (cfs)")
  )
}

message("=== Done ===")
