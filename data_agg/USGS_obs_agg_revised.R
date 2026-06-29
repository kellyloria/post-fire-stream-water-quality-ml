# =============================================================================
# USGS_obs_agg_revised.R
#
# Goals:
#   1. Identify USGS stream gage sites in CA, WA, OR, ID, NV, AZ that have
#      more than 6 months of sub-daily (unit-value / instantaneous) dissolved
#      oxygen (parm_cd 00300) AND water temperature (parm_cd 00010) records.
#
#   2. Download a continuous (unit-value) time-series dataframe of:
#        - Streamflow        parm_cd = 00060  (cfs)
#        - Dissolved oxygen  parm_cd = 00300  (mg/L)
#        - Water temperature parm_cd = 00010  (deg C)
#      for qualifying sites, covering 2011-10-01 through 2024-10-01.
#
# Key design decisions
# ---------------------
#  * Site discovery uses whatNWISsites(service = "uv") — returns only sites
#    with unit-value (sub-daily) data for a given parameter.
#
#  * Data availability / record length uses whatNWISdata() — this function
#    returns a catalogue row per site × parameter × service that includes
#    `begin_date`, `end_date`, and `count_nu` (number of values on record).
#    It avoids readNWISstat() entirely, which (a) caps requests at 10 sites,
#    and (b) does not support statType = "count".
#
#  * "More than 6 months of sub-daily data" is operationalised as the span
#    between begin_date and end_date being > 180 calendar days AND count_nu
#    > 180 (i.e. actual values exist across that span, not just an open
#    period with a gap).
#
# Package requirements:
#   install.packages(c("dataRetrieval", "dplyr", "purrr", "lubridate", "readr"))
# =============================================================================

library(dataRetrieval)
library(dplyr)
library(purrr)
library(lubridate)
library(readr)

# -----------------------------------------------------------------------------
# 0.  Constants
# -----------------------------------------------------------------------------

target_states <- c("CA", "WA", "OR", "ID", "NV", "AZ")

PARM_DO   <- "00300"   # Dissolved oxygen, mg/L
PARM_TEMP <- "00010"   # Water temperature, deg C
PARM_FLOW <- "00060"   # Streamflow / discharge, cfs

START_DATE <- as.Date("2011-10-01")
END_DATE   <- as.Date("2024-10-01")

# >6 months threshold — applied to both span (days) and count of records
MIN_DAYS <- 180


# =============================================================================
# PART 1 — Identify qualifying sites
# =============================================================================
# Step 1a : whatNWISsites()  — find stream sites with uv DO and uv temp
# Step 1b : whatNWISdata()   — get record catalogue (begin/end/count) per site
# Step 1c : filter sites where BOTH DO and temp exceed MIN_DAYS
# =============================================================================

# -----------------------------------------------------------------------------
# Step 1a — Candidate sites (have uv service for DO *and* temp)
# -----------------------------------------------------------------------------

message("=== Step 1a: Querying NWIS for candidate stream sites ===")

fetch_sites <- function(state, parm) {
  message("  state=", state, "  parm=", parm)
  tryCatch(
    whatNWISsites(
      stateCd     = state,
      siteType    = "ST",   # streams only
      parameterCd = parm,
      service     = "uv"    # unit-value (sub-daily) service
    ),
    error = function(e) {
      message("  [WARN] whatNWISsites failed: ", e$message)
      NULL
    }
  )
}

message("  Fetching DO (", PARM_DO, ") sites ...")
do_sites_raw <- map_dfr(target_states, fetch_sites, parm = PARM_DO)

message("  Fetching temperature (", PARM_TEMP, ") sites ...")
temp_sites_raw <- map_dfr(target_states, fetch_sites, parm = PARM_TEMP)

# Intersection: sites that carry BOTH parameters on the uv service
candidate_sites <- intersect(
  unique(do_sites_raw$site_no),
  unique(temp_sites_raw$site_no)
)

message("  Candidate sites (both DO + temp, uv service): ",
        length(candidate_sites))


# -----------------------------------------------------------------------------
# Step 1b — Record availability catalogue via whatNWISdata()
# -----------------------------------------------------------------------------
# whatNWISdata() accepts a vector of site numbers of any length (no hard cap)
# and returns one row per site × parameter × service with:
#   begin_date, end_date, count_nu (total values in the period of record)
#
# We filter to service = "uv" and our two parameters of interest.
# -----------------------------------------------------------------------------

message("=== Step 1b: Fetching data availability catalogue (whatNWISdata) ===")

# Chunk into groups of 100 as a conservative safety measure against
# URL-length limits (whatNWISdata uses GET, not POST).
CHUNK_SIZE <- 100
site_chunks <- split(
  candidate_sites,
  ceiling(seq_along(candidate_sites) / CHUNK_SIZE)
)

fetch_availability <- function(sites_chunk) {
  tryCatch(
    whatNWISdata(
      siteNumber  = sites_chunk,
      service     = "uv",
      parameterCd = c(PARM_DO, PARM_TEMP)
    ),
    error = function(e) {
      message("  [WARN] whatNWISdata chunk failed: ", e$message)
      NULL
    }
  )
}

message("  Fetching availability in ", length(site_chunks),
        " chunk(s) of up to ", CHUNK_SIZE, " sites ...")

avail_raw <- map_dfr(site_chunks, fetch_availability)

message("  Availability catalogue rows returned: ", nrow(avail_raw))


# -----------------------------------------------------------------------------
# Step 1c — Filter: >6 months of sub-daily data for BOTH parameters
# -----------------------------------------------------------------------------
# count_nu   : total number of unit-values on record (each ~15 min)
# begin_date / end_date : period of record
#
# Criterion applied per site × parm:
#   (1) span_days  = end_date − begin_date  > MIN_DAYS
#   (2) count_nu                            > MIN_DAYS
#       (guards against a long-open period that is mostly missing)
# -----------------------------------------------------------------------------

avail_filtered <- avail_raw %>%
  mutate(
    begin_date = as.Date(begin_date),
    end_date   = as.Date(end_date),
    span_days  = as.numeric(end_date - begin_date)
  ) %>%
  filter(
    parm_cd  %in% c(PARM_DO, PARM_TEMP),
    data_type_cd == "uv",
    span_days > MIN_DAYS,
    count_nu  > MIN_DAYS
  )

# Summarise per site × parameter
parm_site_summary <- avail_filtered %>%
  group_by(site_no, parm_cd) %>%
  summarise(
    total_count = sum(count_nu, na.rm = TRUE),
    rec_start   = min(begin_date, na.rm = TRUE),
    rec_end     = max(end_date,   na.rm = TRUE),
    .groups     = "drop"
  )

# Keep only sites that pass the threshold for BOTH DO and temperature
sites_with_do   <- parm_site_summary %>% filter(parm_cd == PARM_DO)   %>% pull(site_no)
sites_with_temp <- parm_site_summary %>% filter(parm_cd == PARM_TEMP) %>% pull(site_no)
qualifying_site_nos <- intersect(sites_with_do, sites_with_temp)

message("  Sites passing >", MIN_DAYS, "-day threshold for BOTH DO and temp: ",
        length(qualifying_site_nos))


# -----------------------------------------------------------------------------
# Step 1d — Attach site metadata (lat/lon, name, state)
# -----------------------------------------------------------------------------

site_meta <- bind_rows(do_sites_raw, temp_sites_raw) %>%
  distinct(site_no, .keep_all = TRUE) %>%
  select(
    site_no,
    station_nm,
    #state_cd,
    dec_lat_va,
    dec_long_va
    #drain_area_va
  ) %>%
  rename(
    site_name     = station_nm,
    lat           = dec_lat_va,
    lon           = dec_long_va
    #drainage_sqmi = drain_area_va
  )

# Wide summary: one row per qualifying site
do_summary   <- sites_with_do   %>% { parm_site_summary[parm_site_summary$parm_cd == PARM_DO,   ] } %>%
  rename(do_count = total_count, do_start = rec_start, do_end = rec_end) %>%
  select(-parm_cd)

temp_summary <- sites_with_temp %>% { parm_site_summary[parm_site_summary$parm_cd == PARM_TEMP, ] } %>%
  rename(temp_count = total_count, temp_start = rec_start, temp_end = rec_end) %>%
  select(-parm_cd)

qualifying_sites <- tibble(site_no = qualifying_site_nos) %>%
  left_join(do_summary,   by = "site_no") %>%
  left_join(temp_summary, by = "site_no") %>%
  left_join(site_meta,    by = "site_no")

message("  Final qualifying site count: ", nrow(qualifying_sites))
print(qualifying_sites %>% select(site_no, site_name, lat, lon,
                                   do_count, do_start, do_end,
                                   temp_count, temp_start, temp_end))

write_csv(qualifying_sites, "USGS_qualifying_sites.csv")
message("  Saved: USGS_qualifying_sites.csv")


# =============================================================================
# PART 2 — Download continuous (unit-value) time-series
# =============================================================================
# readNWISuv() fetches instantaneous (sub-daily) values.
# All three parameters are requested per site in one call to minimise
# round-trips.  renameNWISColumns() replaces X_PARMCD_00000-style column
# names with human-readable equivalents.
# =============================================================================

message("=== Step 2: Downloading unit-value time-series ===")
message("    Parameters : ", PARM_FLOW, " (flow), ",
        PARM_DO, " (DO), ", PARM_TEMP, " (temp)")
message("    Date range : ", START_DATE, " – ", END_DATE)
message("    Sites      : ", nrow(qualifying_sites))

fetch_uv <- function(site) {
  message("  Downloading: ", site)
  tryCatch(
    readNWISuv(
      siteNumbers = site,
      parameterCd = c(PARM_FLOW, PARM_DO, PARM_TEMP),
      startDate   = as.character(START_DATE),
      endDate     = as.character(END_DATE)
    ) %>%
      renameNWISColumns(),
    error = function(e) {
      message("  [WARN] Download failed for site ", site, ": ", e$message)
      NULL
    }
  )
}

uv_list <- map(qualifying_sites$site_no, fetch_uv)
uv_list <- Filter(Negate(is.null), uv_list)
uv_raw  <- bind_rows(uv_list)

message("  Raw rows downloaded: ", nrow(uv_raw))


# -----------------------------------------------------------------------------
# Step 2a — Standardise column names
# -----------------------------------------------------------------------------
# renameNWISColumns() produces names such as:
#   Flow_Inst, Flow_Inst_cd        (00060)
#   DO_Inst,   DO_Inst_cd          (00300)
#   Wtemp_Inst, Wtemp_Inst_cd      (00010)
# The exact names present depend on which parameters are available at each
# site.  any_of() handles sites where one parameter is absent.
# -----------------------------------------------------------------------------

uv_clean <- uv_raw %>%
  select(
    site_no,
    dateTime,
    any_of(c("Flow_Inst",  "Flow_Inst_cd",
             "DO_Inst",    "DO_Inst_cd",
             "Wtemp_Inst", "Wtemp_Inst_cd"))
  ) %>%
  rename(
    datetime_utc  = dateTime,
    any_of(c(
      discharge_cfs  = "Flow_Inst",
      discharge_cd   = "Flow_Inst_cd",
      do_mgl         = "DO_Inst",
      do_cd          = "DO_Inst_cd",
      water_temp_c   = "Wtemp_Inst",
      water_temp_cd  = "Wtemp_Inst_cd"
    ))
  ) %>%
  left_join(
    qualifying_sites %>%
      select(site_no, site_name, state_cd, lat, lon),
    by = "site_no"
  ) %>%
  mutate(date = as.Date(datetime_utc)) %>%
  select(
    site_no, site_name, state_cd, lat, lon,
    datetime_utc, date,
    any_of(c("discharge_cfs", "discharge_cd",
             "do_mgl",        "do_cd",
             "water_temp_c",  "water_temp_cd"))
  )

message("=== Download complete ===")
message("  Rows in uv_clean  : ", nrow(uv_clean))
message("  Sites represented : ", n_distinct(uv_clean$site_no))
message("  Date range        : ",
        min(uv_clean$date, na.rm = TRUE), " to ",
        max(uv_clean$date, na.rm = TRUE))


# -----------------------------------------------------------------------------
# Step 2b — Per-site data-quality summary
# -----------------------------------------------------------------------------

uv_summary <- uv_clean %>%
  group_by(site_no, site_name, state_cd) %>%
  summarise(
    n_records         = n(),
    n_days            = n_distinct(date),
    pct_flow_present  = if ("discharge_cfs" %in% names(.))
                          mean(!is.na(discharge_cfs)) * 100 else NA_real_,
    pct_do_present    = if ("do_mgl"         %in% names(.))
                          mean(!is.na(do_mgl))         * 100 else NA_real_,
    pct_temp_present  = if ("water_temp_c"   %in% names(.))
                          mean(!is.na(water_temp_c))   * 100 else NA_real_,
    date_min          = min(date, na.rm = TRUE),
    date_max          = max(date, na.rm = TRUE),
    .groups           = "drop"
  )

print(uv_summary)

write_csv(uv_clean,   "USGS_uv_timeseries.csv")
write_csv(uv_summary, "USGS_uv_summary.csv")
message("  Saved: USGS_uv_timeseries.csv")
message("  Saved: USGS_uv_summary.csv")


# =============================================================================
# PART 3 — Diagnostic plots (requires ggplot2 + maps)
# =============================================================================

library(ggplot2)
library(maps)

# -- 3a.  Map of qualifying sites --------------------------------------------

western_map <- map_data("state") %>%
  filter(region %in% c("california", "washington", "oregon",
                        "idaho", "nevada", "arizona"))

ggplot() +
  geom_polygon(
    data  = western_map,
    aes(x = long, y = lat, group = group),
    fill  = "grey95", color = "grey60"
  ) +
  geom_point(
    data  = qualifying_sites,
    aes(x = lon, y = lat, color = state_cd),
    size  = 2, alpha = 0.8
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(
    title    = "USGS stream gages — >6 months sub-daily DO & temperature",
    subtitle = paste0("CA, WA, OR, ID, NV, AZ  |  n = ",
                      nrow(qualifying_sites), " sites"),
    x = "Longitude", y = "Latitude", color = "State"
  )

# ggsave("USGS_qualifying_sites_map.png", width = 8, height = 6, dpi = 300)


# -- 3b.  Example time-series — up to 3 randomly selected sites --------------

set.seed(42)
example_sites <- qualifying_sites %>%
  slice_sample(n = min(3, nrow(qualifying_sites))) %>%
  pull(site_no)

plot_data <- uv_clean %>%
  filter(site_no %in% example_sites)

if ("do_mgl" %in% names(plot_data)) {
  print(
    ggplot(plot_data %>% filter(!is.na(do_mgl)),
           aes(x = datetime_utc, y = do_mgl)) +
      geom_line(alpha = 0.5, color = "steelblue", linewidth = 0.3) +
      facet_wrap(~ site_no, ncol = 1, scales = "free_x") +
      theme_minimal() +
      labs(title = "Dissolved oxygen — example sites",
           x = NULL, y = "DO (mg/L)")
  )
}

if ("water_temp_c" %in% names(plot_data)) {
  print(
    ggplot(plot_data %>% filter(!is.na(water_temp_c)),
           aes(x = datetime_utc, y = water_temp_c)) +
      geom_line(alpha = 0.5, color = "firebrick", linewidth = 0.3) +
      facet_wrap(~ site_no, ncol = 1, scales = "free_x") +
      theme_minimal() +
      labs(title = "Water temperature — example sites",
           x = NULL, y = "Temperature (°C)")
  )
}

if (all(c("do_mgl", "discharge_cfs") %in% names(plot_data))) {
  print(
    ggplot(plot_data %>% filter(!is.na(discharge_cfs)),
           aes(x = datetime_utc, y = discharge_cfs)) +
      geom_line(alpha = 0.5, color = "darkgreen", linewidth = 0.3) +
      scale_y_log10() +
      facet_wrap(~ site_no, ncol = 1, scales = "free_x") +
      theme_minimal() +
      labs(title = "Streamflow — example sites (log scale)",
           x = NULL, y = "Discharge (cfs)")
  )
}
