# install.packages(c("dataRetrieval", "dplyr", "purrr", "readr", "stringr"))
library(dataRetrieval)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)

# -----------------------------
# Goal:
# Download surface-water dissolved oxygen and water temperature
# observations for stream sites in Western US.
#
# Notes:
# - statecode = "US:32" is Nevada in WQP queries.
# - siteType = "Stream" limits to streams.
# - sampleMedia = "Water" limits to water samples.
# - Multiple characteristic names are used because WQP naming differs
#   across USGS and EPA/WQX providers.
# -----------------------------

# Nevada FIPS in WQP format
state_code <- "US:32"

# Western states
state_codess <- c(
  "US:06",  # California
  "US:41",  # Oregon
  "US:53",  # Washington
  "US:16",  # Idaho
  "US:49",  # Utah
  "US:04",  # Arizona
  "US:32"   # Nevada
)

# Use several characteristic names to avoid missing records due to naming differences
characteristics <- c(
  # dissolved oxygen variants
  "Dissolved oxygen (DO)",
  "Dissolved oxygen",
  "Oxygen",
  
  # water temperature variants
  "Temperature, water",
  "Water temperature",
  "Temperature"
)

# Pull each characteristic separately.
# This is usually more reliable than one very large multi-characteristic call.


characteristics <- c(
  "Dissolved oxygen (DO)",
  "Dissolved oxygen",
  "Oxygen",
  "Temperature, water",
  "Water temperature",
  "Temperature"
)

wqp_list <- map(
  characteristics,
  \(ch) {
    message("Downloading: ", ch)
    
    out <- tryCatch(
      readWQPdata(
        statecode = state_codes,
        siteType = "Stream",
        sampleMedia = "Water",
        characteristicName = ch
      ),
      error = function(e) {
        message("  Failed for ", ch, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(out) || nrow(out) == 0) return(NULL)
    
    # Force all columns to character so bind_rows() works
    out %>%
      mutate(across(everything(), as.character)) %>%
      as_tibble()
  }
)

wqp_raw <- bind_rows(wqp_list)


sites <- whatWQPsites(statecode = "US:32", siteType = "Stream") %>%
  transmute(
    MonitoringLocationIdentifier,
    site_lat = as.numeric(LatitudeMeasure),
    site_lon = as.numeric(LongitudeMeasure)
  )

wqp_clean <- wqp_raw %>%
  mutate(
    ActivityStartDate = as.Date(ActivityStartDate),
    ResultMeasureValue = suppressWarnings(as.numeric(ResultMeasureValue)),
    characteristic_lower = str_to_lower(CharacteristicName),
    parameter_group = case_when(
      str_detect(characteristic_lower, "dissolved oxygen|oxygen") ~ "dissolved_oxygen",
      str_detect(characteristic_lower, "temperature") ~ "temperature",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(parameter_group)) %>%
  distinct() %>%
  left_join(sites, by = "MonitoringLocationIdentifier")


head(wqp_clean)
names(wqp_clean)

# write.csv(wqp_clean, "./nevada_gage_wq_dat.csv")


do_counts <- wqp_clean %>%
  filter(parameter_group == "dissolved_oxygen") %>%
  filter(!is.na(site_lat), !is.na(site_lon)) %>%
  group_by(MonitoringLocationIdentifier, site_lat, site_lon) %>%
  summarise(n_obs = n(), .groups = "drop")





library(ggplot2)

ggplot(do_counts, aes(x = site_lon, y = site_lat)) +
  geom_point(alpha = 0.7, size = 2) +
  coord_fixed() +
  labs(
    title = "Nevada stream dissolved oxygen observation sites",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()



ggplot(do_counts, aes(x = site_lon, y = site_lat, size = n_obs)) +
  geom_point(alpha = 0.7) +
  coord_fixed() +
  labs(
    title = "Nevada stream dissolved oxygen sites",
    subtitle = "Point size = number of observations",
    x = "Longitude",
    y = "Latitude",
    size = "Observations"
  ) +
  theme_minimal()



library(ggplot2)
library(maps)

nv_map <- map_data("state") %>%
  filter(region == "nevada")

ggplot() +
  geom_polygon(
    data = nv_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey50"
  ) +
  geom_point(
    data = do_counts,
    aes(x = site_lon, y = site_lat, size = n_obs),
    alpha = 0.7
  ) +
  coord_fixed(1.3) +
  labs(
    title = "Nevada stream dissolved oxygen monitoring sites",
    subtitle = "Point size = number of observations",
    x = "Longitude",
    y = "Latitude",
    size = "Observations"
  ) +
  theme_minimal()








site_counts <- wqp_clean %>%
  filter(!is.na(site_lat), !is.na(site_lon)) %>%
  group_by(parameter_group, MonitoringLocationIdentifier, site_lat, site_lon) %>%
  summarise(n_obs = n(), .groups = "drop")

ggplot() +
  geom_polygon(
    data = nv_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey50"
  ) +
  geom_point(
    data = site_counts,
    aes(x = site_lon, y = site_lat, size = n_obs, color = parameter_group),
    alpha = 0.2
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(
    title = "Nevada stream observation sites",
    x = "Longitude",
    y = "Latitude",
    size = "Observations",
    color = "Parameter"
  )



DOplot <- ggplot() +
  geom_polygon(
    data = nv_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey50"
  ) +
  geom_point(
    data = site_counts %>% filter(parameter_group == "dissolved_oxygen"),
    aes(x = site_lon, y = site_lat, size = n_obs),
    alpha = 0.4,
    color = "blue"
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(
    title = "Nevada stream dissolved oxygen sites",
    x = "Longitude",
    y = "Latitude",
    size = "Observations"
  )


# ggsave(plot = DOplot, filename = paste("./NV_stream_DO.png",sep=""),width=5,height=4,dpi=300)


tempplot <- ggplot() +
  geom_polygon(
    data = nv_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey50"
  ) +
  geom_point(
    data = site_counts %>% filter(parameter_group == "temperature"),
    aes(x = site_lon, y = site_lat, size = n_obs),
    alpha = 0.4,
    color = "red"
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(
    title = "Nevada stream temperature sites",
    x = "Longitude",
    y = "Latitude",
    size = "Observations"
  )


# ggsave(plot = tempplot, filename = paste("./NV_stream_temp.png",sep=""),width=5,height=4,dpi=300)




set.seed(131)  # for reproducibility

random_sites <- wqp_clean %>%
  filter(!is.na(parameter_group)) %>%
  distinct(MonitoringLocationIdentifier) %>%
  slice_sample(n = 3) %>%
  pull(MonitoringLocationIdentifier)

library(ggplot2)

ts_data <- wqp_clean %>%
  filter(
    MonitoringLocationIdentifier %in% random_sites,
    parameter_group %in% c("dissolved_oxygen", "temperature"),
    !is.na(ActivityStartDate),
    !is.na(ResultMeasureValue)
  )

ggplot(ts_data, aes(x = ActivityStartDate, y = ResultMeasureValue)) +
  geom_point(alpha = 0.4) +
  geom_line(alpha = 0.5) +
  facet_grid(parameter_group ~ MonitoringLocationIdentifier, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Time series of dissolved oxygen and temperature",
    x = "Date",
    y = "Value",
    subtitle = "Three randomly selected Nevada stream sites"
  )


#USGS-10244950, 119USBR_QXDVO11, 21NEV1_WQX_W2

sites_of_interest <- c(
  "USGS-10244950",
  #"119USBR_QXDVO11",
  #"21NEV1_WQX-W2"
  "21Nev1_WQX-SB5"
)

ts_data <- wqp_clean %>%
  filter(
    MonitoringLocationIdentifier %in% sites_of_interest,
    parameter_group %in% c("dissolved_oxygen", "temperature"),
    !is.na(ActivityStartDate),
    !is.na(ResultMeasureValue)
  )

# compute scaling factor
do_range <- range(ts_data$ResultMeasureValue[ts_data$parameter_group == "dissolved_oxygen"], na.rm = TRUE)
temp_range <- range(ts_data$ResultMeasureValue[ts_data$parameter_group == "temperature"], na.rm = TRUE)

scale_factor <- diff(do_range) / diff(temp_range)


library(ggplot2)

example_dat <- ggplot() +
  # DO (left axis)
  geom_line(
    data = ts_data %>% filter(parameter_group == "dissolved_oxygen"),
    aes(x = ActivityStartDate, y = ResultMeasureValue, group = MonitoringLocationIdentifier),
    color = "blue",
    alpha = 0.7
  ) +
  
  # Temperature (scaled to DO axis)
  geom_line(
    data = ts_data %>% filter(parameter_group == "temperature"),
    aes(x = ActivityStartDate, y = ResultMeasureValue * scale_factor, group = MonitoringLocationIdentifier),
    color = "red",
    alpha = 0.7
  ) +
  
  facet_wrap(~ MonitoringLocationIdentifier, scales = "free_x") +
  
  scale_y_continuous(
    name = "Dissolved Oxygen (%)",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Temperature"
    )
  ) +
  
  theme_minimal() +
  labs(
    title = "Raw Dissolved Oxygen (blue) and Temperature (red)",
    x = "Date"
  )



example_dat2 <- ggplot() +
  # DO (left axis)
  geom_line(
    data = ts_data %>% filter(parameter_group == "dissolved_oxygen"),
    aes(x = ActivityStartDate, y = ResultMeasureValue, group = MonitoringLocationIdentifier),
    color = "blue",
    alpha = 0.7
  ) +
  
  # Temperature (scaled to DO axis)
  geom_line(
    data = ts_data %>% filter(parameter_group == "temperature"),
    aes(x = ActivityStartDate, y = ResultMeasureValue * scale_factor, group = MonitoringLocationIdentifier),
    color = "red",
    alpha = 0.7
  ) +
  
  facet_wrap(~ MonitoringLocationIdentifier, scales = "free_x") +
  
  scale_y_continuous(
    name = "Dissolved Oxygen (%)",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Temperature"
    )
  ) +
  
  theme_minimal() +
  labs(
    title = "Raw Dissolved Oxygen (blue) and Temperature (red)",
    x = "Date"
  )


sites_of_interest <- c(
 # "USGS-10244950",
  "119USBR_QXDVO11",
  "21NEV1_WQX-W2",
  "21NEV1_WQX-SB5"
)

do_data <- wqp_clean %>%
  filter(
    MonitoringLocationIdentifier %in% sites_of_interest,
    parameter_group == "dissolved_oxygen",
    !is.na(ActivityStartDate),
    !is.na(ResultMeasureValue)
  )

library(ggplot2)

example_dat2_DO <- ggplot(do_data%>%filter(ActivityStartDate> as.Date("1981-03-01")), aes(x = ActivityStartDate, y = ResultMeasureValue)) +
  geom_point(alpha = 0.4, color = "blue") +
  geom_line(alpha = 0.6, color = "blue") +
  facet_grid(MonitoringLocationIdentifier ~ ., scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Dissolved Oxygen Time Series (Selected Sites)",
    x = "Date",
    y = "Dissolved Oxygen (mg/L)"
  )

# ggsave(plot =example_dat2_DO , filename = paste("./NV_stream_DO_example.png",sep=""),width=5,height=4,dpi=300)




temp_data <- wqp_clean %>%
  filter(
    MonitoringLocationIdentifier %in% sites_of_interest,
    parameter_group == "temperature",
    !is.na(ActivityStartDate),
    !is.na(ResultMeasureValue)
  )

library(ggplot2)


example_dat2_temp <- ggplot(temp_data%>%filter(ActivityStartDate> as.Date("1981-03-01")), aes(x = ActivityStartDate, y = ResultMeasureValue)) +
  geom_point(alpha = 0.4, color = "red") +
  geom_line(alpha = 0.6, color = "red") +
  facet_grid(MonitoringLocationIdentifier ~ ., scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Date",
    y = "Temperature (C)"
  )

do_data$ActivityStartDate

# ggsave(plot =example_dat2_temp , filename = paste("./NV_stream_temp_example.png",sep=""),width=5,height=4,dpi=300)


library(dplyr)
library(stringr)

## ------------------------------------------------------------------
## 1.  Derive a very simple region field (if you don't already have one)
## ------------------------------------------------------------------
# This example splits Nevada roughly into "North" and "South"
# You can replace the `case_when()` logic with your own rule set
wqp_clean <- wqp_clean %>% 
  mutate(
    # If you already have a region column, comment out the next line
    Region = case_when(
      site_lat >= 39.5 ~ "North",
      site_lat <  39.5 ~ "South",
      TRUE ~ NA_character_
    )
  )

## ------------------------------------------------------------------
## 2.  Summary: streams by collecting organization
## ------------------------------------------------------------------


# Define a character vector of organization names you want to drop
drop_orgs <- c(
  "AZDEQ_SW",
  "Great Lakes Environmental Center (Volunteer)",
  "USGS Arizona Water Science Center",
  "USGS Utah Water Science Center",
  "Utah Department Of Environmental Quality"
)

# Compute the stream summary, then filter out those names
streams_by_org <- wqp_clean %>%
  group_by(OrganizationFormalName) %>%          # or whatever the exact column name is
  summarise(
    n_streams = n_distinct(MonitoringLocationIdentifier),
    .groups = "drop"
  ) %>%
  filter(!OrganizationFormalName %in% drop_orgs)

## ------------------------------------------------------------------
## 3.  Summary: streams by geographic region
## ------------------------------------------------------------------
streams_by_region <- wqp_clean %>%
  group_by(Region) %>%
  summarise(
    n_streams = n_distinct(MonitoringLocationIdentifier),
    .groups = "drop"
  )

## ------------------------------------------------------------------
## 4.  Summary: streams by both region and collecting organization
## ------------------------------------------------------------------
streams_by_region_org <- wqp_clean %>%
  group_by(Region, OrganizationFormalName) %>%
  summarise(
    n_streams = n_distinct(MonitoringLocationIdentifier),
    .groups = "drop"
  )

## ------------------------------------------------------------------
## 5.  View / export the results
## ------------------------------------------------------------------
print(streams_by_org)
print(streams_by_region)
print(streams_by_region_org)

# Optional: write to CSV for reporting
write.csv(streams_by_org, "streams_by_organization.csv", row.names = FALSE)
write.csv(streams_by_region, "streams_by_region.csv", row.names = FALSE)
write.csv(streams_by_region_org, "streams_by_region_and_organization.csv", row.names = FALSE)



