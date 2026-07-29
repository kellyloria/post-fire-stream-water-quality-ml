# post-fire-stream-water-quality-ml

Exploratory project into characterizing stream water quality changes following wildfire using ML

1. Download USGS observations of flow, temperature, and where possible DO in USGS_obs_download.R

2. Process downloaded .parquet files into .rds file for time series analysis in USGS_processing.R

3. Plot coverage and broad dataset summaries in plotting.R

4. Link NHD attributes for landcover classifications from https://www.usgs.gov/national-hydrography/nhd-watershed-tool in COMID_NHD_attributes.R
