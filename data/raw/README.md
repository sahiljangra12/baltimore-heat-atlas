# Raw inputs (not tracked in git)

The raw inputs total about 140 MB. They are attached to the repository's GitHub Release as `raw_data.zip`. Unzip that file here so this folder contains:

| File | Source |
|---|---|
| `landsat/bhs_ls_*.tif` | Landsat 8/9 C2 L2 crops (bands lwir11, qa_pixel, red, nir08), summers 2019–2025, Microsoft Planetary Computer data API |
| `landsat/bhs_worldcover2021.tif` | ESA WorldCover 2021 v200, 10 m |
| `modis.json.gz` | MODIS MYD11A2/MOD11A2 v061 8-day LST crops (base64 GeoTIFFs) |
| `bhs_places_2025.csv` | CDC PLACES 2025 census-tract release (`cwsq-ngmh`), Baltimore City, crude prevalence |
| `bhs_acs2024_tracts.json` | ACS 5-year 2020–2024, 18 tables with margins of error (Census Reporter API) |
| `bhs_svi2022_maryland.csv` | CDC/ATSDR SVI 2022, Maryland |
| `bhs_tracts2020.geojson` | 2020 census tracts (TIGERweb) |
| `bhs_misc_bundle.json` | 2020 blocks with population, ZCTAs, cooling centres 2026, ACIS station data |
| `bhs_HHI_Data.zip` → `hhi/` | CDC/ATSDR Heat & Health Index 2024 |
| `holc_baltimore_tract2020.gpkg` | Mapping Inequality HOLC × 2020-tract crosswalk |
| `csa2020.geojson.gz` | BNIA Community Statistical Areas 2020 |
