# EO4PM
EO4PM (Earth Observation for Phenological Metrics) is an automated, transferable algorithm written in R designed to extract Land Surface Phenology (LSP) metrics from dense satellite Earth Observation (EO) time series

[![R-build](https://img.shields.io/badge/R->=%204.0.0-blue.svg)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL--v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Paper MDPI](https://img.shields.io/badge/DOI-10.3390%2Frs14030721-green.svg)](https://doi.org/10.3390/rs14030721)

**EO4PM (Earth Observation for Phenological Metrics)** is an automated, transferable algorithm written in R designed to extract Land Surface Phenology (LSP) metrics and temporal statistics from dense satellite Earth Observation (EO) time series.

Version 2 introduces major upgrades in command-line interface (CLI) mode, advanced temporal smoothing routines, expanded phenological metrics extraction, fallback fitting routines, C++-optimized multi-cycle cycle cutting, and a modular architecture providing standalone R functions.

---

## 🔑 Key Strengths

- **User-Defined Vegetation Indices & Biophysical Parameters (BIOPAR):** Flexible input compatibility with standard vegetation indices (e.g., NDVI, EVI, NIRv, kNDVI) as well as biophysical parameters like Leaf Area Index (**LAI**), FAPAR, and FCOVER. Using LAI helps avoid signal saturation over dense forest canopies and high-biomass ecosystems.
- **Sensor Agnostic:** Fully independent of satellite platform or sensor type. Works seamlessly with Sentinel-2, Landsat, MODIS, virtual constellations (e.g., Harmonized Landsat-Sentinel - HLS), and Synthetic Aperture Radar (SAR) time series.
- **R Ecosystem & CLI Ready:** Modular code structure supporting both interactive R sessions and automated high-performance CLI/batch execution.

---

## 🚀 What's New in Version 2

### 1. Advanced Temporal Smoothing (`smoothEO4PM`)
- **Noise Removal:** Automated removal of small drops and (optionally) small spikes, with specialized strategies to preserve peak maximum values in phenological cycles.
- **Daily Interpolation & Iterative Weighting:** Daily interpolation using Stineman or Whittaker smoothing algorithms with multiple iterative passes updating weight matrices (`iters`, `wFUN_type = "wTSM"`).
- **Secondary Whittaker Fitting (`fitMax`):** Optional second-pass Whittaker smoother specifically designed to improve fitting and alignment around peak/maximum vegetation values.

### 2. Upgraded Phenological Metrics Extraction (`phenoEO4PM`)
- **Multi-Cycle Phenology (`cut_cycle`):** Integration of the C++ optimized `cut_cycle` algorithm (originally developed in `phenofit`), enabling robust identification and partitioning of multiple vegetative/growing cycles per season (e.g., double cropping or understory dynamics).
- **Logit Regression Fallback System:** Automatic fallback system using logit regression whenever the standard Gu method fails to estimate phenological transition stages, preventing data gaps or script execution failures.
- **Expanded Metrics Suite:** Comprehensive extraction of transition dates, development rates, durations, and time-integrated production metrics.

### 3. Standalone R Functions & Modular Architecture
EO4PM v2 provides clean, individual R functions for each processing stage:
- `smoothEO4PM()`: Temporal smoothing, outlier cleaning, and daily interpolation.
- `phenoEO4PM()`: Phenological metrics extraction and multi-cycle detection.
- `plotEO4PM()`: Custom `ggplot2` visualization and figure export.
- `QFbitDecoder()`: Bitwise Quality Flag decoding and confidence assessment.

---

## 📄 Scientific Reference

If you use EO4PM in your research, please cite the foundational paper:

> **Filipponi, F., Smiraglia, D., & Agrillo, E. (2022).** Earth Observation for Phenological Metrics (EO4PM): Temporal Discriminant to Characterize Forest Ecosystems. *Remote Sensing*, 14(3), 721.  
> DOI: [10.3390/rs14030721](https://doi.org/10.3390/rs14030721)

```bibtex
@article{filipponi2022earth,
  title={Earth Observation for Phenological Metrics (EO4PM): Temporal Discriminant to Characterize Forest Ecosystems},
  author={Filipponi, Federico and Smiraglia, Daniela and Agrillo, Emiliano},
  journal={Remote Sensing},
  volume={14},
  number={3},
  pages={721},
  year={2022},
  publisher={MDPI},
  doi={10.3390/rs14030721}
}
```

---

## 🛠️ Requirements & Installation

### R Dependencies
EO4PM v2 requires **R (>= 4.0.0)** and the following R packages:

```R
install.packages(c("sf", "terra", "phenofit", "ggplot2", "gridExtra", "stinepack", "ptw"))
```

---

## 💻 Usage Examples

The following examples are based on `phenoEO4PM_examples.R`.

### 1. Single-Pixel Processing Pipeline

```R
suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(phenofit)
  library(ggplot2)
  library(gridExtra)
})

terra::setGDALconfig("IGNORE_XY_AXIS_NAME_CHECKS", "YES")

# Source standalone EO4PM functions
source("smoothEO4PM.R")
source("phenoEO4PM.R")
source("plotEO4PM.R")
source("EO4PM_QF_bitDecoder.R")

# Read input NetCDF raster time series
input_file <- "/path/to/S2_L3A_TS_10m_LAI.nc"
r <- terra::rast(input_file)

# Get observation dates and daily target sequence
time_dim <- terra::time(r)
obs_dates <- as.integer(as.Date(time_dim, tz = "UTC"))
tout <- seq(from = obs_dates[1], to = obs_dates[length(obs_dates)], by = 1)

# Target pixel index
p <- 101

# 1. Temporal smoothing
tsg <- smoothEO4PM(y = as.double(r[p]), t = obs_dates, tout = tout)

# 2. Extract phenological metrics
phenometrics <- phenoEO4PM(y = tsg$y, t = tsg$t, c_win = 10)

# 3. Decode Quality Flags
phenometrics <- QFbitDecoder(x = phenometrics)

# 4. Generate ggplot2 visualization
EO4PM_ggplot <- plotEO4PM(ts = tsg, pm = phenometrics)
plot(EO4PM_ggplot)
```

---

### 2. Advanced Parameter Customization

#### Customizing Temporal Smoothing (`smoothEO4PM`)
```R
tsg <- smoothEO4PM(
  y = as.double(r[p]), 
  t = obs_dates, 
  tout = tout, 
  spikes_removal = "none", 
  win_savgol = 35, 
  lambda = 100, 
  minValue = 0.05, 
  iters = 4, 
  wFUN_type = "wTSM", 
  fitMax = TRUE, 
  plot = TRUE
)
```

#### Customizing Phenometrics & Multi-Cycle Detection (`phenoEO4PM`)
```R
phenometrics <- phenoEO4PM(
  y = tsg$y, 
  t = tsg$t, 
  minValue_ylu = 0.1, 
  s_lag = 7, 
  maxExtendMonth = 3, 
  ypeak_min = 0.8, 
  minpeakdistance = 75, 
  length_min = 40, 
  r_min = 0.01, 
  r_max = 0.1, 
  rtrough_max = 0.6, 
  max_season = 2, 
  force = TRUE
)
```

---

### 3. Custom Plotting & Output (`plotEO4PM`)

```R
# Export plot directly to file without legend
plotEO4PM(
  ts = tsg, 
  pm = phenometrics, 
  legend = FALSE, 
  outfile = "output.png"
)

# Custom plot title and restricted date range
EO4PM_ggplot <- plotEO4PM(
  ts = tsg, 
  pm = phenometrics, 
  title = "Sentinel-2 LAI Phenology", 
  ylab = expression(LAI ~ (m^2 / m^2)), 
  start_date = "2021-01-01", 
  end_date = "2023-12-31"
)
plot(EO4PM_ggplot)
```

---

### 4. Batch Vector Processing (GeoPackage Points)

```R
input_points <- "/path/to/vector_points.gpkg"
output_folder <- "/path/to/output_plots"

p_in <- sf::st_read(dsn = input_points, quiet = TRUE)
raw_ts <- terra::extract(r, p_in, ID = FALSE)

dir.create(path = output_folder, showWarnings = FALSE, recursive = TRUE)

for (p in 1:nrow(raw_ts)) {
  # Temporal smoothing & TS plot export
  output_ts_plot <- normalizePath(file.path(output_folder, sprintf("TS_Profile_%05d.png", p)), mustWork = FALSE)
  tsg <- smoothEO4PM(y = as.double(raw_ts[p, ]), t = obs_dates, tout = tout, plot_filename = output_ts_plot)
  
  # Phenometrics extraction
  phenometrics <- phenoEO4PM(y = tsg$y, t = tsg$t)
  
  # Geographic metadata for plot subtitle
  p_coords <- sf::st_coordinates(sf::st_transform(p_in[p, ], 4326))
  p_title <- sprintf("Profile %05d", p)
  p_subtitle <- sprintf("Latitude: %.6f - Longitude: %.6f", p_coords[2], p_coords[1])
  output_plot <- normalizePath(file.path(output_folder, sprintf("Profile_%05d.png", p)), mustWork = FALSE)
  
  # Plot and export
  plotEO4PM(
    ts = tsg, 
    pm = phenometrics, 
    title = p_title, 
    subtitle = p_subtitle, 
    outfile = output_plot, 
    legend = FALSE, 
    ylab = expression(LAI ~ (m^2 / m^2))
  )
}
```

---

### 5. Extraction from Pre-Smoothed Daily Time Series (`TSG`)

```R
input_file_tsg <- "/path/to/S2_L3A_TSG_10m_LAI.nc"
r_tsg <- terra::rast(input_file_tsg)

p <- 1
phenometrics <- phenoEO4PM(
  y = as.double(r_tsg[p]), 
  t = as.integer(as.Date(terra::time(r_tsg), tz = "UTC"))
)

EO4PM_ggplot <- plotEO4PM(
  ts = data.frame(y = as.double(r_tsg[p]), t = as.integer(as.Date(terra::time(r_tsg), tz = "UTC"))), 
  pm = phenometrics
)
```

---

## 📊 Complete Phenological Variables Table

Below is the updated list of variables provided by EO4PM v2 (as detailed in `EO4PM_variables_table.docx`):

| Variable Name | Description / Long Name |
| :--- | :--- |
| `PM_mask` | Phenological Metrics mask |
| `VCx_season` | Seasons number |
| `VCx_SoC_date` | Start of Cycle date |
| `VCx_SoC_doy` | Start of Cycle (DOY) |
| `VCx_SoC_value` | Start of Cycle LAI value |
| `VCx_EoC_date` | End of Cycle date |
| `VCx_EoC_doy` | End of Cycle (DOY) |
| `VCx_EoC_value` | End of Cycle LAI value |
| `VCx_SoS_date` | Start of Season date |
| `VCx_SoS_doy` | Start of Season (DOY) |
| `VCx_SoS_value` | Start of Season LAI value |
| `VCx_SGS_date` | Start of Growing Season date |
| `VCx_SGS_doy` | Start of Growing Season (DOY) |
| `VCx_SGS_value` | Start of Growing Season LAI value |
| `VCx_greenup_date` | Greenup date |
| `VCx_greenup_doy` | Greenup (DOY) |
| `VCx_greenup_value` | Greenup LAI value |
| `VCx_greenup_rate` | LAI Greenup rate |
| `VCx_SMP_date` | Start of Maturity Plateau date |
| `VCx_SMP_doy` | Start of Maturity Plateau (DOY) |
| `VCx_SMP_value` | Start of Maturity Plateau LAI value |
| `VCx_PoS_date` | Peak of Season date |
| `VCx_PoS_doy` | Peak of Season (DOY) |
| `VCx_PoS_value` | Peak of Season LAI value |
| `VCx_EGS_date` | End of Growing Season date |
| `VCx_EGS_doy` | End of Growing Season (DOY) |
| `VCx_EGS_value` | End of Growing Season LAI value |
| `VCx_senescence_date` | Senescence date |
| `VCx_senescence_doy` | Senescence (DOY) |
| `VCx_senescence_value` | Senescence LAI value |
| `VCx_senescence_rate` | LAI Senescence rate |
| `VCx_EoS_date` | End of Season date |
| `VCx_EoS_doy` | End of Season (DOY) |
| `VCx_EoS_value` | End of Season LAI value |
| `VCx_seasonal_amplitude` | LAI seasonal amplitude |
| `VCx_DoS` | Duration of Season |
| `VCx_LMP` | Length of Maturity Plateau |
| `VCx_maturity_plateau_slope` | Maturity Plateau slope LAI value |
| `VCx_STI` | Seasonal Time Integrated LAI value |
| `VCx_NSTI` | Net Seasonal Time Integrated LAI value |
| `VCx_SMV` | Seasonal Mean Value LAI value |
| `VCx_PoS_AG_difference` | Peak of Season Asymmetric Gaussian difference LAI value |
| `VCx_STI_AG_difference` | Seasonal Time Integrated Asymmetric Gaussian difference LAI value |
| `VCx_QF` | Quality Flags |

---

## 🚩 Quality Flags (Bit Decoding)

Quality flags (`VCx_QF`) encode detailed diagnostic information across 15 bits:

| Bit | Quality Flag Meaning |
| :---: | :--- |
| **0** | Processed seasonal cycle |
| **1** | Pixel or date input problem |
| **2** | First derivative error |
| **3** | Identified peak date is before reference year (before 1st January) |
| **4** | No increasing period found |
| **5** | No decreasing period found |
| **6** | Gu stages identification failed |
| **7** | Estimated Gu stages are not chronologically ordered |
| **8** | Estimated Gu stages are not within input temporal range or have missing metrics values |
| **9** | Cannot detect relative minimum (SoS) using inflection point |
| **10** | Identified Gu stages are not within the analyzed time period |
| **11** | Failed to cut cycles |
| **12** | Identified more than specified `max_season` |
| **13** | Incomplete set of Phenological Metrics (e.g. due to out-of-range issues) |
| **14** | Used Asymmetric Gaussian fit after Gu phenometric estimate failure (only applied after successful cut-cycle) |

Bit flags can be decoded automatically using `QFbitDecoder()`:

```R
# Decode quality flags and append to results
phenometrics <- QFbitDecoder(x = phenometrics, append = TRUE)
```

---

## ⚖️ License

This project is licensed under the **GNU General Public License v3.0 (GPL v3)**. See the `LICENSE` file for details.
