# Folder `data`

This folder contains the data needed to run the project scripts.

**Important:**

- The file `RF_data_Carbon_ForestFloor.csv` is an example with a structure similar to the real data.
- The original data is not shared publicly for confidentiality reasons.
- The code was designed to work fully with this simulated version.

## 🔒 Confidentiality Notice:
If you have questions or need access to the real data, you can contact the repository author. Please contact:

**Zaira Rosario Pérez-Vázquez**  
zairpv@gmail.com
---

# 📂 File description

## `RF_data_Carbon_ForestFloor.csv`

### Description
This file includes simulated carbon stock estimates in the **forest floor** divided into two layers:
- `L` (litter)
- `FH` (fermentation-humus)
Measurements were recorded for **918 observations** across three sampling campaigns: **2013, 2018, and 2023**.

### Variable Dictionary

 Column Name        | Description                                                        | Type        |
|--------------------|--------------------------------------------------------------------|-------------|
| `YEAR`             | Year of measurement                                                | Numeric     |
| `PSP`              | Sampling plot                                     | Numeric     |
| `SITE`             | Site identifier within each PSP                                    | Numeric     |
| `LAYER`            | Organic layer: `L` (litter) or `FH` (fermentation-humus)           | Character   |
| `C_STOCKS`         | Simulated carbon stock in the layer (Mg C / ha)                    | Numeric     |
| `UTMX`, `UTMY`     | UTM coordinates of the site (X, Y)                                 | Numeric     |
| `STAND_AGE`        | Stand age at the time of sampling (years)                          | Numeric     |
| `BASAL_AREA`       | Tree basal area (m²/ha)                                            | Numeric     |
| `DOM_HEIGHT`       | Dominant height (m)                                                | Numeric     |
| `SPECIES_RICHNESS` | Number of species recorded at the site                            | Numeric     |
| `SHANNON_INDEX`    | Shannon diversity index                                            | Numeric     |
| `GAP_FRACTION`     | % of canopy openness (gap fraction)                                | Numeric     |
| `CANOPY_COVER`     | % of total canopy cover                                            | Numeric     |
| `ELEVATION`        | Elevation above sea level (m)                                      | Numeric     |
| `SLOPE`            | Slope of the terrain (%)                                           | Numeric     |
| `ASPECT`           | Aspect in degrees (0–360)                                          | Numeric     |
| `CONDITION`        | Forest management condition: `Managed` or `Unmanaged`              | Character   |
| `LAYER_NUM`        | Numerical encoding of the layer (`L` = 1, `FH` = 2)                | Numeric     |
| `CON_NUM`          | Numerical encoding of condition (`Managed` = 1, `Unmanaged` = 2)   | Numeric     |

---
