# Data Codebook - Progresa-Oportunidades RDD Analysis

## Overview
This codebook documents all variables in the prepared datasets for the Progresa-Oportunidades regression discontinuity design analysis.

---

## Main Dataset: `merged_entry_iml_rg.csv`

### Identification Variables

| Variable | Type | Description | Values/Range |
|----------|------|-------------|--------------|
| `cve_edo` | Integer | State code | 1-32 |
| `cve_mun` | Integer | Municipality code | 1-999 |
| `cve_loc` | Integer | Locality code | 1-9999 |
| `loc_id` | Character | Unique locality identifier | 9-digit string (EEMMMMLLLL) |
| `cve_mun_full` | Character | State-municipality code | 5-digit string (EEMMM) |

### Program Entry Variables

| Variable | Type | Description | Values/Range |
|----------|------|-------------|--------------|
| `entry_year` | Integer | Year locality entered program | 1995-2010 |
| `entry_year_30` | Integer | Entry year from source 1 | 1995-2010 |
| `entry_year_enca` | Integer | Entry year from source 2 | 1995-2010 |
| `entry_source` | Character | Data source for entry year | "entry_30", "entry_enca", "both", "none" |
| `post_2002` | Binary | Entered after 2002 (Oportunidades) | 0 = No, 1 = Yes |
| `post_2003` | Binary | Entered after 2003 | 0 = No, 1 = Yes |
| `entry_period` | Character | Categorized entry period | See categories below |
| `program_phase` | Character | Program at entry | "Progresa", "Oportunidades", "Not Enrolled" |

#### Entry Period Categories:
- "1997-1998 (Progresa Phase I)"
- "1999-2000 (Progresa Phase II)"
- "2001-2002 (Transition)"
- "2003-2004 (Oportunidades Early)"
- "2005-2006 (Oportunidades Expansion)"
- "2007+ (Late Entry)"
- "Never Entered"

### Marginalization Index Variables

| Variable | Type | Description | Range |
|----------|------|-------------|-------|
| `iml_95` | Numeric | Marginalization index 1995 | -3 to 3 (approx.) |
| `iml_00` | Numeric | Marginalization index 2000 | -3 to 3 (approx.) |
| `iml_05` | Numeric | Marginalization index 2005 | -3 to 3 (approx.) |

### Demographic Variables

| Variable | Type | Description | Units |
|----------|------|-------------|-------|
| `pobtot_95` | Integer | Total population 1995 | Number of people |
| `pobtot_00` | Integer | Total population 2000 | Number of people |
| `pobtot_05` | Integer | Total population 2005 | Number of people |

### Infrastructure & Service Variables (1995)

| Variable | Type | Description | Units |
|----------|------|-------------|-------|
| `analf_95` | Numeric | Illiteracy rate 1995 | Percentage (0-100) |
| `sinag_95` | Numeric | Households without water 1995 | Percentage (0-100) |
| `sinel_95` | Numeric | Households without electricity 1995 | Percentage (0-100) |
| `pisotier_95` | Numeric | Households with dirt floors 1995 | Percentage (0-100) |

### Infrastructure & Service Variables (2000)

| Variable | Type | Description | Units |
|----------|------|-------------|-------|
| `analf_00` | Numeric | Illiteracy rate 2000 | Percentage (0-100) |
| `sinag_00` | Numeric | Households without water 2000 | Percentage (0-100) |
| `sinel_00` | Numeric | Households without electricity 2000 | Percentage (0-100) |
| `pisotier_00` | Numeric | Households with dirt floors 2000 | Percentage (0-100) |

### Infrastructure & Service Variables (2005)

| Variable | Type | Description | Units |
|----------|------|-------------|-------|
| `analf_05` | Numeric | Illiteracy rate 2005 | Percentage (0-100) |
| `sinag_05` | Numeric | Households without water 2005 | Percentage (0-100) |
| `sinel_05` | Numeric | Households without electricity 2005 | Percentage (0-100) |
| `pisotier_05` | Numeric | Households with dirt floors 2005 | Percentage (0-100) |
| `sinrefri_05` | Numeric | Households without refrigerator 2005 | Percentage (0-100) |

### Regional Classification

| Variable | Type | Description | Values |
|----------|------|-------------|--------|
| `region41` | Integer | Regional classification code | 1-53 (excluding 17-26, 29-30) |

---

## RDD Analysis Dataset: `rdd_analysis_dataset.csv`

All variables from main dataset plus:

| Variable | Type | Description | Range |
|----------|------|-------------|-------|
| `optimal_cutoff` | Numeric | Region-specific optimal IML cutoff | -1.5 to -0.5 (approx.) |
| `iml_95_normalized` | Numeric | Distance from regional cutoff | -3 to 3 (approx.) |
| `above_cutoff` | Binary | Locality above regional cutoff | 0 = No, 1 = Yes |
| `dist_to_cutoff` | Numeric | Absolute distance to cutoff | 0 to 3 (approx.) |

---

## Panel Dataset: `iml_panel_90_10.csv`

| Variable | Type | Description | Values/Range |
|----------|------|-------------|--------------|
| `mge_id` | Character | Locality identifier (alternative coding) | Varies |
| `year` | Integer | Year of observation | 1990, 1995, 2000, 2005, 2010 |
| `entry_year` | Integer | Year locality entered program | 1995-2010 |
| `region41` | Integer | Regional classification | 1-53 |
| `post_2002` | Binary | Entered after 2002 | 0 = No, 1 = Yes |
| `post_2003` | Binary | Entered after 2003 | 0 = No, 1 = Yes |
| `treated` | Binary | Locality treated in given year | 0 = No, 1 = Yes |

Plus various IML component variables that vary by year.

---

## Regional Cutoffs Dataset: `regional_cutoffs.csv`

| Variable | Type | Description | Values/Range |
|----------|------|-------------|--------------|
| `region` | Integer | Region identifier | 1-53 |
| `optimal_cutoff` | Numeric | Optimal IML 1995 cutoff for RDD | -1.5 to -0.5 |
| `n_obs` | Integer | Number of localities in region | Varies |
| `rd_estimate` | Numeric | RDD treatment effect estimate | -1 to 0 |
| `se` | Numeric | Standard error of estimate | 0 to 1 |
| `pvalue` | Numeric | P-value for RDD estimate | 0 to 1 |
| `bandwidth` | Numeric | Optimal bandwidth used | 0.1 to 2 |
| `significance` | Character | Significance level | "", "*", "**", "***" |

---

## Missing Data Codes

- `NA`: Missing value (data not available)
- Blank cells in CSV: Missing value

---

## Data Quality Notes

1. **Encoding**: Original CSV files use Latin-1 encoding for Spanish characters
2. **Municipality Codes**: Standardized to 3-digit format with leading zeros
3. **Locality IDs**: Created as 9-digit concatenation: EEMMMMLLLL (Estado-Municipio-Localidad)
4. **Year Ranges**: 
   - IML scores: 1990, 1995, 2000, 2005, 2010
   - Program entry: 1995-2010
5. **Regional Coverage**: Regions 17-26, 29-30 are not used in current classification

---

## Variable Transformations

### From Raw to Clean

1. **IML 1995 Variables**:
   - `cve_ent` → `cve_edo`
   - `pobtot95` → `pobtot_95`
   - `indice95` → `iml_95`

2. **IML 2000 Variables**:
   - `cve_ent` → `cve_edo`
   - `imlr2000` → `iml_00`
   - `anal_iml` → `analf_00`

3. **IML 2005 Variables**:
   - `cve_ent` → `cve_edo`
   - `im_05` → `iml_05`
   - `anal05` → `analf_05`

### Derived Variables

- `post_2002`: Binary indicator derived from `entry_year > 2002`
- `iml_95_normalized`: Calculated as `iml_95 - optimal_cutoff`
- `above_cutoff`: Binary indicator from `iml_95_normalized > 0`
- `treated` (panel): Binary based on `year >= entry_year`

---

## Usage Guidelines

1. **Primary Key**: Use `loc_id` for cross-sectional analyses, `mge_id` + `year` for panel analyses
2. **Treatment Variable**: Use `post_2002` for main RDD specifications
3. **Running Variable**: Use `iml_95` for raw scores, `iml_95_normalized` for normalized RDD
4. **Sample Selection**: Filter on `!is.na(optimal_cutoff)` for valid RDD sample
5. **Regional Analysis**: Always include `region41` fixed effects when pooling regions

---

