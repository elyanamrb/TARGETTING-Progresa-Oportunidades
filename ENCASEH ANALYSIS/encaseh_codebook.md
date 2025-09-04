# ENCASEH Codebook
## Encuesta de Características Socioeconómicas de los Hogares (2000-2004)

---

## Overview
The ENCASEH (Household Socioeconomic Characteristics Survey) is the primary dataset used to determine household eligibility for Mexico's Progresa/Oportunidades conditional cash transfer program. This codebook documents variables across survey years 2000, 2001, 2002, and 2004.

**Note:** 2003 data is missing from the series.

---

## File Structure

| Year | Filename | Observations | Variables |
|------|----------|--------------|-----------|
| 2000 | fams_00.dta | ~37,000+ | 32 |
| 2001 | Familias_2001.dta | ~1,600,000+ | 24 |
| 2002 | Familias_2002.dta | ~500,000+ | 24 |
| 2004 | Familias_2004.dta | ~1,800,000+ | 24 |

---

## Core Variables (All Years)

### Identification Variables

| Variable | Description | Type | Notes |
|----------|-------------|------|-------|
| **FAM_ID** | Family identifier | String | Unique family ID |
| **FOLIO** | Household folio number | String | Administrative identifier |
| **ENT_FED** | State code | Numeric (1-32) | Federal entity code |
| **MUNICI** | Municipality code | Numeric | Municipality within state |
| **LOCALI** | Locality code | Numeric | Locality within municipality |
| **REGION** | Program region | Numeric (1-41/53) | Administrative regions for program implementation |

### Program Variables

| Variable | Description | Type | Values | Notes |
|----------|-------------|------|--------|-------|
| **PUN_INC** | Poverty score | Continuous | Variable range | See scoring system notes |
| **RESULTADO** | Program acceptance outcome | Categorical | 1, 3 | 1 = Accepted into program, 3 = Not accepted |
| **TIPO_CAP** | Targeting method (2001+) | Categorical | 4, 5, 6 | 4 = Ventanilla (self-selection), 5 = Rural census, 6 = Urban census |

---

## Asset Variables

### Standard Assets (All Years)

| Variable | Description | Values | Coding |
|----------|-------------|--------|--------|
| **LICUADORA** | Blender | 0, 1, 2, 9, -1 | 1 = Has asset, 2 = Doesn't have, 9 = Missing/refused, 0/-1 = Special cases |
| **REFRI** | Refrigerator | Same as above | Same coding scheme |
| **EST_GAS** | Gas stove | Same as above | Same coding scheme |
| **CAL_AGUA** | Water heater | Same as above | Same coding scheme |
| **RADIO** | Radio | Same as above | Same coding scheme |
| **TOCADISC** | Turntable/record player | Same as above | Same coding scheme |
| **TELEVIS** | Television | Same as above | Same coding scheme |
| **VIDEOCAS** | VCR | Same as above | Same coding scheme |
| **LAVADORA** | Washing machine | Same as above | Same coding scheme |
| **VENTILA** | Fan | Same as above | Same coding scheme |
| **AUTO_PROP** | Own car | Same as above | Same coding scheme |
| **CAMION_PRO** | Own truck/bus | Same as above | Same coding scheme |
| **TRACTOR** | Tractor (2000, 2004) | Same as above | Not present in 2001-2002 |

### SIR Variables (2000 Only)
In 2000, each asset variable has a corresponding "SIR_" variable (e.g., SIR_LIC, SIR_REFR). The exact meaning is undocumented but likely relates to asset condition or alternative measurement.

---

## Key Findings and Data Characteristics

### 1. Dual Scoring System (PUN_INC)

The poverty score underwent major changes:

| Year | System 1 (Normal) | System 2 (High) | Notes |
|------|-------------------|-----------------|-------|
| 2000 | 83.1% | 16.9% | Dual system active |
| 2001 | 99.0% | 1.0% | Transition year |
| 2002 | 100% | 0% | Single system |
| 2004 | 100% | 0% | Single system |

**System 1 (Normal):** PUN_INC typically < 10
**System 2 (High):** PUN_INC > 100 (reaching up to 1,294)

### 2. PUN_INC Distribution Statistics

| Year | Mean | Median | SD | Min | Max |
|------|------|--------|----|----|-----|
| 2000 | 98.73 | 0.56 | 227.11 | -27.88 | 1,294.33 |
| 2001 | 7.23 | 1.71 | 61.70 | -640.86 | 1,914.84 |
| 2002 | 1.52 | 1.47 | 1.05 | -2.08 | 7.84 |
| 2004 | 1.16 | 1.06 | 1.09 | -1.97 | 8.94 |

### 3. Targeting Method Evolution (TIPO_CAP)

| Year | Ventanilla | Rural Census | Urban Census |
|------|------------|--------------|--------------|
| 2001 | 0% | 64.8% | 35.2% |
| 2002 | 0% | 35.1% | 64.9% |
| 2004 | 18.6% | 42.0% | 39.4% |

**Key Finding:** Ventanilla (self-selection) mechanism introduced in 2004.

### 4. Acceptance Rates by Targeting Method (2004)

| Method | Acceptance Rate |
|--------|----------------|
| Rural Census | 73.4% |
| Urban Census | 51.2% |
| Ventanilla | 51.1% |

---

## Regional Structure

- Program divided into 41-53 administrative regions
- Regions do not necessarily respect state borders
- Each municipality assigned to one region
- Region assignment available through external file: `regiones41.dta`

### Regional Assignment Fields
- `cve_edo`: State code
- `cve_mun`: Municipality code  
- `region41`: Assigned region number

---

## Important Analysis Notes

### 1. PCA Analysis Results
- **PUN_INC is NOT primarily asset-based**
- Correlation between asset index and PUN_INC: 0.001 to 0.464
- Despite excellent KMO values (0.784-0.975)
- Suggests complex poverty measurement beyond assets

### 2. Missing Data Patterns
- Significant missing regions for Ventanilla cases
- Asset variables have various missing codes (0, -1, 9)
- RESULTADO primarily uses values 1 and 3 (not 0 and 1)

### 3. Sample Size Considerations
- Very large datasets (up to 1.8M observations)
- Sampling may be necessary for computational efficiency
- Regional analysis requires minimum thresholds for validity

### 4. Temporal Inconsistencies
- Variable availability changes across years
- Scoring system undergoes major transformation
- Targeting methods evolve significantly

---

## Data Quality Issues

1. **Inconsistent Asset Coding:** Values vary across years (0, 1, 2, 9, -1)
2. **Dual Scoring System:** Requires careful handling in 2000-2001
3. **Missing 2003 Data:** Gap in temporal analysis
4. **Regional Assignment:** Many observations lack region codes initially
5. **Ventanilla Documentation:** Limited information on self-selection criteria

---

## Recommended Data Preparation

1. **Standardize Asset Coding:**
   - Recode to binary: 1 = has asset, 0 = doesn't have
   - Treat 9, 0, -1 as missing/special cases

2. **Handle Dual Scoring:**
   - Create indicator for scoring system
   - Consider separate analyses or filtering

3. **Assign Regions:**
   - Merge with `regiones41.dta` using municipality keys
   - Create muni_key = paste(ENT_FED, "_", MUNICI)

4. **Create Analysis Variables:**
   - `accepted` = (RESULTADO == 1)
   - `scoring_system` = (PUN_INC > 100)
   - `is_ventanilla` = (TIPO_CAP == 4)

---

## Citation

When using ENCASEH data, cite as:
```
ENCASEH (Encuesta de Características Socioeconómicas de los Hogares) [2000-2004]. 
Mexico: [Appropriate Government Agency].
Used for evaluation of Progresa/Oportunidades program targeting.
```

---

## Contact and Updates

For questions about this codebook or to report data issues:
- Created: 2024
- Last Updated: [Current Date]
- Version: 1.0

---

## Appendix: State Codes

| Code | State Name |
|------|------------|
| 1 | Aguascalientes |
| 2 | Baja California |
| 3 | Baja California Sur |
| 4 | Campeche |
| 5 | Coahuila |
| 6 | Colima |
| 7 | Chiapas |
| 8 | Chihuahua |
| 9 | Ciudad de México |
| 10 | Durango |
| 11 | Guanajuato |
| 12 | Guerrero |
| 13 | Hidalgo |
| 14 | Jalisco |
| 15 | México |
| 16 | Michoacán |
| 17 | Morelos |
| 18 | Nayarit |
| 19 | Nuevo León |
| 20 | Oaxaca |
| 21 | Puebla |
| 22 | Querétaro |
| 23 | Quintana Roo |
| 24 | San Luis Potosí |
| 25 | Sinaloa |
| 26 | Sonora |
| 27 | Tabasco |
| 28 | Tamaulipas |
| 29 | Tlaxcala |
| 30 | Veracruz |
| 31 | Yucatán |
| 32 | Zacatecas |
