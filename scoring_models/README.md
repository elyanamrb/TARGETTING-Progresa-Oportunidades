# Progresa/Oportunidades Scoring Model OCR — Technical Documentation

**Project:** Machine-readable extraction of poverty-scoring SPSS syntax from
`41_REGIONES_Comandos.pdf` (284 pages, scanned)
**Data source:** Mexico's Progresa/Oportunidades conditional cash transfer program
**Pipeline:** PDF → OCR (Claude Vision API) → structured Excel catalogs

---

## 1. Source Files

Three PDF files make up the primary source material for this project.

### 1.1 `41_REGIONES_Comandos.pdf` — Main input (284 pages, 5.9 MB, scanned)

- **Primary input** to the OCR pipeline
- Contains the SPSS syntax implementing the poverty-scoring model for each region
  (variable weights, IF-categorical scores, DISCRIMINANT procedures, thresholds)
- **Covers 21 of 41 regions:** R1, R2, R3, R6, R8, R9, R10, R12, R13, R15, R16,
  R28, R32, R34, R39, R40, R43, R44, R45, R46, R47
- The remaining 20 regions (R4, R5, R7, R11, R14, R27, R31, R33, R35–R38, R41,
  R42, R48–R53) are not present in this file; their SPSS code may exist in a
  separate document

### 1.2 `41_REGIONES_Cutoffs.pdf` — Threshold table (1 page, 12 KB)

- Single-page reference table listing the poverty classification threshold
  (*punto de corte*) for all 41 regions across all phases (Fase 1–9)
- Each cell records the model name (`aa`, `yy`, `bb`, `primero`, `segundo`,
  `s/nombre`, `zz`) and the corresponding threshold value
- Covers all 41 regions, including the 20 regions absent from `Comandos.pdf`
- The six model-pattern classification (Types A–F) in Section 6.2 is based on
  analysis of this file

### 1.3 `41_REGIONES_Modelos_1997_2000.pdf` — Model documentation (15 pages, 3.4 MB)

- Detailed documentation of the scoring models across all 41 regions
- Likely provides broader region coverage than `Comandos.pdf`
- **Not yet processed** by this pipeline; if needed, the same three-step OCR
  workflow can be applied

### Relationship Between the Three Files

| File | Content | Regions | Pipeline status |
|---|---|---|---|
| `Comandos.pdf` | SPSS execution code | 21 of 41 | Processed (this pipeline) |
| `Cutoffs.pdf` | Threshold summary table | All 41 | Reference only (not OCR'd) |
| `Modelos_1997_2000.pdf` | Detailed model descriptions | Likely all 41 | Not yet processed |

---

## 2. Repository Overview

| File | Role |
|---|---|
| `ocr_progresa.py` | **Step 1** — Convert PDF pages to text via Claude Vision API |
| `parse_to_excel.py` | **Step 2** — Parse OCR text into a first-pass Excel summary |
| `build_catalog.py` | **Step 3** — Split text into block files and build the final Excel catalog |
| `output/full_ocr_text.txt` | Concatenated OCR text for all 284 pages |
| `output/pages_NNN_NNN.txt` | Per-batch OCR text (intermediate output) |
| `output/progress.json` | Checkpoint file (allows resuming interrupted OCR runs) |
| `blocks/*.txt` | One text file per region × fase × modelo block |
| `progresa_models.xlsx` | Intermediate Excel (Step 2 output) |
| `progresa_scoring_catalog.xlsx` | Final Excel catalog (Step 3 output) |

---

## 3. Pipeline Scripts

### 3.1 `ocr_progresa.py` — PDF OCR (Step 1)

**Purpose:** Sends each page of the scanned PDF to the Claude Vision API and saves the
recognized text.

**Usage:**
```
python ocr_progresa.py --pdf 41_REGIONES_Comandos.pdf --output-dir ./output
```

**Key options:**

| Option | Default | Description |
|---|---|---|
| `--pdf` | *(required)* | Path to input PDF |
| `--output-dir` | `./output` | Directory for OCR output files |
| `--batch-size` | `4` | Pages per API call |
| `--start-page` | `1` | Resume from this page |
| `--model` | `claude-sonnet-4-...` | Claude model ID |

**How it works:**
1. Converts PDF pages to images using `pdf2image` (requires Poppler)
2. Sends batches of page images to the Claude Vision API with an OCR prompt
3. Appends recognized text to `full_ocr_text.txt`, prefixing each page with a
   marker: `--- Page N/284 ---`
4. Saves `progress.json` after each batch so runs can be resumed with `--start-page`

**Cost / runtime estimate:** ~$5–6 and ~30 min for all 284 pages at batch size 4.

---

### 3.2 `parse_to_excel.py` — First-pass Parser (Step 2)

**Purpose:** Splits the full OCR text into region/model blocks by detecting SPSS
header patterns, extracts key parameters, and writes an intermediate Excel file.

**Usage:**
```
python parse_to_excel.py --input ./output/full_ocr_text.txt --output progresa_models.xlsx
```

**How it works:**
1. Detects `--- Page N/284 ---` markers to build a character-position → page-number
   index
2. Finds region headers matching patterns like `REGIÓN 1. LOS ALTOS`, `REGION 3`,
   or `PUNTAJE PARA LA REGION 3`
3. Within each header context, looks for `FASE N` and `MODELO yy/aa/bb`
4. Splits the text at header positions to form blocks; discards blocks without
   `COMPUTE` or `DISCRIMINANT` (i.e., non-SPSS content)
5. For each block, extracts: continuous weights, IF-categorical scores, threshold,
   SYSMIS defaults, DISCRIMINANT variables, PRIORS, and the CALIF formula
6. Writes two Excel sheets: **Summary** (one row per block) and
   **Variables_Weights** (one row per variable per block)

---

### 3.3 `build_catalog.py` — Final Catalog Builder (Step 3)

**Purpose:** Uses the page-range information in the Step 2 Excel to slice the OCR
text precisely into per-block text files, re-parses them with an extended parser,
and produces a richly formatted four-sheet Excel catalog.

**Usage:**
```
python build_catalog.py \
    --ocr-text ./output/full_ocr_text.txt \
    --initial-excel ./progresa_models.xlsx \
    --blocks-dir ./blocks \
    --output ./progresa_scoring_catalog.xlsx
```

**How it works:**
1. Reads `page_start` / `page_end` from the Step 2 Excel Summary sheet
2. Slices `full_ocr_text.txt` at `--- Page N/284 ---` markers to extract the
   exact text for each block; writes one `.txt` file per block to `blocks/`
3. Re-parses each block with an extended parser that also captures composite
   variables, second-stage CALIF formulas, viejitos adjustments, PERCEP
   preprocessing, and vehicle filters
4. Writes `progresa_scoring_catalog.xlsx` with four sheets (see Section 5)

---

## 4. Output: `progresa_models.xlsx` (Step 2)

### Sheet: Summary

One row per detected region × fase × modelo block (104 rows).

| Column | Type | Description |
|---|---|---|
| `region` | integer | Region number (1–47; not all 41 values present — see note) |
| `region_name` | string | Region name as read from OCR (may include noise) |
| `fase` | string | Phase number (`1`–`8`); empty if not detected |
| `modelo` | string | Model identifier: `yy`, `aa`, or `bb`; empty if not detected |
| `page_start` | integer | First PDF page of this block |
| `page_end` | integer | Last PDF page of this block |
| `n_weight_vars` | integer | Number of continuous variables with COMPUTE weights |
| `n_if_vars` | integer | Number of categorical variables assigned via IF statements |
| `threshold` | float | Poverty classification threshold (value N in `xxcalif < N`) |
| `priors` | string | DISCRIMINANT priors (e.g., `EQUAL` or `.42 .58`) |
| `discriminant_vars` | string | Space-separated variable list from `/VARIABLES=` |
| `calif_formula` | string | Right-hand side of `COMPUTE xxcalif = ...` |

> **Note on region coverage:** The PDF contains 41 regions but only 21 are detected
> by the header regex. The remaining 20 regions use non-standard header formats
> (e.g., headers in Spanish only, or part of multi-page layouts) and are present
> in the raw OCR text but not parsed into separate rows. This is a known limitation
> of the current regex-based parser and does not reflect a problem with OCR quality.

### Sheet: Variables_Weights

One row per variable per block (~1,354 rows).

| Column | Type | Description |
|---|---|---|
| `region` | integer | Region number |
| `region_name` | string | Region name |
| `fase` | string | Phase number |
| `modelo` | string | Model identifier |
| `page_start` | integer | First PDF page of this block |
| `page_end` | integer | Last PDF page of this block |
| `variable` | string | SPSS variable name (uppercase) |
| `weight` | float | Linear weight from `COMPUTE AAvar = weight * variable` |
| `sysmis_default` | float | System-missing replacement value from `RECODE` |
| `score_when_N` | float | IF-categorical score when `variable = N` (one column per value of N) |

---

## 5. Output: `progresa_scoring_catalog.xlsx` (Step 3)

### Sheet: Catalog

One row per block. The primary reference sheet.

| Column | Description |
|---|---|
| `block_id` | Sequential identifier (B001–B104) |
| `file` | Filename of the corresponding block text file in `blocks/` |
| `region` | Region number |
| `fase` | Phase number |
| `modelo` | Model identifier |
| `prefix` | Detected variable prefix (`yy`, `aa`, `bb`, `xz`, `zz`, or `?`) |
| `pages` | Page range string, e.g., `p001-p004` |
| `n_continuous` | Number of continuous (linearly weighted) variables |
| `n_binary_if` | Number of binary (2-category) IF variables |
| `n_multi_if` | Number of multi-category (3+ category) IF variables |
| `n_composites` | Number of composite dummy variables (constructed via DO IF / RECODE) |
| `threshold` | Poverty classification threshold |
| `threshold_on` | Name of the CALIF variable the threshold applies to |
| `has_viejitos` | `Y` if a viejitos (elderly/small-household) score adjustment is present |
| `calif1_extra_adj` | Additional adjustments to the CALIF1 variable beyond viejitos |
| `calif2_additions` | Second-stage additions to CALIF (e.g., `+ viejitos + enseres`) |
| `vehi_filter` | `Y` / rule text if a vehicle ownership filter overrides poverty status |
| `has_percep_preproc` | `Y` if PERCEP=0 / DEPEN=99 preprocessing is applied |
| `priors` | DISCRIMINANT prior probabilities |

### Sheet: Formulas

One row per block. Focus on the scoring formula structure.

| Column | Description |
|---|---|
| `block_id` | Block identifier |
| `file` | Block text filename |
| `region`, `fase`, `modelo` | Block identifiers |
| `calif_formula` | Full text of `COMPUTE xxcalif = ...` (first stage) |
| `calif2_additions` | Second-stage CALIF update (the additive terms after `xxcalif + ...`) |
| `viejitos_rule` | IF condition that triggers the viejitos score adjustment |
| `composites` | Composite variable names and their construction conditions |

### Sheet: Variables

One row per variable per block (~1,834 rows). Extended version of the Step 2
Variables_Weights sheet.

| Column | Description |
|---|---|
| `block_id` | Block identifier |
| `file` | Block text filename |
| `region`, `fase`, `modelo` | Block identifiers |
| `variable` | SPSS variable name |
| `type` | `continuous` (COMPUTE weight) or `categorical` (IF scores) |
| `weight` | Linear weight coefficient (continuous variables only) |
| `sysmis_default` | RECODE SYSMIS replacement value |
| `if_N` | Score assigned when variable = N (categorical variables; one column per category) |

### Sheet: Discriminant

One row per block. Contains the DISCRIMINANT procedure inputs.

| Column | Description |
|---|---|
| `block_id` | Block identifier |
| `file` | Block text filename |
| `region`, `fase`, `modelo` | Block identifiers |
| `priors` | Prior probabilities (`EQUAL` or numeric pair) |
| `disc_vars` | Variables listed in `/VARIABLES=` of the DISCRIMINANT command |

---

## 6. Model Structure Analysis

### 6.1 What the Model Names Mean

Each region's scoring algorithm classifies households into poor / non-poor through a
sequential discriminant analysis pipeline. The model name (yy / aa / bb / primero /
segundo / s-nombre / zz) identifies the role of each stage:

#### `yy` (= `primero` = Model 1 — Coarse screen)

- **Applied first** to all households in the region
- Uses dummy variable groups (e.g., HH1–HH6, ESCJ1–ESCJ4, OCUPJ0–OCUPJ5\_9)
  with categorical IF-score assignments
- `PRIORS = EQUAL` (prior probability of poor = non-poor = 50%)
- Threshold is **higher** (e.g., `yyCALIF < 558`), casting a wide net for poverty
  candidates
- Households classified poor by `yy` are flagged `yypobre = 1` and **exit the
  pipeline**; they do not proceed to the `aa` model

#### `aa` (= `segundo` = Model 2 — Main classification)

- Applied to households **not classified poor by `yy`**
- Uses **continuous linear weights**: `COMPUTE AAVAR = coefficient * variable`
  summed into `AACALIF`
- `PRIORS ≈ .42 .58` (non-poor 42%, poor 58% — asymmetric)
- Threshold is **lower** (e.g., `AACALIF < 368`), setting a stricter standard
- Includes optional **viejitos adjustment**: small or elderly households receive a
  score boost
- Households classified poor by `aa` (and not already poor by `yy`) become `nvopob`
  (new poor)

#### `bb` (= Model 3 — Alternative/supplementary)

- Uses the **same variables as `aa`** but with different coefficients and threshold
- Identifies households poor under `bb` but not under `aa` (`nvo_pob`)
- Households already classified poor by `aa` are excluded:
  `IF (aapobre = 1) nvo_pob = 0`
- Functions as a second sweep to capture households that `aa` missed

#### `s/nombre` (unnamed)

- A single scoring model with no named stage
- Applied to regions that do not use the two-stage yy → aa pipeline

#### `zz` (Region 2 only)

- Appears to serve the same first-stage role as `yy`, used in one region only

---

### 6.2 Six Model Pattern Types Across 41 Regions

Analysis of the SPSS blocks identifies six structural patterns:

| Type | Pattern | Regions (n) | Description |
|---|---|---|---|
| **A** | aa / bb only | ~17 | Most common. `aa` is the main model; `bb` is an alternate model for some phases. No `yy` first-stage screen — households go directly to discriminant analysis. |
| **B** | yy → aa → bb | ~4 | R1, R5, R12, R15. Full three-model pipeline: coarse screen (`yy`) → main classification (`aa`) → supplementary sweep (`bb`). |
| **C** | primero / segundo | ~8 | R4, R7, R11, R27, R33–R37. Functionally identical to Type B but uses Spanish stage names instead of `yy` / `aa`. |
| **D** | s/nombre only | ~4 | R16, R39, R50, R52. Single unnamed model for 1–3 phases; no multi-stage pipeline. |
| **E** | Mixed s/nombre + aa/bb | ~7 | R3, R6, R8, R9, R10, R13, R40. Early or late phases use unnamed models; middle phases use `aa`/`bb`. |
| **F** | zz variant | 1 | R2 only. Uses `zz` prefix in place of `yy` for the first-stage screen. |

**Functional equivalences across naming conventions:**

| Role | Name variants |
|---|---|
| First-stage coarse screen | `yy`, `primero`, `zz` |
| Main classification | `aa`, `segundo` |
| Supplementary/alternative | `bb` |
| Single-stage (unnamed) | `s/nombre` |

---

## 7. OCR Quality Notes

- OCR was performed by Claude Vision API (`claude-sonnet-4`), one batch of 4 pages
  per API call
- Numeric values (coefficients, thresholds) are generally well-recognized
- Variable names in uppercase are reliably extracted
- Some header lines contain OCR noise (extra spaces, accents dropped), which is why
  the parser uses flexible regex patterns
- The `blocks/` text files are the authoritative OCR output for each block and can
  be inspected directly to verify any value in the Excel outputs
- 21 of 41 regions are detected by the automated header parser; the remaining 20
  are present in the raw OCR text but require manual identification of their header
  format

---

*Documentation generated for the Progresa/Oportunidades OCR extraction pipeline.
Data source: `41_REGIONES_Comandos.pdf` (internal Progresa technical document).*
