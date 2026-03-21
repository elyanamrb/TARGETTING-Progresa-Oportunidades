# Progresa SPSS Syntax OCR Project

## Purpose
OCR and structure `41_REGIONES_Comandos.pdf` (284 pages, scanned PDF).
This PDF contains the poverty-scoring models for each region of Mexico's
Progresa/Oportunidades program, written as SPSS syntax.

## Pipeline Overview

```
[Step 1] scripts/ocr_progresa.py    : PDF → output/full_ocr_text.txt (Claude API vision OCR)
[Step 2] scripts/parse_to_excel.py  : output/full_ocr_text.txt → progresa_models.xlsx (block boundary detection)
[Step 3] scripts/build_catalog.py   : output/full_ocr_text.txt + progresa_models.xlsx → blocks/ + catalog.xlsx
```

## What To Do

### If full_ocr_text.txt and progresa_models.xlsx already exist
Steps 1 and 2 are already complete. Only run Step 3:

```bash
pip install pandas openpyxl
python scripts/build_catalog.py \
    --ocr-text output/full_ocr_text.txt \
    --initial-excel progresa_models.xlsx \
    --blocks-dir ./blocks \
    --output progresa_scoring_catalog.xlsx
```

### To redo everything from the PDF
```bash
pip install anthropic pdf2image Pillow openpyxl pandas
# macOS: brew install poppler / Linux: apt install poppler-utils

# Step 1: OCR (requires ANTHROPIC_API_KEY)
python scripts/ocr_progresa.py --pdf source_pdf/41_REGIONES_Comandos.pdf --output-dir ./output

# Step 2: Block boundary detection
python scripts/parse_to_excel.py --input ./output/full_ocr_text.txt --output progresa_models.xlsx

# Step 3: Block splitting + catalog creation
python scripts/build_catalog.py \
    --ocr-text ./output/full_ocr_text.txt \
    --initial-excel progresa_models.xlsx \
    --blocks-dir ./blocks \
    --output progresa_scoring_catalog.xlsx
```

## Step Details

### Step 1: ocr_progresa.py
- Sends 284 pages to the Claude API vision in batches of `--batch-size` (default 4)
- Progress is tracked in `output/progress.json`; interrupted runs resume automatically
- Batch size affects only API efficiency, not OCR accuracy

### Step 2: parse_to_excel.py
- Parses the OCR text with regex to identify block boundaries (page_start/page_end) per region × model
- Headers have significant spelling variation; if parse results are sparse, the regex patterns may need adjustment

### Step 3: build_catalog.py
- Uses page_start/page_end from Step 2 to split the OCR text into 104 blocks
- Extracts parameters from each block via regex and writes a 4-sheet Excel catalog

## Outputs

### blocks/ directory (104 text files)
- Naming convention: `B001_R01_F2_yy_p001-004.txt`
  - `B001` = sequential ID (matches the ID column in Excel)
  - `R01` = region number
  - `F2` = phase (FX = unknown)
  - `yy` = model name (XX = unknown)
  - `p001-004` = PDF page range
- Each file contains the complete SPSS syntax for that model

### progresa_scoring_catalog.xlsx (4 sheets)
- **Catalog**: Classification overview. Variable composition (continuous weights / binary IF / multi-category IF / composite variables), thresholds, viejitos adjustment flag, two-stage calif correction, vehicle filter, PERCEP preprocessing, etc.
- **Formulas**: Calif calculation formulas (stage 1 + stage 2 additions), viejitos rules, composite variable definitions
- **Variables**: All variable weights, IF-conditional scores, SYSMIS default values (~1,800 rows)
- **Discriminant**: DISCRIMINANT analysis variable lists and PRIORs
- All sheets have an `ID` column (B001–B104) matching the block filenames

## PDF Structure (Reference)
- Each region spans 5–15 pages. Region 1 (Los Altos) is pp. 1–13, Region 2 (Sierra) is pp. 14–29, etc.
- Each region contains multiple models (MODELO yy, aa, bb, zz)
- Header spelling varies: "REGION 1. LOS ALTOS", "PUNTAJE PARA LA REGION 3", etc.
- SPSS commands used: COMPUTE, IF, RECODE, DISCRIMINANT, FREQUENCIES, MEANS, FILTER, etc.

## File Structure
- `CLAUDE.md` — This file (English)
- `CLAUDE_JP.md` — Japanese version of this file
- `README.md` — Detailed English documentation (for co-authors)
- `scripts/ocr_progresa.py` — Step 1: Claude API vision OCR
- `scripts/parse_to_excel.py` — Step 2: Block boundary detection
- `scripts/build_catalog.py` — Step 3: Block splitting + catalog Excel creation
- `source_pdf/` — Source PDF files (3 files)
- `output/full_ocr_text.txt` — OCR text for all 284 pages
- `progresa_models.xlsx` — Step 2 output
- `progresa_scoring_catalog.xlsx` — Step 3 output (4-sheet catalog)
- `blocks/` — Block text files (104 files)
