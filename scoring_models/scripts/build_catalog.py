#!/usr/bin/env python3
"""
Step 3: Split OCR text into blocks and create structured catalog Excel
======================================================================
Takes the OCR output (full_ocr_text.txt) and the initial parsed Excel
(progresa_models.xlsx) to produce:
  1. Individual block text files (one per region×fase×modelo)
  2. A structured catalog Excel with 4 sheets

Usage:
    python build_catalog.py \
        --ocr-text ./output/full_ocr_text.txt \
        --initial-excel ./progresa_models.xlsx \
        --blocks-dir ./blocks \
        --output ./progresa_scoring_catalog.xlsx

Requirements:
    pip install pandas openpyxl
"""

import os
import re
import argparse
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter


# ============================================================
# Part 1: Split OCR text into block files
# ============================================================

def split_into_blocks(ocr_text_path: str, initial_excel_path: str, blocks_dir: str) -> pd.DataFrame:
    """Split OCR text into per-block files using page boundaries from the initial Excel."""
    df = pd.read_excel(initial_excel_path, sheet_name="Summary")

    with open(ocr_text_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    total_lines = len(lines)

    # Build page number -> line number mapping
    page_to_line = {}
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("--- Page ") and "/284" in stripped:
            try:
                pnum = int(stripped.split("Page ")[1].split("/")[0])
                page_to_line[pnum] = i
            except ValueError:
                pass

    max_page = max(page_to_line.keys()) if page_to_line else 284
    os.makedirs(blocks_dir, exist_ok=True)

    filenames = []
    block_ids = []
    for idx, (_, row) in enumerate(df.iterrows(), 1):
        region = int(row["region"])
        fase = row["fase"]
        modelo = row["modelo"]
        ps = int(row["page_start"])
        pe = int(row["page_end"])

        start_line = page_to_line.get(ps, 0)
        end_line = page_to_line.get(pe + 1, total_lines) if pe < max_page else total_lines

        block_text = "".join(lines[start_line:end_line])

        fase_str = f"F{int(fase)}" if pd.notna(fase) else "FX"
        modelo_str = modelo if pd.notna(modelo) else "XX"
        block_id = f"B{idx:03d}"
        fname = f"{block_id}_R{region:02d}_{fase_str}_{modelo_str}_p{ps:03d}-{pe:03d}.txt"

        with open(os.path.join(blocks_dir, fname), "w", encoding="utf-8") as f:
            f.write(block_text)

        block_ids.append(block_id)
        filenames.append(fname)

    df["block_id"] = block_ids
    df["file"] = filenames
    print(f"  Split into {len(df)} block files in {blocks_dir}")
    return df


# ============================================================
# Part 2: Parse each block for structured info
# ============================================================

def parse_block(text: str) -> dict:
    """Extract all structured information from a single SPSS block."""

    # --- Continuous weights: COMPUTE AAvar = coeff * VAR ---
    # Pattern 1: COMPUTE AAVAR = coeff * VAR
    wt = re.findall(
        r"COMPUTE\s+(?:AA|YY|BB|XZ|ZZ)\w+\s*=\s*(-?\d+\.?\d*)\s*\*\s*(\w+)",
        text, re.I,
    )
    weight_vars = {v.upper(): float(c) for c, v in wt}
    # Pattern 2: IF (VAR >= 0) xxVAR = VAR * coeff (used in early phases)
    wt2 = re.findall(
        r"IF\s*\(\s*(\w+)\s*>=\s*0\s*\)\s*(?:aa|yy|bb|xz|zz)\w+\s*=\s*\1\s*\*\s*(-?\d+\.?\d*)",
        text, re.I,
    )
    for v, c in wt2:
        if v.upper() not in weight_vars:
            weight_vars[v.upper()] = float(c)

    # --- IF categorical scores: IF (VAR=val) AAvar=score ---
    ifs = re.findall(
        r"IF\s*\(\s*(\w+)\s*=\s*(\d+)\s*\)\s*(?:AA|YY|BB|XZ|ZZ)\w+\s*=\s*(-?\d+\.?\d*)",
        text, re.I,
    )
    if_vars: dict[str, dict[int, float]] = {}
    for var, cond, score in ifs:
        if_vars.setdefault(var.upper(), {})[int(cond)] = float(score)

    # --- Composite variables (tianle, indisp, patron, enseres) ---
    composites = []
    for m in re.finditer(
        r"COMPUTE\s+(\w+)\s*=\s*0\s*\..*?DO\s+IF\s*\((.+?)\)\s*\..*?RECODE\s+\1",
        text, re.I | re.DOTALL,
    ):
        composites.append({
            "name": m.group(1).lower(),
            "condition": m.group(2).strip().replace("\r", "")[:100],
        })

    # --- Calif formula (1st stage) ---
    calif_formula = ""
    m = re.search(
        r"COMPUTE\s+(?:aa|yy|bb|xz|zz)calif\s*=\s*(.+?)(?:\.\s*\r?\n|EXECUTE)",
        text, re.I | re.DOTALL,
    )
    if m:
        calif_formula = re.sub(r"\s+", " ", m.group(1).replace("\r", "").replace("\n", " ")).strip()

    # --- Second-stage calif update (FILTRO additions) ---
    calif2_formula = ""
    m = re.search(
        r"COMPUTE\s+(?:aa|yy|bb|xz|zz)calif\s*=\s*(?:aa|yy|bb|xz|zz)calif\s*\+\s*(.+?)(?:\.\s*\r?\n|EXECUTE)",
        text, re.I | re.DOTALL,
    )
    if m:
        calif2_formula = re.sub(r"\s+", " ", m.group(1).replace("\r", "").replace("\n", " ")).strip()

    # --- Viejitos rule ---
    viejitos_rule = ""
    m = re.search(r"IF\s*\(\s*(.+?)\)\s*(?:viejitos|VIEJITOS)\s*=\s*(-?\d+)", text, re.I)
    if m:
        viejitos_rule = m.group(1).strip().replace("\r", "")[:120]

    # --- AACALIF1 extra adjustments beyond viejitos ---
    calif1_extra = ""
    m = re.search(
        r"(?:aa|yy|bb)calif1\s*=\s*(?:aa|yy|bb)calif1\s*(.+?)(?:\.\s*\r?\n|EXECUTE)",
        text, re.I | re.DOTALL,
    )
    if m:
        calif1_extra = re.sub(r"\s+", " ", m.group(1).replace("\r", "").replace("\n", " ")).strip()[:120]

    # --- Threshold ---
    threshold = None
    thresh_var = ""
    # Pattern 1: DO IF (xxcalif < N) or (xxcalif1 < N)
    m = re.search(r"((?:aa|yy|bb|xz|zz)calif\w*)\s*<\s*=?\s*(-?\d+\.?\d*)", text, re.I)
    if m:
        thresh_var = m.group(1).lower()
        threshold = float(m.group(2))
    # Pattern 2: if (xxcalif<=N) pobre=1  (older style with <=)
    if threshold is None:
        m = re.search(r"((?:aa|yy|bb|xz|zz)calif\w*)\s*<=\s*(-?\d+\.?\d*)", text, re.I)
        if m:
            thresh_var = m.group(1).lower()
            threshold = float(m.group(2))

    # --- Post-threshold vehicle filter ---
    vehi_filter = ""
    if re.search(r"vehi\s*=\s*0.*(?:pobre|POBRE)", text, re.I):
        vehi_filter = "vehi=0 → remove pobre"

    # --- PERCEP preprocessing ---
    has_percep = bool(re.search(r"PERCEP\s*=\s*0.*DEPEN\s*=\s*99", text, re.I))

    # --- DISCRIMINANT vars ---
    disc_vars = ""
    m = re.search(r"/VARIABLES\s*=\s*([\w\s_]+?)(?:\r?\n\s*/|\r?\n\s*\r?\n)", text, re.I)
    if m:
        disc_vars = " ".join(m.group(1).split()).upper()

    # --- PRIORS ---
    priors = ""
    m = re.search(r"/PRIORS\s*=?\s*([\d\.\s]+|EQUAL)", text, re.I)
    if m:
        priors = m.group(1).strip()

    # --- SYSMIS defaults ---
    sysmis = {
        v.upper(): float(d)
        for v, d in re.findall(
            r"RECODE\s+(?:aa|yy|bb|xz|zz)(\w+)\s*\(\s*SYSMIS\s*=\s*(-?\d+\.?\d*)\s*\)",
            text, re.I,
        )
    }

    # --- Variable prefix ---
    has_yy = bool(re.search(r"(?:COMPUTE|IF\s*\()\s*.*\byy", text, re.I))
    has_aa = bool(re.search(r"(?:COMPUTE|IF\s*\()\s*.*\baa", text, re.I))
    has_bb = bool(re.search(r"(?:COMPUTE|IF\s*\()\s*.*\bbb", text, re.I))
    has_xz = bool(re.search(r"(?:COMPUTE|IF\s*\()\s*.*\bxz", text, re.I))
    has_zz = bool(re.search(r"(?:COMPUTE|IF\s*\()\s*.*\bzz", text, re.I))
    prefix = "yy" if has_yy else ("aa" if has_aa else ("bb" if has_bb else ("xz" if has_xz else ("zz" if has_zz else "?"))))

    n_binary_if = sum(1 for mappings in if_vars.values() if len(mappings) == 2)
    n_multi_if = sum(1 for mappings in if_vars.values() if len(mappings) > 2)

    return {
        "prefix": prefix,
        "n_continuous": len(weight_vars),
        "n_binary_if": n_binary_if,
        "n_multi_if": n_multi_if,
        "n_composites": len(composites),
        "composites": "; ".join(f"{c['name']}({c['condition'][:50]})" for c in composites),
        "calif_formula": calif_formula[:200],
        "calif2_additions": calif2_formula[:150],
        "has_viejitos": bool(viejitos_rule),
        "viejitos_rule": viejitos_rule,
        "calif1_extra_adj": calif1_extra,
        "threshold": threshold,
        "threshold_on": thresh_var,
        "vehi_filter": vehi_filter,
        "has_percep_preproc": has_percep,
        "priors": priors,
        "disc_vars": disc_vars[:200],
        "weight_vars": weight_vars,
        "if_vars": if_vars,
        "sysmis": sysmis,
    }


def parse_all_blocks(df: pd.DataFrame, blocks_dir: str) -> pd.DataFrame:
    """Parse all block files and add columns to the dataframe."""
    parsed_rows = []
    for _, row in df.iterrows():
        fpath = os.path.join(blocks_dir, row["file"])
        with open(fpath, "r", encoding="utf-8") as f:
            text = f.read()
        info = parse_block(text)
        info["chars"] = len(text)
        parsed_rows.append(info)

    pdf = pd.DataFrame(parsed_rows)
    # Merge back
    for col in pdf.columns:
        df[col] = pdf[col].values
    return df


# ============================================================
# Part 3: Write Excel catalog
# ============================================================

def write_excel(df: pd.DataFrame, output_path: str):
    """Write the 4-sheet structured Excel."""
    wb = Workbook()
    hdr_font = Font(bold=True, size=10, color="FFFFFF")
    hdr_fill = PatternFill("solid", fgColor="2F5496")
    wrap = Alignment(wrap_text=True, vertical="top")
    border = Border(bottom=Side(style="thin", color="D0D0D0"))

    def write_sheet(ws, col_defs, data_df):
        for ci, (key, width, label) in enumerate(col_defs, 1):
            c = ws.cell(row=1, column=ci, value=label)
            c.font = hdr_font
            c.fill = hdr_fill
            c.alignment = Alignment(wrap_text=True, vertical="bottom", horizontal="center")
            ws.column_dimensions[get_column_letter(ci)].width = width
        for ri, (_, row) in enumerate(data_df.iterrows(), 2):
            for ci, (key, _, _) in enumerate(col_defs, 1):
                val = row.get(key, "")
                if isinstance(val, bool):
                    val = "Y" if val else ""
                if pd.isna(val) if not isinstance(val, (str, bool)) else False:
                    val = ""
                c = ws.cell(row=ri, column=ci, value=val)
                c.alignment = wrap
                c.border = border
        ws.auto_filter.ref = f"A1:{get_column_letter(len(col_defs))}{len(data_df)+1}"
        ws.freeze_panes = "A2"

    # Sheet 1: Catalog
    ws1 = wb.active
    ws1.title = "Catalog"
    write_sheet(ws1, [
        ("block_id", 7, "ID"),
        ("file", 28, "Block File"),
        ("region", 6, "Region"),
        ("fase", 6, "Fase"),
        ("modelo", 7, "Modelo"),
        ("prefix", 6, "Prefix"),
        ("pages", 12, "Pages"),
        ("n_continuous", 10, "Cont.\nWeights"),
        ("n_binary_if", 10, "Binary\nIF Vars"),
        ("n_multi_if", 10, "Multi-cat\nIF Vars"),
        ("n_composites", 10, "Composite\nVars"),
        ("threshold", 10, "Threshold"),
        ("threshold_on", 12, "Thresh.\nVariable"),
        ("has_viejitos", 9, "Viejitos\nAdj."),
        ("calif1_extra_adj", 20, "CALIF1\nExtra Adj."),
        ("calif2_additions", 25, "CALIF\n2nd Stage"),
        ("vehi_filter", 18, "Vehicle\nFilter"),
        ("has_percep_preproc", 9, "PERCEP\nPreproc"),
        ("priors", 10, "Priors"),
    ], df)

    # Sheet 2: Formulas
    ws2 = wb.create_sheet("Formulas")
    write_sheet(ws2, [
        ("block_id", 7, "ID"),
        ("file", 28, "Block File"),
        ("region", 6, "Region"),
        ("fase", 6, "Fase"),
        ("modelo", 7, "Modelo"),
        ("calif_formula", 80, "CALIF Formula (1st stage)"),
        ("calif2_additions", 40, "CALIF 2nd Stage Additions"),
        ("viejitos_rule", 50, "Viejitos Rule"),
        ("composites", 50, "Composite Variables"),
    ], df)

    # Sheet 3: Variables detail
    ws3 = wb.create_sheet("Variables")
    var_rows = []
    for _, row in df.iterrows():
        weights = row.get("weight_vars", {})
        if_v = row.get("if_vars", {})
        sysmis = row.get("sysmis", {})
        if not isinstance(weights, dict):
            weights = {}
        if not isinstance(if_v, dict):
            if_v = {}
        if not isinstance(sysmis, dict):
            sysmis = {}

        all_vars = sorted(set(list(weights.keys()) + list(if_v.keys())))
        for var in all_vars:
            vr = {
                "block_id": row["block_id"],
                "file": row["file"], "region": row["region"],
                "fase": row["fase"], "modelo": row["modelo"],
                "variable": var,
                "type": "continuous" if var in weights else "categorical",
                "weight": weights.get(var, ""),
                "sysmis_default": sysmis.get(var, ""),
            }
            if var in if_v:
                for cond, score in sorted(if_v[var].items()):
                    vr[f"if_{cond}"] = score
            var_rows.append(vr)

    df_vars = pd.DataFrame(var_rows)
    if_cols = sorted(
        [c for c in df_vars.columns if c.startswith("if_")],
        key=lambda x: int(x.split("_")[1]),
    )
    var_col_defs = [
        ("block_id", 7, "ID"),
        ("file", 28, "Block File"),
        ("region", 6, "Region"),
        ("fase", 6, "Fase"),
        ("modelo", 7, "Modelo"),
        ("variable", 14, "Variable"),
        ("type", 10, "Type"),
        ("weight", 10, "Weight"),
        ("sysmis_default", 12, "SYSMIS\nDefault"),
    ] + [(c, 10, f"Score\n({c})") for c in if_cols]
    write_sheet(ws3, var_col_defs, df_vars)

    # Sheet 4: Discriminant
    ws4 = wb.create_sheet("Discriminant")
    write_sheet(ws4, [
        ("block_id", 7, "ID"),
        ("file", 28, "Block File"),
        ("region", 6, "Region"),
        ("fase", 6, "Fase"),
        ("modelo", 7, "Modelo"),
        ("priors", 10, "Priors"),
        ("disc_vars", 120, "DISCRIMINANT Variables"),
    ], df)

    wb.save(output_path)
    print(f"  Saved Excel: {output_path}")
    print(f"    Catalog:      {len(df)} blocks")
    print(f"    Variables:    {len(df_vars)} rows")
    print(f"    Discriminant: {len(df)} blocks")


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Build structured catalog from OCR output")
    parser.add_argument("--ocr-text", required=True, help="Path to full_ocr_text.txt")
    parser.add_argument("--initial-excel", required=True, help="Path to progresa_models.xlsx (from parse_to_excel.py)")
    parser.add_argument("--blocks-dir", default="./blocks", help="Output dir for block text files")
    parser.add_argument("--output", default="./progresa_scoring_catalog.xlsx", help="Output Excel path")
    args = parser.parse_args()

    print("Step 1: Splitting OCR text into blocks...")
    df = split_into_blocks(args.ocr_text, args.initial_excel, args.blocks_dir)

    # Add page range string
    df["pages"] = df.apply(lambda r: f"p{int(r['page_start']):03d}-{int(r['page_end']):03d}", axis=1)

    print("Step 2: Parsing block contents...")
    df = parse_all_blocks(df, args.blocks_dir)

    print("Step 3: Writing catalog Excel...")
    write_excel(df, args.output)

    print(f"\nDone!")
    print(f"  Block files: {args.blocks_dir}/ ({len(df)} files)")
    print(f"  Excel:       {args.output}")


if __name__ == "__main__":
    main()
