#!/usr/bin/env python3
"""
Step 2: Parse OCR text → structured Excel
==========================================
Reads the concatenated OCR output and extracts model parameters
(variables, weights, thresholds, IF-scores, SYSMIS defaults)
per region/fase/modelo.

Usage:
    python parse_to_excel.py --input ./output/full_ocr_text.txt --output progresa_models.xlsx

Requirements:
    pip install openpyxl pandas
"""

import re
import argparse
import pandas as pd


def build_page_index(text: str) -> list[tuple[int, int]]:
    """Return list of (char_pos, page_num) for every '--- Page N/NNN ---' marker."""
    page_marker = re.compile(r'--- Page (\d+)/\d+ ---')
    return [(m.start(), int(m.group(1))) for m in page_marker.finditer(text)]


def pos_to_page(char_pos: int, page_index: list[tuple[int, int]]) -> int:
    """Return the PDF page number containing the given character position."""
    page = 1
    for marker_pos, page_num in page_index:
        if marker_pos <= char_pos:
            page = page_num
        else:
            break
    return page


def parse_region_blocks(text: str) -> list[dict]:
    """
    Split the full OCR text into blocks by region/model header.
    Handles header variations:
      - REGIÓN 1. LOS ALTOS / FASE 2 / MODELO yy
      - PUNTAJE PARA LA REGION 3
      - REGION 3 / FASE 8 / MODELO AA
    """
    page_index = build_page_index(text)

    # Find all headers with their positions
    header_pattern = re.compile(
        r'(?:REGI[ÓO]N|PUNTAJE\s+(?:PARA|DE)\s+LA\s+REGI[ÓO]N)\s*[:\s]*(\d+)'
        r'[\.\s]*([A-ZÁÉÍÓÚÑ\s]*)',
        re.IGNORECASE
    )
    fase_pattern = re.compile(r'FASE\s*(\d+)', re.IGNORECASE)
    modelo_pattern = re.compile(r'MODELO\s*(yy|aa|bb|YY|AA|BB)', re.IGNORECASE)

    # Find all header positions
    headers = []
    for match in header_pattern.finditer(text):
        pos = match.start()
        region_num = int(match.group(1))
        region_name = (match.group(2) or '').strip()

        # Look for FASE and MODELO in the next 200 chars
        context = text[pos:pos + 300]
        fase_match = fase_pattern.search(context)
        modelo_match = modelo_pattern.search(context)

        headers.append({
            'pos': pos,
            'region': region_num,
            'region_name': region_name,
            'fase': fase_match.group(1) if fase_match else '',
            'modelo': modelo_match.group(1).lower() if modelo_match else '',
        })

    # Split text into blocks between headers
    blocks = []
    for i, h in enumerate(headers):
        start = h['pos']
        end = headers[i + 1]['pos'] if i + 1 < len(headers) else len(text)
        block_text = text[start:end]

        # Only keep blocks that have actual SPSS content
        if 'COMPUTE' in block_text.upper() or 'DISCRIMINANT' in block_text.upper():
            page_start = pos_to_page(start, page_index)
            page_end = pos_to_page(end - 1, page_index)
            blocks.append({**h, 'text': block_text,
                           'page_start': page_start, 'page_end': page_end})

    return blocks


def extract_weights(text: str) -> dict:
    """Extract COMPUTE XX_VAR = coeff * VAR patterns."""
    weights = {}
    pattern = re.compile(
        r'COMPUTE\s+(?:AA|YY|BB)(\w+)\s*=\s*(-?\d+\.?\d*)\s*\*\s*(\w+)',
        re.IGNORECASE
    )
    for m in pattern.finditer(text):
        var = m.group(3).upper()
        w = float(m.group(2))
        weights[var] = w
    return weights


def extract_if_scores(text: str) -> dict:
    """Extract IF (VAR=val) XXVAR=score patterns."""
    scores = {}
    pattern = re.compile(
        r'IF\s*\(\s*(\w+)\s*=\s*(\d+)\s*\)\s*(?:AA|YY|BB)\w+\s*=\s*(-?\d+\.?\d*)',
        re.IGNORECASE
    )
    for m in pattern.finditer(text):
        var = m.group(1).upper()
        cond = int(m.group(2))
        score = float(m.group(3))
        scores.setdefault(var, {})[cond] = score
    return scores


def extract_threshold(text: str) -> float | None:
    """Extract DO IF (xxcalif < threshold)."""
    m = re.search(
        r'(?:aa|yy|bb)calif\w*\s*<\s*(-?\d+\.?\d*)',
        text, re.IGNORECASE
    )
    return float(m.group(1)) if m else None


def extract_sysmis(text: str) -> dict:
    """Extract RECODE xxVAR (SYSMIS=default)."""
    defaults = {}
    pattern = re.compile(
        r'RECODE\s+(?:aa|yy|bb)(\w+)\s*\(\s*SYSMIS\s*=\s*(-?\d+\.?\d*)\s*\)',
        re.IGNORECASE
    )
    for m in pattern.finditer(text):
        defaults[m.group(1).upper()] = float(m.group(2))
    return defaults


def extract_discriminant_vars(text: str) -> list[str]:
    """Extract /VARIABLES= ... from DISCRIMINANT."""
    m = re.search(
        r'/VARIABLES\s*=\s*([\w\s_]+?)(?:\n\s*/|\n\s*\n)',
        text, re.IGNORECASE | re.MULTILINE
    )
    if m:
        return [v.strip().upper() for v in m.group(1).split() if v.strip()]
    return []


def extract_priors(text: str) -> str:
    m = re.search(r'/PRIORS\s*=?\s*([\d\.\s]+|EQUAL)', text, re.IGNORECASE)
    return m.group(1).strip() if m else ''


def extract_calif_formula(text: str) -> str:
    m = re.search(
        r'COMPUTE\s+(?:aa|yy|bb)calif\w*\s*=\s*(.+?)(?:\.\s*\n|EXECUTE)',
        text, re.IGNORECASE | re.DOTALL
    )
    if m:
        return re.sub(r'\s+', ' ', m.group(1).replace('\n', ' ')).strip()
    return ''


def build_dataframes(blocks: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    detail_rows = []

    for b in blocks:
        text = b['text']
        weights = extract_weights(text)
        if_scores = extract_if_scores(text)
        threshold = extract_threshold(text)
        sysmis = extract_sysmis(text)
        disc_vars = extract_discriminant_vars(text)
        priors = extract_priors(text)
        formula = extract_calif_formula(text)

        base = {
            'region': b['region'],
            'region_name': b['region_name'],
            'fase': b['fase'],
            'modelo': b['modelo'],
            'page_start': b.get('page_start', ''),
            'page_end': b.get('page_end', ''),
        }

        summary_rows.append({
            **base,
            'n_weight_vars': len(weights),
            'n_if_vars': len(if_scores),
            'threshold': threshold,
            'priors': priors,
            'discriminant_vars': ', '.join(disc_vars),
            'calif_formula': formula,
        })

        all_vars = sorted(set(list(weights.keys()) + list(if_scores.keys())))
        for var in all_vars:
            row = {**base, 'variable': var, 'weight': weights.get(var, ''),
                   'sysmis_default': sysmis.get(var, '')}
            if var in if_scores:
                for cond, score in sorted(if_scores[var].items()):
                    row[f'score_when_{cond}'] = score
            detail_rows.append(row)

    return pd.DataFrame(summary_rows), pd.DataFrame(detail_rows)


def main():
    parser = argparse.ArgumentParser(description="Parse OCR text → Excel")
    parser.add_argument("--input", required=True, help="Path to full_ocr_text.txt")
    parser.add_argument("--output", default="progresa_models.xlsx", help="Output Excel")
    args = parser.parse_args()

    print(f"Reading: {args.input}")
    with open(args.input, encoding='utf-8') as f:
        text = f.read()
    print(f"  {len(text):,} characters")

    print("Parsing region blocks...")
    blocks = parse_region_blocks(text)
    print(f"  {len(blocks)} blocks found:")
    for b in blocks:
        print(f"    R{b['region']:02d} F{b['fase']:>2s} M{b['modelo']:>2s} "
              f"{b['region_name'][:25]:<25s} ({len(b['text']):,} chars)")

    print("\nExtracting parameters...")
    df_summary, df_detail = build_dataframes(blocks)

    print(f"\nWriting: {args.output}")
    with pd.ExcelWriter(args.output, engine='openpyxl') as w:
        df_summary.to_excel(w, sheet_name='Summary', index=False)
        df_detail.to_excel(w, sheet_name='Variables_Weights', index=False)

    print(f"  Summary: {len(df_summary)} rows")
    print(f"  Variables: {len(df_detail)} rows")
    print("Done!")


if __name__ == "__main__":
    main()
