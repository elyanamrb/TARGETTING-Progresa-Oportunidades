#!/usr/bin/env python3
"""
Step 1: OCR scanned PDF using Claude API vision
================================================
Sends pages as images to Claude API and transcribes SPSS syntax.

Each page is sent as an independent image within a batch. The batch size
only affects API call efficiency — it does NOT affect OCR accuracy because:
  - Each page is a separate image in the request
  - The final output (full_ocr_text.txt) concatenates all pages in order
  - Step 2 parsing operates on the concatenated text, not individual batches

Usage:
    python ocr_progresa.py --pdf 41_REGIONES_Comandos.pdf --output-dir ./output
    python ocr_progresa.py --pdf 41_REGIONES_Comandos.pdf --output-dir ./output --batch-size 2  # slower but safer
    python ocr_progresa.py --pdf 41_REGIONES_Comandos.pdf --output-dir ./output --start-page 101  # resume

Requirements:
    pip install anthropic pdf2image Pillow
    System: poppler-utils (brew install poppler / apt install poppler-utils)
"""

import anthropic
import base64
import json
import os
import sys
import argparse
import time
from pathlib import Path
from io import BytesIO

try:
    from pdf2image import convert_from_path
    from pdf2image.pdf2image import pdfinfo_from_path
except ImportError:
    print("ERROR: pdf2image not installed. Run: pip install pdf2image")
    sys.exit(1)

try:
    from PIL import Image
except ImportError:
    print("ERROR: Pillow not installed. Run: pip install Pillow")
    sys.exit(1)


SYSTEM_PROMPT = """\
You are an expert OCR system specializing in SPSS syntax transcription.
You will receive scanned images of SPSS command files from Mexico's \
Progresa/Oportunidades poverty targeting program.

Rules:
1. Transcribe ALL text exactly as written — do not interpret, summarize, or correct.
2. Preserve SPSS command structure: COMPUTE, IF, RECODE, DISCRIMINANT, DO IF, END IF, EXECUTE, FREQUENCIES, etc.
3. Preserve variable names exactly (DEPEN, NINOS, HACINA, OCUPJ_1, etc.)
4. Preserve ALL numeric values (coefficients, thresholds, SYSMIS defaults) with exact decimal places.
5. Preserve region headers (REGIÓN N. NAME / FASE N / MODELO yy/aa).
6. Each image is one page. Output them in order, separated by a blank line and a page marker like: --- Page N/284 ---
7. If a character is genuinely unreadable, mark it as [?]. Do NOT guess.
8. Pages may start or end mid-command — that is fine, just transcribe what is visible.
"""


def pdf_page_to_base64(pdf_path: str, page_num: int, dpi: int = 150) -> str:
    """Convert a single PDF page to a base64-encoded PNG."""
    images = convert_from_path(pdf_path, first_page=page_num, last_page=page_num, dpi=dpi)
    buf = BytesIO()
    images[0].save(buf, format="PNG")
    return base64.standard_b64encode(buf.getvalue()).decode("utf-8")


def ocr_batch(
    client: anthropic.Anthropic,
    pdf_path: str,
    page_start: int,
    page_end: int,
    dpi: int = 150,
    model: str = "claude-sonnet-4-20250514",
) -> str:
    """Send a batch of page images to Claude API for transcription."""
    content = []

    for page_num in range(page_start, page_end + 1):
        b64 = pdf_page_to_base64(pdf_path, page_num, dpi=dpi)
        content.append({
            "type": "text",
            "text": f"--- Image for page {page_num}/284 ---",
        })
        content.append({
            "type": "image",
            "source": {
                "type": "base64",
                "media_type": "image/png",
                "data": b64,
            },
        })

    content.append({
        "type": "text",
        "text": "Transcribe all pages above following the rules in the system prompt.",
    })

    response = client.messages.create(
        model=model,
        max_tokens=16000,
        system=SYSTEM_PROMPT,
        messages=[{"role": "user", "content": content}],
    )

    return response.content[0].text


def get_total_pages(pdf_path: str) -> int:
    info = pdfinfo_from_path(pdf_path)
    return info["Pages"]


def rebuild_full_text(output_dir: Path) -> None:
    """Rebuild full_ocr_text.txt from individual batch files in page order."""
    batch_files = sorted(output_dir.glob("pages_*.txt"))
    full_path = output_dir / "full_ocr_text.txt"
    with open(full_path, "w", encoding="utf-8") as out:
        for bf in batch_files:
            with open(bf, "r", encoding="utf-8") as f:
                out.write(f.read())
    print(f"  Rebuilt {full_path} from {len(batch_files)} batch files")


def main():
    parser = argparse.ArgumentParser(
        description="OCR Progresa SPSS syntax PDF using Claude API vision"
    )
    parser.add_argument("--pdf", required=True, help="Path to scanned PDF")
    parser.add_argument("--output-dir", default="./output", help="Output directory (default: ./output)")
    parser.add_argument("--start-page", type=int, default=None, help="Resume from this page")
    parser.add_argument("--end-page", type=int, default=None, help="Stop at this page")
    parser.add_argument(
        "--batch-size", type=int, default=4,
        help="Pages per API call (default: 4). Lower = safer, higher = faster. "
             "Does not affect accuracy — only API efficiency."
    )
    parser.add_argument("--dpi", type=int, default=150, help="Image DPI (default: 150)")
    parser.add_argument("--model", default="claude-sonnet-4-20250514", help="Claude model")
    parser.add_argument("--delay", type=float, default=1.0, help="Seconds between API calls")
    parser.add_argument("--rebuild-only", action="store_true",
                        help="Only rebuild full_ocr_text.txt from existing batch files")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.rebuild_only:
        rebuild_full_text(output_dir)
        return

    if not os.path.exists(args.pdf):
        print(f"ERROR: PDF not found: {args.pdf}")
        sys.exit(1)

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        print("ERROR: Set ANTHROPIC_API_KEY environment variable")
        sys.exit(1)

    client = anthropic.Anthropic(api_key=api_key)

    total_pages = get_total_pages(args.pdf)
    end_page = args.end_page or total_pages

    # Load progress
    progress_file = output_dir / "progress.json"
    if progress_file.exists():
        with open(progress_file) as f:
            progress = json.load(f)
    else:
        progress = {"completed_batches": [], "last_page_completed": 0}

    # Determine start page
    if args.start_page is not None:
        start_page = args.start_page
    elif progress["last_page_completed"] > 0:
        start_page = progress["last_page_completed"] + 1
        print(f"Resuming from page {start_page} (last completed: {progress['last_page_completed']})")
    else:
        start_page = 1

    n_batches = -(-((end_page - start_page + 1)) // args.batch_size)  # ceiling division
    print(f"PDF: {args.pdf} ({total_pages} pages)")
    print(f"Pages {start_page}–{end_page}, batch size {args.batch_size} → {n_batches} API calls")
    print(f"Output: {output_dir}")
    print(f"Model: {args.model}")
    print()

    batch_num = 0
    for batch_start in range(start_page, end_page + 1, args.batch_size):
        batch_end = min(batch_start + args.batch_size - 1, end_page)
        batch_key = f"{batch_start:03d}_{batch_end:03d}"
        batch_file = output_dir / f"pages_{batch_key}.txt"

        if batch_key in progress["completed_batches"]:
            print(f"  Skip pages {batch_start}–{batch_end} (done)")
            continue

        batch_num += 1
        print(f"  [{batch_num}/{n_batches}] Pages {batch_start}–{batch_end} ...", end=" ", flush=True)

        retries = 0
        while retries < 3:
            try:
                t0 = time.time()
                text = ocr_batch(client, args.pdf, batch_start, batch_end,
                                 dpi=args.dpi, model=args.model)
                elapsed = time.time() - t0

                with open(batch_file, "w", encoding="utf-8") as f:
                    f.write(text)
                    f.write("\n\n")

                progress["completed_batches"].append(batch_key)
                progress["last_page_completed"] = batch_end
                with open(progress_file, "w") as f:
                    json.dump(progress, f, indent=2)

                print(f"OK ({elapsed:.1f}s, {len(text)} chars)")
                time.sleep(args.delay)
                break

            except anthropic.RateLimitError:
                retries += 1
                wait = 30 * retries
                print(f"RATE LIMITED, waiting {wait}s (retry {retries}/3)")
                time.sleep(wait)
            except Exception as e:
                retries += 1
                if retries >= 3:
                    print(f"FAILED after 3 retries: {e}")
                    print(f"Resume with: --start-page {batch_start}")
                    # Still rebuild what we have so far
                    rebuild_full_text(output_dir)
                    sys.exit(1)
                print(f"Error: {e}, retrying ({retries}/3)...")
                time.sleep(5)

    # Rebuild the combined file from all batch files
    rebuild_full_text(output_dir)
    print(f"\nDone! {len(progress['completed_batches'])} batches completed.")
    print(f"Next step: python parse_to_excel.py --input {output_dir}/full_ocr_text.txt")


if __name__ == "__main__":
    main()
