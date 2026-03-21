# Progresa SPSS Syntax OCR Project

## 目的
`41_REGIONES_Comandos.pdf`（284ページ、スキャンPDF）をOCRで読み取り、構造化データにする。
このPDFにはメキシコのProgresa/Oportunidadesプログラムの各地域の貧困スコアリングモデルが
SPSSシンタックスとして記述されている。

## パイプライン概要

```
[Step 1] scripts/ocr_progresa.py    : PDF → output/full_ocr_text.txt (Claude API vision OCR)
[Step 2] scripts/parse_to_excel.py  : output/full_ocr_text.txt → progresa_models.xlsx (ブロック境界特定)
[Step 3] scripts/build_catalog.py   : output/full_ocr_text.txt + progresa_models.xlsx → blocks/ + catalog.xlsx
```

## 何をすべきか

### full_ocr_text.txt と progresa_models.xlsx が既にある場合
Step 1, 2 は完了済み。Step 3 だけ実行すればよい：

```bash
pip install pandas openpyxl
python scripts/build_catalog.py \
    --ocr-text output/full_ocr_text.txt \
    --initial-excel progresa_models.xlsx \
    --blocks-dir ./blocks \
    --output progresa_scoring_catalog.xlsx
```

### PDFから全部やり直す場合
```bash
pip install anthropic pdf2image Pillow openpyxl pandas
# macOS: brew install poppler / Linux: apt install poppler-utils

# Step 1: OCR (ANTHROPIC_API_KEY が必要)
python scripts/ocr_progresa.py --pdf source_pdf/41_REGIONES_Comandos.pdf --output-dir ./output

# Step 2: ブロック境界特定
python scripts/parse_to_excel.py --input ./output/full_ocr_text.txt --output progresa_models.xlsx

# Step 3: ブロック分割 + カタログ作成
python scripts/build_catalog.py \
    --ocr-text ./output/full_ocr_text.txt \
    --initial-excel progresa_models.xlsx \
    --blocks-dir ./blocks \
    --output progresa_scoring_catalog.xlsx
```

## 各ステップの詳細

### Step 1: ocr_progresa.py
- 284ページを `--batch-size`（デフォルト4）ページずつClaude API visionに送ってOCR
- `output/progress.json` で進捗管理。途中で止まっても自動再開
- batch-sizeはAPI効率のみに影響。OCR精度には無関係

### Step 2: parse_to_excel.py
- OCRテキストを正規表現で解析し、地域×モデルごとのブロック境界（page_start/page_end）を特定
- ヘッダーの表記揺れが多いため、パース結果が少ない場合は正規表現の調整が必要

### Step 3: build_catalog.py
- Step 2 の page_start/page_end を使って OCR テキストを 104 ブロックに分割
- 各ブロックから以下を正規表現で抽出し、4シートの Excel にまとめる

## 出力物

### blocks/ ディレクトリ（104個のテキストファイル）
- 命名規則: `B001_R01_F2_yy_p001-004.txt`
  - `B001` = 連番（ExcelのID列と対応）
  - `R01` = 地域番号
  - `F2` = フェーズ（FX = 不明）
  - `yy` = モデル名（XX = 不明）
  - `p001-004` = PDFのページ範囲
- 各ファイルには該当モデルの完全なSPSSシンタックスが入っている

### progresa_scoring_catalog.xlsx（4シート）
- **Catalog**: 分類一覧。変数構成（連続重み/バイナリIF/多カテゴリIF/合成変数の数）、
  閾値、viejitos調整有無、2段階calif修正、車フィルター、PERCEP前処理等
- **Formulas**: calif計算式（第1段階＋第2段階追加分）、viejitosルール、合成変数定義
- **Variables**: 全変数の重み係数、IF条件スコア、SYSMISデフォルト値（約1,800行）
- **Discriminant**: DISCRIMINANT分析の変数リストとPRIORs
- 全シートに `ID` 列（B001〜B104）があり、blocks/ のファイル名と対応

## PDFの構造（参考）
- 各地域は5-15ページ。地域1(Los Altos)はp.1-13, 地域2(Sierra)はp.14-29くらい
- 各地域内に複数モデル（MODELO yy, aa, bb, zz）がある
- ヘッダーの表記揺れあり: 「REGIÓN 1. LOS ALTOS」「PUNTAJE PARA LA REGION 3」等
- SPSSコマンド: COMPUTE, IF, RECODE, DISCRIMINANT, FREQUENCIES, MEANS, FILTER等

## ファイル構成
- `CLAUDE.md` — このファイル
- `README.md` — 英語ドキュメント（共著者向け）
- `scripts/ocr_progresa.py` — Step 1: Claude API vision OCR
- `scripts/parse_to_excel.py` — Step 2: ブロック境界特定
- `scripts/build_catalog.py` — Step 3: ブロック分割 + カタログExcel作成
- `source_pdf/` — 元PDFファイル（3つ）
- `output/full_ocr_text.txt` — 全284ページのOCRテキスト
- `progresa_models.xlsx` — Step 2出力
- `progresa_scoring_catalog.xlsx` — Step 3出力（4シートカタログ）
- `blocks/` — ブロック分割テキストファイル（104個）
