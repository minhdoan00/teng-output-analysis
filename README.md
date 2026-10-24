# TENG Output Analysis (`teng-output-analysis`)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A comprehensive Python suite and agent skill for analyzing and visualizing **Triboelectric Nanogenerator (TENG)** electrical performance—including:
- Open-circuit voltage ($V_{\text{OC}}$)
- Short-circuit current ($I_{\text{SC}}$)
- Transferred charge integration ($I_{\text{SC}} \to Q_{\text{SC}}, Q_{\text{SC,max}}$)
- Triboelectric surface charge density ($\sigma$)

Includes publication-grade **Origin-style** waveform plotting ($300\text{ DPI}$) with embedded peak-to-peak statistics and **4-panel charge-integration diagnostic figures**.

---

## Key Features

1. **Probabilistic Signal Recognition (`Voc` vs `Isc`)**:
   - Automatic classification using a 3-tier hierarchy based on oscilloscope vertical scale (`#voltbase`), paired channel context, and physical amplitude, assigning confidence levels (`EXPLICIT`, `HIGH`, `MODERATE`, `AMBIGUOUS`).
   - Human-in-the-loop guidelines to prevent dogmatic misclassification.
2. **Resting-State DC Baseline Correction & Windowing**:
   - Automatic time-window extraction (default: $[1.0\text{ s}, 7.0\text{ s}]$) discarding startup edge transients.
   - 5-point rolling-variance resting-state DC baseline estimation preventing peak asymmetry bias.
3. **Zero-Gap Midpoint Windowing (`Isc -> Qsc`)**:
   - Detects dominant polarity and contact-separation sequence (`[+ -> -]` vs `[- -> +]`).
   - Resolves closely spaced bipolar pulses ($50\text{--}60\text{ ms}$) without truncation.
4. **Origin-Style Scientific Figures**:
   - Bold sans-serif typography, inward ticks, 2.0 pt axis spines.
   - Embedded peak-to-peak amplitude annotation in the top-left corner ($V_{pp}$, $I_{pp}$, or $Q_{sc}$).

---

## Directory Structure

```text
teng-output-analysis/
├── SKILL.md                 # Agent skill specification (for Antigravity / Claude / Cursor)
├── README.md                # Project documentation
├── scripts/
│   └── teng_analyzer.py    # Standalone Python CLI & analysis engine
└── requirements.txt         # Python dependencies
```

---

## Installation

```bash
git clone https://github.com/minhdoan00/teng-output-analysis.git
cd teng-output-analysis
pip install -r requirements.txt
```

---

## Quick Usage

### 1. Analyze a Single CSV File
```bash
python scripts/teng_analyzer.py single --input path/to/waveform.csv --outdir ./output --lang vi
```

### 2. Batch-Process a Folder of Oscilloscope CSV Files
```bash
python scripts/teng_analyzer.py folder --input path/to/folder --t-start 2.0 --t-end 7.0 --lang vi
```

### 3. Batch-Process a Multi-Column Excel File
```bash
python scripts/teng_analyzer.py excel --input Data/Book1.xlsx --outdir Data/NIPS --lang vi
```

---

## Author & Citation

- **Author**: Minh Doan ([@minhdoan00](https://github.com/minhdoan00))
- **Institution**: VNU University of Engineering and Technology (VNU-UET)
- **License**: MIT

