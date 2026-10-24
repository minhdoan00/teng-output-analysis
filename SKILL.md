---
name: teng-output-analysis
description: >-
  Analyzes and visualizes Triboelectric Nanogenerator (TENG) electrical output
  performance—including open-circuit voltage (Voc), short-circuit current (Isc),
  short-circuit transferred charge integration (Isc -> Qsc, Qsc,max), and surface
  charge density (sigma)—from digital oscilloscope CSV files or multi-column Excel
  workbooks. Generates publication-ready Origin-style plots with embedded
  peak-to-peak annotations and 4-panel charge-integration diagnostic figures.
  Use whenever a user asks to process TENG electrical measurements, calculate
  Voc/Isc/Qsc/sigma, integrate current pulses to charge, or plot TENG waveforms.
---

# TENG Electrical Output Analysis & Visualization (`Voc`, `Isc`, `Isc -> Qsc`)

## Overview

This skill provides a complete, reproducible scientific workflow for analyzing and visualizing the electrical output performance of **Triboelectric Nanogenerators (TENGs)** from raw digital oscilloscope recordings (`.csv`) and multi-column experimental spreadsheets (`.xlsx`). It unifies two core analytical pipelines:

1. **Origin-Style Waveform & Peak-to-Peak Analysis (`Voc`, `Isc`, `Qsc`)**:
   - Automated time-window slicing (default: discarding the first $1.0\text{ s}$ transient and extracting the subsequent $6.0\text{ s}$ steady-state window from $1.0\text{ s} \to 7.0\text{ s}$, shifting the plotted time axis to $0.0 \to 6.0\text{ s}$).
   - Resting-state DC baseline offset removal, dual-polarity peak detection, cycle pairing, and statistical calculation ($\text{Mean} \pm \text{SD}$ of positive peaks, negative peaks, and peak-to-peak amplitude $V_{\text{pp}}$ or $I_{\text{pp}}$).
   - Publication-grade **Origin-style** single-panel plots ($300\text{ DPI}$, bold sans-serif typography, inward ticks) with the **peak-to-peak value embedded in the top-left corner** and the **sample identifier/legend in the top-right corner**.
2. **Short-Circuit Current-to-Charge Integration (`Isc -> Qsc,max` & Surface Charge Density $\sigma$)**:
   - Automatic contact-separation phase-order recognition (`[+ -> -]` vs `[- -> +]`).
   - **Zero-Gap Midpoint Windowing** to partition contiguous cycle boundaries at the exact midpoint of inter-pulse resting intervals, eliminating blind spots for closely spaced bipolar pulses ($50\text{--}60\text{ ms}$ separation).
   - Bounded zero-crossing trapezoidal integration of positive ($Q_+$) and negative ($|Q_-|$) half-cycles, determining per-cycle maximum transferred charge $Q_{\text{SC,max}} = \max(Q_+, |Q_-|)$ and surface charge density $\sigma = Q_{\text{SC,max}} / A$.
   - Comprehensive **4-panel diagnostic figures** ($300\text{ DPI}$) for full-waveform inspection, peak verification, multi-cycle cumulative charge trajectories, and cycle-by-cycle statistical distribution.

---

## Language & Localization Rule (Mandatory)

> **CRITICAL INSTRUCTION**: Although this skill specification is written in English for international distribution, **you MUST always respond, present statistical tables, and generate explanations in the exact language used by the user when invoking the skill** (e.g., Vietnamese, English, Japanese, Korean, Chinese, French, etc.).
>
> - **User Prompt in Vietnamese** $\rightarrow$ Respond in Vietnamese; use Vietnamese plot axis labels (`LABELS["vi"]`) unless the user explicitly asks for English figure labels.
> - **User Prompt in English** $\rightarrow$ Respond in English; use English plot axis labels (`LABELS["en"]`).
> - **User Prompt in any other language** $\rightarrow$ Respond in that language and localize plot labels accordingly (or use standard English scientific labels if requested for international journal submission).

---

## Physical & Measurement Circuit Model

### 1. Transimpedance Current Scaling (`Isc`)
When short-circuit current $I_{\text{SC}}(t)$ is measured via a shunt resistor and preamplifier connected to an oscilloscope voltage channel (Column 2 of the CSV/Excel data):

$$R_{\text{trans}} = R_{\text{shunt}} \times G_{\text{amp}} \times \left(\frac{\text{OSC\_PROBE\_SETTING}}{\text{PROBE\_ATTEN}}\right)$$

$$I_{\text{SC}}(\mu\text{A}) = \frac{V_{\text{raw}} \times \text{UNIT\_V}}{R_{\text{trans}}} \times 10^6, \qquad I_{\text{SC}}(\text{nA}) = I_{\text{SC}}(\mu\text{A}) \times 10^3$$

**Default Hardware Parameters** (configurable per experiment):
| Parameter | Symbol | Default Value | Description |
| :--- | :---: | :---: | :--- |
| `CSV_VOLTAGE_UNIT` | $\text{UNIT\_V}$ | `"mV"` ($10^{-3}\text{ V}$) | Unit of raw Column 2 values (`"V"`, `"mV"`, `"uV"`) |
| `R_SHUNT_OHM` | $R_{\text{shunt}}$ | $100\times 10^3\ \Omega$ | Shunt resistance ($100\text{ k}\Omega$ for standard, $10\text{ k}\Omega$ for repeat/high-current) |
| `AMP_GAIN` | $G_{\text{amp}}$ | $10.0$ | Current preamplifier gain ($\times 10$) |
| `PROBE_ATTEN` | $P_{\text{phys}}$ | $10.0$ | Physical probe attenuation ratio ($\times 10$) |
| `OSC_PROBE_SETTING` | $P_{\text{osc}}$ | $10.0$ | Oscilloscope probe ratio setting ($\times 10$) |
| `VOC_MULTIPLIER` | $M_{\text{Voc}}$ | $1.0$ | Additional multiplier for $V_{\text{OC}}$ probe |
| `AREA_CM2` | $A$ | $25.0\text{ cm}^2$ | Active triboelectric contact area ($5\times 5\text{ cm}$) |

**Shunt Resistor Selection & Conversion Scaling**:
- **Standard ($R_{\text{shunt}} = 100\text{ k}\Omega$)**: $R_{\text{trans}} = 10^6\text{ V/A} = 1\text{ V}/\mu\text{A} = 1000\text{ mV}/\mu\text{A} \implies 1\text{ mV}_{\text{raw}} \leftrightarrow 1\text{ nA} = 0.001\ \mu\text{A}$.
- **High-Current / Repeat ($R_{\text{shunt}} = 10\text{ k}\Omega$)**: $R_{\text{trans}} = 10^5\text{ V/A} = 0.1\text{ V}/\mu\text{A} = 100\text{ mV}/\mu\text{A} \implies 1\text{ mV}_{\text{raw}} \leftrightarrow 10\text{ nA} = 0.01\ \mu\text{A}$ ($10\times$ current conversion).
- **Auto-Correction Rule for Repeat Measurements & Sub-$5\ \mu\text{A}$ Traces**: In multi-run experiments (such as electrospun composite TENG batches), repeat measurements often switch to $R_{\text{shunt}} = 10\text{ k}\Omega$. If an $I_{\text{SC}}$ file has `"repeat"` in its name or if the calculated $I_{\text{pp}} < 5.0\ \mu\text{A}$ under $100\text{ k}\Omega$, the pipeline automatically re-scales with $R_{\text{shunt}} = 10\text{ k}\Omega$ and logs the parameter in the output summary.

### 2. Open-Circuit Voltage Scaling (`Voc`)
$$V_{\text{OC}}(\text{V}) = V_{\text{raw}} \times \text{UNIT\_V} \times \text{VOC\_MULTIPLIER}$$

### 3. Transferred Charge (`Qsc`) & Surface Charge Density ($\sigma$)
For each contact-separation cycle $k$ bounded by $[t_{\text{start},k}, t_{\text{end},k}]$:

$$Q_{+,k} = \int_{t_{+,1}}^{t_{+,2}} \max\big(I_{\text{corr}}(t), 0\big)\,dt \quad [\text{nC}], \qquad |Q_{-,k}| = \int_{t_{-,1}}^{t_{-,2}} \max\big(-I_{\text{corr}}(t), 0\big)\,dt \quad [\text{nC}]$$

$$Q_{\text{sym},k} = \frac{Q_{+,k} + |Q_{-,k}|}{2}, \qquad Q_{\text{SC,max},k} = \max\big(Q_{+,k},\, |Q_{-,k}|\big)$$

$$\sigma_k \left[\mu\text{C/m}^2\right] = \frac{Q_{\text{SC,max},k} \times 10^{-9}\text{ C}}{A \times 10^{-4}\text{ m}^2} \times 10^6 = \frac{10 \times Q_{\text{SC,max},k}\ [\text{nC}]}{A\ [\text{cm}^2]}$$

---

## Standard Analytical Workflow

### Step 1: Data Ingestion, Metadata Parsing & Probabilistic Signal Recognition (`Voc` vs `Isc`)
- **Oscilloscope CSVs**: Automatically skip metadata header lines containing `#` (e.g., `#timebase`, `,#voltbase`, `#size`), extract `#voltbase` vertical scale metadata, and parse numeric `(t_s, raw)` pairs. Sort files in natural numerical order (`WaveData1239.csv` before `WaveData12310.csv`).
- **Probabilistic Signal-Type Classification (`auto_detect_signal_kind`)**:
  When filenames lack explicit `Vol_` or `Cur_` prefixes (e.g., `WaveData1042.csv`), classify waveforms using a tiered physical model and assign an explicit **Confidence Level** (`EXPLICIT`, `HIGH`, `MODERATE`, or `AMBIGUOUS`):
  1. **Tier 1 (Explicit Filename / Header Keyword - `EXPLICIT`)**:
     Matches `vol`, `voc`, `voltage` $\rightarrow$ `Voc`; `cur`, `isc`, `current` $\rightarrow$ `Isc`; `qsc`, `charge` $\rightarrow$ `Qsc`.
  2. **Tier 2 (Oscilloscope `#voltbase` Scale & Paired Channel Context)**:
     - `#voltbase >= 5 V/div` ($5,000,000\ \mu\text{V}$):
       - If $V_{\text{pp,raw}} \ge 4000\text{ mV}$ $\rightarrow$ **`HIGH CONFIDENCE Voc`** (standard open-circuit voltage setup).
       - If $V_{\text{pp,raw}} < 4000\text{ mV}$ $\rightarrow$ **`MODERATE CONFIDENCE Voc`** (likely weak sample or custom setting; prompt user if in doubt).
     - `#voltbase <= 1 V/div` ($1,000,000\ \mu\text{V}$):
       - If $V_{\text{pp,raw}} \le 6000\text{ mV}$ $\rightarrow$ **`HIGH CONFIDENCE Isc`** (standard transimpedance preamplifier range).
       - If $V_{\text{pp,raw}} > 6000\text{ mV}$ $\rightarrow$ **`MODERATE CONFIDENCE Isc`** (verify preamplifier gain).
     - `#voltbase == 2 V/div` ($2,000,000\ \mu\text{V}$):
       - With peer recording: `#voltbase_peer >= 5 V/div` $\rightarrow$ **`HIGH CONFIDENCE Isc`**; `#voltbase_peer <= 0.5 V/div` $\rightarrow$ **`HIGH CONFIDENCE Voc`**.
       - Standalone without peer recording $\rightarrow$ **`AMBIGUOUS`** (intermediate scale; ask user for confirmation).
  3. **Tier 3 (Raw Amplitude Fallback when `#voltbase` is Absent)**:
     - $V_{\text{pp,raw}} \ge 12000\text{ mV}$ $\rightarrow$ **`MODERATE CONFIDENCE Voc`**.
     - $V_{\text{pp,raw}} \le 3500\text{ mV}$ $\rightarrow$ **`MODERATE CONFIDENCE Isc`**.
     - $3500 < V_{\text{pp,raw}} < 12000\text{ mV}$ $\rightarrow$ **`AMBIGUOUS`** (intermediate range; must confirm with user).
- **Human-in-the-Loop Clarification Rule (Crucial)**:
  - **Never present automated predictions as dogmatic certainty** when ground-truth labels are absent. State the physical evidence and the estimated likelihood (e.g., *"There is a high probability (~85-90%) that this file represents Voc based on the 10 V/div scale and 13.6 V peak-to-peak amplitude"*).
  - Whenever confidence is `AMBIGUOUS` or when working interactively on critical publications, **ask the user to confirm the signal type** rather than making unverified assumptions.
- **Multi-Column Excel Workbooks (`.xlsx`)**: Iterate over column pairs `(Col 2j-1 = Time, Col 2j = Signal)`, reading signal type (`Voltage` vs `Current`) from Row 1, units (`V` vs `uA`) from Row 2, sample name from Row 3, and time-series data from Row 4 downward.

### Step 2: Time-Window Slicing & Resting-State DC Baseline Removal
1. **Default Steady-State Window (`DEFAULT_T_START = 1.0`, `DEFAULT_T_END = 7.0`)**:
   - Discard the first $1.0\text{ s}$ ($0 \to 1.0\text{ s}$) to eliminate startup mechanical/electrical transients.
   - Retain the next $6.0\text{ s}$ ($1.0\text{ s} \to 7.0\text{ s}$) for analysis.
   - Discard the trailing segment ($t > 7.0\text{ s}$).
   - When plotting with `X_FROM_ZERO = True`, shift the time axis by $x = t - t_{\text{start}}$ so the horizontal axis spans $0.0\text{ s} \to 6.0\text{ s}$.
2. **Resting-State DC Baseline Estimation**:
   - Compute a 5-point centered rolling standard deviation $s_{\text{local}}(t)$ on the windowed signal.
   - Select flat resting points where $s_{\text{local}}(t) \le \max(\varepsilon, P_{55}(s_{\text{local}}))$.
   - Estimate the DC offset $y_{\text{base}}$ from the median/histogram mode of the resting points (preventing asymmetric positive/negative TENG peaks from biasing the zero line) and subtract it: $y_{\text{corr}}(t) = y(t) - y_{\text{base}}$.

### Step 3: Peak Detection, Phase-Order Recognition & Zero-Gap Midpoint Windowing
1. **Adaptive Peak Detection**:
   - Positive threshold: $\theta_+ = \max\big(\text{floor},\, 0.28 \times P_{98}(\max(y_{\text{corr}}, 0))\big)$.
   - Negative threshold: $\theta_- = \max\big(\text{floor},\, 0.28 \times P_{98}(\max(-y_{\text{corr}}, 0))\big)$.
   - Minimum inter-peak distance: `MIN_DIST_S = 0.12` s.
2. **Peak-to-Peak Pairing (`Voc` and `Isc`)**:
   - Estimate fundamental period $T = \text{median}(\Delta t_{\text{pos}})$.
   - Pair each positive peak $t_p$ with the deepest negative peak $t_q$ within $\pm 0.5\,T$, yielding per-cycle $y_{\max} = y_{\text{corr}}(t_p)$, $y_{\min} = y_{\text{corr}}(t_q)$, and peak-to-peak amplitude $y_{\text{pp}} = y_{\max} - y_{\min}$.
   - Exclude weak outlier cycles where $y_{\text{pp}} < 0.5 \times \text{median}(y_{\text{pp}})$.
3. **Zero-Gap Midpoint Windowing (`Isc -> Qsc` Integration)**:
   - Select the more stable peak train as `anchor_peaks` and the opposite polarity as `partner_peaks`.
   - Measure median separation before ($\Delta t_{\text{before}}$) and after ($\Delta t_{\text{after}}$) each anchor peak to detect whether the cycle phase order is `[- -> +]` or `[+ -> -]`.
   - Place cycle boundary $B_{k+1}$ at the exact midpoint of the quiet resting interval between cycle $k$ and cycle $k+1$ ($c_{\text{end}}[k] = c_{\text{start}}[k+1]$).
   - Integrate the largest positive lobe and largest negative lobe within $[c_{\text{start}}[k], c_{\text{end}}[k]]$ up to zero-crossings ($0.5\%$ foot tolerance) using `scipy.integrate.trapezoid`.

---

## Visualization Standards

### 1. Origin-Style Single/Multi-Trace Plot (`origin_plot`)
- **Figure Dimensions & Resolution**: `figsize = (6.6, 5.4)` inches, exported at `300 DPI`.
- **Typography & Frame**: Bold sans-serif (`Arial`, `Liberation Sans`, `DejaVu Sans`), axis spine width `2.0 pt`, axis labels `24 pt bold`, tick labels `20 pt bold`.
- **Ticks**: Inward-pointing major (`length=8, width=2`) and minor (`length=4, width=1.5`) ticks on bottom and left axes; `MaxNLocator(nbins=5, steps=[1, 2, 5, 10])`.
- **Headroom**: Add $22\%$ vertical headroom above $y_{\max}$ (`head = 0.22` for 1 trace, `0.18 + 0.09 * n` for $n$ traces) and $8\%$ below $y_{\min}$.
- **Top-Right Annotation**:
  - Single trace ($n = 1$): Sample name at `(0.97, 0.97)` in axes coordinates (`ha="right", va="top", fontsize=18-19, fontweight="bold"`).
  - Multiple traces ($n \ge 2$): Frameless legend at `loc="upper right"` with text color matched to line color.
- **Top-Left Peak-to-Peak Annotation (`SHOW_PEAK_PEAK_TEXT = True`)**:
  - Placed at `(0.04, 0.97 - k*0.075)` in axes coordinates (`ha="left", va="top", fontsize=18, fontweight="bold"`).
  - Formatted with mathematical symbols and uncertainty:
    - `Voc`: `$V_{pp}$ = 12.19 ± 2.14 V`
    - `Isc`: `$I_{pp}$ = 1.56 ± 0.43 µA`
    - `Qsc`: `$Q_{sc}$ = 18.45 ± 1.20 nC`

### 2. Four-Panel Charge Integration Diagnostic Figure (`plot_file_analysis`)
- **Figure Dimensions**: `2x2` grid, `figsize = (14.4, 9.4)` inches, `300 DPI`.
- **Panel (a)**: Full raw oscilloscope waveform $I_{\text{raw}}(t)$, horizontal dashed DC baseline, and shaded active analysis window.
- **Panel (b)**: Baseline-corrected short-circuit current $I_{\text{SC}}(t)$ with green circle markers at positive peaks (`t_pos_s`) and red downward-triangle markers at negative peaks (`t_neg_s`).
- **Panel (c)**: Zoomed view of 3 central cycles showing shaded positive lobe area ($Q_+$) and negative lobe area ($|Q_-|$) on the primary axis, overlaid with dual-axis cumulative charge curves $Q_+(t)$ and $|Q_-(t)|$ in $\text{nC}$.
- **Panel (d)**: Grouped bar chart of $Q_+$ and $|Q_-|$ per cycle, overlaid with the $Q_{\text{SC,max}}$ line trajectory and horizontal dashed mean $\pm$ SD reference line (annotated with $\sigma$ in $\mu\text{C/m}^2$).

---

## Reference Implementation

A complete, standalone Python implementation combining both the **Origin-Style Peak-to-Peak Suite** and the **4-Panel $Q_{\text{SC}}$ Integration Engine** is provided in [`scripts/teng_analyzer.py`](./scripts/teng_analyzer.py).

### Quick CLI / Script Usage

```bash
# 1. Analyze a single CSV file (Origin plot + 4-panel Qsc diagnostic)
python scripts/teng_analyzer.py single --input Data/Vol-Cur/Cur_B3-E5.csv --kind Isc --lang vi --qsc-diag

# 2. Batch process a folder of oscilloscope CSV files (window 2.0s -> 7.0s, with custom R_shunt options)
python scripts/teng_analyzer.py folder --input Data/Vol-Cur --t-start 2.0 --t-end 7.0 --lang vi --r-shunt 100000 --r-shunt-repeat 10000

# 3. Batch process a multi-column Excel workbook (e.g., Book1.xlsx -> Data/NIPS)
python scripts/teng_analyzer.py excel --input Data/Book1.xlsx --outdir Data/NIPS --t-start 1.0 --t-end 7.0 --lang vi
```

### Core Python API Usage

```python
from scripts.teng_analyzer import (
    analyze_file,
    origin_plot,
    analyze_teng_qsc_cycles,
    plot_qsc_4panel,
    batch_process_folder,
    batch_process_excel,
)

# Batch-process an entire folder of Vol_*.csv and Cur_*.csv files
df_summary = batch_process_folder(
    folder_path="Data/Vol-Cur",
    window=(1.0, 7.0),
    lang="vi"  # or "en" matching user prompt language
)

# Batch-process a multi-column Excel file into a target directory
df_excel_summary = batch_process_excel(
    excel_path="Data/Book1.xlsx",
    out_dir="Data/NIPS",
    window=(1.0, 7.0),
    lang="vi"
)
```

---

## Verification & Quality Checklist

Before delivering results to the user, always verify:
1. **Language Match**: Confirm that your response text, summary tables, and figure axis labels match the language of the user's prompt.
2. **Unit Integrity**: Check whether raw values are in `mV` or `V`. Under the standard $100\text{ k}\Omega \times 10$ setup with `mV` raw units, TENG $V_{\text{OC}}$ is typically $1\text{--}100\text{ V}$ and $I_{\text{SC}}$ is typically $0.5\text{--}20\ \mu\text{A}$. If values are off by $10^3$, check `CSV_VOLTAGE_UNIT`.
3. **Windowing Consistency**: Verify that the $1.0\text{ s} \to 7.0\text{ s}$ window (or user-specified window) was applied *before* baseline estimation and peak detection, and that the plotted time axis starts cleanly at $0.0\text{ s}$.
4. **Top-Left Annotation Visibility**: Ensure the `$V_{pp}$`, `$I_{pp}$`, or `$Q_{sc}$` annotation in the top-left corner does not collide with the waveform peaks or the top-right sample label.
5. **Charge Conservation Check (`Isc -> Qsc`)**: Compare $Q_+$ and $|Q_-|$ across steady cycles. Large systematic discrepancies indicate either uncompensated DC baseline drift or truncated cycle boundaries.
6. **Probabilistic Confidence & Human-in-the-Loop Clarification**: If the signal type was inferred automatically without ground-truth metadata, explicitly state the confidence level (`HIGH`, `MODERATE`, or `AMBIGUOUS`) and physical evidence. Never present inferences as absolute certainty; for ambiguous or borderline cases, prompt the user to confirm the channel assignment.

