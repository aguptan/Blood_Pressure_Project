# BP Oscillometry (MATLAB)

MATLAB tools for analyzing oscillometric blood pressure (BP) signals.
This repository includes signal segmentation, filtering, envelope analysis, and validation routines for estimating Mean Arterial Pressure (MAP) from cuff oscillograms.

---

## Features

* Segment oscillogram signals into inflation/deflation cycles
* Design and apply high-pass Butterworth filters
* Extract oscillation envelopes and locate MAP peak
* Compare estimates against device-reported MAP
* Summarize trial results in clean tables
* Validation scripts for regression and group delay analysis

Function Role:

* Load test oscillogram data
* Segment cycles
* Filter & compute envelopes
* Estimate MAP
* Compare with device values
* Print and plot results

---

## Repository Structure

```
bp-oscillometry/
│
├─ src/+bp_osc/           # Core functions
│   ├─ analyze_oscillogram.m
│   ├─ compare_with_device.m
│   ├─ design_hp_filter.m
│   ├─ load_bp_data.m
│   ├─ segment_cuff_cycles.m
│   ├─ summarize_results.m
│
├─ AlgorithmTesting/      # Example script + demo figures
│   └─ run_bp_oscillometry.m
│
├─ validation/            # Validation scripts & result figures
│   ├─ AlgorithmValidation.m
│   ├─ GroupDelayAnalysis.m
│
├─ README.md              # This file
```

---

## Data

**No raw patient/device data is included.**


## Technologies

* **Language**: MATLAB (R2021a or newer)
* **Toolboxes**: Signal Processing Toolbox
* **Algorithms**: FIR/IIR filtering, Butterworth design, oscillation envelope detection, regression validation

---




