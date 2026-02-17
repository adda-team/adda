# ADDA Simulation GUI - Hemisphere Integration (hsph)

Version: 5.1  
Author: Mattia Andrini (mattia.andrini@unicatt.it)

## Description

Computes integrated Forward/Backward scattering efficiencies over user-defined hemispherical domains, accounting for full azimuthal (phi) dependence. The tool also automates ADDA simulations over a wavelength spectrum and features Free Space \& Surface modes, custom shape support, and Mie theory comparison. Read the guide (`hsph_Csca_guide.pdf`) for a detailed description.

## Requirements

* Python 3.x (Libs: `numpy`, `matplotlib`, `scipy`, `tkinter`). `miepython` is optional.
* Executable: `adda.exe` (CPU) or `adda_ocl.exe` (GPU acceleration).
  The script auto-detects these in the local folder, `win64/` subfolder, or standard source directories.

## Quick Start

1. Run the script: `python hsph_Csca.py`
2. GUI: Click "Browse" to select Input File.
3. Set parameters (size, grid, etc.) and click "Run, Upload and Plot".

## Input File Format (`.txt`)
Provide 4 space-separated columns (no headers). One row per wavelength. Format: 
```
<Wavelength(nm)> <Particle Real Refractive Index n> <Particle Imaginary Refractive Index k> <Substrate n>
```
Example: `400 1.50 0.01 1.0`

See included example file: `lambda_nk_test.txt`

## Outputs
Saved in `SimulationResults/` folder:

* `Simulation_result_Q_*.txt`: Efficiency Factors (Qext, Qsca, Qabs) \& Mie data.
* `Simulation_result_C_*.txt`: Scattering Cross Sections (micrometers^2).

# Citation & Contact

If used for research, please cite: https://doi.org/10.1021/acs.analchem.5c00168

For questions: mattia.andrini@unicatt.it

