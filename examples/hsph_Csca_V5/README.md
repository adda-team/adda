ADDA Simulation GUI - Hemisphere Integration (hsph) - Version 5.0
Author: Mattia Andrini (mattia.andrini@unicatt.it)



DESCRIPTION
Computes integrated Forward/Backward scattering efficiencies over user-defined hemispherical domains, accounting for full azimuthal (phi) dependence. The tool also automates ADDA simulations over a wavelength spectrum and features Free Space \& Surface modes, custom shape support, and Mie theory comparison. Read the guide for a detailed description.



REQUIREMENTS

* Python 3.x (Libs: numpy, matplotlib, scipy, tkinter). miepython is optional.
* Executable: 'adda.exe' (CPU) or 'adda\_ocl.exe' (GPU acceleration).
  The script auto-detects these in the local folder, 'win64' subfolder, or standard source directories.



QUICK START

1. Run the script: python hsph\_V5.py
2. GUI: Click "Browse" to select Input File.
3. Set parameters (size, grid, etc.) and click "Run, Upload and Plot".



INPUT FILE FORMAT (.txt)
Provide 4 space-separated columns (no headers). One row per wavelength.
Format: \[Wavelength(nm)] \[Particle Real Refractive Index n] \[Particle Imaginary Refractive Index k] \[Substrate n]
Example: 400 1.50 0.01 1.0
See included example file: lambda\_nk\_test.txt



OUTPUTS
Saved in "SimulationResults" folder:

* Simulation\_result\_Q\_\*.txt: Efficiency Factors (Qext, Qsca, Qabs) \& Mie data.
* Simulation\_result\_C\_\*.txt: Scattering Cross Sections (micrometers^2).



CITATION \& CONTACT
If used for research, please cite: https://doi.org/10.1021/acs.analchem.5c00168
For questions: mattia.andrini@unicatt.it

