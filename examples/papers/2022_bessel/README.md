The following Python scripts reproduce the data plotted in four figures of the paper: Glukhova S.A. and Yurkin M.A. Vector Bessel beams: General classification and scattering simulations, [_Phys. Rev. A_ **106**, 033508](https://doi.org/10.1103/PhysRevA.106.033508) (2022). The data was originally produced with ADDA v.1.5.0-alpha (443f00a).

The module `bb_module.py` includes common functions, while the following python scripts run ADDA and produce corresponding figures:
* `fig12_sphere.py`
* `fig13_coated_sphere.py`
* `fig14_extrapolation.py` (requires 8.2 GB of RAM)
* `fig15_cube.py`

Other folders (some are produced by scripts):
* `dda/` - contains raw ADDA output in subfolders `/option_#/` (0 - Fig.12; 1 - Fig.13; 2-4 - Fig.15). `/extrapolation/` - contains several `run_...` folders with raw ADDA output for Fig.14.
* `ref/glmt/` - contains raw data for the scattering intensities of Bessel beams obtained with the generalized Lorenz-Mie theory (GLMT) for Figs. 12 and 13 (`/option_0/` and `/option_1/` subfolders, respectively).
* `particles/` - contains images of scattering particles presented in Figs. 12,13, and 15.
* `saved/` - final figures in PDF format.
