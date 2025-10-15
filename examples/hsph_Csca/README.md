The `hsph_Csca_ADDA_V2` python code allows the user to perform controlled simulations using ADDA 1.4.0, read and save the results, while integrating the scattering cross section in two angular hemispheres.
The code is featuring also a user-friendly GUI for Windows environments (10 and 11 have been tested). 
Please, read carefully the guide and the comments in the code as it is still an experimental version.

To start, put `hsph_Csca_ADDA_V2.py`, `lambda_nk_test.txt` and `all_dir_params.dat` in the same folder where the ADDA executable is (usually `win64/`), then run the python script.
An user-friendly and self-explanatory GUI should appear, with some preset parameters. 
The preset values simulate a sphere of 50nm radius, 32x32x32 scattering grid, in free space, with a theta grid going from 0° to 180° with steps of 0.5°,
a phi grid going from 0° to 360° with steps of 2°, and with internal ADDA integration to check the consistency.
To run the simulation click on the "Run" button. To perform the integrations and save the results, click on the "Upload" button.
To plot the result, click the "Plot" button. These three can be combined in one go by clicking the first button "Run, upload and plot".

A more detailed explanation of the code functioning and parameters meaning is included in the pdf `hsph_Csca_guide_V2`. If you use this script, please cite the following paper (where the core functionality was developed):
Andrini M., Federici S., and Gavioli L. Refractive index of benchmark polystyrene nanoplastics by optical modeling of UV–vis spectra, [_Anal. Chem._ **97**, 19419–19426](http://doi.org/10.1021/acs.analchem.5c00168) (2025).

The directory structure is the following:
- `alldir_params.dat` - default file with parameters of standard Csca calculation (to be used as a reference)
- `hsph_Csca_ADDA_V2.py` - main script (executable)
- `hsph_Csca_guide_V2.pdf` - user guide
- `lambda_nk_test.txt` - sample input file with refractive-index spectrum
- `miepython.py` - [external](https://miepython.readthedocs.io) Lorenz-Mie code

