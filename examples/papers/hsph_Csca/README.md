The hsph_Csca_ADDA_V2.0 python code allows the user to perform controlled simulations using ADDA 1.4.0, read and save the results, while integrating the scattering cross section in two angular hemispheres.
The code is featuring also a user-friendly GUI for Windows environment (10 and 11 have been tested).
Please, read carefully the guide and the comments in the code as it is still an experimental version.
To start, put hsph_Csca_ADDA_V2.0.py, lambda_nk_test.txt and all_dir_params.dat in the same folder where the adda executable is (usually win64), then run the python script.
An user-friendly and self-explanatory GUI should appear, with an preset parameters that should run out of the box.
The preset values simulate a sphere of 50nm radius, 32x32x32 scattering grid, in free space, with a theta grid going from 0° to 180° with steps of 0.5°,
a phi grid going from 0° to 360° with steps of 2°, with internal ADDA integration to check the consistency.
To run the simulation click on the "Run" box. To perform the integrations and save the results, clik on the "Upload" box.
To plot the result, click the "Plot" box.