The included Python3 scripts reproduce Figure 8 and parts of Figure 11 of the paper: 
Virkki, A. and Yurkin, M. A. (2025), Microwave scattering by rough polyhedral particles on a surface, [JQSRT, In revision](https://arxiv.org/abs/arXiv:2501.10019). 
The computations were originally conducted using ADDA v.1.4.0, but also tested with v.1.5.0-alpha3.

To run ADDA, you can use the Python script `makerunfiles_surf.py` to generate a run file. By default, it generates the file `run_sphere` for 10 commands: 8 for a surface case (4 incidence angles for 2 size parameters) and 2 for spheres of corresponding sizes in free space. The file can be executed immediately if ADDA binaries are available in one of the standard locations. 

Folders produced by the script `makerunfiles_surf.py`:
* `SurfSphere/` - For raw ADDA output for spheres on (and in touch with) a surface
* `NoSurfSphere/` - For raw ADDA output for spheres in free space

The following python scripts produce corresponding figures `fig8.pdf` and `fig11A.pdf` (or `fig11C.pdf`):
* `fig8.py` – Figure 8: a mosaic of the scattering-matrix elements (1,1), (1,2), (2,2), (3,3), and (4,4) in the scattering plane for a sphere with a size parameter 6, with a comparison between the surface case at normal incidence and the back- and forward-scattering in the free-space case in the first column.
* `fig11.py` – Figure 11A (default; sp = 1.0) or Figure 11C (change sp = 6.0 on line 19): log10-scale intensity (1,1-element) in polar plots for the zenith and nadir views and one for both hemispheres in the scattering plane (along the horizontal axis of the polar plots).

The scripts compute and plot by default the incidence angles 0°, 20°, 40°, and 60°.

The steps to run the (unchanged) scripts on the command line to generate `fig8.pdf` and `fig11A.pdf`:
```bash
./makerunfiles_surf.py
sh run_sphere
./fig8.py
./fig11.py 
```
(On Windows you may need to run Python3 scripts by prepending them with `python` instead of `./`, depending on the Python installation).

The data for polyhedrons from the same paper takes a lot of time to recompute. Ensemble-averaged scattering matrices for polyhedrons are stored at https://doi.org/10.5281/zenodo.15040283, and separate scripts are available at https://github.com/a-virkki/ADDA-grid-plot-codes to plot some of the paper's figures based on these data.
