The included Python3 scripts reproduce Figure 8 and parts of Figure 11 of the paper: 
Virkki, A. and Yurkin, M. (2025), Microwave scattering by rough polyhedral particles on a surface, [JQSRT, In revision](https://arxiv.org/abs/arXiv:2501.10019). 
The computations were originally conducted using ADDA v.1.4.0.

To run ADDA, you can use the Python script `makerunfiles_surf.py` to generate a run file. By default, it generates the file `run_sphere` for 10 commands: 8 for a surface case (4 incidence angles for 2 size parameters) and 2 for spheres of corresponding sizes in free space. The file has to be in the same location as ADDA and has to have the permission to run. 

Folders produced by the script `makerunfiles_surf.py`:
* `SurfSphere/` - For raw ADDA output for spheres on (and in touch with) a surface
* `NoSurfSphere/` - For raw ADDA output for spheres in free space

The following python scripts produce corresponding figures `fig8.pdf` and `fig11A.pdf` (or `fig11C.pdf`):
* `fig8.py` – Figure 8: a mosaic of the scattering-matrix elements (1,1), (1,2), (2,2), (3,3), and (4,4) in the scattering plane for a sphere with a size parameter 6, with a comparison between the surface case at normal incidence and the back- and forward-scattering in the free-space case in the first column.
* `fig11.py` – Figure 11A (default; sp = 1.0) or Figure 11C (change sp = 6.0 on line 19): log10-scale intensity (1,1-element) in polar plots for the zenith and nadir views and one for both hemispheres in the scattering plane (along the horizontal axis of the polar plots).

The scripts compute and plot by default the incidence angles 0°, 20°, 40°, and 60°.

The steps to run the (unchanged) scripts on the command line to generate `fig8.pdf` and `fig11A.pdf`:
* `python3 makerunfiles_surf.py`
* `chmod u+x run_sphere`
* `./run_sphere`
* `python3 fig8.py`
* `python3 fig11.py` 

