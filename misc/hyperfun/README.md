The following describes creation of shape files for ADDA using [Hyperfun](http://hyperfun.org) - a high-level language for description of particle shape through implicit functional relations.

### 1) Create Hyperfun input file 

Or use example `example.hf`, which describes a deformed prolate spheroid, as used in [this paper]( http://doi.org/10.1111/j.1600-0889.2011.00559.x).

### 2) Convert Hyperfun into WRL file

Under Windows this can be done using the [HyperFun Polygonizer](http://hyperfun.org):
* Start command prompt under Windows
* Change to directory where the Hyperfun input file are
* Execute Hyperfun using 
  ```
  PATH\hfp example.hf -wrl example.wrl -g 180
  ```
  This creates a file in the WRL format, containing a surface description of the particle; the option `-g` sets the grid density of this file.

The HyperFun Polygonizer is [open-source](https://sourceforge.net/projects/hyperfun/), so it should be possible to compile it for other OS as well. But it is not trivial, since no general makefile is provided (only the project file for Visual Studio). 

!!! There is important issue with the current version of Polygonizer (2.03). It seems that it chooses arbitrary order of vertices when defining the faces. I.e. it is not consistent with vertex normals that are also written to WRL file. The problem is that software that reads these files (see below) discards explicit normals and deducts face alignment from the vertex order. The same (arbitrary) alignment is then saved to OBJ file, which breaks down the 'pip' program.
  
There are two solutions to this issue:
* Either use additional command line option `-usedc 0` to use the older mesh generation routines in Polygonizer. This may be inferior in terms of speed and size of mesh, but should be free of the described issue,
* or align normals during the next step.

### 3) Convert WRL file to OBJ file

Under Windows this can be done using [Accutrans](http://www.micromouse.ca), which is trialware, but trial period is not enforced. Open WRL file and save it as a Wavefront OBJ file - you will also see the 3D image of the model - which is nice for extra verification. In between, you may align the face normals of the model, if needed. Click "Align->Flip Polygons - Start/Stop" in the menu. The faces will be colored either white (correct) or red (wrong). If you see red faces click "Align->Preset All Layers". It is also possible to use Accutrans from the command line, like `PATH\at3d_2-13-0 -noshow example.wrl` but it needs preliminary setting of the output format inside the program. It is also not possible to align normals in this (batch) mode.

Another GUI application (GPL, available for all OS) is [MeshLab](http://meshlab.sourceforge.net/). It can also align normals using "Filters->Normals, Curvatures and Orientation->Re-Orient all faces coherently" in the menu.

There is also a simple command-line tool [meshconv](https://www.patrickmin.com/meshconv/), available for any OS. It should be used like `meshconv -c obj example.wrl` but it does not allow alignment of normals.

### 4) Convert OBJ file into DDSCAT shape format (that is readable by ADDA)

For that you should use "Point in Polyhedron" tool, provided in [misc/pip](../pip) in ADDA repository.

Alternatively, you may use the following script to generate several discretized shape files corresponding to different size parameters (and fixed dpl). For that you will need
* Bash. Under Windows you may use [MSYS](http://www.mingw.org/wiki/MSYS)
* Compiled executable of `calc_pip_arg.f90`. For that run `make` in current directory. 64-bit Windows executable is available in the corresponding [ADDA package](https://github.com/adda-team/adda/releases).

To execute the tool chain, first, adjust `make_datfiles`:
* variable `dpl` defines the number of dipoles per wavelength, 
* array `size_array` defines the size parameters, for which the shape files should be created.

Second, run `./make_datfiles` with the name of the obj-file, without extension, as argument, for example `./make_datfiles example`. This will create the files for ADDA. In particular, the shape file for size parameter 12.0 will be `example_12.0.dat`. The computation time can be very long (several days).

### 5) Simple test using only command line tools

The following can also be performed also with GUI tools described above. Running the following three commands should produce `shape.dat` similar (if not identical) to the provided `shape_test.dat`.
```
PATH\hfp example.hf -wrl example.wrl -g 180 -usedc 0
PATH\meshconv -c obj -o shape example.wrl
PATH\pip 20
```

Suggestions are welcome at josef.gasteiger@lmu.de
