"Point in Polyhedron" transforms `.obj` files into DDSCAT7 shape format (readable by ADDA).

Authors: Roman Schuh and Thomas Wriedt, using routines by John Burkardt.

To compile run `make` in current directory. 64-bit Windows executable is available in the corresponding [ADDA package](https://github.com/adda-team/adda/releases).

Usage: 
```
pip [<grid> [<filename>]]
```

Executable is named `pip`, it accepts up tp two command line parameters: 
- `<grid>` - maximum shape size (number of dipoles) along the largest dimension that is determined automatically. If omitted, the default value of 80 is used (and further arguments cannot be used).
- `<filename>`. Input shape is read from `<filename>.obj` and DDSCAT7 shape is saved into `<filename>.dat`. If omitted, PIP will use the default filename `shape` (read `shape.obj` and save `shape.dat`). If `<filename>` includes extension, then it will be used for input file (instead of `.obj`).

It should be possible to read other 3D formats, which are supported by routines in `ivread_wr.f90` - see comments in the source files. However, only .obj format is sufficiently tested. Moreover, the algorithm requires consistent alignment of face normals - they should all point outward. Should work for very large number of dipoles (limited only by memory and computational time).

To test executable run `pip 20` in the current folder. It will use provided `shape.obj` that defines a star-shaped object and produce `shape.dat` which should be identical to provided `shape_test.dat`.

The reference for this tool is: R. Schuh, "Arbitrary particle shape modeling in DDSCAT and validation of simulation results", in T. Wriedt & A. G. Hoekstra (Editors): [Proceedings of the DDA - Workshop](http://diogenes.iwt.uni-bremen.de/vt/laser/papers/DDA-Workshop-final-proceedings.pdf), Institut fur Werkstofftechnik, Bremen (23 March 2007), p.22-24.