This folder contains ADDA executables for 64-bit Windows and all required DLLs. Hence, these executables are ready to run, except for:
* MPI executables (including `misc/neafield/nearfield.exe`) require Microsoft MPI [to be installed](https://github.com/adda-team/adda/wiki/InstallingMPI#windows)
* OpenCL executables require GPU with [OpenCL support](https://github.com/adda-team/adda/wiki/OpenCL#windows).
* All auxilliary files, like the [manual](https://github.com/adda-team/adda/blob/master/doc/manual.pdf) or [sample input files](https://github.com/adda-team/adda/tree/master/input), are supplied with the main (source) ADDA package. It should be [obtained separately](https://github.com/adda-team/adda/releases). Note that those text files will most probably contain Unix-style line endings. This is OK for ADDA but may cause inconveniences if you try to edit them in Windows Notepad. Use [more universal editor](https://notepad-plus-plus.org/) instead.

All executables are console applications, which should be run from a terminal (not by clicking the mouse). On Windows either use built-in `cmd.exe` or install one of Linux-type shells (like MSYS Bash or Git Bash) or FAR File Manager. You can also run ADDA through a system call from any programming/scripting language.

We also provide executables for `misc/` packages. They are located in separate folders for convenience, but require DLLs inside `misc` folder. For instructions on using this executables see the corresponding `misc/...` folders in the main ADDA package. 

In general, one either need to add path to required DLLs to the PATH environmental variable or to keep executables and DLLs in the same folder.

### Directory structure

* `misc/` - executables for miscellaneous packages and three DLLs (provided by gfortran 8.1 in a MinGW64 suite). 
* `adda.exe` - executable of sequential version of ADDA
* `adda_mpi.exe` - executable of MPI version of ADDA
* `adda_ocl.exe` - executable of OpenCL version of ADDA
* `adda_spa.exe` - executable of sequential sparse version of ADDA
* `adda_spa_mpi.exe` - executable of MPI sparse version of ADDA
* `clFFT.dll` - DLL for clFFT 2.12.2
* `libfftw3-3.dll` - DLL for FFTW 3.3.5
