Tools used by ADDA developers for building, testing, and packaging of ADDA:
* `build` - script to build ADDA with various flags and optionally copy the obtained executables
* `build_degug` - compiles ADDA (in debug mode) with various compilation flags
* `build_misc` - builds misc tools and optionally copies (installs) the obtained executables
* `commit_docs` - script to commit `docs/` and `src/const.h` for release
* `release_sequence.txt` - sequence of tasks to be made during release
* `svg_images.txt` - typical workflow (guidelines) for producing `.svg` images
* `test_new` - simple wrapper to build and test ADDA
* `versions.txt` - note on files that need to be checked prior to release (e.g. contain version number)
* `win_all.sh` - script to perform all Windows-specific tasks during release
* `win_commit.sh` - script to run a sample simulation and commit `sample/` and  `win64/`
* `zip_packages` - script to export (archive) source and binary packages for a given version


