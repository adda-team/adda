Utils for automatic testing of ADDA by running equivalent command line combinations and comparing the output. Currently, it contains only a single test for superellipsoids, but should be further developed into a general framework (similar to `tests/comp2exec`).

### Directory structure

* `superellipsoid.sh` - script to test superellipsoid shape

* `bb_equiv.py` - the comparison of equivalent cases (command lines) for the scattering of Bessel beams calculated in ADDA. Differences are shown when |diff| > 1e-08. See results in bb_results.txt