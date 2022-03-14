Utils for automatic testing of ADDA by running equivalent command line combinations and comparing the output. Currently, it contains only a few tests, but should be further developed into a general framework (similar to `tests/comp2exec`).

### Directory structure

* `bb_equiv.py` - test for Bessel beams. Differences are shown when |diff| > 1e-08. See results in bb_results.txt
* `superellipsoid.sh` - script to test superellipsoid shape

