Utils for automatic testing of ADDA by running equivalent command line combinations and comparing the output. Currently, it contains only a few tests, but should be further developed into a general framework (similar to `tests/comp2exec`).

### Directory structure

* `bb_equiv.py` - tests for Bessel beams. Relative differences are shown when > 1e-08. See results in bb_results.txt
* `superellipsoid.sh` - test for superellipsoid shape
