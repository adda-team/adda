Utils for automatic cross-testing different ADDA executables (by comparing the simulation results with some tolerance). The main part is a bash script `comp2exec`. Under Windows one may use MSYS to run it. This script contains comments in the head and throughout the file to explain its usage. In particular, one may specify the particular executables to cross test (a few typical pairs can be selected by command line options to this script).

The particular tests to perform are specified by files `suite*`, by a list of command lines with certain additions and abbreviations (see comments in the file). Currently, these files provide an extensive list of tests for testing future ADDA versions against a stable release. Some of this tests may normally produce errors, however they can be silenced by tuning the ignore patterns inside `comp2exec`. A different suite file can also be used, if specified by command line option to 'comp2exec'.

### Directory structure

* `comp2exec` - main script to run the test
* `diff_numeric.awk` - script to compare two files with given numerical tolerance
* `suite` - the default test suite
* `suite_rd` - shortened suite for testing rectangular-dipoles mode
* `suite_sparse` - shortened suite for testing sparce mode
* `suite_surf` - shortened suite for testing sufrace mode
* `test_all` - a wrapper to run comp2exec with several predefined sets of parameters to perform extensive testing
* ...  - many sample input files to be used in tests
