Various examples of using ADDA. Mostly these are tutorial examples, which illustrate certain features of ADDA and should not consume too much computational time. The only exception is subfolder `papers/`, which additionally aims to reproduce specific published simulations (thus, may require a large cluster). All examples are designed to be run out of the box assuming that compiled binaries are available (either by default compilation or pre-compiled for Windows). Some only run ADDA with proper command lines, others - also plot or post-process the obtained results. See further details inside corresponding folders:
* `bessel/` - examples involving excitation by Bessel beams
* `eels/` - examples involving excitation by a fast electron (will be added after #155)
* `papers/` - examples reproducing published results
* `find_adda` - common script to locate ADDA binaries and set corresponding environmental variables
