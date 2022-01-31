The following Bash scripts reproduce the data plotted in two figures of the paper: Inzhevatkin K.G. and Yurkin M.A. Uniform-over-size approximation of the internal fields for scatterers with low refractive-index contrast, [_J. Quant. Spectrosc. Radiat. Transfer_ **277**, 107965](http://doi.org/10.1016/j.jqsrt.2021.107965) (2022). Also contains Python scripts to calculate approximations of the internal field using various versions of the Wentzel-Kramers-Brillouin (WKB) approximation and interface to exact Mie solution for a sphere (in separate folders).

The data was originally produced with ADDA v.1.4.0 and the scripts were last tested with ADDA v.1.5.0-alpha (6c0b401). During final testing, we uncovered that they produce slightly different convergence results for "WKBr I" (simpler version). That is because a preliminary version of the algorithm of "WKBr I" was erroneously used for producing Figs. 15 and 16 (this error does not affect other results in the paper). We decided not to publish a corrigendum, since the convergence curves are qualitatively the same and all conclusions remain valid.

Most scripts allow manual configuration of the parameters inside them, but can also be run out of the box. The only necessary extra effort is providing external tools for calculation of reference Mie solution (see the corresponding folder). The directory structure is the following:

* `cluster/` - realistic scripts for working at the compute cluster, require very large simulation time
* `Mie_solution/` - interface for calculation of the exact electric field in a sphere
* `WKBr/` – scripts for calculation of the WKBr approximation (of internal field) for a sphere
* `fig15.sh`, `fig16.sh` – simplified calculation scripts (wrappers) for the corresponding figures of the paper
* `calc.sh` - main calculation script, contains a lot of descriptive comments
