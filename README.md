### Important note
This project has been moved from the Google Code platform. Most functionality is already here, but the migration has not been completed yet - the progress can be tracked at [#210](https://github.com/adda-team/adda/issues/210). In particular, the wiki pages are located on a [separate branch](../../tree/wiki) and some of their content, related to the Google Code or Subversion, is outdated. If you encounter any artefacts, you may also turn to the [archived repository](https://code.google.com/archive/p/a-dda). All links to the previous repository should be automatically redirected to the archived one, but the formatting of some wiki pages there is broken.

### Description

`ADDA` is a C software package to calculate scattering and absorption of electromagnetic waves by particles of arbitrary shape and composition using the discrete dipole approximation (DDA). Its main feature is the ability to run on a multiprocessor system or multicore processors (parallelizing a _single_ DDA simulation). It can also employ modern GPUs to accelerate computations. `ADDA` is intended to be a [versatile tool](../wiki/Features.md), suitable for a wide variety of applications ranging from interstellar dust and atmospheric aerosols to metallic nanoparticles and biological cells. Its applicability is limited only by [available computer resources](../wiki/LargestSimulations.md).

`ADDA` [originated](../wiki/EarlyHistory.md) at the University of Amsterdam but has then evolved into an open-source international project. We recommend to start with one of the following:
* [Getting Started](../wiki/GettingStarted.md)
* [Manual](doc/manual.pdf)
* [FAQ](../wiki/FAQ.md)
* [Features](../wiki/Features.md)
* [Release Notes](../wiki/ReleaseNotes.md)
* [Downloads](../wiki/Downloads.md)
* [Compiling ADDA](../wiki/CompilingADDA.md)

If you choose to use `ADDA`, please [subscribe](mailto:adda-announce+subscribe@googlegroups.com) to [announcement mailing list](http://groups.google.com/group/adda-announce); "registered" users of `ADDA` will be notified when updates to the code are made. If you publish results obtained using `ADDA`, you should acknowledge the source of the code. The general reference is - M. A. Yurkin and A. G. Hoekstra, “The discrete-dipole-approximation code ADDA: capabilities and known limitations,” [J. Quant. Spectrosc. Radiat.  112, 2234-2247 (2011)](http://dx.doi.org/10.1016/j.jqsrt.2011.01.031).
One may also look at [a list of more specific references](../wiki/References.md).

We encourage users to provide feedback in any possible way. If you have any questions, please [write](mailto:adda-discuss@googlegroups.com) to the [discussion group](http://groups.google.com/group/adda-discuss). If you have any suggestions or a bug report, submit it directly to the [issue tracker](https://github.com/adda-team/adda/issues), taking advantage of the open development process. Please also provide feedback on the existing issues - you may star issues that you find important (Google account is required) and/or provide some meaningful comments.

Finally, we maintain a [list of publications that use ADDA](../wiki/Publications.md). Please let us know of any papers that should be included in this list. If you are doing very large simulations with `ADDA`, consider them to be included in the corresponding [wiki page](../wiki/LargestSimulations.md). There is also a list of [ADDA awards](../wiki/Awards.md).
