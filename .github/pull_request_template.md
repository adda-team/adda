_Thank you for contributing to ADDA. Please remove the Italic text (`_..._`) throughout this template after following the corresponding instructions_

### Description
_Please explain the changes you made here. Answer both why and what_

### Related issues
_Any nontrivial pull request should first be discussed in the issue_ 

Fixes #...  or Related to #...

### Types of changes
_What types of changes does your code introduce to ADDA? Put an `x` in the boxes that apply_

- [ ] Bugfix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)

### Checklist
_Put an `x` in the boxes that apply. You can also fill these out after creating the PR. If you're unsure about any of them, don't hesitate to ask. This is simply a reminder of what we are going to look for before merging your code._

- [ ] I have read the [contributing guidelines](https://github.com/adda-team/adda/wiki/InstructionCommitters)
- [ ] Code compiles correctly in all relevant regimes (at least, `make seq`)
- [ ] New code does not rely on any Fortran or C++ sources or is disabled by `NO_FORTRAN` or `NO_CPP` preprocessor macros, respectively. 
- [ ] No warnings appear during debug compilation (at least, `make seq OPTIONS=DEBUG`, but better `devtools/build_debug`) or they are discussed below
- [ ] Tests pass locally with my changes (at least, `sh comp2exec seq` in `tests/2exec`, but better `devtools/test_new [seq]`). If any errors appear, they are discussed below.
- [ ] I have added tests that prove my fix is effective or that my feature works. And these tests pass. This includes new command line in suite files in `tests/2exec` (and potentially new ignore patterns). In some cases, it is desirable to add new tests to `tests/equiv`. 
- [ ] I have added/extended necessary documentation (if appropriate). If suggesting changes to the manual, I have used "Track changes" in the doc file.

### Further comments
_If this is a relatively large or complex change, explain why you chose the solution you did and what alternatives you considered. 
If that was discussed in issue or somewhere else, refer to it here_
