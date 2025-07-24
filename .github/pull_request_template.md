<!--Thank you for contributing to ADDA. The following contains the instructions in the comments. You may remove or leave them.-->

<!--The following is a general template for changes to the main source code of ADDA. There are a few simpler ones for other cases, listed below. Click them in the Preview mode. If proceeding with a general template, please remove this list.-->

### _Simpler specialized templates (click or remove)_
- [Miscellaneous tools](?quick_pull=1&template=misc.md&title=[MISC]+_Replace+With+Suitable+Title_&labels=comp-Misc)
- [Examples](?quick_pull=1&template=example.md&title=[EXAMPLE]+_Replace+With+Suitable+Title_&labels=comp-Example)

### Summary

<!--Briefly describe the new feature(s), enhancement(s), or bugfix(es) included in this pull request.-->

### Related issues

<!--Any nontrivial pull request should first be discussed in the issue. If such issue exists, please mention the issue number here as `related to #...`. Also refer to issues, which discuss possible implementation options, if you chose one of them. Use the phrases `fixes #...` or `closes #...`, when you want an issue to be automatically closed when the pull request is merged.-->

### Types of changes

<!--What types of changes does your code introduce to ADDA? Put an `x` in the boxes that apply (replacing the space between square brackets).-->

- [ ] Bugfix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds or improves functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)

### Implementation notes

<!--Provide any relevant details about how the changes are implemented, how correctness was verified, how other features - if any - are affected. If this is a relatively large or complex change, explain why you chose the solution you did and what alternatives you considered. If that was discussed in issue or somewhere else, refer to it here.-->

### Checklist

<!--Put an `x` in the boxes that apply (replacing the space between square brackets). Typically, all boxes need to be checked before the final merge, but you can also fill these out after creating the PR. If you're unsure about any of them, don't hesitate to ask. If you think that some of them are not relevant, discuss this above.-->

- [ ] I have read the [contributing guidelines](https://github.com/adda-team/adda/wiki/InstructionCommitters).
- [ ] The new code complies with the existing [code style](https://github.com/adda-team/adda/wiki/CodeStyleGuide).
- [ ] The code compiles correctly in all relevant regimes (at least, `make seq`).
- [ ] The new code does not rely on any Fortran or C++ sources or is disabled by `NO_FORTRAN` or `NO_CPP` preprocessor macros, respectively.
- [ ] The change neither adds or removes files; otherwise, these changes are reflected in `README.md` in corresponding folders.  
- [ ] No warnings appear during debug compilation (at least, `make seq OPTIONS=DEBUG`, but better `devtools/build_debug`) or they are discussed above.
- [ ] Tests pass locally with my changes (at least, `sh comp2exec seq` in `tests/2exec`, but better `devtools/test_new [seq]`). If any errors appear, they are discussed above.
- [ ] I have added tests that prove my fix is effective or that my feature works. And these tests pass. This includes new command line in suite files in `tests/2exec` (and potentially new ignore patterns). In some cases, it is desirable to add new tests to `tests/equiv`. 
- [ ] I have added/extended necessary documentation (if appropriate). If suggesting changes to the manual, I have used "Track changes" in the doc file.
- [ ] I have looked through all changes introduced by this pull request (line by line), using `git diff` or, better, some GUI tool, to ensure that no unexpected changes are introduced.

