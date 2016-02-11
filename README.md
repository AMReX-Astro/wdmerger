# wdmerger
*A software package for simulating white dwarf mergers with CASTRO*

## Getting Started

There is a user's guide in `documentation/` that will guide you 
through installing the required software and running simulations. 
Enter `make` at the command line to build the user's guide. You 
can also find the user's guide (updated nightly) at 
http://astro.sunysb.edu/mkatz/files/wdmergeruserguide.pdf.

## Branches

The `master` branch is the stable branch, which should in general
work correctly. Active development of the code is done on the
`development` branch; this has the latest features but is not
guaranteed to be error-free. This mirrors the structure of the
other BoxLib codes. Consequently, for whichever branch you choose,
you should choose the same branch in `CASTRO`, `BoxLib`, and
`Microphysics` (if you are using that) to ensure that all codes
are in sync with each other. When you do a `git pull` on any one
of these codes, you should do a `git pull' on all the others, as
we often submit simultaneous updates to multiple codes that are
intended to all work together.

## Contact

Please contact Max Katz (maximilian.katz@stonybrook.edu) at
Stony Brook University if you need help with `wdmerger`.
General `CASTRO` and `BoxLib` questions should be directed
to the `CASTRO` user's mailing list/Google group,
castro-help@googlegroups.com.
