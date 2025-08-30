# Holmes et al. (2019)

This code is copied from [ContinuousMIEstimation](https://github.com/EmoryUniversityTheoreticalBiophysics/ContinuousMIEstimation/) with the following changes:

+ matlab `.m` files in the root directory have been deleted
+ In `MIxnyn.C:
    - remove inclusion of `mex.h`
    - comment out function `mexFunction`
    - included `extern "C"` around MIxnyn
    - mention above changes in preamble comment
+ compiled the shared library `libMIxnyn.so` via the command `gcc -shared -fPIC -Wl,--export-dynamic -o libMIxnyn.so MIxnyn.C`
