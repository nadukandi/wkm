# Computing the Wave-Kernel Matrix Functions

## About ###

This repository contains two MATLAB functions to compute (`wkm.m`) and test 
(`test_wkm.m`) the Wave-Kernel Matrix Functions `cosh(sqrt(A))` and 
`sinhc(sqrt(A))`, where `A` is **any** `n x n` matrix, `sinhc(z) = sinh(z)/z` 
for any nonzero scalar `z`, and `sinhc(0) = 1`.
Details on the underlying algorithms can be found in the
[MIMS EPrint 2018.4](http://eprints.maths.manchester.ac.uk/2621/).

The function `wkm.m` has the following signature.
```
[C, S] = wkm (A)
```
where `C = cosh(sqrt(A))` and `S = sinhc(sqrt(A))`.
Type `help wkm` at the MATLAB command prompt for more information.

The function `test_wkm.m` runs some tests for a suite of matrices chosen from
the [Matrix Computation Toolbox](http://www.maths.manchester.ac.uk/~higham/mctoolbox)
and the matrix function literature.
```
test_wkm('mct_testmats');
test_wkm('expm_testmats');
test_wkm('logm_testmats');
```
Type `help test_wkm` at the MATLAB command prompt for more information.
The Symbolic Math Toolbox (to use variable precision arithmetic) is assumed to 
be present for the tests.
The raw output in the files `results_test_wkm.txt` and
`results_multiplier_test_wkm.txt` were generated using MATLAB 2017b.
The MATLAB script `plot_error.m` contains the processed data and the code used 
to generate the figures in the
[MIMS EPrint 2018.4](http://eprints.maths.manchester.ac.uk/2621/).
The source code of the MATLAB function `funm_condest1.m` is included in
`test_wkm.m` to compute the condition number estimates of matrix functions.
`funm_condest1.m` is a part of Higham's
[Matrix Function Toolbox](http://www.maths.manchester.ac.uk/~higham/mftoolbox).

#### Reference
P. Nadukandi and N. J. Higham, "Computing the Wave-Kernel Matrix Functions",
*Manchester Institute for the Mathematical Sciences*,
[MIMS EPrint 2018.4](http://eprints.maths.manchester.ac.uk/2621/), 2018.

#### Funding
Nadukandi was supported by an individual fellowship from the European Union's 
Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie
grant [702138](https://cordis.europa.eu/project/rcn/200435_en.html).  
Higham was supported by Engineering and Physical Sciences Research Council 
grant [EP/P020720/1](http://gow.epsrc.ac.uk/NGBOViewGrant.aspx?GrantRef=EP/P020720/1).

## Clone this repository ###
```
git clone https://github.com/nadukandi/wkm.git
cd wkm
```
