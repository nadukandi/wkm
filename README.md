# Computing the Wave-Kernel Matrix Functions

## About ###

This repository contains two MATLAB functions to compute (`wkm.m`) and test 
(`test\_wkm.m`) the Wave-Kernel Matrix Functions `cosh(sqrt(A))` and 
`sinhc(sqrt(A))`, where `A` is **any** `n x n` matrix, `sinhc(z) = sinh(z)/z` 
for any nonzero scalar `z`, and `sinhz(0) = 1`.

Function `wkm.m` has the following signature.
```
[C, S] = wkm (A)
```
where `C = cosh(sqrt(A))` and `S = sinhc(sqrt(A))`.
Type `help wkm` at the MATLAB command prompt for more information.

Function `test\_wkm.m` runs some tests for a suite of matrices.
```
test_wkm('mct\_testmats');
test_wkm('expm\_testmats');
test_wkm('logm\_testmats');
```
Type `help test\_wkm` at the MATLAB command prompt for more information.

#### Reference
P. Nadukandi and N. J. Higham, Computing the Wave-Kernel Matrix Functions,
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
