# MQDT
An implementation of Multichannel Quantum Defect Theory for atomic Rydberg and autoionization states.
This toolkit contains mqdt and mqdt-post to solve the MQDT equation and postprocess the results.
The mqdt-genEigen is a test version that utilizes the geneal eigenvalue solver for the continuum energy region to significantly speed up the solution of the determinant equation. 
[However, the implementation now relies on the lapack package, which is only limited to double precission, the eigen-wavefunction mixing coefficients can be problematic.]
## install
make
make -f Makefile-ifort

before running the make, you make need to modify the Makefile a bit. 
