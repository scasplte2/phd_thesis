========================= README.TXT ===============================

ARTICLE INFORMATION:

Journal:    J. Chem. Phys.

Authors:   Tijs Karman, Liebeth M. C. Janssen, and Gerrit C. Groenenboom

Title: A renormalized potential following propagation algorithm for
		solving the coupled-channels equations

DEPOSIT INFORMATION

Description:
MATLAB routines for the developed potential following algorithms.
Test calculations with constant and linear potentials,
as well as the Secrest and Johnson test problem.[J. Chem. Phys., 45, 4556 (1966)]

Total No. of Files: 29

Filenames: README.TXT + the 28 files:

	Qsin.m				Matlab functions for the propagators
	Qairy.m				--
	numerov.m			Renormalized Numerov
	seig.m				Diagonalization wrapper
	airexp1.m			Helper functions for the asymptotic expansion of the Airy functions
	airexp1e.m			--
	airexp1o.m 			--
        airexp2.m			--
        airexp2e.m			--
        airexp2o.m			--
	eenxk.m				--
	eenxkm.m			--
	Xairy.m				--
	Yairy.m				--
	matchBoundaryConditions.m       Contains function to match to scattering (S-matrix) boundary conditions.
	Vy.m				Helper function for Secrest and Johnson test problem.
	sinc_dvr.m			Sinc function discrete variable representation of Secrest and Johnson test problem.
	wmat_const.m			returns W-matrix for constant potential
	wmat_lin.m			returns W-matrix for linear potential
	wmat_sj.m			returns W-matrix for Secrest and Johnson problem
	test1.sh			test jobs with output
	test1.out			--
        test2.sh			--
        test2.out			--
        test3.sh			--
        test3.out			--
        test4.sh			--
        test4.out			--

File types: all ASCII.
Special Instructions: See Below
Contact Information:
        Tijs Karman
        Theoretical Chemistry, Institute of Molecules and Materials
        Radboud University Nijmegen
	Heyendaalseweg 135
	6525 AJ Nijmegen, The Netherlands 
         Phone: +31-(0)24-3653037 Email: t.karman@science.ru.nl
-----------------------
Special instructions:

(1) Usage, input and output of propagators

(2) Comments on implementation

(3) Running the test calculations


======================================================
(1) Usage, input and output of propagators
======================================================

numerov.m contains an implementation of the renormalized Numerov algorithm.
Usage is: Q=numerov(Rgrid,E,mu,wmat_fun,wmat_vars)
	Rgrid is a vector of radial grid points. This should be equidistant, or the results will be meaningless.
	E is scattering energy (scalar)
	mu is reduced mass (scalar)
	wmat_fun is the handle of a function that is called with arguments R, wmat_vars, 
		and should return the wmat, except for the term with the scattering energy, at separation R.
	wmat_vars is passed to wmat_fun.

	Returns final Q matrix.

Qsin.m contains an implementation of the renormalized constant reference potential algorithm.
Usage is: Q=Qsin(Rgrid,E,mu,Q,wmat_fun,wmat_vars)
	Rgrid is the radial grid.
	E is the scattering energy. This can be a vector of energies, for which the transformation to the adiabatic is performed once.
	mu is the reduced mass.
	Q is a list, where the ith entry corresponds to the first Q matrix at the ith energy. This Q matrix relates the solution in the first two radial grid points supplied.
	wmat_fun is the handle of a function that is called with arguments R, wmat_vars,
		and should return the wmat, except for the term with the scattering energy, at separation R.
	wmat_vars is passed to wmat_fun.

	Returns Q is the list of Q matrices relating the solutions in the last two grid points.

Qairy.m contains an implemenation of the renormalized linear reference potential algorithm.
Usage is exactly the same as for Qsin: Q=Qairy(Rgrid,E,mu,Q,wmat_fun,wmat_vars)

======================================================
(2) Comments on implementation
======================================================

The present implementation uses notation which is consistent with the Appendix of our paper.
The algorithm loops over the three-point intervals, labeled with iR.
On each interval we determine the reference potential, by calling the function wmat_fun.
The transformation to the locally adiabatic basis is determined by diagonalizing a W-matrix.
This is done by calling seig, which is a wrapper function, that at present calls the MATLAB built-in
function eig(), but this can be replaced by, for example, a call to dsyev.
We transform all W-matrices with this transformation, and extract the diagonal
(the off-diagonal elements are ignored, which is the locally adiabatic approximation).

Then, we loop over all scattering energies. These shift the diagonal of the W-matrix.
We calculate \tile{A} and \tilde{C}, in terms of X^\pm and Y^\pm, as is done in the appendix.
Then we transform back to the primitive basis, and propagate Q using
Q_n+1 = -(A Q_n + 1 )^-1 C.

For the Airy functions, the evaluation of X^\pm and Y^\pm is written such that it is conceptually clear (hopefully),
but perhaps not so efficient.
Specifically, the code loops over the different adiabats, instead of calculating the asymptotic expansions,
and the factors X^\pm and Y^\pm in a vectorized fashion.
Moreover, some sums which are required in the asymptotic expansion repeatedly are computed multiple times,
and these sums are evaluated to 2 times the machine epsilon, which might be more than required.

======================================================
(3) Running the test calculations
======================================================

Included test jobs are named job1.sh through job4.sh
Example output is included in the corresponding job*.out

job1.sh         
Runs test calculation with constant potential, and compare to the exact solution, which is a (hyperbolic) sine.
This calculation should be exact for the constant and linear reference potential algorithms (Qsin and Qairy).
The Numerov is not exact, but if the step size is small compared to wave length of the solution, it should still be accurate.
If you take larger steps, the Numerov results become meaningless, whereas the Qsin and Qairy results are still exact.

You can change the radial grid by changing this line:
        Rgrid=0:2*pi/(npoint-1):2*pi;
One can take the grid non-equidistant, but then the Numerov results become meaningless.
One should take the first gridpoint equal to zero, or the exact result to compare with is wrong (or you should calculate this from a linear combination of sines and cosines, instead)
alternatively you can adjust npoint to choose a different number of equidistant grid points on the same range.

The potential can be changed by changing this line:
        pars.k2=-5;

job2.sh
Runs test calculation with linear potential, and compare to the exact solution, which is an Airy function.
For this calculation, only Qairy should be exact.

Again, the grid can be changed by changing npoint or the line
	Rgrid=0:2*pi/(npoint-1):2*pi;

The potential is changed by altering the lines
        pars.W0=-5;
        pars.W1=0.1;
Obviously, when choosing the slope W1=0, the Qsin algorithm is exact again.
The "exact" result is calculated numerically using the routines of Amos,
so this result may become meaningless for small slopes.

job3.sh
Runs test calculation for the Secrest and Johnson test problem.
Results compared to Numerov calculations with 10000 points.

job4.sh
Runs test calculation again for the Secrest and Johnson test problem,
now using Numerov in the short-range, and the Qsin and Qairy in the long-range.
Also, calculations are run for three energies. The diagonalization (in Qsin and Qairy) is not repeated for subsequent energies.


Tijs Karman,
Nijmegen,
May 2014

