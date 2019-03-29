The Hermite DVR is suited for an expansion in eigenfunctions of the harmonic
oscillator (Hermite polynomials times an exponential). In contrast to other DVR
methods, the grid points are not natively bounded to some interval, but can lie
anywhere (of course, they are always located somewhere around the minimum of the
harmonic oscillator).

Furthermore, the Gauss-Hermite quadrature is conveniently defined in scaled
coordinates (corresponding to m * omega = 1, r_e = 0).  To use them for real
problems, you have to supply the shift and the properties of the harmonic
potential you want to use.


The following attributes can be set:

mass
    the mass of the particle. Usually not used independently, but needed
    together with omega to scale the points we get from Gaussian quadrature.

omega
    the angular frequency of the potential. The force constant is related to
    this by k = m * omega^2. Needed for scaling of the quadrature points.

r_e
    the equilibrium position of the harmonic oscillator. Needed for shifting the
    quadrature points.

n_pts
    the number of points (== number of basis functions)


Further internally used attributes are:

momentum
    the grid representation of the momentum operator in DVR(!) basis

kin
    the grid representation of the kinetic energy Hamiltonian
    in DVR basis.

kinexpo
    the grid representation of the kinetic energy short-time
    propagator in DVR basis.

trafo_expand
    The transformation matrix for the DVR=>FBR transform

trafo_reconstruct
    The transformation matrix for the FBR=>DVR transform
