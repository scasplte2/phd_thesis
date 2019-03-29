The fft DVR accepts several additional attributes:

mass
    the mass thet enters the kinetic energy (a.u. as everything).
    The kinetic energy operator is assumed to be
    -1/(2*mass) d^2/dx^2

n_pts
    the number of grid points (equal in position and momentum grid)

x_min
    the lower bound of the position grid

x_max
    the upper bound of the position grid

periodic
    if set to true, assume periodic boundary conditions for the Hamiltonian,
    otherwise use kinetic energy corrections for non-periodicity. Only used
    for the TISE.


Internally used further attributes are:

p_min
    the lower bound of the momentum grid

p_max
    the upper bound of the momentum grid

kin
    the matrix of the kinetic energy 
    (conveniently a matrix in the momentum grid multiplied elementwise to
    perform the action of the kinetic energy operator)

intern_kin
    an fftshifted version of kin. Is never handed out.

intern_kinexpo
    like intern_kin, but exponentiated, so that it can be used with the split
    operator method.
