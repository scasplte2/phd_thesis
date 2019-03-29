The Legendre DVR is used for problems where the minor quantum number (most often
called m) is conserved while the major quantum number (typically named l) can
vary. Note that for homonuclear diatomic molecules, symmetry arguments require
that l is always either even or odd. If your potential does not fulfill this
automatically (linearly polarized lasers do), you have to add the required
symmetry by hand (not yet implemented).

The Legendre DVR considers the quantum number l as momentum.

The Legendre DVR accepts several additional attributes:

mass
    The mass thet enters the kinetic energy (a.u. as everything).
    The kinetic energy operator is assumed to be
	L^2/(2mR^2)

R_dof
	Index of the degree of freedom that is used for the R in the kinetic energy
	operator L^2/(2mR^2). May not be set together with R_0 (raises an error)

R_0
	for a rigid rotator the fixed value of R that enters the kinetic energy operator.
	May not be set together with R_dof

l_max
	The maximum angular momentum of the spectral basis, which is also the number
	of grid points in position space.

m_0
	the (constant) value of the minor quantum number. The spectral basis goes from
    l = m_0 to (and including) l_max, therefore it has (l_max - m_0 + 1) elements.


Internally used further attributes are:

kin
    the matrix of the kinetic energy. Multiply elementwise with the wave
    function in the spectral basis to imitate the action of the kinetic energy
    operator. Note that there are some pecularities here, because we commonly
    cast the problem internally to a pseudo-2D form, so this should not be used
    from outside.

kinexpo
    like kin, but exponentiated, so that it can be used with the split
    operator method.

trafo_expand
    The matrix that you multiply a (reshaped) object (wave function) in
    pseudospectral representation with to get the spectral basis representation.

trafo_reconstruct
	The matrix that you multiply an object in spectral representation with
	to get the pseudospectral representation.
