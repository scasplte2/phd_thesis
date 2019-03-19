========================= README.TXT ===============================

ARTICLE INFORMATION:

Journal:    Phys. Rev. A

Authors:   Tijs Karman and Gerrit C. Groenenboom

Title: Cold magnetically-trapped ^2D scandium atoms: II Scattering dynamics


DEPOSIT INFORMATION

Description:
Routines for evaluating scattering cross sections for collisions of
arbitrary state Russel-Saunders coupled atoms due to magnetic dipole-dipole
coupling in the Born approximation.

Total No. of Files: 5

Filenames: README.TXT + the 4 files:

	born_field.m	Contains a function to evaluate scattering cross sections with magnetic fields present
	born_nofield.m	Contains a script to evaluate scattering cross sections for the field-free case
	ff_3jm.F*	Fortran function for evaluating 3-jm symbols
	ff_3jm.mexa64*	mex file for interfacing matlab and the above fortran function

Filetypes: all ASCII.
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

born_nofield.m is a script that can be executed directly in matlab.
	constants set some natural constants and the reduced mass in mu
	at present, this corresponds to the case of Sc-Sc, but could be altered.
	parameters are set under the comment %input:
		LMAX is the highest partial wave included
		ja and jb are the total electronic angular momenta
		LA and SA are the orbital and spin angular momenta of the atom (L=2 and S=1/2 for Sc-Sc)
		ma and mb are the initial projections
		parity determined (-1)^l
		al, as and gj calculate the lande g factor from LS coupling.
			if this does not apply (for heavy elements) you could simply set gj to the experimental value, for instance.

	final states (m_A and m_B quantum numbers) are contained in a list called mps.
		This is generated to include all possible values, except for the initial state.
		Alternatively, you could set it mps=[ma mb]; in order to calculate elastic cross sections (odd parity only)

	Then, the actual calculation is performed,
	and the resulting cross section is returned in the variable sigma.

born_field. contains a function that calculates the cross section including magnetic fields for a range of field strengths and energies.
	input is a list of energies (Egrid) and magnetic field strength (Bgrid). Input in atomic units.
	Parameters for Sc-Sc are hard-coded as in born_nofield.m,
	and can be adjusted for other atom-atom systems or different final states as discussed above.

	Function returns a matrix with cross sections.
	The first index loop over the energies, the second over the magnetic field strength.


