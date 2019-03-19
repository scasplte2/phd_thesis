#!/bin/csh -f
# $1 lmax $2 M $3 parity $4 permutation symmetry $5 magnetic field.
# field included.

set job = $0:r
set out = $job'_lmax'$1'_M'$2'_par'$3'_per'$4'_B'$5'_Emin'$6'_Emax'$7'_Esteps'$8

mkdir -p $job

matlab -nosplash -nojvm -singleCompThread >& ./$job/$out.out << EOF

generate_18j

EOF

