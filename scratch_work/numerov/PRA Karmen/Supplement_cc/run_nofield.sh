#!/bin/csh -f
# $1 lmax $2 parity $3 permutation symmetry.
# fieldfree.

set job = $0:r
set out = $job'_B0_lmax'$1'_par'$2'_per'$3'_Emin'$4'_Emax'$5'_Esteps'$6'_Rmax'$7

mkdir -p $job

matlab -nosplash -nojvm -singleCompThread >& ./$job/$out.out << EOF

Rmax=$7;
[totalElastic,totalInelastic,Egrid,lmax,permsym,parity,Ms]=ScSc_nofield($1,$2,$3,1e$4,1e$5,$6,Rmax);

save ./$job/$out.mat totalElastic totalInelastic Egrid parity permsym lmax Ms Rmax

EOF
