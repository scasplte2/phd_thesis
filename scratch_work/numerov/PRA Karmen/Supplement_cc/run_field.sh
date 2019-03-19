#!/bin/csh -f
# $1 lmax $2 M $3 parity $4 permutation symmetry $5 magnetic field.
# field included.

set job = $0:r
set out = $job'_lmax'$1'_M'$2'_par'$3'_per'$4'_B'$5'_Emin'$6'_Emax'$7'_Esteps'$8

mkdir -p $job

matlab -singleCompThread >& ./$job/$out.out << EOF

lmax=$1;
M=$2;
parity=$3;
permsym=$4;
B=1e$5;	%in kelvin

[totalElastic,totalInelastic,Egrid,B]=ScSc_field(lmax,M,parity,permsym,B,1e$6,1e$7,$8);
gauss=4.254382549059841e-10;

save ./$job/$out.mat totalElastic totalInelastic Egrid B gauss lmax M parity permsym

EOF

