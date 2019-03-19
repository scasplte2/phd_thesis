function c_g=c_gm(J);
%C_GM      Clebsch-Gordan coefficient
%     C_GM( [j1 j2 j3 ; m1 m2 m3 ] ) returns <j1 m1 j2 m2 | j3 m3>

fac=(-1)^(J(1,1)-J(1,2)+J(2,3))*sqrt(2*J(1,3)+1);
J(2,3)=-J(2,3);
c_g=fac*ff_3jm(J);
