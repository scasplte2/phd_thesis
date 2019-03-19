function [y]=morse_eval(xgrid,alpha,a,b,c)
%geeft A+Bexp(-alphaR)+Cexp(-2alphaR) op het zelfde grid als xgrid..

z=exp(-alpha*xgrid);
y=a+b*z+c*z.^2;

end
