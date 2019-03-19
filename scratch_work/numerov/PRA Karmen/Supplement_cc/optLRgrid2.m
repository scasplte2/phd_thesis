function Rgrid = optLRgrid2(Rgrid,n,drmax,Rend)
%tkarman 9/11 2012.
%
% function Rgrid = optLRgrid2(Rgrid,n,drmax,Rend)
%
%"opt grid", the spacing doubles every n steps,
%and there is a maximum dr.

a=2^(1/n);
dr=Rgrid(2)-Rgrid(1);
while(Rgrid(end)<Rend)
dr=dr*a;
dr=min(dr,drmax);
Rgrid=[Rgrid Rgrid(end)+dr];
end

end
