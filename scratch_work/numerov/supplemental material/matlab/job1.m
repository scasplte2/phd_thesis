% Constant potential test Qsin and Qairy should be exact.
format long
mu    = 2/3;

npoint=5;
Rgrid=0:2*pi/(npoint-1):2*pi;
npoint=length(Rgrid);
fprintf('Grid ranges from %f to %f and contains %d points.\n',Rgrid(1),Rgrid(npoint),npoint);

pars.k2=-5;

tic();
Qnnum = numerov(Rgrid, 0, mu, @wmat_const,pars);
toc()
[nchan,~] = size(Qnnum);
Ls    = zeros(nchan, 1);
tic();
Qnsin{1}=zeros(nchan,nchan);
Qnsin = Qsin(Rgrid,0,mu,Qnsin,@wmat_const,pars);
toc()
tic();
Qnair{1}=zeros(nchan,nchan);
Qnair = Qairy(Rgrid,0,mu,Qnair,@wmat_const,pars);
toc()

%match numerov
[S1,Kopen] = matchBoundaryConditions(Rgrid,Ls,mu,Qnnum,0,0);
I     = eye(size(Kopen));
Snum     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);

%match sine
[S1,Kopen] = matchBoundaryConditions(Rgrid,Ls,mu,Qnsin{1},0,0);
I     = eye(size(Kopen));
Ssin     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);

%match airy
[S1,Kopen] = matchBoundaryConditions(Rgrid,Ls,mu,Qnair{1},0,0);
I     = eye(size(Kopen));
Sair     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);

Qnnum
Qnsin=Qnsin{1}
Qnair=Qnair{1}
if(pars.k2<0)
	Qexact=sin(sqrt(abs(pars.k2))*Rgrid(npoint-1))/sin(sqrt(abs(pars.k2))*Rgrid(npoint))
else
    Qexact=sinh(sqrt(abs(pars.k2))*Rgrid(npoint-1))/sinh(sqrt(abs(pars.k2))*Rgrid(npoint))
end




Rspan = linspace(0,10,500);
[Qexact,psi] = deal(zeros(1,length(Rspan)));
for i = 1:length(Rspan)
    Rgrid=Rspan(i):2*pi/(npoint-1):(2*pi+Rspan(i));
    Qexact(i)=sin(sqrt(abs(pars.k2))*Rgrid(npoint-1))/sin(sqrt(abs(pars.k2))*Rgrid(npoint));
    
    if i > 1
        psi(i) = psi(i - 1)/Qexact(i);
    else
        psi(i) = 1e-6;
    end
end
plot(Rspan,psi)