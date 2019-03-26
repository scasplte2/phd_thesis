function renormalizedNumerovAdiabatic
% From
% http://www2.chem.umd.edu/groups/alexander/teaching/inelastic_scattering.pdf
% Avilable at
% http://www2.chem.umd.edu/groups/alexander/teaching/matlab_files/renormalized_numerov_adiab_final.m

% to carry out renormalized numerov calculation in the locally adiabatic basis
% for cl+h2
tocm1=219474.6;
% muconv converts from atomic mass units (C=12) to atomic units (m-electron=1)
muconv = 5.48579903e-04;
% reduced mass is 2.85 atomic mass units
mu = 2.85/muconv;
% spin-orbit constant A in cm-1
a = 293.3;

% integration parameters (N is the number of sectors)
rmin=3.5;
rmax=22;
N=200;
ind=0;
for iN=[1000:1000];
    ind=ind+1;
    NN(ind)=iN;
    
    N=iN;
    
    h=(rmax-rmin)/(N);
    rr=rmin:h:rmax;
    
    % add one more point, for final propagation to r_{N+1}
    rr=[rr rmax+h];
    
    
    disp(['renormalized numerov adiabatic basis, ' num2str(rmin) ' < R < ' num2str(rmax) ...
        ', N = ' int2str(N) ', h = ' num2str(h)])
    eetot=1000;
    lne=length(eetot);
    clear transprob
    
    % unit matrix
    unit=diag([1 1]);
    sq2=sqrt(2);
    for ie=1:lne
        etot=eetot(ie);
        if etot < 2*a
            disp(['etot = ' num2str(etot) ' .lt. 2*a = ' num2str(2*a) '; abort'])
            break
        end
        % debroglie wavelength in atomic units
        debrog=2*pi/sqrt(2*mu*etot/tocm1);
        disp(['Etot = ', num2str(etot) ' cm^{-1}; deBroglie wavelength = ', num2str(debrog)])
        tp=0d0;
        iturnpi=0;
        iturnsig=0;
        
        % calculate T matrix at r0
        vv=vsigpiclh2(rr(1));
        vsig=vv(1);
        vpi=vv(2);
        
        % check that we are starting inside the classical turning point
        if (vpi<etot)
            disp(['vpi = ' num2str(vpi) ' < etot = ' num2str(etot) '; ie = ' int2str(ie) ';  abort'])
            return
        end
        if (vsig<etot)
            disp(['vsig = ' num2str(vsig) ' < etot = ' num2str(etot) '; ie = ' int2str(ie) ';  abort'])
            return
        end
        
        % construct V matrix in asymptotic basis in atomic units
        % this is step a of the operational outline
        vmat=[ vsig*(2/3)+vpi*(1/3)-a   sq2*(vpi-vsig)/3
               sq2*(vpi-vsig)/3         vsig/3+2*vpi/3+2*a ]./tocm1;
        
        % diagonalize the V matrix
        [XMm1, vMm1]=eig(vmat);
        % this is step b of the operational outline
        wmat=2*mu*(unit*etot/tocm1-vmat);
        % construct (diagonal) t_tilde matrix in atomic units [Eq. 126]
        ttMm1=-h*h*2*mu*(etot/tocm1-diag(vMm1))/12;
        % this is tilde T0
        % construct tilde u0 [eq. 127]
        uMm1=(2+10*ttMm1)./(1-ttMm1);
        uMm1=diag(uMm1);
        ttMm1=diag(ttMm1);
        
        
        for M=1:N+1
            % M is the sector index
            % in this loop, we propagate the ratio matrix in the locally-adiabatic
            % basis from r_{M-1} to r_M
            
            % on entry we know the X and the \tilde u matrix at the left-hand side of the sector
            % (r=rM-1)
            
            % determine ratio matrix at the right hand side of the sector (step b
            % and b')
            if M==1
                % first step:  \tilde R_1 = \tilde u_0
                RtM=uMm1;
            else
                % all other steps, solve Eq. 128 for next ratio matrix (use linear eqn followed by
                % matrix multiply)
                % first linear equation, using Matlab's backlash function
                B=RtMm1\OMm2Mm1;
                % then matrix multiply
                C=OMm2Mm1'*B;
                RtM=uMm1-C;
            end
            % calculate V at the right hand side of this sector (r=r_M)
            % NOTE r_0 = rr(1), so that r_M=rr(M+1)
            vv=vsigpiclh2(rr(M+1));
            vsig=vv(1);vpi=vv(2);
            % construct V matrix in atomic units
            vmat=[vsig*(2/3)+vpi*(1/3)-a sq2*(vpi-vsig)/3;
                sq2*(vpi-vsig)/3 vsig/3+2*vpi/3+2*a]/tocm1;
            [XM vM]=eig(vmat);
            % construct overlap matrix O_{M-1,M} [Eq. 121]
            OMm1M=XMm1'*XM;
            % this is step b of the operational outline
            % construct (diagonal) t_tilde matrix in atomic units [Eq. 126]
            wmat=2*mu*(unit*etot/tocm1-vmat);
            % next line defines t_tilde at right hand side of sector
            ttM=-h*h*2*mu*(etot/tocm1-diag(vM))/12;
            % construct tilde u0 [eq. 127]
            uM=(2+10*ttM)./(1-ttM);
            uM=diag(uM);
            ttM=diag(ttM);
            % go back for new step
            if M<N+1
                ttMm2=ttMm1;
                XMm2=XMm1;
                ttMm1=ttM;
                XMm1=XM;
                RtMm1=RtM;
                uMm1=uM;
                OMm2Mm1=OMm1M;
            end
        end
        % at the end of the last step:
        % RtM contains R-tilde_{N+1}
        % RtMm1 contains R-tilde_N
        % XM contains X_{N+1}
        % XMm1 contains X_N
        % XMm2 contains X_{N-1}
        % ttM contain T_{N+1}
        % ttMm1 contains T_{N}
        % ttMm2 contains T_{N-1}
        % OMm1M contains O(N,N+1)
        % OMm2Mm1 contains O(N-1,N)
        
        % now determine log-derivative matrix [Eq. (131)]   
        YN=(XM*(0.5*unit-ttM)*inv(unit-ttM)*OMm1M'*RtM ...
            -XMm2*(0.5*unit-ttMm2)*inv(unit-ttMm2)*inv(RtMm1)*OMm2Mm1) ...
            *(unit-ttMm1)*XMm1'/h;
        disp(['log-derivative matrix at r = ',num2str(rr(N+1))]);
        disp(YN)
        % now match to boundary conditions
        % get wavevectors (this is sqrt of wmat)
        rN=rr(N+1);
        k=sqrt(-ttM*12/(h*h));
        % collect just the diagonal elements as a vector
        k=diag(k);
        i=sqrt(-1);
        % now expand to a full diagonal matrix
        h2=diag(i*exp(-k*i*rN)./sqrt(k));
        h2p=-i*diag(k).*h2;
        h1=conj(h2);
        h1p=i*diag(k).*h1;
        
        % solve for S matrix [Eq. 109]
        disp('S matrix')
        S=(h2p-YN*h2)\(h1p-YN*h1);
        disp(S)
        disp('transition probabilities')
        disp(S.*conj(S));
        transprob(ie)=abs(S(1,2))^2;
        
        % to check symmetry of S matrix
        disp([' N = ', int2str(N) '; [S(1,2)+S(2,1)]/2 = ', num2str((S(1,2)+S(2,1))/2) ...
            '; [S(1,2)-S(2,1)]/2 = ', num2str((S(1,2)-S(2,1))/2)])
        
        
        
    end
    tproba(ind)=transprob(ie);
end
end

function v = vsigpiclh2(r)

lam1 = [  0.813      0.677    ];
lam2 = [  1.2014     2.2061   ];
lam3 = [  5.5701     6.2093   ];
C1   = [  3.7457e3  -7.7718e3 ];
C2   = [  6.728e5   -7.3837e7 ];
C3   = [ -1.2773e5   3.2149e7 ];
C4   = [  3.4733e6   2.9743e6 ];

v = [ C1(1).*exp(-lam1(1).*r) + (C2(1) + C3(1).*r).*exp(-lam2(1).*r) - C4(1)/2.*(1 + tanh(1.2.*(r - lam3(1)))).*r.^(-6)
      C1(2).*exp(-lam1(2).*r) + (C2(2) + C3(2).*r).*exp(-lam2(2).*r) - C4(2)/2.*(1 + tanh(1.2.*(r - lam3(2)))).*r.^(-6) ];
end
