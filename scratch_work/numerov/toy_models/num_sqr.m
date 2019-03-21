function en = num_sqr
    
figure;
en = zeros(1,10); numPeaks = -1;
for j = 1:length(en);
    N = 1e3;
    Psi = zeros(1,N);

    g2 = 200;
    v = -1*ones(N);

    eps = -1 + (j-1)*(1/length(en));
    deps = 0.05;

    k2 = g2*(eps - v);
    Psi = wavefunction(Psi, k2, N);
    P1 = Psi(N-1);

    while abs(deps) > 1e-12;
        k2 = g2*(eps - v);
        Psi = wavefunction(Psi, k2, N);
        P2 = Psi(N-1);

        if P1*P2<0
            deps = -deps/2; 
        elseif P1*P2>0 && length(findpeaks(Psi)) ~= numPeaks
            en(j) = eps;
            numPeaks = length(findpeaks(Psi));
        end
        
        eps = eps + deps;
        P1 = P2;
    end
    
    plot(1:N,Psi); hold on
end

function Psi = wavefunction(Psi, k2, N)
    d = (1/(N-1));
    Psi(1) = 0;
    Psi(2) = 1e-4;

    for i = 3:N
        Psi(i) = ( (2 - (10/12)*d^2*k2(i-1))*Psi(i-1) - (1 + (d^2/12)*k2(i-2))*Psi(i-2) )/( 1 + (d^2/12)*k2(i) );
    end
end
end