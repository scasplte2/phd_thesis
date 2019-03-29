function pseudo (S,T)
% Check whether T is the pseudo-inverse of S 
% where 3 special cases may or may not arise:
% If S has linearly independent    rows, then S*S' is invertible and S*T=I
% If S has linearly independent columns, then S'*S is invertible and T*S=I
% If S has linearly independent rows AND columns, then S is invertible 
% and the pseudoinverse is identical to the inverse, T = S^{-1} 
s = size (S);
t = size (T);

if s(1)~=t(2)
    log.error ('Wrong dimensions s1,t2 in testing pseudo-inverse')
end
if s(2)~=t(1)
    log.error ('Wrong dimensions s2,t1 in testing pseudo-inverse')
end

I1=eye(s(1));
I2=eye(s(2));
STI =   S*T - I1;
TSI =   T*S - I2;
STS = S*T*S - S;
TST = T*S*T - T;
log.disp ('Transformation matrices S and T (pseudo-inverse)')
log.disp (['Dimension of S : ' int2str(s(1)) ' x ' int2str(s(2))])
log.disp (['Dimension of T : ' int2str(t(1)) ' x ' int2str(t(2))])
log.disp (['||( ST-I)|| / ||I|| = ' num2str(norm(STI,'fro')/norm(I1,'fro'))])
log.disp (['||( TS-I)|| / ||I|| = ' num2str(norm(TSI,'fro')/norm(I2,'fro'))])
log.disp (['||(STS-S)|| / ||S|| = ' num2str(norm(STS,'fro')/norm(S, 'fro'))])
log.disp (['||(TST-T)|| / ||T|| = ' num2str(norm(TST,'fro')/norm(T, 'fro'))])
log.disp ('   ')

figure (42)
thisplot = vis.styles;
show_logo (thisplot)
az=-20; el=10;
%
subplot(2,2,1)
mesh(real(S*T))
view(az,el)
axis ([1 s(1) 1 s(1) -inf inf])
title('real ( S T )')
set ( gca, ...
    'LineWidth',     thisplot.l_thick, ...
    'FontName',      thisplot.f_name,  ...
    'FontSize',      thisplot.f_large, ...
    'FontWeight',    thisplot.f_heavy )

%
subplot(2,2,2)
mesh(imag(S*T))
view(az,el)
axis ([1 s(1) 1 s(1) -inf inf])
title('imag ( S T )')
set ( gca, ...
    'LineWidth',     thisplot.l_thick, ...
    'FontName',      thisplot.f_name,  ...
    'FontSize',      thisplot.f_large, ...
    'FontWeight',    thisplot.f_heavy )

%
subplot(2,2,3)
mesh(real(T*S))
view(az,el)
axis ([1 t(1) 1 t(1) -inf inf])
title('real ( T S )')
set ( gca, ...
    'LineWidth',     thisplot.l_thick, ...
    'FontName',      thisplot.f_name,  ...
    'FontSize',      thisplot.f_large, ...
    'FontWeight',    thisplot.f_heavy )

%
subplot(2,2,4)
mesh(imag(T*S))
view(az,el)
axis ([1 t(1) 1 t(1) -inf inf])
title('imag ( T S )')
set ( gca, ...
    'LineWidth',     thisplot.l_thick, ...
    'FontName',      thisplot.f_name,  ...
    'FontSize',      thisplot.f_large, ...
    'FontWeight',    thisplot.f_heavy )


