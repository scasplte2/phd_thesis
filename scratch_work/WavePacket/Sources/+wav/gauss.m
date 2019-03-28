%------------------------------------------------------------------------------
%
% This function creates a coherent superposition of complex-valued Gaussian
% wavepackets of given width, centered at r_0 or k_0 in position or momentum
% space. Wavefunction not yet normalized. Note that it comes in two flavors.
% If the input DOF is not supplied, we assume correlated description,
% otherwise only one degree of freedom is set up.
%
% The function parameter tells us which degree of freedom we are supposed to
% fill up. The resulting wave function is returned by this function.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% see the README file for license details.

function gauss = gauss(dir)

global psi space

%% Correlated Gaussian
if nargin == 0
    % Number of Gaussian wavepackets (number of rows in parameter matrix)
    n = size ( psi.init.corr.pos_0, 1 );

    % Set default values
    if ~isfield(psi.init.corr, 'mom_0')
        psi.init.corr.mom_0 = zeros(size(psi.init.corr.mom_0));
    end
    if ~isfield(psi.init.corr, 'factor')
        psi.init.corr.factor = ones(n, 1);
    end

    % Output
    util.disp (' ')
    util.disp ('*******************************************************')
    util.disp ( ['Coherent superposition of ' int2str(n)] )
    util.disp ( 'Gaussian minimum uncertainty wave packet(s)          ' )
    util.disp ( '                                                     ' )
    util.disp ( '           N     DOF      [                    ( Rk - R0gk )^2 ] ' )
    util.disp ( 'psi(R) =  Sum C  Prd  exp [ I K0gk*(Rk-R0gk) - (-----------)   ] ' )
    util.disp ( '          g=1  g k=1      [                    (   2*W_gk  )   ] ' )
    util.disp ( '                                                     ' )
    util.disp ( ' Note that this wavepacket corresponds to a          ' )
    util.disp ( ' coherent (Glauber) state of a harmonic oscillator   ' )
    util.disp ( ' V(R) = k/2 (R-R0)^2 with mass m for the following   ' )
    util.disp ( ' choice of the width parameter W = (4*k*m)^(-1/4)    ' )
    util.disp ( '                                                     ' )
    util.disp ('*******************************************************')


    psi.dvr.grid_ND{1} = zeros(size(space.dvr.grid_ND{1}));

    % Loop over Gaussian packets and sum them up
    for g=1:n
        if n>1
            util.disp (' ')
            util.disp ( [ 'Gaussian packet #' int2str(g)])
        end
        util.disp ( [ 'Mean value of position(s) R0   : ' num2str(     psi.init.corr.pos_0(g,:)) ] )
        util.disp ( [ 'Mean value of momentum(s) K0   : ' num2str(     psi.init.corr.mom_0(g,:)) ] )
        util.disp ( [ 'Prefactor                 C    : ' num2str(     psi.init.corr.factor(g))  ] )
        util.disp ( [ 'Position uncertainties    W    : ' num2str(     psi.init.corr.width(g,:)) ] )
        util.disp ( [ 'Momentum uncertainties 1/(2*W) : ' num2str(0.5./psi.init.corr.width(g,:)) ] )

        gauss = ones(size(space.dvr.grid_ND{1}));

        for k = 1:space.size.n_dim
            gauss = gauss .* exp (  1i*(space.dvr.grid_ND{k}-psi.init.corr.pos_0(g,k)) *  psi.init.corr.mom_0(g,k)...
                                     -((space.dvr.grid_ND{k}-psi.init.corr.pos_0(g,k)) / (psi.init.corr.width(g,k)*2)).^2  );
        end

        psi.dvr.grid_ND{1} = psi.dvr.grid_ND{1} + psi.init.corr.factor(g)*gauss;
    end

%% Only along one dimension
else
    % Number of Gaussian wavepackets
    n = length ( psi.init.dof{dir}.pos_0 );

    % Set default values
    if ~isfield ( psi.init.dof{dir}, 'mom_0' )
         psi.init.dof{dir}.mom_0 = zeros(n,1);
    end
    if ~isfield(psi.init.dof{dir}, 'factor')
        psi.init.dof{dir}.factor = ones(n,1);
    end

    util.disp (' ')
    util.disp ( '****************************************************' )
    util.disp ( ['Initial wavefunction for DOF :' int2str(dir)] )
    if n>1
        util.disp ( ['Coherent superposition of ' int2str(n)] )
    end
    util.disp ( 'Gaussian minimum uncertainty wave packet(s)         ' )
    util.disp ( '                                                    ' )
    util.disp ( '           N         [                 ( R - R0g )^2 ] ' )
    util.disp ( 'psi(R) =  Sum C  exp [ I K0g*(R-R0g) - (---------)   ] ' )
    util.disp ( '          g=1  g     [                 (  2*Wg   )   ] ' )
    util.disp ( '                                                    ' )
    util.disp ( ' Note that this wavepacket corresponds to a         ' )
    util.disp ( ' coherent (Glauber) state of a harmonic oscillator  ' )
    util.disp ( ' V(R) = k/2 (R-R0)^2 with mass m for the following  ' )
    util.disp ( ' choice of the width parameter W = (4*k*m)^(-1/4)   ' )
    util.disp ( '                                                    ' )
    util.disp ( '****************************************************' )


    gauss = zeros(size(space.dvr.grid_ND{1}));

    % Loop over Gaussian packets and sum them up
    for g=1:n
        if n>1
            util.disp (' ')
            util.disp ( [ 'Gaussian packet #' int2str(g)])
        end
        util.disp ( [ 'Mean value position       R0 : ' num2str(     psi.init.dof{dir}.pos_0(g)) ] )
        util.disp ( [ 'Mean value of momentum    K0 : ' num2str(     psi.init.dof{dir}.mom_0(g)) ] )
        util.disp ( [ 'Prefactor                 C  : ' num2str(     psi.init.dof{dir}.factor(g))] )
        util.disp ( [ 'Position uncertainty      W  : ' num2str(     psi.init.dof{dir}.width(g)) ] )
        util.disp ( [ 'Momentum uncertainty 1/(2*W) : ' num2str(0.5./psi.init.dof{dir}.width(g)) ] )

        gauss = gauss + psi.init.dof{dir}.factor(g) * ...
            exp (  1i*(space.dvr.grid_ND{dir}-psi.init.dof{dir}.pos_0(g)) *  psi.init.dof{dir}.mom_0(g)...
            -((space.dvr.grid_ND{dir}-psi.init.dof{dir}.pos_0(g)) / (psi.init.dof{dir}.width(g)*2)).^2  );
    end
end
