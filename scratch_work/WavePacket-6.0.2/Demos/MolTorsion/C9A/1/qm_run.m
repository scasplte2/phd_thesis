% Copyright (C) 2009,2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

% Our final goal is to simulate the absorption spectrum of C9A.
% This calculation here only serves for calculating the electronic
% and vibrational groundstate that we can lift to the excited state
% (/calculate overlap with excited states) later on.

% Run the propagation and extract the nearly degenerate two vibrational states
% in a manner suitable for later use. First the state with gerade symmetry,
% then the other one
qm_setup();
state=wave;
qm_init(state,'g');
qm_bound(state);
eigen(state, 1);
dlmwrite('state_g.dat', state.dvr{1}, 'precision', 10, 'delimiter', ' ');
qm_cleanup

qm_setup
state=wave;
qm_init(state,'u');
qm_bound(state);
eigen(state, 1);
dlmwrite('state_u.dat', state.dvr{1}, 'precision', 10, 'delimiter', ' ');
qm_cleanup();
