% Copyright (C) 2009,2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

% Run the relaxation, and copy the resulting relaxed wavefunction to the
% other run directories so that it can be used as input for the interpolating
% initial wave function
qm_setup();
state = wave;
qm_init(state);
qm_propa(state);

global space

output = cat(2, space.dvr{1}, state.dvr{1});
dlmwrite('../2/wav_1.dat', output, ' ');

qm_cleanup();
