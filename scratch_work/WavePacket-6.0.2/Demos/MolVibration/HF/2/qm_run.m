% Copyright (C) 2009-2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

qm_setup();
state = wave;
qm_init(state);
qm_propa(state);
qm_cleanup();