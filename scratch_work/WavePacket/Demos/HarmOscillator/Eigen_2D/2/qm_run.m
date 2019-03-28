% Copyright (C) 2008,2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.


% qm_propa does the real calculation
qm_setup();
qm_init();
qm_propa();

% Now plot the autocorrelation minus 1 and save it to a file
global time

plot(abs(time.acf.grid).^2 - 1);
saveas(gcf, 'acf.jpg');

qm_cleanup();
