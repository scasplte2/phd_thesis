% Copyright (C) 2009-2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

% Run the propagation and calculate the ratio OD/OH by hand.

qm_setup();
state=wave;
qm_init(state);
qm_propa(state);

% The expectation value is set such that it gives the 
% propability of forming OH [+ D] in the simple picture,
% i.e., R_{OD} > R_{OH}. What we want, however, is the
% ratio OD / OH. OD + OH gives the norm Ne of the excited
% surface (absorption can be neglected). Simple algebra
% gives then 
% OD/OH = (Ne - OH) / OH
global expect time

ratio = expect.amo{2}.cha{2} ./ expect.amo{1}.cha{2};

figure(3);
plot(time.steps.m_grid, ratio, 'LineWidth', 2);
xlabel('t (a.u.)');
ylabel('OD/OH');
saveas(gcf, 'ratio.jpg');

qm_cleanup();
