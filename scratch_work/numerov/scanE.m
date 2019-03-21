numPnts = 5;
[x, R] = deal(cell(1,numPnts));
[S_el, del_el] = deal(zeros(1,numPnts));

T = logspace(-1,-8,5);
for i = 1:length(T)
    tic;
    [x{i}, ~, R{i}, ~, S_el(i), del_el(i)] = sr_1S0_scatter(T(i), 0, [5 100e3], [84 84]);
    toc
end