function ctabs=generate_ctabs
%function ctabs=generate_ctabs
% generates labels for body fixed effective potential in ctabs

ctabs=[];

for(S=0:1)
for(lambda=-4:4)
for(c=1:5-abs(lambda))

ctabs=[ctabs; S c lambda];

end
end
end

end
