%% Creating all the Wells et al plots
for slnum = 1:length(SLNames)
    dat(slnum).by_SL = by_SL(dat(slnum).SL, dat(slnum).crunched);
end
