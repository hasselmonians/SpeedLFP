function res = nova(vehicle, drug)
% Perform a multivariate codependent nova explosion.
% [intercept slope] = f(group) 

%% anocova

freq = cell2mat(arrayfun(@(x) x.bands.freq.meanSig, drug.crunched.rat,'UniformOutput',0));
freq_d = freq(:,1) - freq(:,2);
freq = cell2mat(arrayfun(@(x) x.bands.freq.meanSig, vehicle.crunched.rat,'UniformOutput',0));
freq_v = freq(:,1) - freq(:,2);

freq = [freq_d; freq_v];
drug_bin = [ones(size(freq_d)) ; zeros(size(freq_v))];
rs = 5:2.5:27.5;
rs = repmat(rs(:), (length(freq_d)+length(freq_v))/10,1); 


[~, res.atab, res.ctab, res.stats] = aoctool(rs, freq, drug_bin,[],[],[],[],'off',5);

% atab --> field --> p that different
% ctab --> field --> p that non-zero
%                    p that non-zero for group

end
