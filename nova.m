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


%[~, res(1).atab, res(1).ctab, res(1).stats] = aoctool(rs, freq, drug_bin,[],[],[],[],'off',1);
%[~, res(2).atab, res(2).ctab, res(2).stats] = aoctool(rs, freq, drug_bin,[],[],[],[],'off',2);
%[~, res(3).atab, res(3).ctab, res(3).stats] = aoctool(rs, freq, drug_bin,[],[],[],[],'off',3);
%[~, res(4).atab, res(4).ctab, res(4).stats] = aoctool(rs, freq, drug_bin,[],[],[],[],'off',4);

[~, res.atab, res.ctab, res.stats] = aoctool(rs, freq, drug_bin,[],[],[],[],'off',5);

% atab --> field --> p that different
% ctab --> field --> p that non-zero
%                    p that non-zero for group

end