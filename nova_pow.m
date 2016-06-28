function res = nova_pow(vehicle, drug)
% Perform a multivariate codependent nova explosion.
% [intercept slope] = f(group) 

%% anocova

pow = cell2mat(arrayfun(@(x) x.bands.pow.meanSig, drug.crunched.rat,'UniformOutput',0));
pow_d = pow(:,1) - pow(:,2);
pow = cell2mat(arrayfun(@(x) x.bands.pow.meanSig, vehicle.crunched.rat,'UniformOutput',0));
pow_v = pow(:,1) - pow(:,2);

pow = [pow_d; pow_v];
drug_bin = [ones(size(pow_d)) ; zeros(size(pow_v))];
rs = 5:2.5:27.5;
rs = repmat(rs(:), (length(pow_d)+length(pow_v))/10,1); 


%[~, res(1).atab, res(1).ctab, res(1).stats] = aoctool(rs, pow, drug_bin,[],[],[],[],'off',1);
%[~, res(2).atab, res(2).ctab, res(2).stats] = aoctool(rs, pow, drug_bin,[],[],[],[],'off',2);
%[~, res(3).atab, res(3).ctab, res(3).stats] = aoctool(rs, pow, drug_bin,[],[],[],[],'off',3);
%[~, res(4).atab, res(4).ctab, res(4).stats] = aoctool(rs, pow, drug_bin,[],[],[],[],'off',4);
[~, res.atab, res.ctab, res.stats] = aoctool(rs, pow, drug_bin,[],[],[],[],'off',5);

end