function res = nova_mag(vehicle, drug)
% Perform a multivariate codependent nova explosion.
% [intercept slope] = f(group) 

%% anocova

mag = cell2mat(arrayfun(@(x) x.bands.mag.meanSig, drug.crunched.rat,'UniformOutput',0));
mag_d = mag(:,1) - mag(:,2);
mag = cell2mat(arrayfun(@(x) x.bands.mag.meanSig, vehicle.crunched.rat,'UniformOutput',0));
mag_v = mag(:,1) - mag(:,2);

mag = [mag_d; mag_v];
drug_bin = [ones(size(mag_d)) ; zeros(size(mag_v))];
rs = 5:2.5:27.5;
rs = repmat(rs(:), (length(mag_d)+length(mag_v))/10,1); 


%[~, res(1).atab, res(1).ctab, res(1).stats] = aoctool(rs, mag, drug_bin,[],[],[],[],'off',1);
%[~, res(2).atab, res(2).ctab, res(2).stats] = aoctool(rs, mag, drug_bin,[],[],[],[],'off',2);
%[~, res(3).atab, res(3).ctab, res(3).stats] = aoctool(rs, mag, drug_bin,[],[],[],[],'off',3);
%[~, res(4).atab, res(4).ctab, res(4).stats] = aoctool(rs, mag, drug_bin,[],[],[],[],'off',4);
[~, res.atab, res.ctab, res.stats] = aoctool(rs, mag, drug_bin,[],[],[],[],'off',5);


end