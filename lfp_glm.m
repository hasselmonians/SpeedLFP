function res = lfp_glm(vehicle, drug)

% Version History:
% 0.01 - wchapman 20160121: Intial writing (no code, just theory)

%% LFP as a General Linear Model
% #general=/=generalized
%
% for the sake of sanity, we assume no interaction of drug*exposure. Also,
% we work on the assumption that Well's normalization is okay to do
% beforehand ... for now. 

% drug = [0,1]         (logistic)
% exposure = [0,1,2,3] (class)
% speed = \real        (covariate)

% fr = drug + exposure + speed + (speed*drug) + (speed*exposure)

speedDim = drug.crunched.rat(1).bands.freq.speedDim; 

freq_drug = arrayfun(@(x) x.bands.freq.meanSig, drug.crunched.rat,'UniformOutput',0);
freq_vehi = arrayfun(@(x) x.bands.freq.meanSig, vehicle.crunched.rat,'UniformOutput',0);

trial_vehi = repmat([1 2 3 4],size(freq_vehi,2),1);
trial_drug = repmat([1 2 3 4],size(freq_drug,2),1);

manip_vehi = repmat([0 0 0 0], size(freq_vehi,2),1);
manip_drug = repmat([0 1 1 1],size(freq_drug,2),1);

% Create data structure

drugdat = [cell2mat(arrayfun(@(x) [1*ones(size(speedDim)) speedDim x{1}], freq_drug(:,1),'UniformOutput',0));
cell2mat(arrayfun(@(x) [2*ones(size(speedDim)) speedDim x{1}], freq_drug(:,2),'UniformOutput',0));
cell2mat(arrayfun(@(x) [3*ones(size(speedDim)) speedDim x{1}], freq_drug(:,3),'UniformOutput',0));
cell2mat(arrayfun(@(x) [4*ones(size(speedDim)) speedDim x{1}], freq_drug(:,4),'UniformOutput',0))];
drugdat = [zeros(size(drugdat,1),1), drugdat];
drugdat((size(drugdat,1)/4 : size(drugdat,1)/2),1) = 1;
drugdat = [drugdat drugdat(:,1).*drugdat(:,end) drugdat(:,2)];

vehidat = [cell2mat(arrayfun(@(x) [1*ones(size(speedDim)) speedDim x{1}], freq_vehi(:,1),'UniformOutput',0));
cell2mat(arrayfun(@(x) [2*ones(size(speedDim)) speedDim x{1}], freq_vehi(:,2),'UniformOutput',0));
cell2mat(arrayfun(@(x) [3*ones(size(speedDim)) speedDim x{1}], freq_vehi(:,3),'UniformOutput',0));
cell2mat(arrayfun(@(x) [4*ones(size(speedDim)) speedDim x{1}], freq_vehi(:,4),'UniformOutput',0))];
vehidat = [zeros(size(vehidat,1),1), vehidat];
vehidat((size(vehidat,1)/4 : size(vehidat,1)/2),1) = 1;
vehidat = [vehidat vehidat(:,1).*vehidat(:,end) vehidat(:,2)];


dat = [drugdat; vehidat];

% [intercept, manip, trial, speed, mani*speed, trial*speed]


[res.b, res.dev, res.stats] = glmfit(dat(:,1:end-1),dat(:,end));

%yfit = glmval(b, x, 'probit', 'size', n);
%        plot(x, y./n, 'o', x, yfit./n, '-')
end