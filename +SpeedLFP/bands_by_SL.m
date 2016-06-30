function [res] = bands_by_SL(crunched)
% ???
%
% Runs and generates figures for a single SL (See Monaghan SFN 2015
% poster). As in line with this original analysis, this script does NOT
% account for vehicle data, etc. 
%
% AKA: Assumes no habituation etc.
%
% version history:
% 0.01 - ??? Taken from poster code
% 0.10 - wchapman 20160121: Converted to work with "crunched" input.


%% frequency: 
rat = crunched.rat;
bands = arrayfun(@(x) x.bands, crunched.rat);
freq = arrayfun(@(x) x.freq, bands);
meanSig = arrayfun(@(x) x.meanSig, freq,'UniformOutput',0);
tao = arrayfun(@(x) x.tao, freq,'UniformOutput',0);
speedDim = bands(1).freq.speedDim;

[normFreqMat,normFreqMean, FreqMean] = normalize(meanSig, tao); 

% FreqMean = average frequency for each animal for each session averaged across running speeds
% normFreqMean = average of baseline session values, normalizes other sessions to baseline
% normFreqMat = for each speed dim, the mean signal is normalized to
% mean baseline across rats
          

% testing for effects of drug / lights / gravity / netherVsNormal Realm / A vs B:
for i = 1:size(normFreqMat,1)
    for k = 1:size(normFreqMat,2)
        [~, ~, ~, bsig(i,k), yintsig(i,k)] = CMBHOME.Utils.LinearRegression(speedDim,normFreqMat{i,k});
    end
end
[h_b p_b] = ttest(bsig(:,1),bsig(:,2),0.05);
[h_yint p_yint] = ttest(yintsig(:,1),yintsig(:,2),0.05);

% Fitting:
numCond = size(rat,2);
R = NaN(numCond,1);
P = R;
b = R;
y_int = R;

for k = 1:numCond
    mu = nanmean(cell2mat(normFreqMat(:,k)')'); %mean freq for each speed bin for each condition (speed bin x condition matrix)
    sigma = nanstd(cell2mat(normFreqMat(:,k)')') / size(normFreqMat,1); %standard error for mu
    
    f = mu;
    sd = speedDim;
    bads = isnan(f);
    f = f(~bads); sd = sd(~bads);
    
    [y_mt, R(k), P(k), b(k), y_int(k)] = CMBHOME.Utils.LinearRegression(sd',f);
    
    bads = find(bads);
        
    %for i = 1:length(bads)
    %    y_mt = [y_mt(1:bads(i)-1) NaN y_mt(bads(i):end)];
    %end
    
    mu_k(:,k) = mu;
    sigma_k(:,k) = sigma;
    y_m(:,k) = y_mt(:);
end

% Consolidate
res.freq.mu = mu_k;
res.freq.sigma = sigma_k;
res.freq.R = R;
res.freq.P = P;
res.freq.b = b;
res.freq.y_int = y_int;
res.freq.p_b = p_b;
res.freq.p_yint = p_yint;

%%%%%%%%%%%%%% plot
markerColor = custom_map(4);
figure
hold on

numCond = numel(res.freq.b);
for k = 1:numCond
   tv{k} = ['Condition_' num2str(k)];
end
legend(tv,'Location','Best') %added location
hold on
for k = 1:4
    h=errorbar(speedDim,res.freq.mu(:,k),res.freq.sigma(:,k),'.');
    set(h,'Color',markerColor(k,:));
    set(h,'LineWidth',3);
    %h=errorbar(0, (freq.y_int),  nanstd(freq.y_int) / sqrt(length(freq.y_int)), '.');
    %set(h,'Color',markerColor(k,:));
    %set(h,'LineWidth',3);
end

for k = 1:4
   h = plot(speedDim,res.freq.y_int(k)+speedDim*res.freq.b(k),'Color',markerColor(k,:));
   h = plot([0 speedDim(1)], [res.freq.y_int(k) res.freq.y_int(k)+speedDim(1)*res.freq.b(k)], 'LineStyle', ':', 'Color',markerColor(k,:));
   h = plot([0], [res.freq.y_int(k)], 'square', 'Color',markerColor(k,:),'MarkerSize',13,'LineWidth',2);
   %set(h,'LineWidth',0.1);
   xl = get(gca,'xlim');
   xlim([-1 xl(2)])
end

for k = 1:4
   plot(speedDim,res.freq.mu(:,k),'square','MarkerSize',13,'MarkerFaceColor',markerColor(k,:),'MarkerEdgeColor',markerColor(k,:));
end

xlabel('Running Speed (cm/s)','FontSize',24)
ylabel('Frequency (Hz)','FontSize',24)
set(gca,'fontsize',18)
axis square



end


