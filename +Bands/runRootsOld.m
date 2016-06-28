function [freq, mag, pow, speedDim] = runRoots(roots, varargin)

%% Parse!

p = inputParser;
p.addParamValue('ifPlot',   1,  @(x) (x==0)||(x==1));
p.addParamValue('speedDim', 5:5:30, @(x) isvector(x));
p.addParamValue('freqBand', [6 10], @(c) numel(c) == 2);
p.addParamValue('markerColor', customMap(size(roots,2)), (@(c) (iscellstr(c) && numel(c)==size(roots,2)) || (size(c,1)==size(roots,2) && size(c,2)==3)));
p.addParamValue('figNums',  [NaN NaN NaN], @(x) numel(x)==3);
p.parse(varargin{:});

fn = fieldnames(p.Results);
for i = 1:length(fn)
    eval([fn{i} '=' 'p.Results.' fn{i} ';']);
end

%% Setup
if ~matlabpool('size') % open the pool, because this is slow.
    matlabpool open
end

numRats = size(roots,1);
numCond = size(roots,2);

%% frequency:

for i = 1:numel(roots)
    [meanSig{i}, sigmaSig{i},speedDim,~,~,~,~,~, tao{i}] = Bands.speedVSfrequency(roots(i),varargin{:},'ifPlot',0);
end

meanSig = reshape(meanSig,size(roots,1),size(roots,2));
sigmaSig = reshape(sigmaSig,size(roots,1),size(roots,2));
tao = reshape(tao,size(roots,1),size(roots,2));

[normFreqMat,normFreqMean, FreqMean] = normalize(meanSig, tao);

% Fitting:
R = NaN(numCond,1);
P = R;
b = R;
y_int = R;

for k = 1:numCond
    mu = nanmean(cell2mat(normFreqMat(:,k)')');
    sigma = nanstd(cell2mat(normFreqMat(:,k)')') / size(normFreqMat,1);
    
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

% Consolidte
freq.mu = mu_k;
freq.sigma = sigma_k;
freq.R = R;
freq.P = P;
freq.b = b;
freq.y_int = y_int;

%% magnitude:

for i = 1:numel(roots)
    [meanSig{i}, sigmaSig{i},speedDim,~,~,~,~,~, tao{i}] = Bands.speedVSmagnitude(roots(i),varargin{:},'ifPlot',0);
end

meanSig = reshape(meanSig,size(roots,1),size(roots,2));
sigmaSig = reshape(sigmaSig,size(roots,1),size(roots,2));
tao = reshape(tao,size(roots,1),size(roots,2));

[normFreqMat,normFreqMean, FreqMean] = normalize(meanSig, tao);

% Fitting:
R = NaN(numCond,1);
P = R;
b = R;
y_int = R;

for k = 1:numCond
    mu = nanmean(cell2mat(normFreqMat(:,k)')');
    sigma = nanstd(cell2mat(normFreqMat(:,k)')') / size(normFreqMat,1);
    
    f = mu;
    sd = speedDim;
    bads = isnan(f);
    %f = f(~bads); sd = sd(~bads);
    
    [y_mt, R(k), P(k), b(k), y_int(k)] = CMBHOME.Utils.LinearRegression(sd',f);
    
    bads = find(bads);
    
    %for i = 1:length(bads)
    %    y_mt = [y_mt(1:bads(i)-1) NaN y_mt(bads(i):end)];
    %end
    
    mu_k(:,k) = mu;
    sigma_k(:,k) = sigma;
    y_m(:,k) = y_mt(:);
end

% Consolidte
mag.mu = mu_k;
mag.sigma = sigma_k;
mag.R = R;
mag.P = P;
mag.b = b;
mag.y_int = y_int;

%% power:

for i = 1:numel(roots)
    [meanSig{i}, ~, ~, ~, ~, ~, ~, tao{i}] = Bands.speedVSpower(roots(i),varargin{:},'ifPlot',0);
end

meanSig = reshape(meanSig,size(roots,1),size(roots,2));
tao = reshape(tao,size(roots,1),size(roots,2));

[normFreqMat,normFreqMean, FreqMean] = normalize(meanSig, tao);

% Fitting:
R = NaN(numCond,1);
P = R;
b = R;
y_int = R;

for k = 1:numCond
    mu = nanmean(cell2mat(normFreqMat(:,k)')');
    sigma = nanstd(cell2mat(normFreqMat(:,k)')') / size(normFreqMat,1);
    
    f = mu;
    sd = speedDim;
    bads = isnan(f);
    %f = f(~bads); sd = sd(~bads);
    
    [y_mt, R(k), P(k), b(k), y_int(k)] = CMBHOME.Utils.LinearRegression(sd',f);
    
    bads = find(bads);
    
    %for i = 1:length(bads)
    %    y_mt = [y_mt(1:bads(i)-1) NaN y_mt(bads(i):end)];
    %end
    
    mu_k(:,k) = mu;
    sigma_k(:,k) = sigma;
    y_m(:,k) = y_mt(:);
end

% Consolidte
pow.mu = mu_k;
pow.sigma = sigma_k;
pow.R = R;
pow.P = P;
pow.b = b;
pow.y_int = y_int;

%% Plotting
if ifPlot == 1
   
   % Frequency:
   if isnan(figNums(1)); figNums(1) = figure;end
   figure(figNums(1));
   hold on
   for k = 1:4
       plot(speedDim,freq.mu(:,k),'o','MarkerSize',3,'MarkerFaceColor',markerColor(k,:),'MarkerEdgeColor',markerColor(k,:));
   end
   
   for k = 1:numCond
       tv{k} = ['Condition_' num2str(k)];
   end
   legend(tv)
   
   for k = 1:4
        h=errorbar(speedDim,freq.mu(:,k),freq.sigma(:,k),'.');
        set(h,'Color',markerColor(k,:));
   end
   
   for k = 1:4
       h = plot(speedDim,freq.y_int(k)+speedDim*freq.b(k),'Color',markerColor(k,:));
   end
   
   xlabel('Running Speed (cm/s)')
   ylabel('Band Frequency')
   
   % Magnitude:
   if isnan(figNums(2)); figNums(2) = figure;end
   figure(figNums(2));
   hold on
   for k = 1:4
       plot(speedDim,mag.mu(:,k),'o','MarkerSize',3,'MarkerFaceColor',markerColor(k,:),'MarkerEdgeColor',markerColor(k,:));
   end
   
   for k = 1:numCond
       tv{k} = ['Condition_' num2str(k)];
   end
   legend(tv)
   
   for k = 1:4
        h=errorbar(speedDim,mag.mu(:,k),mag.sigma(:,k),'.');
        set(h,'Color',markerColor(k,:));
   end
   
   for k = 1:4
       h = plot(speedDim,mag.y_int(k)+speedDim*mag.b(k),'Color',markerColor(k,:));
   end
   
   xlabel('Running Speed (cm/s)')
   ylabel('Band Magnitude')
   
   % Power:
   if isnan(figNums(3)); figNums(3) = figure;end
   figure(figNums(3));
   hold on
   for k = 1
       plot(speedDim,pow.mu(:,k),'o','MarkerSize',3,'MarkerFaceColor',markerColor(k,:),'MarkerEdgeColor',markerColor(k,:));
   end
   
   for k = 1:numCond
       tv{k} = ['Condition_' num2str(k)];
   end
   legend(tv)
   
   for k = 1:4
        h=errorbar(speedDim,pow.mu(:,k),pow.sigma(:,k),'.');
        set(h,'Color',markerColor(k,:));
   end
   
   for k = 1:4
       h = plot(speedDim,pow.y_int(k)+speedDim*pow.b(k),'Color',markerColor(k,:));
   end
   
   xlabel('Running Speed (cm/s)')
   ylabel('Band Power')
   
end

end

function [normExMat, normMeanVal, meanVal] = normalize(exMat, tao)
    % Normalization across each rat, as in Wells 2013
    normVal = NaN(size(exMat,1),1);
    for i = 1:size(exMat,1)
        normVal(i) = nansum(((tao{i,1}./nansum(tao{i,1})).*exMat{i,1}));
    end
    
    normVal = normVal - nanmean(normVal);
    
    %normExMat = NaN(size(exMat));
    normExMat = cell(size(exMat));
    for i = 1:size(exMat,1)
        for k = 1:size(exMat,2)
            meanVal(i,k) = nansum(((tao{i,k}./nansum(tao{i,k})).*exMat{i,k}));
            normMeanVal(i,k) = meanVal(i,k) - normVal(i);
            normExMat{i,k} = exMat{i,k} - normVal(i); 
        end
    end
end

function vals = customMap(num)

    vals = [0 0 0; ...
            0 0 1; ...
            0 1 0; ...
            1 0 0; ...
            0 1 1; ...
            1 0 0; ...
            1 0 1];
        
        vals = vals(1:num,:);
            
end