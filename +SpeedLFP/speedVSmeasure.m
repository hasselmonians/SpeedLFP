function [meanFeature, sigmaFeature, speedDim, ft, tao, fet, speed_us, hVar] = speedVSmeasure(root,varargin)
%
% NOTE: root.active_lfp and root.epoch must be set before use!
%
% Inputs:
%   root: singleton or vector of CMBHOME.Session objects. If a vector, then
%   takes running speed & frequency samples from each session and bins them
%   all together.
%
%
% Parameters: 
%   active_lfp: Which lfp to use? (root.active_lfp)
%   markerColor: 3-element row vector of color to use for plotting ([0 0 1])
%   ifPlot: To plot, or just return? (1)
%   figureVar: Which figure to use? (new figure)
%   plotFit: Plot the fit? (1)
%   plotErr: Plot error bars? (1)
%   speedRange: What speeds (cm/s) to incorporate? ([5 30])
%   speedBinSize: 
%
% Outputs:
%   meanFeature: Mean theta frequency in each bin 
%   sigmaFeature: Standard Deviation in each bin
%   speedDim: Edges of the speed bins used
%   y_m: Linear fit positions
%   R: Linear fit R value
%   P: Linear fit P value
%   b: Linear fit slope
%   y_int: Linear fit y intercept
% 
% [meanFeature, sigmaFeature, speedDim, ft, tao, fet, speed_us, hVar = speedVSmeasure(root, varargin)

%% Parse & setup
p = inputParser;
p.addParamValue('markerColor', [0 0 1], (@(c) ischar(c) || @(c) numel(c) ==3 )); % Can specify as 'b' or [0 0 1]
p.addParamValue('ifPlot', 1, @(c) (c==1 || c==0));
p.addParamValue('figureVar', [], @isnumeric);
p.addParamValue('plotFit', 1, @(c) (c==1 || c==0));
p.addParamValue('plotErr', 1, @(c) (c==1 || c==0));
p.addParamValue('speedDim', (5:1:30)');
p.addParamValue('feature', 'frequency', @(x) ischar(x)); 
p.addParamValue('freqBand', [6 12], @(c) numel(c) == 2)

p.parse(varargin{:});

markerColor = p.Results.markerColor;
ifPlot = p.Results.ifPlot;
figureVar = p.Results.figureVar;
speedDim = p.Results.speedDim;
plotFit = p.Results.plotFit;
plotErr = p.Results.plotErr;
freqBand = p.Results.freqBand;
feature = p.Results.feature;
speedDim = speedDim(:);

%% Get the Variables
% feature from LFP
switch lower(feature)
    case 'frequency'
        % Calculate instantaneous frequency, based on Hilbert transform
        phase = CMBHOME.LFP.InstPhase(CMBHOME.LFP.BandpassFilter(root.b_lfp(root.active_lfp).signal, root.b_lfp(root.active_lfp).fs, freqBand));
        [~, iv] = findpeaks(phase);
        freq = root.lfp.fs ./ diff(iv);
        freq = [NaN;freq];
        freq = interp1(root.b_lfp(root.active_lfp).ts(iv),freq,root.b_lfp(root.active_lfp).ts);
        freq(freq >  freqBand(2)*1.5) = NaN;       % problem cycles discounted
        root.b_lfp(root.active_lfp).b_myvar = freq;
        fet = CMBHOME.Utils.ContinuizeEpochs(root.lfp.myvar);
    case 'amplitude'
        % Instantaneous amplitude from hilbert transform
        root.b_lfp(root.active_lfp).b_myvar = ...
            CMBHOME.LFP.InstAmplitude(CMBHOME.LFP.BandpassFilter(root.b_lfp(root.active_lfp).signal, root.b_lfp(root.active_lfp).fs, freqBand));
        fet = CMBHOME.Utils.ContinuizeEpochs(root.lfp.myvar);
        
    otherwise
        error('Unknown feature')
end

% Upsample the speed:
speed_us = interp1(root.b_ts, root.b_vel, root.b_lfp(root.active_lfp).ts);
root.b_lfp(root.active_lfp).b_myvar = speed_us * root.spatial_scale;
speed_us = CMBHOME.Utils.ContinuizeEpochs(root.lfp.myvar);


%% SpeedMod
[counts, wb] = histc(speed_us,speedDim);
meanFeature = NaN(length(counts),1);
sigmaFeature = NaN(length(counts),1); %standard error in each bin
tao = sigmaFeature;

for i = 1:length(counts)-1
   meanFeature(i) = nanmean(fet(wb==i)); 
   sigmaFeature(i) = nanstd(fet(wb==i)) / sqrt(nansum((wb==i)));
   tao(i) = nansum((wb==i)) / nanmean(root.lfp.fs);
end

bads = isnan(meanFeature) | isnan(speedDim);
meanFeature(bads) = [];
sigmaFeature(bads) = [];
tao(bads) = [];
speedDim(bads) = [];

%% Fit
[ft.y_m, ft.R, ft.P, ft.b, ft.y_int] = CMBHOME.Utils.LinearRegression(speedDim,meanFeature);


%% Plot

if ifPlot
    
    if isempty(figureVar), figureVar = figure;end
    
        hVar(1) = plot(speedDim,meanFeature,'o','MarkerSize',3,'MarkerFaceColor',markerColor,'MarkerEdgeColor',markerColor);            
        hold on
        
        if plotFit == 1
            hVar(2) = plot(speedDim,ft.y_m,'Color',markerColor);
        end
        
        if plotErr == 1
            t = errorbar(speedDim,meanFeature,sigmaFeature,'.');
            set(t,'Color',markerColor);
            hVar(3) = t;
        end

        xlim([speedDim(1) speedDim(end)])
        xlabel('Running Speed (cm/s)');
        ylabel(['Mean Band ' feature]);
end

end