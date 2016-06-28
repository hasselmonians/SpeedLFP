function [meanMag, sigmaMag, speedDim, y_m, R, P, b, y_int, tao, hVar] = speedVSmagnitude(root,varargin)
% [meanFreq, sigmaFreq, speedDim, y_m, R, P, b, y_int ] = speedVSfrequency(root, varargin)
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
%   meanFreq: Mean theta frequency in each bin 
%   sigmaFreq: Standard Deviation in each bin
%   speedDim: Edges of the speed bins used
%   y_m: Linear fit positions
%   R: Linear fit R value
%   P: Linear fit P value
%   b: Linear fit slope
%   y_int: Linear fit y intercept

% wchapman 20130821

%% Parse & setup
p = inputParser;
p.addParamValue('active_lfp', [], @isnumerc)
p.addParamValue('markerColor', [0 0 1], (@(c) ischar(c) || @(c) numel(c) ==3 )); % Can specify as 'b' or [0 0 1]
p.addParamValue('ifPlot', 1, @(c) (c==1 || c==0));
p.addParamValue('figureVar', [], @isnumeric);
p.addParamValue('plotFit', 1, @(c) (c==1 || c==0));
p.addParamValue('plotErr', 1, @(c) (c==1 || c==0));
p.addParamValue('speedDim', (5:1:30)');
p.addParamValue('freqBand', [6 12], @(c) numel(c) == 2)

p.parse(varargin{:});

markerColor = p.Results.markerColor;
ifPlot = p.Results.ifPlot;
figureVar = p.Results.figureVar;
speedDim = p.Results.speedDim;
plotFit = p.Results.plotFit;
plotErr = p.Results.plotErr;
freqBand = p.Results.freqBand;
speedDim = speedDim(:);
%% Upsample the speed:
speed_us = interp1(root.ts, root.vel, root.lfp.ts);
speed_us = speed_us * root.spatial_scale;

%% Extract Magnitude
mag = CMBHOME.LFP.InstAmplitude(CMBHOME.LFP.BandpassFilter(root.lfp.signal, root.lfp.fs, freqBand));

%% SpeedMod
[counts, wb] = histc(speed_us,speedDim);
meanMag = NaN(length(counts),1);
sigmaMag = NaN(length(counts),1); %standard error in each bin
tao = sigmaMag;

for i = 1:length(counts)-1
   meanMag(i) = nanmean(mag(wb==i)); 
   sigmaMag(i) = nanstd(mag(wb==i)) / sqrt(nansum((wb==i)));
   tao(i) = nansum((wb==i)) / nanmean(root.lfp.fs);
end

bads = isnan(meanMag) | isnan(speedDim);
meanMag(bads) = [];
sigmaMag(bads) = [];
tao(bads) = [];
speedDim(bads) = [];

%% Fit
[y_m, R, P, b, y_int] = CMBHOME.Utils.LinearRegression(speedDim,meanMag);


%% Plot

if ifPlot
    
    if isempty(figureVar), figureVar = figure;end
    
        hVar(1) = plot(speedDim,meanMag,'o','MarkerSize',3,'MarkerFaceColor',markerColor,'MarkerEdgeColor',markerColor);            
        hold on
        
        if plotFit == 1
            hVar(2) = plot(speedDim,y_m,'Color',markerColor);
        end
        
        if plotErr == 1
            t = errorbar(speedDim,meanMag,sigmaMag,'.');
            set(t,'Color',markerColor);
            hVar(3) = t;
        end

        xlim([speedDim(1) speedDim(end)])
        xlabel('Running Speed (cm/s)');
        ylabel('Mean Band Amplitude (mV)');
end

end