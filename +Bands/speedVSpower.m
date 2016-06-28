function [meanPower, speedDim, y_m, R, P, b, y_int, tao, hVar] = speedVSpower(root,varargin)
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

%% Get the signal:
sig = root.lfp.signal;

%% SpeedMod
[counts, wb] = histc(speed_us,speedDim);
meanPower = NaN(length(counts),1);
tao = meanPower;

for i = 1:length(counts)-1
   if sum(wb==i) > 0
       [meanPower(i), ~,S{i}] = spectrum(sig(wb==i),root(1).lfp.fs,10,9,[0 120], freqBand);
       tao(i) = nansum((wb==i)) / nanmean(root.lfp.fs);
   else
       meanPower(i) = NaN;
       S{i} = NaN;
       tao(i) = NaN;
   end
end

bads = isnan(meanPower) | isnan(speedDim);
meanPower(bads) = [];
tao(bads) = [];
speedDim(bads) = [];

%% Fit
[y_m, R, P, b, y_int] = CMBHOME.Utils.LinearRegression(speedDim,meanPower);


%% Plot

if ifPlot
    
    if isempty(figureVar), figureVar = figure;end
    
        hVar(1) = plot(speedDim,meanPower,'o','MarkerSize',3,'MarkerFaceColor',markerColor,'MarkerEdgeColor',markerColor);            
        hold on
        
        if plotFit == 1
            hVar(2) = plot(speedDim,y_m,'Color',markerColor);
        end
        
        xlim([speedDim(1) speedDim(end)])
        xlabel('Running Speed (cm/s)');
        ylabel('Band Power');
end

end


function [totalSpectralPower, f,S,signalPower] = spectrum(signal, Fs, TBP, n_tapers, f_range, band)
% root.plot_lfp_spectrogram(lfp_ind, windowsize, windowsin, bandwidth,
% f_range);
%
% Plots a spectrogram in jet colormap for root.b_lfp(lfp_ind). Params for
% the moving FFT analysis above.
%
% Uses the Chronux toolbox.
%
% andrew bogaard 3 april 2010

    import CMBHOME.Utils.*
    
    %AddChronuxPackage;
    CMBHOME.Utils.AddChronuxPackage;
    if ~exist('TBP', 'var')
        TBP = 10;
    end
    
    if ~exist('n_tapers', 'var')
        n_tapers = 9;
    end
    
    if iscell(signal), signal = vertcat(signal{:}); end
        
    params.tapers = [TBP, n_tapers];
    params.fpass = f_range;
    params.Fs = Fs;
    params.pad = 2;

    [S,f]=mtspectrumc(signal,params);

    stdf = .2/mean(diff(f)); % elements in a .2 hz std
    
    kernel = pdf('Normal', -3*stdf:3*stdf, 0, stdf); % smooth with gaussian kernel
    
    S = conv(S, kernel, 'same'); % plot(f, S)
    
    [xmax,imax] = extrema(S); % find local maxima and see if any lie within theta

    ftmp = f(imax);
    stmp = S(imax);

    [a_peak,ind] = max(stmp(ftmp > band(1) & ftmp < band(2))); % intrinsic frequency

    ftmp = ftmp(ftmp > band(1) & ftmp < band(2));

    if ~isempty(ind)                                    % if theta peak was found
        f_peak = ftmp(ind)';     
        a_peak_av = mean(S(f>f_peak-1 & f<f_peak+1));
        peak = S(f==f_peak);
    else                                                % if theta peak was not found
        f_peak = 0;
        peak = 0;
        a_peak_av = 0;
    end

    power_ratio = a_peak_av / mean(S);
        
    
    a_peak_av = mean(S(f>f_peak-1 & f<f_peak+1));

    inds = f>band(1) & f<band(2);
    totalSpectralPower = sum(S(inds));
    
    signal_power = mean(signal.^2);

    if a_peak_av<1.5*signal_power
        f_peak = [];
    end

end
