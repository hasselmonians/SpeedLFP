function SL = ReRunBands_alt(SL)
%
% Does stuff

% Version History:
%{
wchapman 20150202: Initial writing. Need to rerun bands function on all SLs
                   to account for new selection criteria (see email from Caitlin).
                   Algorithm is:
                   1) Delta/Theta ration > 2 for atleast 5% of session
                   2) If multiple tetrodes with 1, then one with peak thetaness from spectrogram
%}


    for i = 1:size(SL,1)
        disp(i)
        c = load(SL(i,1).fname);
        root = c.root;
        clear pr
        for k = 1:length(root.b_lfp)
            pr(k) = calc_spectrum(root.b_lfp(k));
        end

        [~, ActiveInd] = max(pr);

        for k = 1:4
            try
                SL(i,k).active_lfp = ActiveInd;
                c = load(SL(i,k).fname);
                root = c.root;
                root.active_lfp = ActiveInd;

                root = root.FixDir;
                root = root.FixPos;
                root.b_vel = [];
                root = root.AppendKalmanVel;
                root.epoch = SL(i,k).epoch;

                [SL(i,k).bands2.freq.meanSig, SL(i,k).bands2.freq.sigmaSig,SL(i,k).bands2.freq.speedDim,~,~,~,~,SL(i,k).bands2.freq.intercepts, SL(i,k).band2s.freq.tao] = Bands.speedVSfrequency(root,'ifPlot',0,'speedDim', (5:2.5:30)','freqBand', [6 12]);
                [SL(i,k).bands2.mag.meanSig, SL(i,k).bands2.mag.sigmaSig,SL(i,k).bands2.mag.speedDim,~,~,~,~,SL(i,k).bands2.mag.intercepts, SL(i,k).bands2.mag.tao] = Bands.speedVSmagnitude(root,'ifPlot',0,'speedDim', (5:2.5:30)','freqBand', [6 12]);
                [SL(i,k).bands2.pow.meanSig, ~, ~, ~, ~, ~, SL(i,k).bands2.pow.intercepts, SL(i,k).bands2.pow.tao] = Bands.speedVSpower(root,'ifPlot',0,'speedDim', (5:2.5:30)','freqBand', [6 12]);
            end
        end

    end

end


function power_ratio = calc_spectrum(lfp)

    import CMBHOME.Utils.*
    AddChronuxPackage;


    TBP = 10;
    n_tapers = 9;
    signal = lfp.signal;

    f_low = 0; % Hz
    f_high = 20; % Hz
    f_range = [f_low, f_high];
    bandwidth = .25; %Hz
    windowsize = 10; % seconds
    windowinc = 2; % seconds
    TBP = 10;
    tapers = 9;

    if iscell(signal), signal = vertcat(signal{:}); end

    params.tapers = [TBP, n_tapers];
    params.fpass = f_range;
    params.Fs = lfp.fs;
    params.pad = 2;

    [S,f]=mtspectrumc(signal,params);

    stdf = .2/mean(diff(f)); % elements in a .2 hz std

    kernel = pdf('Normal', -3*stdf:3*stdf, 0, stdf); % smooth with gaussian kernel

    S = conv(S, kernel, 'same'); % plot(f, S)

    [xmax,imax] = extrema(S); % find local maxima and see if any lie within theta

    ftmp = f(imax);
    stmp = S(imax);

    [a_peak,ind] = max(stmp(ftmp > 7 & ftmp < 11)); % intrinsic frequency

    ftmp = ftmp(ftmp > 7 & ftmp < 11);

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


end



%{
close all; clear; clc;
addpath(genpath('/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/CMBs'))
addpath(genpath('/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/SLs'))

addpath /projectnb/hasselmogrp/wchapman/CMBHOME_stable/
addpath(genpath('mle_rhythmicity'))
addpath(genpath('/projectnb/hasselmogrp/wchapman/mle_rhythmicity'))
addpath(genpath('/media/wchapman/RatBrains/Dropbox (hasselmonians)'' /Bill/Projects/SpeedModulation/'))
addpath(genpath('/home/wchapman/Downloads'))

addpath /media/wchapman/RatBrains/Dropbox' (hasselmonians)'/Bill/Projects/SpeedModulation/


SLNames = {'SL_8-oh-dpat_infusion_output',...    %1
           'SL_8-oh-dpat_systemic_output',...    %2
           'SL_diazepam_infusion_output',...     %3
           'SL_diazepam_systemic_output',...     %4
           'SL_DMSO_systemic_output',...         %5
           'SL_PBS_infusion_output',...          %6
           'SL_PBS_systemic_output'};            %7

clc

for i = 1:length(SLNames)
    clc
    disp(i)
    disp('------------------------')
    load(SLNames{i});
    SL = ReRunBands_alt(SL);
    save(['/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/SLs/', SLNames{i} '.mat'], 'SL')

end

%}
