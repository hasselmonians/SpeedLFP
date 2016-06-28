function SLe = AppendAll(SLname,rowvar)


%{
addpath(genpath('/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/CMBs'))
addpath(genpath('/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/SLs'))

addpath /projectnb/hasselmogrp/wchapman/CMBHOME_stable/
addpath(genpath('mle_rhythmicity'))
addpath(genpath('/projectnb/hasselmogrp/wchapman/mle_rhythmicity'))
addpath(genpath('/media/wchapman/RatBrains/Dropbox (hasselmonians)'' /Bill/Projects/SpeedModulation/'))
addpath(genpath('/home/wchapman/Downloads'))

addpath /media/wchapman/RatBrains/Dropbox' (hasselmonians)'/Bill/Projects/SpeedModulation/
load(SLName)
for i = 179:numel(SL)
try
        tic;AppendAll(SLName,i);
        disp([num2str(i) ' good'])
        toc
    catch
        disp([num2str(i) ' bad'])
    end
end
%}

    %% Setup
    c=load(SLname);
    SLe = c.SL(rowvar);
    clear c

    if isempty(SLe.cel)
        SLe.cel = [NaN NaN];
    end

    if ~isnan(SLe.cel(1))

        c = load(SLe.fname);
        root = c.root;
        root.cel = SLe.cel;
        if numel(root.b_lfp) > 1
            root.active_lfp = SLe.active_lfp;
        else
            root.active_lfp = 1;
        end

        root = root.FixDir;
        root = root.FixPos;
        root.b_vel = [];
        root = root.AppendKalmanVel;
        root.epoch = SLe.epoch;

        %% Behavioral information:
        % median running speed
        SLe.behavior.rs_median = median(root.vel);
        % distribution of running speed
        SLe.behavior.rs_bins = 0:root.spatial_scale^-1:prctile(root.vel,95);
        SLe.behavior.rs_vals = histc(root.vel, SLe.behavior.rs_bins);
        SLe.behavior.spatial_scale = root.spatial_scale;
        SLe.behavior.fs_video = root.fs_video;

        %% Basics:

        %nspk
        SLe.nspk = sum(cellfun(@(x) length(x), root.cel_i));

        %f
        SLe.f = SLe.nspk / sum(root.epoch(:,2) - root.epoch(:,1));

        %Autocorr
        [SLe.f_intrins, SLe.theta_index, ms] = root.IntrinsicFrequency2(root.cel,'theta_skipping',1);
        SLe.acorr_cnts = ms.cor;
        SLe.acorr_lags = ms.lag;
        SLe.thetaskip = ms.theta_skipping;

        %ratemap
        SLe.ratemap = root.RateMap(root.cel);

        %ratemapac
        SLe.ratemapac = root.AutoCorr(root.cel);

        % Occupancy:
        SLe.occupancy = root.Occupancy;

        %polar things
        [SLe.ratemap_polar, ...
         SLe.mra, SLe.mrl, SLe.polar_peak]=root.plot_polar_rate_map(root.cel);

        SLe.u2 = root.HDWatsonsU2(root.cel);

        %gridness1
        SLe.gridness1 = root.Gridness(root.cel);

        %gridness2
        SLe.gridness2 = root.Gridness2(root.cel);

        %gridness3
        SLe.gridness3 = root.Gridness(root.cel,'grid3',1);

        %si
        SLe.si = root.SpatialInformation(root.cel);

        % running speed:
        [SLe.speedMod.f, SLe.speedMod.v] = root.VelocityRate(root.cel);

         %% Phase precession:
         try
            SLe.phase_precession = root.plot_phase_precession();
         end

        %% Thetamod gamma:
        %[SLe.tmg.modindex, SLe.tmg.thetarange, SLe.tmg.gammarange, ...
        %    SLe.tmg.powPhsDists, SLe.tmg.bincenters, ...
        %    SLe.tmg.partial] = thetaModGamma(root);

        %% Rhythmicity analyses:
        try
            SLe.rhythmicity.stats_all = runMLE(root);
        end

        %% Non-binned speed modulation:

        %epoch it
        root.epoch = [-inf inf];
        ee = root.b_ts(find(root.b_vel*root.spatial_scale > 100, 1));

        cts = CMBHOME.Utils.ContinuizeEpochs(root.cel_ts) - root.ts(1);
        ts = root.ts - root.ts(1);
        sp = CMBHOME.Utils.ContinuizeEpochs(root.vel) * root.spatial_scale;
        cps = CMBHOME.Utils.ContinuizeEpochs(root.cel_vel) * root.spatial_scale;

        try
            SLe.poi_sfr.all = ana.fit_speed(cts, root.fs_video, ts, sp, cps, 0);
        end

        %% Incorporate runRoots things:
        [SLe.bands.freq.meanSig, SLe.bands.freq.sigmaSig,SLe.bands.freq.speedDim,~,~,~,~,SLe.bands.freq.intercepts, SLe.bands.freq.tao] = Bands.speedVSfrequency(root,'ifPlot',0,'speedDim', (5:2.5:30)');
        [SLe.bands.mag.meanSig, SLe.bands.mag.sigmaSig,SLe.bands.mag.speedDim,~,~,~,~,SLe.bands.mag.intercepts, SLe.bands.mag.tao] = Bands.speedVSmagnitude(root,'ifPlot',0,'speedDim', (5:2.5:30)');
        [SLe.bands.pow.meanSig, ~, ~, ~, ~, ~, SLe.bands.pow.intercepts, SLe.bands.pow.tao] = Bands.speedVSpower(root,'ifPlot',0,'speedDim', (5:2.5:30)');
    else

    end

    %% Cleanup:
    close all
    mkdir(['preprocess/' SLname])
    save(['preprocess/' SLname '/SLe' num2str(rowvar) '.mat'],'SLe');

end
