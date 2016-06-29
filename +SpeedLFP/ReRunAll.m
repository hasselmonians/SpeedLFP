function SL = ReRunAll(SL)

    % Reruns all fields only for cells that have been added. AKA, looks for
    % SL.cel set, and empty SL.f

    for  i = 1:numel(SL)

        if ~isnan(SL(i).cel(1)) && isempty(SL(i).f)
            SLe = SL(i);
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

            %% Clean up
            SLe.tmg = [];

            try
                pp.rho = SLe.phase_precession.rho;
                pp.p = SLe.phase_precession.p;
                pp.s = SLe.phase_precession.s;
                pp.b = SLe.phase_precession.b;
                pp.is_precessing = SLe.phase_precession.is_precessing;
                SLe.phase_precession = pp;
            catch
                SLe.phase_precession = NaN;
            end

            SLe.rhythmicity.stats_all.stats0.x = [];
            SLe.rhythmicity.stats_all.stats0.x_ = [];
            SLe.rhythmicity.stats_all.stats0.inds = [];

            SLe.rhythmicity.stats_all.stats0f.TS = [];
            SLe.rhythmicity.stats_all.stats0f.YAXIS = [];
            SLe.rhythmicity.stats_all.stats0f.stats0 = [];
            SLe.rhythmicity.stats_all.stats0f.x = [];
            SLe.rhythmicity.stats_all.stats0f.x_ = [];
            SLe.rhythmicity.stats_all.stats0f.y = [];
            SLe.rhythmicity.stats_all.stats0f.y_ = [];

            SLe.rhythmicity.stats_all.statsA.stats0 = [];
            SLe.rhythmicity.stats_all.statsA.x = [];
            SLe.rhythmicity.stats_all.statsA.x_ = [];
            SLe.rhythmicity.stats_all.statsA.y = [];
            SLe.rhythmicity.stats_all.statsA.y_ =[];

            SLe.rhythmicity.stats_all.statsf.stats0 = [];
            SLe.rhythmicity.stats_all.statsf.x = [];
            SLe.rhythmicity.stats_all.statsf.x_ = [];
            SLe.rhythmicity.stats_all.statsf.y = [];
            SLe.rhythmicity.stats_all.statsf.y_ = [];

            SLe.poi_sfr.all.cnt = [];
            SLe.poi_sfr.all.speed = [];
            SLe.poi_sfr.all.ts = [];
            SLe.poi_sfr.all.stats.cnt = [];
            SLe.poi_sfr.all.stats.speed = [];
            SLe.poi_sfr.all.stats.spk_speed = [];
            SLe.poi_sfr.all.stats.spk_ts = [];
            SLe.poi_sfr.all.stats.c = [];
            SLe.poi_sfr.all.stats.fitobject = [];
            SLe.poi_sfr.all.stats.fitobject2 = [];

            SLe.poi_sfr.all.c = [];
            SLe.poi_sfr.all.fitobject = [];
            SLe.poi_sfr.all.fitobject2 = [];
            SLe.poi_sfr.all.spk_speed = [];
            SLe.poi_sfr.all.spk_ts = [];

            %% Output
            fn = fieldnames(SLe);
            for k = 1:length(fn)
                eval(['SL(i).' fn{k} '=SLe.' fn{k} ';'])
            end
            

        end

    end


end


%{
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
    SL = ReRunAll(SL);
    save(['/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/SLs/', SLNames{i} '.mat'], 'SL')
end
%}