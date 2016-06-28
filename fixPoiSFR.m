SLNames = {'SL_8-oh-dpat_infusion_output',...    %1
           'SL_8-oh-dpat_systemic_output',...    %2
           'SL_diazepam_infusion_output',...     %3
           'SL_diazepam_systemic_output',...     %4
           'SL_DMSO_systemic_output',...         %5
           'SL_PBS_infusion_output',...          %6
           'SL_PBS_systemic_output'};            %7
       
for i = 1:length(SLNames)
    load(SLNames{i});
    
    for k = 1:numel(SL)
        if isstruct(SL(k).poi_sfr)
            if isnan(SL(k).poi_sfr.all.pF_lin)
                
                load(SL(k).fname);
                root = root.FixPos;
                root.b_vel = [];
                root = root.AppendKalmanVel;
                root.cel = SL(k).cel;
                root.epoch = [-inf inf];

                [root, success] = CMBHOME.Utils.RunningEpochs(root, 1, 30, 0.5);

                cts = CMBHOME.Utils.ContinuizeEpochs(root.cel_ts) - root.b_ts(1);
                ts = CMBHOME.Utils.ContinuizeEpochs(root.ts) - root.b_ts(1);
                sp = CMBHOME.Utils.ContinuizeEpochs(root.vel) * root.spatial_scale;
                cps = CMBHOME.Utils.ContinuizeEpochs(root.cel_vel) * root.spatial_scale;
                
                SL(k).poi_sfr.all = ana.fit_speed(cts, root.fs_video, ts, sp, cps, 0);
                
            end
        end        
    end
    
    save(['/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/SLs/' SLNames{i} '.mat'], 'SL')

end
