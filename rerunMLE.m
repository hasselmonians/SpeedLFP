function rerunMLE()


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
            disp(['[' num2str(i) ', ' num2str(k) ']'])
            try
                load(SL(k).fname);
                root.cel = SL(k).cel;
                stats_all = runMLE(root);
                
                stats_all.stats0.x = [];
                stats_all.stats0.x_ = [];
                stats_all.stats0.inds = [];

                stats_all.stats0f.TS = [];
                stats_all.stats0f.YAXIS = [];
                stats_all.stats0f.stats0 = [];
                stats_all.stats0f.x = [];
                stats_all.stats0f.x_ = [];
                stats_all.stats0f.y = [];
                stats_all.stats0f.y_ = [];

                stats_all.statsA.stats0 = [];
                stats_all.statsA.x = [];
                stats_all.statsA.x_ = [];
                stats_all.statsA.y = [];
                stats_all.statsA.y_ =[];

                stats_all.statsf.stats0 = [];
                stats_all.statsf.x = [];
                stats_all.statsf.x_ = [];
                stats_all.statsf.y = [];
                stats_all.statsf.y_ = [];
                
                SL(k).rhythmicity.stats_all = stats_all;
                
            catch
                disp(['[' num2str(i) ', ' num2str(k) ']: BAD'])
            end
            close all
        end 
        save(['/media/wchapman/RatBrains/Dropbox/UnitRecordingData/Caitlin''s Data/Analyses/SLs/' SLNames{i}])
    end
end