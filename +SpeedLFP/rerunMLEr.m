function rerunMLEr()

    warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
    warning('off','stats:mle:IterLimit');
    
    SLNames = {'SL_8-oh-dpat_infusion_output',...    %1
               'SL_8-oh-dpat_systemic_output',...    %2
               'SL_diazepam_infusion_output',...     %3
               %'SL_diazepam_systemic_output',...     %4
               %'SL_DMSO_systemic_output',...         %5
               %'SL_PBS_infusion_output',...          %6
               %'SL_PBS_systemic_output'};            %7
               };

           
    parfor i = 1:length(SLNames)
        disp(['starting: ' num2str(i)])
        subfunc(SLNames{i})
    end
end


function subfunc(SLName)

    load(SLName);clc;
    for k = 1:numel(SL)
        try
            try
                temp = SL(k).rhythmicity.stats_all.stats0.phat(4);
                if temp < 4 || temp > 12
                    temp = 1;
                else
                    temp = 0;
                end
            catch
                temp = 1;
            end

            if temp==1
                disp(['starting: ' SLName ', ' num2str(k)])
                load(SL(k).fname);
                root.cel = SL(k).cel;
                root.b_vel =[];
                root = root.AppendKalmanVel;

                stats_all = runMLE(root);

                SL(k).rhythmicity.stats_all = stats_all;
            end  

            close all
        catch
            
        end
    end 
    try
        save(['/media/wchapman/RatBrains/Dropbox/UnitRecordingData/Caitlin''s Data/Analyses/SLs/' SLName])
    catch
        save(['/home/wchapman/Desktop/' SLName])
    end
end
