function AppendTheta(SLName)
    %{
    for i = 1:length(SLNames)
        AppendTheta(SLNames{i});
    end
    %}
    load(SLName)

    %{
    for i = 1:numel(SL)
        try
            load(SL(i).fname);root.cel=SL(i).cel;root.active_lfp=SL(i).active_lfp;
            root.epoch=[-inf inf]; t=root.cel_thetaphase{1};
            SL(i).thetaphase = t;
        catch
            SL(i).thetaphase = NaN;
        end
        
        try
            SL(i).theta_mrl = CMBHOME.Utils.circ.circ_r(SL(i).thetaphase);
            SL(i).theta_mra = rad2deg(CMBHOME.Utils.circ.circ_mean(SL(i).thetaphase));
        catch
            SL(i).theta_mrl = NaN;
            SL(i).theta_mra = NaN;
        end
        
    end
    %}
    
    for i = 1:size(SL,1)
        for k = 1:4
            try
                [SL(i,k).theta_p, SL(i,k).theta_tbl] = CMBHOME.Utils.circ.circ_wwtest(SL(i,1).thetaphase, SL(i,k).thetaphase);
            catch
                SL(i,k).theta_p = NaN;
                SL(i,k).theta_tbl = {};
            end
        end
    end
    
    for i = 1:numel(SL)
        try
            [SL(i).theta_pr, SL(i).theta_zr] = CMBHOME.Utils.circ.circ_rtest(SL(i).thetaphase);
        catch
            SL(i).theta_pr = NaN;
            SL(i).theta_zr = NaN;
        end
    end
    
    save(['/media/wchapman/RatBrains/Dropbox/UnitRecordingData/Caitlin''s Data/Analyses/SLs/' SLName],'SL');
    
end

