function [threshs] = shuffle_threshs(SL)
    % [thresh] = shuffle_grid(SL);
    %
    % Shuffles each session, generates gridness threshold for the
    % population
    
    %http://www.nature.com/neuro/journal/v16/n3/pdf/nn.3311.pdf
    %{ 
    For each permutation, the entire sequence of spikes fired by the 
    cell was time shifted along the rat’s path by a random interval between
    a lower limit of 20 s and an upper limit of 20 s less than the length 
    of the trial, with the end of the trial wrapped to the beginning.
    %}
    
    %{
    load /media/wchapman/RatBrains/Dropbox' (hasselmonians)'/UnitRecordingData/Caitlin''''s' Data'/Analyses/SLs/SL_8-oh-dpat_infusion

    threshs = shuffle_threshs(SL);
    load /media/wchapman/RatBrains/Dropbox' (hasselmonians)'/UnitRecordingData/Caitlin''''s' Data'/Analyses/SLs/SL_diazepam_infusion
    threshs(2) = shuffle_threshs(SL);
    save('threshs.mat','threshs');
    %}
    
    addpath(genpath('/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/CMBs'))
    addpath /projectnb/hasselmogrp/wchapman/CMBHOME_stable/

    for i = 1:size(SL,1)
        
        disp(i/size(SL,1))
        
        load(SL(i).fname);
        root.cel = SL(i).cel;
        spk = root.spike(root.cel(1),root.cel(2));
        spk.ts = spk.ts-root.b_ts(1);
        root.b_ts = root.b_ts - root.b_ts(1);
        spkn = spk;
        
        for k = 1:10
            ts = spk.ts;
            
            a = 20; b = root.b_ts(end)-20;
            r = a + (b-a).*rand(length(ts),1);
            ts = mod(ts + r, root.b_ts(end));
            spkn = CMBHOME.Spike('ts',ts,'vid_ts',root.b_ts);
            
            root.spike(50,k) = spkn;
            root = root.AlignSpike2Session;
            root.cel = [50,k];
            
            gr(i,k) = root.Gridness(root.cel,'grid3',1);
            mrl(i,k) = CMBHOME.Utils.circ.circ_r(clean(deg2rad(root.cel_headdir{1})));
            si(i,k) = root.SpatialInformation(root.cel);
            
            
        end
        

    end

    threshs.grid3 = prctile(gr(:),95);
    threshs.mrl = prctile(mrl(:), 95);
    threshs.si = prctile(si(:), 95);

end