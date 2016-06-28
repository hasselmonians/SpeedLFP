function ctps = assign_cells(SL, thresh)
% ctps = assign_cells(SL, threshStruct)
% 
% Creates an assignment of cell types, following the thresholds passed in
% through threshStruct. threshStruct can be generated by shuffle_grids
    
    if ~exist('thresh','var')
        thresh.grid3 = 0.34;
        thresh.mrl = .3;
        thresh.p_rhythm = 0.05;
        thresh.si = 2;
    else
        thresh.p_rhythm = 0.05;
    end
    
    thresh.fr = 10;

    for i = 1:size(SL,1)
        try
            rhythmic(i) = SL(i).rhythmicity.stats_all.stats0.p_rhyth;
        catch
            rhythmic(i) = 1;
        end
    end
    rhythmic = rhythmic(:);
    
    gridness = arrayfun(@(x) x.gridness3, SL(:,1));
    si = arrayfun(@(x) x.si, SL(:,1));
    dropout = arrayfun(@(x) isnan(x.cel(1)), SL(:,3));
    mrl = arrayfun(@(x) x.mrl, SL(:,1));
    fr = arrayfun(@(x) x.f, SL(:,1));
    
    ctps{1,1} = 'all';
    ctps{1,2} = ~dropout & (fr<=thresh.fr);
    
    ctps{2,1} = 'grids';
    ctps{2,2} = (gridness>=thresh.grid3) & (fr<=thresh.fr) & ~dropout;
    
    ctps{3,1} = 'conjunctive';
    ctps{3,2} = (gridness>=thresh.grid3) & (mrl>=thresh.mrl) & (fr<=thresh.fr) &~dropout;

    ctps{4,1} = 'nonconjunctive';
    ctps{4,2} = (gridness>=thresh.grid3) & ~(mrl>=thresh.mrl) & (fr<=thresh.fr) &~dropout;
    
    ctps{5,1} = 'hd';
    ctps{5,2} = (mrl>=thresh.mrl) & ~(gridness>=thresh.grid3) &(fr<=thresh.fr) & ~dropout;

    ctps{6,1} = 'spatial';
    ctps{6,2} = (si>=thresh.si) & ~(mrl>=thresh.mrl) & ~(gridness>=thresh.grid3) &(fr<=thresh.fr) & ~dropout;

    ctps{7,1} = 'theta';
    ctps{7,2} = (rhythmic<thresh.p_rhythm) & ~(si>=thresh.si) & ...
                ~(mrl>=thresh.mrl) & ~(gridness>=thresh.grid3) & (fr<=thresh.fr) &~dropout;
            
    ctps{8,1} = 'interneuron';
    ctps{8,2} = ~dropout & ~(fr<=thresh.fr);

end