%% History
% version history:
% 0.01 - wchapman 20160120. Initial writing based on conversation
%
% 0.02 - wchapman 20160121. Broke in to multiple scripts. Below history is
%        overall flow information.
%      - Crunching -> by_SL is immediate concern
%
% 0.10 - wchapman 20160127. 
%      - Crunching complete
%      - Control Section complete: Switched to using ANOCOVA for LFP and
%        ANOVA for single unit measures. Outputs are located in
%        controlled.FIELDNAME
%
%      - Example plots for Caitlin in plots_SingleUnitBars and
%        plots_WellsStyle
%
%      - p order for controlled: [time, group, subject matched, groupTime]
%        ** aka, look at #4

%% Master Flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MASTER FLOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

0) ImplantRat
1) Collect some data
2) Arrange in to SL format
3) AppendAll.m on each SL
4) unpack.m -> save as "SL_..._output.mat"
5) Manually do "nearest_six"
5) shuffle_threshs.m
    - dpat_infusion, dpat_systemic, & ?
    - diaz_infusion, diaz_infusion, & ?
5) crunching.m 
6) by_SL.m
7) ???
8) Thesis

%}

%% Setup
%%{
cd /media/wchapman/RatBrains/Dropbox' (hasselmonians)'/UnitRecordingData/Caitlin''''s' Data'/Analyses/BillLivesHere/

addpath(genpath('/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/CMBs'))
addpath(genpath('/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/SLs'))
addpath /projectnb/hasselmogrp/wchapman/CMBHOME_stable/
addpath /media/wchapman/RatBrains/Dropbox/Bill/Projects/SpeedModulation/mle_rhythmicity
%startupJason
%}

clear all; close all; clc
%load output.mat
load threshs.mat

SLNames = {'SL_8-oh-dpat_infusion_output',...    %1
           'SL_8-oh-dpat_systemic_output',...    %2
           'SL_diazepam_infusion_output',...     %3
           'SL_diazepam_systemic_output',...     %4
           'SL_DMSO_systemic_output',...         %5
           'SL_PBS_infusion_output',...          %6
           'SL_PBS_systemic_output'};            %7

threshmat = [2 2 1 1 1 2 2];  % which set of thresholds to use
controlmat = [6 7 6 5];       % For the experiments, which control to use
%controlmat(controlmat==7) = 6; %tempoarry, currently rerunning things

%%
for slnum = 1:length(SLNames)
    load(SLNames{slnum});
    
    
    SL = NaNSL(SL);
    
    dat(slnum).SL = SL;
    dat(slnum).ctps = assign_cells(dat(slnum).SL, threshs(threshmat(slnum)));
    dat(slnum).crunched = crunching(dat(slnum).SL, dat(slnum).ctps);
    dat(slnum).by_SL = by_SL(dat(slnum).SL, dat(slnum).crunched);
    
    temp = arrayfun(@(x) x.medianspeed, dat(slnum).crunched.session);
    [h, p] = ttest(temp(:,2),temp(:,1));
    dat(slnum).p_median_speed = p; 
end

%arrayfun(@(x) x.p_median_speed, dat) < 0.05/7

%% Things needing control
for slnum = 1:length(controlmat)
    controlled(slnum).lfp_freq = nova(dat(controlmat(slnum)), dat(slnum));
    controlled(slnum).lfp_mag = nova_mag(dat(controlmat(slnum)), dat(slnum));
    controlled(slnum).lfp_pow = nova_pow(dat(controlmat(slnum)), dat(slnum));
    
    controlled(slnum).singleunit = anovas(dat(controlmat(slnum)), dat(slnum));
    controlled(slnum).singleunit_filtered = anovas_filtered(dat(controlmat(slnum)), dat(slnum));
    
end
% Single unit rules: compare must be significant, then can check drug/veh

controlled = controlled(:);


%-------- remove theta phase of spikes before saving:
for i = 1:length(dat)
    SL = dat(i).SL;
    %SL = rmfield(SL,'thetaphase');
    
    for k = 1:numel(SL)
        try
            SL(k).phase_precession.ts = [];
            SL(k).phase_precession.cs = [];
            SL(k).phase_precession.filtered_lfp = [];
            SL(k).phase_precession.filtered_lfp_phase = [];
            SL(k).phase_precession.theta_cycle = [];
            SL(k).phase_precession.lfpts = [];
            SL(k).phase_precession.b_raw = [];

        catch
        end
    end
    
    dat(i).SL = SL;
end

% ------- save
save('output.mat','controlled','dat')




%% Fix infusion thingy
load(SLNames{1});
for i = 1:40
    load(['/media/wchapman/RatBrains/Dropbox/Bill/Projects/other/cmonagha/preprocess/SL_8-oh-dpat_infusion/SLe' num2str(i) '.mat'])
    fns = fieldnames(SLe);
    for k = 1:length(fns)
        eval(['SL(' num2str(i) ').' fns{k} '=SLe.' fns{k} ';']);
    end
    
end
save(SLNames{1},'SL');
%%
%Need from Bill:

%-Constrain Jason’s rhythmicity measures to be within theta (6-12 Hz), rerun all
%-Rerun stats on Jason’s rhythmicity measures, make sure the  p3 values not all zeros (ensure it’s working)
    % -If no fit/value for one of first three sessions, shouldn’t plot










