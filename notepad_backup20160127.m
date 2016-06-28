%% History
% version history:
% 0.01 - wchapman 20160120. Initial writing based on conversation
%
% 0.02 - wchapman 20160121. Broke in to multiple scripts. Below history is
%        overall flow information.
%      - Crunching -> by_SL is immediate concern
%
% 0.03 - 

%% backup = 20160127

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
%{
cd /media/wchapman/RatBrains/Dropbox' (hasselmonians)'/UnitRecordingData/Caitlin''''s' Data'/Analyses/BillLivesHere/

addpath(genpath('/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/CMBs'))
addpath(genpath('/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin''s Data/Analyses/SLs'))
addpath /projectnb/hasselmogrp/wchapman/CMBHOME_stable/
%}

clear all; close all; clc
load output.mat
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

%% Intial Crunching
for slnum = 1:length(SLNames)
    load(SLNames{slnum});
    
    if slnum == 3
        SL(1:15,:) = []; % no short term recovery
    end
    
    dat(slnum).SL = SL;
    dat(slnum).ctps = assign_cells(dat(slnum).SL);
    dat(slnum).crunched = crunching(dat(slnum).SL, dat(slnum).ctps);
    dat(slnum).by_SL = by_SL(dat(slnum).SL, dat(slnum).crunched);
    
    temp = arrayfun(@(x) x.medianspeed, dat(slnum).crunched.session);
    [h, p] = ttest(temp(:,2),temp(:,1));
    dat(slnum).p_median_speed = p; 
end

arrayfun(@(x) x.p_median_speed, dat) < 0.05/7

%% Things needing control
for slnum = 1:length(controlmat)
    controlled(slnum).lfp_freq = nova(dat(slnum), dat(controlmat(slnum)));
    controlled(slnum).lfp_mag = nova_mag(dat(slnum), dat(controlmat(slnum)));
    controlled(slnum).lfp_pow = nova_pow(dat(slnum), dat(controlmat(slnum)));
    
    controlled(slnum).singleunit = anovas(dat(slnum), dat(controlmat(slnum)));
end
% Single unit rules: compare must be significant, then can check drug/veh

controlled = controlled(:);
save('output.mat','controlled','dat')


