clear all; close all; clc
load output.mat
load threshs.mat

%addpath H:\Dropbox' (hasselmonians)'\UnitRecordingData\Caitlin''''s' Data'\Analyses\BillLivesHere\

SLNames = {'SL_8-oh-dpat_infusion_output',...    %1
           'SL_8-oh-dpat_systemic_output',...    %2
           'SL_diazepam_infusion_output',...     %3
           'SL_diazepam_systemic_output',...     %4
           'SL_DMSO_systemic_output',...         %5
           'SL_PBS_infusion_output',...          %6
           'SL_PBS_systemic_output'};            %7

threshmat = [2 2 1 1 1 2 2];  % which set of thresholds to use
controlmat = [6 7 6 5];       % For the experiments, which control to use

