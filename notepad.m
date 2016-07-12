
%% Setup
addpath /media/wchapman/RatBrains/git/CMBHOME

sc(1,1).fname = '/media/wchapman/RatBrains/Dropbox/UnitRecordingData/Caitlin''s Data/Analyses/CMBs/diazepam_infusion/CMBobject_clamps-69.mat';
sc(1,2).fname = '/media/wchapman/RatBrains/Dropbox/UnitRecordingData/Caitlin''s Data/Analyses/CMBs/diazepam_infusion/CMBobject_clamps-70.mat';
sc(1,3).fname = '/media/wchapman/RatBrains/Dropbox/UnitRecordingData/Caitlin''s Data/Analyses/CMBs/diazepam_infusion/CMBobject_clamps-71.mat';
sc(1,4).fname = '/media/wchapman/RatBrains/Dropbox/UnitRecordingData/Caitlin''s Data/Analyses/CMBs/diazepam_infusion/CMBobject_clamps-72.mat';

for i = 1:4
    sc(1,i).active_lfp = 5;
end

features = {'frequency','amplitude'};

%% Run on each session
for i = 1:size(sc,1)
    for k = 1:size(sc,2)
        load(sc(i,k).fname);
        root.active_lfp = sc(i,k).active_lfp;
        root = root.FixTime;
        root.b_vel = [];
        root = root.AppendKalmanVel;
        
        for m = 1:numel(features)
            [sc(i,k).(features{m}).mean, sc(i,k).(features{m}).sigma, sc(i,k).(features{m}).speedDim,...
                sc(i,k).(features{m}).ft, sc(i,k).(features{m}).tao, sc(i,k).(features{m}).fet, ...
                sc(i,k).(features{m}).speed_us] = SpeedLFP.speedVSmeasure(root,'feature',features{m});
        end
    end
end

%% Wells et al normalization
sc = SpeedLFP.normalize(sc,features);

%% Setup the vectors for regression
rs = [];
c = [];
s = [];
d = [];
si = [];

for i = 1:size(sc,1)
    for k = 1:size(sc,2)
        rs = [rs; sc(i,k).amplitude.speedDim];
        c  = [c ; ones(size(sc(i,k).amplitude.speedDim))*i];
        s  = [s ; ones(size(sc(i,k).amplitude.speedDim))*k];
        d  = [d ; ones(size(sc(i,k).amplitude.speedDim))*(k==2)];
        si = [si; sc(i,k).frequency.mean];
    end
end


[~, res.atab, res.ctab, res.stats] = aoctool(rs, si, s);