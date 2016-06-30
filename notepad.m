%% run for a single channel/session
features = {'frequency','amplitude'};
% sc = "session*channel"
for i = 1:numel(features)
    [sc(1).(features{i}).mean, sc(1).(features{i}).sigma, sc(1).(features{i}).speedDim,...
        sc(1).(features{i}).ft, sc(1).(features{i}).tao, sc(1).(features{i}).fet, ...
        sc(1).(features{i}).speed_us] = SpeedLFP.speedVSmeasure(root,'feature',features{i});
end

%% Compare a single channel across sessions:
% transform sc to nXm where n = number of channels, m=number of conditions
%sc(2) = sc; sc = (sc(:))';

rs1 = sc(1).amplitude.speed_us;
s1 = sc(1).amplitude.fet;
rs2 = sc(2).amplitude.speed_us;
s2 = sc(2).amplitude.fet;

ident = [zeros(size(rs1)); ones(size(rs2))];

[~, res.atab, res.ctab, res.stats] = aoctool([rs1;rs2], [s1;s2], ident,[],[],[],[],'off',5);

%%

rs1 = sc(1).amplitude.speedDim;
s1 = sc(1).amplitude.mean;
rs2 = sc(2).amplitude.speedDim;
s2 = sc(2).amplitude.mean;

ident = [zeros(size(rs1)); ones(size(rs2))];

[~, res.atab, res.ctab, res.stats] = aoctool([rs1;rs2], [s1;s2], ident);