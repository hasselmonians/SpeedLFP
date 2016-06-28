fields = {
                'f'
        'f_intrins'
      'theta_index'
        'thetaskip'
              'mrl'
               'u2'
        'gridness2'
        'gridness3'
               'si'
        'coherence'
    'poi_sfr_inter'
    'poi_sfr_slope'
       'nearestSix'
       };
   
ctps = {
             'all'
             'grids'
       'conjunctive'
    'nonconjunctive'
                'hd'
           'spatial'
             'theta'
       'interneuron'
    
};

% controlled(1).singleunit.cmp.(fields{i}).(ctps{k}).table{1,2}
p1 = {2,6};  %--> Time
p2 = {3,6};  %--> Group
p4 = {4,6};  %--> Interaction
p3 = {5,6};  %--> Matching

clear p

for i = 1:length(fields)
    for k = 1:length(ctps)
        p(:,i,k) = arrayfun(@(x) x.singleunit.cmp.(fields{i}).(ctps{k}).table{1,2}{4,6}, controlled);
    end
end


clc
for exp = 1:4
    [a, b] = ind2sub(size(squeeze(p(exp,:,:))), find(squeeze(p(exp,:,:))<(0.05/32)));
    
    %printit
    disp('----------------------------------')
    disp(['|         Experiment' num2str(exp) '            |'])
    disp('----------------------------------')
    disp([ctps(b) fields(a)]) 
end


%% filtered things:
fields = {
    'svfr_R2'
    'svfr_slope'
    'svfr_inter'
    'poi_sfr_R2'
    'poi_sfr_inter'
    'poi_sfr_slope'
    'rhythmicity_base_tau'
    'rhythmicity_base_b'
    'rhythmicity_base_c'
    'rhythmicity_base_f'
    'rhythmicity_base_s'
    'rhythmicity_base_r'
    'rhythmicity_skip'
    'rhythmicity_frequency0'
    'rhythmicity_frequency_intercept'
    'rhythmicity_frequency_slope'
    'rhythmicity_amplitude0'
    'rhythmicity_amplitude_intercept'
    'rhythmicity_amplitude_slope'
       };
   
ctps = {
             'all'
             'grids'
       'conjunctive'
    'nonconjunctive'
                'hd'
           'spatial'
             'theta'
       'interneuron'
};

% controlled(1).singleunit.cmp.(fields{i}).(ctps{k}).table{1,2}
p1 = {2,6};  %--> Time
p2 = {3,6};  %--> Group
p4 = {4,6};  %--> Interaction
p3 = {5,6};  %--> Matching

clear p

for i = 1:length(fields)
    for k = 1:length(ctps)
        p(:,i,k) = arrayfun(@(x) x.singleunit_filtered.cmp.(fields{i}).(ctps{k}).table{1,2}{4,6}, controlled);
    end
end


disp('-------------------------------------------------------------------------')
disp('|                            Filtered                                   |')
disp('-------------------------------------------------------------------------')

for exp = 1:4
    [a, b] = ind2sub(size(squeeze(p(exp,:,:))), find(squeeze(p(exp,:,:))<(0.05/1)));
    
    %printit
    disp('----------------------------------')
    disp(['|         Experiment' num2str(exp) '            |'])
    disp('----------------------------------')
    disp([ctps(b) fields(a)]) 
end