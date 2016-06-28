%% dpat--------------------------------------------------------------------
% SL_8-oh-dpat_infusion_output
load SL_8-oh-dpat_infusion_output
rows = {1:13 14:18 19:53 54:65 66:74 75:80 81:88 89:95 96:101 102:112}';
dose =  [1 1 1 1 1 1 1 2 2 2];
rat = [1 1 2 3 2 4 4 5 6 7];

% Cell
dpat(1).SL = SL(cell2mat(rows(find(dose==1))'),:);
dpat(2).SL = SL(cell2mat(rows(find(dose==2))'),:);
dpat=dpat(:);

for slnum = 1:length(dpat)
    dpat(slnum).ctps = assign_cells(dpat(slnum).SL, threshs(threshmat(slnum)));
    dpat(slnum).crunched = crunching(dpat(slnum).SL, dpat(slnum).ctps);
    dpat(slnum).by_SL = by_SL(dpat(slnum).SL, dpat(slnum).crunched);
end

dpat_controlled(1).lfp_freq = nova(dpat(1), dpat(2));
dpat_controlled(1).lfp_mag = nova_mag(dpat(1), dpat(2));
dpat_controlled(1).lfp_pow = nova_pow(dpat(1), dpat(2));

dpat_controlled(1).singleunit = anovas(dpat(1), dpat(2));
dpat_controlled(1).singleunit_filtered = anovas_filtered(dpat(1), dpat(2));



%% diaz--------------------------------------------------------------------
% diaz
load SL_diazepam_infusion_output
rows = {1:19 20:28 29:39 40:59 60:74 75:84 85:92};
dose =  [1 1 1 2 2 2 2];
rat  =  [1 2 3 4 5 5 5];

% Cell
diaz(1).SL = SL(cell2mat(rows(find(dose==1)))',:);
diaz(2).SL = SL(cell2mat(rows(find(dose==2)))',:);
diaz=diaz(:);

for slnum = 1:length(diaz)
    diaz(slnum).ctps = assign_cells(diaz(slnum).SL, threshs(threshmat(slnum)));
    diaz(slnum).crunched = crunching(diaz(slnum).SL, diaz(slnum).ctps);
    diaz(slnum).by_SL = by_SL(diaz(slnum).SL, diaz(slnum).crunched);
end

diaz_controlled(1).lfp_freq = nova(diaz(1), diaz(2));
diaz_controlled(1).lfp_mag = nova_mag(diaz(1), diaz(2));
diaz_controlled(1).lfp_pow = nova_pow(diaz(1), diaz(2));

diaz_controlled(1).singleunit = anovas(diaz(1), diaz(2));
diaz_controlled(1).singleunit_filtered = anovas_filtered(diaz(1), diaz(2));

%% output
save(['/media/wchapman/RatBrains/Dropbox/UnitRecordingData/Caitlin''s Data/Analyses/BillLivesHere/drugs.mat'],'dpat','dpat_controlled','diaz','diaz_controlled');
