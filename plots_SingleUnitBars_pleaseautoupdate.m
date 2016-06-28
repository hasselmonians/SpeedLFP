%% Single Unit Plots
 
for slnum = 4;
    %1:length(dat)
    %
    % To change these figure:
    % 1) Set SLNum to 1-7 (see line 46 of notepad)
    % 2) Set line 8 x.f to x.FIELDNAME 
    % 3) Change line 100 to cell type desired (see assign_cells.m)
    clear f;
    f = arrayfun(@(x) x.gridness3, dat(slnum).crunched.cell,'UniformOutput',0);
    ctp = 4;
    
    inds = find(~(sum(cell2mat(cellfun(@(x) isnan(x), f(ctp,1:3),'UniformOutput',0)),2)));
    
    figure
    hold on
    datum = cell2mat(f(ctp,:));
    %datum = datum(inds,:)/(.2750/2);
    bar(nanmean(datum,1),'w','LineWidth',3);
    %plot(cell2mat(f(ctp,:))');
    plot(datum')
    ss = nanstd(datum,0,1) ./ arrayfun(@(x) sqrt(numel(x)), f(ctp,:));
    errorbar(1:4, nanmean(datum,1), ss,'k.','LineWidth',3, 'MarkerSize',0.1);
    title([SLNames(slnum) dat(1).ctps(ctp)],'FontSize',24,'Fontname','Times','FontWeight','Bold')
    ylabel('Spatial Information','FontSize',18,'Fontname','Times','FontWeight','Bold')
    %ylim([0 4.5])
    set(gca,'fontsize',18, 'Fontname', 'Times','FontWeight','Bold')
    set(gcf,'renderer','painters');
    set(gcf,'PaperPositionMode','auto')
    %set(gca,'box','on')
    %set(gca,'XTick',0)
    %set(gca,'YAxisLocation','right')
    %set(gca,'YColor','w')
    axis square
end


% Single unit rules: compare must be significant, then can check drug/veh

%{
Cell Types:
ctps{1,1} = 'all';   
ctps{2,1} = 'grids';   
ctps{3,1} = 'conjunctive';
ctps{4,1} = 'nonconjunctive'
ctps{5,1} = 'hd';
ctps{6,1} = 'spatial';
ctps{7,1} = 'theta';
ctps{8,1} = 'interneuron';
%}

%% Filtered Bar Plots

clear filt; clear val; clear inds; clear ss;

filter_col = [1 1 1 1];  % Filter each of the 4 columns based on field val in which column
                         % [1 1 1 1] to filter based on baseline
                         % [1 2 3 4] to filter based on self

field = 'poi_sfr_slope';
field_filt = 'poi_sfr_p';
thresh = 0.05;
ctp = 2;

for slnum = 1:length(dat)
    
    filt = eval(['arrayfun(@(x) x.' field_filt '(:),dat(slnum).crunched.cell,''UniformOutput'',0);']);
    val =  eval(['arrayfun(@(x) x.' field '(:),dat(slnum).crunched.cell,''UniformOutput'',0);']);
    
    for k = 1:4
        inds(:,k) = cellfun(@(x) x<=thresh, filt(:,filter_col(k)),'UniformOutput',0);
    end
    
    val = cellfun(@(x,y) x(y), val, inds,'UniformOutput',0);
    
    figure
    hold on
    bar(cellfun(@(x) nanmean(x), val(ctp,:)),'w','LineWidth',3);
    plot(cell2mat(cellfun(@(x) (x), val(ctp,:),'UniformOutput',0))');
    ss = cellfun(@(x) nanstd(x), val(ctp,:)) ./ arrayfun(@(x) sqrt(numel(x)), val(ctp,:));
    errorbar(1:4, cellfun(@(x) nanmean(x), val(ctp,:)),ss,'k.','LineWidth',3, 'MarkerSize',0.1);
    title([SLNames(slnum) dat(1).ctps(ctp)],'FontSize',24)
    set(gca,'fontsize',18)
    axis square
end




%% Do them all!
%/media/wchapman/RatBrains/Dropbox (hasselmonians)/UnitRecordingData/Caitlin's Data/Dissertation/Figures/diazepam_systemic

fields = {'f', ...
        'f_intrins', ...
        'theta_index', ...
        'thetaskip', ...
        'mrl', ...
        'u2', ...
        'gridness2', ...
        'gridness3', ...
        'si', ...
        'nearestSix', ...
        'coherence', ...
        'svfr_R2', ...
        'svfr_slope', ...
        'svfr_inter', ...
        'poi_sfr_R2', ...
        'poi_sfr_inter', ...
        'poi_sfr_slope', ...
        'poi_sfr_F', ...
        'rhythmicity_base_tau', ...
        'rhythmicity_base_b', ...
        'rhythmicity_base_c', ...
        'rhythmicity_base_f', ...
        'rhythmicity_base_s', ...
        'rhythmicity_base_r', ...
        'rhythmicity_skip', ...
        'rhythmicity_frequency0', ...
        'rhythmicity_frequency_intercept', ...
        'rhythmicity_frequency_slope', ...
        'rhythmicity_amplitude0', ...
        'rhythmicity_amplitude_intercept', ...
        'rhythmicity_amplitude_slope'};

filters = {};

for slnum = 1:length(dat)    
    for field = 1:length(fields)
        val =  eval(['arrayfun(@(x) x.' fields{field} '(:),dat(slnum).crunched.cell,''UniformOutput'',0);']);
        for ctp = 1:size(dat(1).crunched.cell,1)
            figure
            hold on
            bar(cellfun(@(x) nanmean(x), val(ctp,:)),'w','LineWidth',3);
            plot(cell2mat(cellfun(@(x) (x), val(ctp,:),'UniformOutput',0))');
            ss = cellfun(@(x) nanstd(x), val(ctp,:)) ./ arrayfun(@(x) sqrt(numel(x)), val(ctp,:));
            errorbar(1:4, cellfun(@(x) nanmean(x), val(ctp,:)),ss,'k.','LineWidth',3, 'MarkerSize',0.1);
            title([SLNames(slnum) dat(1).ctps(ctp)],'FontSize',24)
            set(gca,'fontsize',18)
            axis square
            ylabel(fields{field})
            
            mkdir(['autofigures' filesep SLNames{slnum} filesep dat(slnum).ctps{ctp,1}])
            fname = ['autofigures' filesep SLNames{slnum} filesep dat(slnum).ctps{ctp,1} '_' fields{field}];
            
            savefig_fname(fname)
            close all
        end
    end  
end


%% Do them all filtered!

fields = {
    'svfr_R2', ...
    'svfr_slope', ...
    'svfr_inter', ...
    'poi_sfr_R2', ...
    'poi_sfr_inter', ...
    'poi_sfr_slope', ...
    'rhythmicity_base_tau', ...
    'rhythmicity_base_b', ...
    'rhythmicity_base_c', ...
    'rhythmicity_base_f', ...
    'rhythmicity_base_s', ...
    'rhythmicity_base_r', ...
    'rhythmicity_skip', ...
    'rhythmicity_frequency0', ...
    'rhythmicity_frequency_intercept', ...
    'rhythmicity_frequency_slope', ...
    'rhythmicity_amplitude0', ...
    'rhythmicity_amplitude_intercept', ...
    'rhythmicity_amplitude_slope'};

filters = {
    'svfr_p', ...
    'svfr_p', ...
    'svfr_p', ...
    'poi_sfr_p', ...
    'poi_sfr_p', ...
    'poi_sfr_p', ...
    'rhythmicity_p_rhythmic', ...
    'rhythmicity_p_rhythmic', ...
    'rhythmicity_p_rhythmic', ...
    'rhythmicity_p_rhythmic', ...
    'rhythmicity_p_rhythmic', ...
    'rhythmicity_p_rhythmic', ...
    'rhythmicity_p_skip', ...
    'rhythmicity_p_freq', ...
    'rhythmicity_p_freq', ...
    'rhythmicity_p_freq', ...
    'rhythmicity_p_amp', ...
    'rhythmicity_p_amp', ...
    'rhythmicity_p_amp', ...
    };

thresh = 0.05;
filter_col = [1 1 1 1];  % Filter each of the 4 columns based on field val in which column
                         % [1 1 1 1] to filter based on baseline
                         % [1 2 3 4] to filter based on self
                         
for slnum = 1:length(dat)    
    for field = 1:length(fields)
        filt =  eval(['arrayfun(@(x) x.' filters{field} '(:),dat(slnum).crunched.cell,''UniformOutput'',0);']);
        val =  eval(['arrayfun(@(x) x.' fields{field} '(:),dat(slnum).crunched.cell,''UniformOutput'',0);']);
        
        for k = 1:4
            inds(:,k) = cellfun(@(x) x<=thresh, filt(:,filter_col(k)),'UniformOutput',0);
        end
        val = cellfun(@(x,y) x(y), val, inds,'UniformOutput',0);
        
        for ctp = 1:size(dat(1).crunched.cell,1)
            figure
            hold on
            bar(cellfun(@(x) nanmean(x), val(ctp,:)),'w','LineWidth',3);
            plot(cell2mat(cellfun(@(x) (x), val(ctp,:),'UniformOutput',0))');
            ss = cellfun(@(x) nanstd(x), val(ctp,:)) ./ arrayfun(@(x) sqrt(numel(x)), val(ctp,:));
            errorbar(1:4, cellfun(@(x) nanmean(x), val(ctp,:)),ss,'k.','LineWidth',3, 'MarkerSize',0.1);
            title([SLNames(slnum) dat(1).ctps(ctp)],'FontSize',24)
            set(gca,'fontsize',18)
            axis square
            ylabel(fields{field})
            
            mkdir(['autofigures_filtered' filesep SLNames{slnum} filesep dat(slnum).ctps{ctp,1}])
            fname = ['autofigures_filtered' filesep SLNames{slnum} filesep dat(slnum).ctps{ctp,1} '_' fields{field}];
            
            savefig_fname(fname)
            close all
        end
    end  
end

