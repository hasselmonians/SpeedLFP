function results = crunching(SL, ctps)
    % results = population_analyses(SL);
    %
    % Runs a series of population analyses on an SL that has already been run
    % through the AppendAll/unpack process. 
    %
    %

    % version history:
    % 0.01 - wchapman 20160120. Initial writing based on conversation
    % 0.02 - wchapman 20160121. Broke in to multiple scripts. Below history
    % is for crunching only.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% Prep Structures %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Single Unit Crunching
    % For each of the cell types, compare the following:

    for i = 1:size(ctps,1)
        for k = 1:size(SL,2)
            inds = find(ctps{i,2});
            
            % Firing Rate
            %try
                celres(i,k).f = arrayfun(@(x) x.f, SL(ctps{i,2},k),'UniformOutput',1);

                % Instrinsic Frequency
                celres(i,k).f_intrins = arrayfun(@(x) x.f_intrins, SL(ctps{i,2},k));

                % Theta_index
                celres(i,k).theta_index = arrayfun(@(x) x.theta_index, SL(ctps{i,2},k));

                % thetaskip
                celres(i,k).thetaskip = arrayfun(@(x) x.thetaskip, SL(ctps{i,2},k));

                % hd mrl
                celres(i,k).mrl = arrayfun(@(x) x.mrl, SL(ctps{i,2},k));

                % hd u2
                celres(i,k).u2 = arrayfun(@(x) x.u2, SL(ctps{i,2},k));

                % gridness 2
                celres(i,k).gridness2 = arrayfun(@(x) x.gridness2, SL(ctps{i,2},k));

                % gridness 3
                celres(i,k).gridness3 = arrayfun(@(x) x.gridness3, SL(ctps{i,2},k));

                % spatial coherence across sessions
                %res(i,k).spatial = arrayfun(@(x) x.f, SL(ctps{i,2},k));

                % spatial information
                celres(i,k).si = arrayfun(@(x) x.si, SL(ctps{i,2},k));
                
                %nearest 6
                celres(i,k).nearestSix = arrayfun(@(x) x.nearestSix, SL(ctps{i,2},k));
                
                % theta stuff
                celres(i,k).theta_mrl = arrayfun(@(x) x.theta_mrl, SL(ctps{i,2},k));
                celres(i,k).theta_mra = arrayfun(@(x) x.theta_mra, SL(ctps{i,2},k));
                celres(i,k).theta_p =   arrayfun(@(x) x.theta_p,   SL(ctps{i,2},k));
                %celres(i,k).theta_tbl = arrayfun(@(x) x.theta_tbl, SL(ctps{i,2},k));
                celres(i,k).theta_pr = arrayfun(@(x) x.theta_pr,   SL(ctps{i,2},k));
                celres(i,k).theta_zr = arrayfun(@(x) x.theta_zr,   SL(ctps{i,2},k));
                
                ac = arrayfun(@(x,y) CMBHOME.Utils.moserac(x.ratemap,y.ratemap), SL(ctps{i,2},k), SL(ctps{i,2},1),'UniformOutput',0);
                celres(i,k).coherence = NaN(size(celres(i,k).nearestSix));
                for m = 1:length(ac)
                    [r, c] = size(ac{m});
                    cent = floor([r/2 c/2]);
                    d = floor([r c] - cent)/5;
                    d = sqrt(d(1)*d(2));
                    
                    [rr cc] = meshgrid(1:c,1:r);
                    C = sqrt((rr-cent(2)).^2+(cc-cent(1)).^2)<=d;
                    C = C.*ac{m};
                    C = nanmax(C(:));
                    celres(i,k).coherence(m) = C;
                end
            %{
            catch
                celres(i,k).f = NaN;
                celres(i,k).f_intrins = NaN;
                celres(i,k).theta_index = NaN;
                celres(i,k).thetaskip = NaN;
                celres(i,k).mrl = NaN;
                celres(i,k).u2 = NaN;
                celres(i,k).gridness2 = NaN;
                celres(i,k).gridness3 = NaN;
                celres(i,k).nearestSix = NaN;
                
            end
            %}
            % speed mod
            spdind = 32;
            for m = 1:length(inds)
                ind = inds(m);
                
                try
                    freq = SL(ind,k).speedMod.f(1:spdind);
                    vel = SL(ind,k).speedMod.v(1:spdind);
                    [~, R, P, b, y_int] = CMBHOME.Utils.LinearRegression(freq,vel);

                    celres(i,k).svfr_R2(m) = (R);
                    celres(i,k).svfr_p(m) = (P);
                    celres(i,k).svfr_slope(m) = (b);
                    celres(i,k).svfr_inter(m) = (y_int);
                catch
                    
                    celres(i,k).svfr_R2(m) = NaN;
                    celres(i,k).svfr_p(m) = NaN;
                    celres(i,k).svfr_slope(m) = NaN;
                    celres(i,k).svfr_inter(m) = NaN;
                    
                end
            end

            
            % poisson speed mod
            for m = 1:length(inds)
                ind = inds(m);
                
                try
                    celres(i,k).poi_sfr_R2(m) = SL(ind,k).poi_sfr.all.R2_lin;
                    celres(i,k).poi_sfr_inter(m) = SL(ind,k).poi_sfr.all.mle_linear(1);
                    celres(i,k).poi_sfr_slope(m) = SL(ind,k).poi_sfr.all.mle_linear(2);
                    celres(i,k).poi_sfr_p(m) = SL(ind,k).poi_sfr.all.pF_lin;
                    celres(i,k).poi_sfr_F(m) = SL(ind,k).poi_sfr.all.F_lin;
                catch
                    celres(i,k).poi_sfr_R2(m) = NaN;
                    celres(i,k).poi_sfr_inter(m) =  NaN;
                    celres(i,k).poi_sfr_slope(m) =  NaN;
                    celres(i,k).poi_sfr_p(m) =  NaN;
                    celres(i,k).poi_sfr_F(m) =  NaN;
                end
            end
            
            % Rhythcmicity Measures
            for m = 1:length(inds)
                ind = inds(m);
                
                try
                    celres(i,k).rhythmicity_p_rhythmic(m) = SL(ind,k).rhythmicity.stats_all.stats0.p_rhyth;
                    
                    celres(i,k).rhythmicity_base_tau(m) =  SL(ind,k).rhythmicity.stats_all.stats0.phat(1);
                    celres(i,k).rhythmicity_base_b(m) =  SL(ind,k).rhythmicity.stats_all.stats0.phat(2);
                    celres(i,k).rhythmicity_base_c(m) =  SL(ind,k).rhythmicity.stats_all.stats0.phat(3);
                    celres(i,k).rhythmicity_base_f(m) =  SL(ind,k).rhythmicity.stats_all.stats0.phat(4);
                    celres(i,k).rhythmicity_base_s(m) =  SL(ind,k).rhythmicity.stats_all.stats0.phat(5);
                    celres(i,k).rhythmicity_base_r(m) =  SL(ind,k).rhythmicity.stats_all.stats0.phat(6);
                    
                    celres(i,k).rhythmicity_p_skip(m) = SL(ind,k).rhythmicity.stats_all.stats0.p_sk;
                    celres(i,k).rhythmicity_skip(m) = SL(ind,k).rhythmicity.stats_all.stats0.phat(5);
       
                    celres(i,k).rhythmicity_frequency0(m) = SL(ind,k).rhythmicity.stats_all.stats0.phat(4);
                    celres(i,k).rhythmicity_p_freq(m) = SL(ind,k).rhythmicity.stats_all.p_freq;
                    celres(i,k).rhythmicity_frequency_intercept(m) = SL(ind,k).rhythmicity.stats_all.statsA.phat(5);
                    celres(i,k).rhythmicity_frequency_slope(m) = SL(ind,k).rhythmicity.stats_all.statsA.phat(6);
                    
                    celres(i,k).rhythmicity_amplitude0(m) = SL(ind,k).rhythmicity.stats_all.stats0.a;
                    celres(i,k).rhythmicity_p_amp(m) = SL(ind,k).rhythmicity.stats_all.p_a;
                    celres(i,k).rhythmicity_amplitude_intercept(m) = SL(ind,k).rhythmicity.stats_all.statsA.phat(8);
                    celres(i,k).rhythmicity_amplitude_slope(m) = SL(ind,k).rhythmicity.stats_all.statsA.phat(9);
                    
                catch
                    celres(i,k).rhythmicity_p_rhythmic(m) = NaN;
                    
                    celres(i,k).rhythmicity_base_tau(m) = NaN;
                    celres(i,k).rhythmicity_base_b(m) =  NaN;
                    celres(i,k).rhythmicity_base_c(m) =  NaN;
                    celres(i,k).rhythmicity_base_f(m) =  NaN;
                    celres(i,k).rhythmicity_base_s(m) =  NaN;
                    celres(i,k).rhythmicity_base_r(m) =  NaN;
                    
                    celres(i,k).rhythmicity_p_skip(m) = NaN;
                    celres(i,k).rhythmicity_skip(m) = NaN;
       
                    celres(i,k).rhythmicity_frequency0(m) = NaN;
                    celres(i,k).rhythmicity_p_freq(m) = NaN;
                    celres(i,k).rhythmicity_frequency_intercept(m) = NaN;
                    celres(i,k).rhythmicity_frequency_slope(m) = NaN;
                    
                    celres(i,k).rhythmicity_amplitude0(m) = NaN;
                    celres(i,k).rhythmicity_p_amp(m) = NaN;
                    celres(i,k).rhythmicity_amplitude_intercept(m) = NaN;
                    celres(i,k).rhythmicity_amplitude_slope(m) = NaN;
                    
                end
                
            end
        end
    end

    %% Session Crunching:
    %median running speed (one sample for each session)
    fnames = arrayfun(@(x) x.fname, SL(1,:),'UniformOutput',0);
    for i = 2:size(SL,1)
        if ~strcmp(SL(i,1).fname,fnames{end,1})
           fnames = [fnames; arrayfun(@(x) x.fname, SL(i,:),'UniformOutput',0)]; 
        end
    end
    
    for i = 1:numel(fnames)
        try
            c = load(fnames{i});
            sesres(i).medianspeed = c.root.spatial_scale * nanmedian(c.root.vel);
            sesres(i).fname = fnames{i};
        catch %missing session
            sesres(i).fname = 'NULL';
            sesres(i).medianspeed = NaN;
        end
    end
    
    sesres = reshape(sesres, size(fnames));
    
    %% LFP (Rat) Crunching:
    % should combine to make one data point for each rat. 
    
    for i = 1:numel(SL)
        try
            fname = SL(i).fname;
            under = strfind(fname,'_');
            hyper = strfind(fname, '-');
            hyper = hyper(end);
            ratname{i} = fname(under+1:hyper-1);
        catch
            ratname{i} = 'NULL';
        end
    end
    ratname = reshape(ratname, size(SL));
    bads = arrayfun(@(x) ~isstruct(x.bands),SL);
    ratnames = unique(ratname);
    

    for i = 1:length(ratnames)
        for k = 1:size(SL,2)
            rows = find(cellfun(@(x) ~isempty(x), strfind(ratname(:,k), ratnames{i})));
            rows = setdiff(rows, find(bads(:,k)));
            
            if ~isempty(rows)
                rows = rows(1);

                rat(i,k).name = ratnames{i};
                rat(i,k).ind = rows;
                rat(i,k).fname = SL(rows, k).fname;
                rat(i,k).bands = SL(rows, k).bands;
            else
                rat(i,k).name = 'NULL';
                rat(i,k).ind = NaN;
                
            end
        end
    end
   
    %% Output
    results.rat = rat;
    results.cell = celres;
    results.session = sesres;
    
end
