function res = anovas(vehicle, drug)

    % dat(1).crunched.cell.FIELD    
    
    %% Single fields
    flds = {'f','f_intrins','theta_index','thetaskip','mrl','u2','gridness2','gridness3','si','coherence','poi_sfr_inter','poi_sfr_slope','nearestSix','theta_mra','theta_mrl'};
    for k = 1:length(flds)
        
        % Get the data out
        v = ['x_v = arrayfun(@(x) x.' flds{k} ', vehicle.crunched.cell, ''UniformOutput'', 0);'];
        d = ['x_d = arrayfun(@(x) x.' flds{k} ', drug.crunched.cell, ''UniformOutput'', 0);'];
        eval(v); eval(d);
        
        for ctp = 1:8 % For each cell grouping
            for bl = 1:4
                for xp = 1:4
                    vbl = x_v{ctp,bl}(:);
                    vxp = x_v{ctp,xp}(:);
                    dbl = x_d{ctp,bl}(:);
                    dxp = x_d{ctp,xp}(:);
                    
                    bads=isnan(vbl) | isnan(vxp);
                    vbl=vbl(~bads); vxp=vxp(~bads);
                    
                    bads=isnan(dbl) | isnan(dxp);
                    dbl=dbl(~bads); dxp=dxp(~bads);
                    
  
                    try %~(isempty(vbl) + isempty(vxp) + isempty(dbl) + isempty(dxp))
                        % compare
                        [p, table] = anova_rm({[vbl vxp], [dbl dxp]},0);
                        
                        eval(['cmp.' flds{k} '.' mapping(ctp) '.p1(bl,xp) = p(1);']);
                        eval(['cmp.' flds{k} '.' mapping(ctp) '.p2(bl,xp) = p(2);']);
                        eval(['cmp.' flds{k} '.' mapping(ctp) '.p3(bl,xp) = p(3);']);
                        eval(['cmp.' flds{k} '.' mapping(ctp) '.p4(bl,xp) = p(4);']);
    
                        % drug
                        p1 = '[~, ';
                        p2 = ['drg.' flds{k} '.' mapping(ctp) '.p(bl,xp),'];
                        p3 = ['drg.' flds{k} '.' mapping(ctp) '.est{bl,xp},'];
                        p4 = 'temp] = ttest(dbl,dxp);';
                        eval([p1 p2 p3 p4])

                        flds2 = fieldnames(temp);
                        for m = 1:length(flds2)
                            p1 = ['drg.' flds{k} '.' mapping(ctp) '.stats_' flds2{m} '(bl,xp)='];
                            p2 = ['temp.' flds2{m} ';'];
                            eval([p1 p2])
                        end

                        % veh
                        p1 = '[~, ';
                        p2 = ['veh.' flds{k} '.' mapping(ctp) '.p(bl,xp),'];
                        p3 = ['veh.' flds{k} '.' mapping(ctp) '.est{bl,xp},'];
                        p4 = 'temp] = ttest(vbl,vxp);';
                        eval([p1 p2 p3 p4])

                        flds2 = fieldnames(temp);
                        for m = 1:length(flds2)
                            p1 = ['veh.' flds{k} '.' mapping(ctp) '.stats_' flds2{m} '(bl,xp)='];
                            p2 = ['temp.' flds2{m} ';'];
                            eval([p1 p2])
                        end
                    catch
                        eval(['cmp.' flds{k} '.' mapping(ctp) '.p1(bl,xp) = NaN;']);
                        eval(['cmp.' flds{k} '.' mapping(ctp) '.p2(bl,xp) = NaN;']);
                        eval(['cmp.' flds{k} '.' mapping(ctp) '.p3(bl,xp) = NaN;']);
                        eval(['cmp.' flds{k} '.' mapping(ctp) '.p4(bl,xp) = NaN;']);
                        
                        eval(['drg.' flds{k} '.' mapping(ctp) '.p(bl,xp)=NaN;']);
                        eval(['drg.' flds{k} '.' mapping(ctp) '.est{bl,xp}=NaN;']);

                        %flds2 = fieldnames(temp);
                        %for m = 1:length(flds2)
                        %    p1 = ['drg.' flds{k} '.' mapping(ctp) '.stats_' flds2{m} '(bl,xp)='];
                        %    p2 = ['NaN;'];
                        %    eval([p1 p2])
                        %end

                        % veh
                        eval(['veh.' flds{k} '.' mapping(ctp) '.p(bl,xp)=NaN;']);
                        eval(['veh.' flds{k} '.' mapping(ctp) '.est{bl,xp}=NaN;']);
                    end
                    
                end
            end
        end  
    end
    
    res.drug = drg;
    res.vehi = veh;
    res.cmp = cmp;

end