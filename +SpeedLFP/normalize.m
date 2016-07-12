function [sc] = normalize(sc, features)
    % Normalization across each rat, as in Wells 2013
    %
    % [sc] = SpeedLFP.normalize(sc);
    
    for i = 1:length(features)
        normVal = NaN(size(sc,1),1);
        for k = 1:size(sc,1)
            normVal(k) = nansum(...
            (sc(k,1).(features{i}).tao./nansum(sc(k,1).(features{i}).tao))...
            .*sc(k,1).(features{i}).mean);
        end
        
        normVal = normVal - nanmean(normVal);
        
        for k = 1:size(sc,1)
            for m = 1:size(sc,2)
                meanVal(k,m) = nansum(...
                    (sc(k,m).(features{i}).tao./nansum(sc(k,m).(features{i}).tao))...
                    .*sc(k,m).(features{i}).mean);
                sc(k,m).(features{i}).mean_norm = ...
                    sc(k,m).(features{i}).mean - meanVal(k,m);
            end
        end
    end
end
