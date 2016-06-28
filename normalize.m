function [normExMat, normMeanVal, meanVal] = normalize(exMat, tao)
    % Normalization across each rat, as in Wells 2013
    normVal = NaN(size(exMat,1),1);
    for i = 1:size(exMat,1)
        normVal(i) = nansum(((tao{i,1}./nansum(tao{i,1})).*exMat{i,1}));
    end
    
    normVal = normVal - nanmean(normVal);
    
    %normExMat = NaN(size(exMat));
    normExMat = cell(size(exMat));
    for i = 1:size(exMat,1)
        for k = 1:size(exMat,2)
            meanVal(i,k) = nansum(((tao{i,k}./nansum(tao{i,k})).*exMat{i,k}));
            normMeanVal(i,k) = meanVal(i,k) - normVal(i);
            normExMat{i,k} = exMat{i,k} - normVal(i); 
        end
    end
end
