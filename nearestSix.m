function [d] = nearestSix(root,Pd,ifplot)
% Calculates and returns the distance from the center Acorr peak to the
% nearest 6 peaks.
%
% Inputs: root -- CMBObject, function assumes root.cel is set
%         Pd   -- "PeakDistance", the minimum distance between peaks.
%         Defaults to 7, which works for most cells. Change if issues.
%         ifplot -- For debugging mostly. Defaults to 0
%
% Outputs: d -- The 6 distances (pixels).

if ~exist('Pd','var');Pd = 7;end

if isa(root,'CMBHOME.Session')
    rm = root.RateMap(root.cel);
else
    rm = root;
end
rmA = CMBHOME.Utils.moserac(rm);
%%
[~,maxInds] = CMBHOME.Utils.extrema2(rmA); %linear index
[rowInd colInd] = ind2sub(size(rmA),maxInds);

rr = rowInd - rowInd(1);
cr = colInd - colInd(1);

d = sqrt(rr.^2 + rr.^2);

%% Get rid of false peaks
locs = sqrt(rowInd.^2+colInd.^2);

peakHeight = NaN(length(rowInd),1);
for i = 1:length(rowInd)
    peakHeight(i) = rmA(rowInd(i),colInd(i));
end

[peakHeight inds] = sort(peakHeight,'descend'); % May be unnecessary, seems like they come out in order
rowInd = rowInd(inds);
colInd = colInd(inds);


idelete = zeros(size(peakHeight))<0; 
for i = 1:length(peakHeight)
    if ~idelete(i)
        if i >1
            for k = 1:i-1 %check to see if too close to any of the previous (larger) peaks
                d(i,k) = sqrt((rowInd(i) - rowInd(k)).^2 + (colInd(i) -colInd(k)).^2);
                if d(i,k) < Pd;
                    idelete(i) = idelete(i)+1;
                end
            end
        end
    end
end
inds(idelete) = []; %indeces into rowInds and colInd 

%% Filter
rowInd = rowInd(inds);
colInd = colInd(inds);
d = sqrt((rowInd - rowInd(1)).^2 + (colInd - colInd(1)).^2);
[d inds] = sort(d,'ascend');

try
    d = d(2:7);
catch err
    if strcmp(err.identifier,'MATLAB:badsubscript')
        warning('Totes bad cell');
        d = NaN;
        return
    end
end
    

colInd = colInd(inds(2:7));
rowInd = rowInd(inds(2:7));


%}
%% ifplot
if ~exist('ifplot','var');ifplot=1;end

if ifplot==1
    figure
    imagesc(rmA)
    hold on
    plot(colInd,rowInd,'ko','MarkerFaceColor',[0 0 0])
end

end