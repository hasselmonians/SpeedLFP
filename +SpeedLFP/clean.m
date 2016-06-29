function dat = clean(varargin)

    n = nargin;
    dims = cellfun(@(x) numel(x),varargin);

    if ~all(dims == dims(1))
        error('All inputs must have same dimensions');
    end

    dat = cell2mat(varargin);

    bads = isnan(dat) | ~isfinite(dat);
    bads = find(sum(bads,2));

    dat(bads,:) = [];

end