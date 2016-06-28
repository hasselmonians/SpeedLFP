function SL = NaNSL(SL)

    inds = find(arrayfun(@(x) isempty(x.f), SL));
    
    for i = 1:length(inds)
        sle1 = SL(1);
        sle = SL(inds(i));
        if isempty(sle.cel)
            sle.cel = [NaN NaN];
        end
        
        if ~isnan(sle.cel(1));
            disp(i)
        end
        
        fn = fieldnames(sle);
        for k = 1:length(fn)
            eval(['sle.' fn{k} '=NaN' ';']);
        end
        sle.fname = SL(inds(i)).fname;
        %SLe.cel = SLo(i).cel;
        SL(inds(i)) = sle;
    end



end