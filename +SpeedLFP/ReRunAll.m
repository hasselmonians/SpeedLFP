function SL = ReRunAll(SL)

    % Reruns all fields only for cells that have been added. AKA, looks for
    % SL.cel set, and empty SL.f

    for  i = 1:numel(SL)

        if ~isnan(SL(i).cel(1)) && isempty(SL(i).f)




             %% Incorporate runRoots things:


