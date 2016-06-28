for i = 1:length(dat)
    figure; xlim([-0.7 0.7]); hold on
    for m = 1:2
        clear model
        cnts = arrayfun(@(x) x.acorr_cnts, dat(i).SL(dat(i).ctps{2,2},m),'UniformOutput',0);
        lags = arrayfun(@(x) x.acorr_lags, dat(i).SL(dat(i).ctps{2,2},m),'UniformOutput',0);

        f = fittype('[a*(cos(w*x)+1)+a2*(cos(w*x/2)+1)+b]*exp(-abs(x)/tau)+c*exp(-x^2/tau2^2)',...
                'independent', 'x', ...
                'coefficients', {'a' 'a2' 'b' 'c' 'tau' 'tau2' 'w'}); % make model

        for k = 1:length(cnts)
            options = fitoptions(f); 

            options.Lower = [0 0 0 -inf 0 0  5*pi*2];
            options.Upper = [1 1 1 inf 5 .050 9*pi*2];
            options.MaxIter = 10^3;
            options.MaxFunEvals = 10^4;

            cor = cnts{k}; lag=lags{k};
            options.StartPoint = [  max(cor(lag>.1 & lag<.150))-min(cor(lag>.1 & lag<.150)),... % a
                                    0,... % a2
                                    mean(cor(lag>-.1&lag<.1)),... %b
                                    cor(round(end/2)),... %c
                                    .1,... % tau
                                    .015,... % tau2
                                    2*pi*8]; % w

            [model{k}] = fit(lag(lag~=0), cor(lag~=0), f, options);
        end
        % mean it
        m2 = model{1};
        fns = fieldnames(m2);
        for p = 1:length(fns)
            m2.(fns{p}) = mean(cellfun(@(x) x.(fns{p}), model));
        end
        if m2.c < 0
            m2.c=median(cellfun(@(x) x.c, model));
        end
        
        if m==1
            ca = plot(m2, 'k-');
        else
            ca = plot(m2,'r-');
        end
        set(ca, 'LineWidth', 2.5)
    end
    ylim([0 0.5])
end
addpath /media/wchapman/RatBrains/Dropbox' (hasselmonians)'/Bill/utils/bp
savefigs('means2',1,1)