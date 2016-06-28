function results = by_SL(SL, crunched)


    %% Single Unit Changes as function of drug
    
    
    %% Session Changes as function of drug
    
    
    %% LFP 
    
    % Speed Vs Frequency for all Rats
    figure;
    for i = 1:size(crunched.rat,1)
        for k = 1:4
            try
                subplot(size(crunched.rat,1),4,(i-1)*4+k)
                hold on
                plot(crunched.rat(i,k).bands.freq.meanSig);
                ylim([7.5 8.5])
            end
        end
    end
    
    % Run poster figures
    bands_simple = bands_by_SL(crunched);
    
    results = bands_simple;
    
    %% 
    %keyboard

