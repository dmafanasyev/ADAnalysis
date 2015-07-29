function [ ddata, trend, imf, period, trendidx, residue, filtered, filteredidx ] = emd_processing( data, H, alpha, method, plots, noiseStd)
    %EMD_PROCESSING Process data using empirical mode decomposition
    %   Compute:
    %       - empirical mode decomposition (IMFs) using given method (EMD, EEMD, CEEMDAN)
    %       - detrended signal and trend-cyclical component using 4 criteria of trend index (energy, RZCN, statistical significance, low-frequency)
    %       - estimated periods of IMFs using zero-crossing method
    %       - filtered signal, that obtained after removing of statistically insignificant IMFs
    %
    %   Input:
    %
    %   Output:
    %
    %   Reference(s):
    %       Flandrin, P., Goncalves, P., Rilling, G., 2004. Detrending and denoising with empirical mode decomposition. EUSIPCO.
    %       Moghtader, A., Borgnat, P., Flandrin, P., 2011. Trend filtering: empirical mode decomposition versus l1 and Hodrick-Prescott. Advances in Adaptive Data Analysis 3 (1 and 2), 41–61.
    %       Colominas, M., Schlotthauer, G., Torres, M., Flandrin, P., 2012. Noise-assisted EMD methods in action. Advances in Adaptive Data Analysis 4 (4). 
    %       Afanasyev, D., Fedorova, E., Popov, V., 2014. Fine structure of the price-demand relationship in the electricity market: multi-scale correlation analysis. URL: http://mpra.ub.uni-muenchen.de/58827/.
    %
    %   Copyright (c) 2014 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0 2014.04.22: initial version
    %       1.1 2015.07.30: method renamed from 'deseasonalize_emd' to 'emd_processing'
    %
    
    if (nargin < 2)
        H = 0.5;
    end
    if (nargin < 3)
        alpha = 0.05;
    end
    if (nargin < 4)
        method = 'ceemdan';
    end
    if (nargin < 5)
        plots = 0;
    end
    if (nargin < 6)
        noiseStd = 0.2;
    end
    
    alphaErr = 'Use only alpha equal 0.05 or 0.01';
    hErr = 'Use only H equal 0.2, 0.5 or 0.8';
    methodErr = 'Use only method "emd", "eemd" or "ceemdan"';
    
    % see Flandrin et al., 2004
    if (H == 0.2)
        beta = 0.49;
        if(alpha == 0.05)
            a = 0.46; b = -2.44;
        elseif(alpha == 0.01)
            a = 0.45; b = -1.95;
        else
            error(alphaErr);
        end
    elseif(H == 0.5)
        beta = 0.72;
        if(alpha == 0.05)
            a = 0.47; b = -2.45;
        elseif(alpha == 0.01)
            a = 0.46; b = -1.92;
        else
            error(alphaErr);
        end
    elseif(H == 0.8)
        beta = 1.03;
        if(alpha == 0.05)
            a = 0.45; b = -2.33;
        elseif(alpha == 0.01)
            a = 0.45; b = -1.83;
        else
            error(alphaErr);
        end
    else
        error(hErr);
    end
    
    if(~strcmp(method, 'emd') && ~strcmp(method, 'eemd') && ~strcmp(method, 'ceemdan'))
        error(methodErr);
    end
    
    % see Table 1 in Moghtader et al., 2011 for spline interpolation
    if(alpha == 0.05)
         rzcntr = 2.7030; rzcntl = 1.7232;
    elseif(alpha == 0.01)
         rzcntr = 3.5317; rzcntl = 1.5073;
    end
    
    po = 2.01 + 0.2*(H-0.5) + 0.12*(H-0.5)^2;
    
    % decompose by IMFs
    if(strcmp(method, 'emd'))
        imf = emd(data, 'MAXMODES', 100, 'MAXITERATIONS', 5000)';
    elseif(strcmp(method, 'eemd'))
        imf = eemd(data, noiseStd, 300, 5000)';
    elseif(strcmp(method, 'ceemdan'))
        % see description of CEEMDAN in Colominas et al., 2012
        imf = ceemdan_par(data, noiseStd, 300, 5000)';
    end
    
    nImf = size(imf, 2);
    nObs = size(data, 1);
    dataMean = mean(data);
    
    means = zeros(nImf, 1);
    energy = zeros(nImf, 1);
    noise = zeros(nImf, 1);
    confidence = zeros(nImf, 1);
    period =  zeros(nImf, 1);
    rzcn = zeros(nImf-1, 1);
    numzerlast = 0;
    
    % see energy and RZCN criterias in Moghtader et al., 2011, statistical
    % significance criteria in Flandrin et al., 2004 and low-frequency
    % criteria in Afanasyev et al., 2014
    for i=1:nImf
        [indmin, indmax, indzer] = extr(imf(:, i)');
        numzercur = size(indzer, 2);
        
        if(i == 1)
            noise(1, 1) = sum(imf(:, i).^2, 1)/nObs;
        else
            noise(i, 1) = (noise(1, 1) / beta) * po^(-2*(1-H)*i);
            rzcn(i-1, 1) = numzerlast / numzercur;
            
            if (isinf(rzcn(i-1, 1)) || isnan(rzcn(i-1, 1)))
                rzcn(i-1, 1) = rzcntr;
            end
        end
        
        means(i, 1) = mean(sum(imf(:, 1:i), 2))/dataMean;
        energy(i, 1) = sum(imf(:, i).^2, 1)/nObs;
        confidence(i, 1) = 2^(log2(noise(i, 1)) + 2^(a*i + b));
        period(i, 1) = 4 * (nObs / (size(indmin, 2) + size(indmax, 2) + numzercur));
        numzerlast = numzercur;
    end
    
    trendidx = 0;
    filteredidx = [];
    filtered = zeros(nObs, 1);
    
    for i=nImf:-1:1
        if(i > (nImf/2+1) && energy(i, 1) > energy(i-1, 1) && energy(i, 1) >= confidence(i, 1) && (rzcn(i-1, 1) <= rzcntl || rzcn(i-1, 1) >= rzcntr))
            trendidx = i;
        end
        
        if(energy(i, 1) >= confidence(i, 1))
            filteredidx(1, end+1) = i;
            filtered = filtered + imf(:, i);
        end
    end
    
    % calculate trend, long-term season component & deseasonalized data
    residue = data - sum(imf(:, 1:nImf), 2);
    filtered = filtered + residue;
    trend = sum(imf(:, trendidx:nImf), 2) + residue;
    ddata = data - trend;
    
    % align min values of the raw and deseasonalized data
    ddata = ddata + min(data) - min(ddata);
    
    if (plots)
        figure;
        for i=1:nImf;
            subplot(nImf, 1, i);
            plot(imf(:, i));
            
            if (i == 1)
                title('IMFs');
            end
            
            ylabel(['M', num2str(i)]);
            xlim([1 nObs]);
            
            if (i < nImf);
                set(gca, 'xticklabel', '');
            end
        end
        
        figure;
        subplot(2, 3, 1);
        scatter(1:nImf, means, 'filled');
        xlim([0 nImf+1]);
        ylabel('Mean (standardized)');
        grid on;
        
        subplot(2, 3, 2);
        hold on;
        plot(log2(noise), '-k+');
        plot(log2(confidence), '--rx');
        scatter(1:nImf, log2(energy), 'filled');
        xlim([0 nImf+1]);
        legend('"Noise only" model', [num2str((1-alpha)*100), '% confidence interval'], 'Estimated energies');
        ylabel('log_2(Energy)');
        grid on;
        hold off;
        
        subplot(2, 3, 3);
        hold on;
        scatter(2:nImf, log2(rzcn), 'filled');
        line([0 nImf+1], [log2(rzcntl) log2(rzcntl)], 'LineStyle', '--', 'Color', 'r');
        line([0 nImf+1], [log2(rzcntr) log2(rzcntr)], 'LineStyle', '--', 'Color', 'r');
        xlim([0 nImf+1]);
        ylabel('log_2(RZCN)');
        legend('Ratio of the zero-crossing numbers', [num2str((1-alpha)*100), '% confidence interval']);
        grid on;
        hold off;
        
        subplot(2, 3, 4:6);
        plot([data trend]);
        xlim([1 nObs]);
        legend('Data', 'Trend');
    end
    
    if (nargin == 6)
        save(savePath);
    end
    
end