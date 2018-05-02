function [ ddata, trend, imf, period, trendidx, residue, filtered, filteredidx ] = emd_processing( data, method, plots, alpha, noiseStd, numR, SNRFlag, H)
    %EMD_PROCESSING Process data using empirical mode decomposition (spline interpolation)
    %   Compute:
    %       - empirical mode decomposition (IMFs) using given method (EMD, EEMD, CEEMDAN, iCEEMDAN)
    %       - detrended signal and trend-cyclical component using 4 criteria of trend index (energy, RZCN, statistical significance, low-frequency)
    %       - estimated periods of IMFs using zero-crossing method
    %       - filtered signal, that obtained after removing of statistically insignificant IMFs
    %
    %   Input:
    %       data - time-series to processing (row is time)
    %       method - method of decomposition: 'emd', 'eemd', 'ceemdan' or 'iceemdan'
    %       plots - boolean flag (0 or 1) for graphics plot, default is 0
    %       alpha - significance level, default is 0.05
    %       noiseStd - noise standard deviation, only for ensembles methods
    %       numR - noise ensemble size for EEMD/CEEMDAN/iCEEMDAN
    %       SNRFlag - noise amplitude strategy (only for iCEEMDAN): 1 - SNR increases for every stage, 2 - SNR is constant.
    %       H - assumed Hurst index of first "noise" IMF (for statistical significance test); only values 0.2, 0.5 and 0.8 supported
    %       
    %   Output:
    %       ddata - detrended time-series
    %       trend - trend-cyclic component
    %       imf - IMFs of time-series (row is time, column is IMF)
    %       period - IMFs mean periods (estimated using zero-crossing method)
    %       trendidx - minimum index of IMF for trend-cylic component
    %       residue - residue of decomposition (~0 for complete methods)
    %       filtered - noise filtered time-series
    %       filteredidx - indexes of IMFs that pass statistical significance test
    %
    %   References:
    %       Flandrin, P., Goncalves, P., Rilling, G., 2004. Detrending and denoising with empirical mode decomposition. Proceedings of Eusipco, Wien (Austria), 1581-1584.
    %       Moghtader, A., Borgnat, P., Flandrin, P., 2011. Trend filtering: empirical mode decomposition versus l1 and Hodrick-Prescott. Advances in Adaptive Data Analysis 3 (1 and 2), 41-61.
    %       Colominas, M., Schlotthauer, G., Torres, M., Flandrin, P., 2012. Noise-assisted EMD methods in action. Advances in Adaptive Data Analysis 4 (4).
    %       Afanasyev, D., Fedorova, E., Popov, V., 2015. Fine structure of the price-demand relationship in the electricity market: multi-scale correlation analysis. Energy Economics 51, 215-226.
    %       Colominas, M., Schlotthauer, G., Torres, M., 2014. Improve complete ensemble EMD: A suitable tool for biomedical signal processing. Biomedical Signal Processing and Control, Vol. 14, 19-29
    %
    %   Copyright (c) 2014-2015 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0  2014.04.22: initial version
    %       1.1  2015.07.30: method renamed to 'emd_processing'
    %       1.2  2015.08.04: added improved CEEMDAN (iCEEMDAN)
    %       1.21 2015.08.13: added 2 optional params 'numR' and 'SNRFlag'
    %                        set optimal value of noise ensemble size as default for numR.
    %       1.3  2015.09.02: explicit estimation of empirical Hurst exponent for first IMF
    %                        use of iCEEMDAN routine optimized for parallel computations
    %                        reordered input parameters
    %                 
    
    if(nargin == 0 || ~isvector(data))
        error('Input data must be non empty vector');
    end
    
    if (nargin < 2)
        method = 'ceemdan';
    end
    if (nargin < 3)
        plots = 0;
    end
    if (nargin < 4)
        alpha = 0.05;
    end
    if (nargin < 5)
        noiseStd = 0.2;
    end
    if (nargin < 6)
        % the optimal value for the CEEMDAN and iCEEMDAN (with constant noise amplitude) in sense of the mean relative range (from the run to run) and computional time minimization
        numR = 1000;
    end
    if (nargin < 7)
        SNRFlag = 2;
    end
    if (nargin < 8)
        H = 0;
    end
    
    if(~strcmp(method, 'emd') && ~strcmp(method, 'eemd') && ~strcmp(method, 'ceemdan') && ~strcmp(method, 'iceemdan'))
        error('Use only method "emd", "eemd", "ceemdan" or "iceemdan"');
    end
    
    if(strcmp(method, 'iceemdan') && SNRFlag == 1 && numR < 3000)
        warning('For the results stability of iCEEMDAN with SNR increases for every stage (SNRFlag = 1) use ensemble size N >= 3000');
    end
    
    if size(data,2) > 1
        data = data.';
    end
    
    % decomposition into IMFs
    if(strcmp(method, 'emd'))
        imf = emd(data, 'MAXITERATIONS', Inf)';
        %imf = emdc([], data)';
        residue = imf(:,end);
        imf = imf(:,1:end-1);
    else
        if(strcmp(method, 'eemd'))
            imf = eemd(data, noiseStd, numR, Inf)';
        elseif(strcmp(method, 'ceemdan'))
            imf = ceemdan_par(data, noiseStd, numR, Inf)';
            %imf = ceemdan_fast(data, noiseStd, 1000, Inf)';
        elseif(strcmp(method, 'iceemdan'))
            imf = iceemdan_par(data, noiseStd, numR, Inf, SNRFlag)';
        end
        
        residue = data - sum(imf, 2);
    end
    
    alphaErr = 'Use only alpha equal to 0.1, 0.08, 0.05, 0.03 or 0.01';
    hErr = 'Use only H equal to 0.2, 0.5 or 0.8';
    if(H == 0)
        % calculation of corrected empirical Hurst exponent (R/S analysis based) for first IMF
        He = round(hurst(imf(:,1)), 2);
        % mapping of He to H with available coefficients for the confidence interval computation: 0.45 <= He <= 0.55 mapped to H = 0.5 (random process)
        % TODO: to make request from authors Flandrin et al., 2004 the coefficients for H = 0.1, 0.2, ..., 0.9, 1.0
        if(He < 0.45)
            H = 0.2;
        elseif(He > 0.55)
            H = 0.8;
        else
            H = 0.5;
        end
    end
    
    % see Flandrin et al., 2004
    % a and b available only for 1% and 5% significance levels, but for the compability with RZCN I have grouped it with some other values
    if (H == 0.2)
        beta = 0.49;
        if(alpha == 0.05 || alpha == 0.08 || alpha == 0.1)
            a = 0.46; b = -2.44;
        elseif(alpha == 0.01 || alpha == 0.03)
            a = 0.45; b = -1.95;
        else
            error(alphaErr);
        end
    elseif(H == 0.5)
        beta = 0.72;
        if(alpha == 0.05 || alpha == 0.08 || alpha == 0.1)
            a = 0.47; b = -2.45;
        elseif(alpha == 0.01 || alpha == 0.03)
            a = 0.46; b = -1.92;
        else
            error(alphaErr);
        end
    elseif(H == 0.8)
        beta = 1.03;
        if(alpha == 0.05 || alpha == 0.08 || alpha == 0.1)
            a = 0.45; b = -2.33;
        elseif(alpha == 0.01 || alpha == 0.03)
            a = 0.45; b = -1.83;
        else
            error(alphaErr);
        end
    else
        error(hErr);
    end
    
    % see Table 1 in Moghtader et al., 2011 for spline interpolation
    if(alpha == 0.1)
         rzcntr = 2.4117; rzcntl = 1.7941;
    elseif(alpha == 0.08)
         rzcntr = 2.5000; rzcntl = 1.7708;
    elseif(alpha == 0.05)
         rzcntr = 2.7030; rzcntl = 1.7232;
    elseif(alpha == 0.03)
         rzcntr = 3.0238; rzcntl = 1.6647;
    elseif(alpha == 0.01)
         rzcntr = 3.5317; rzcntl = 1.5073;
    end
    
    % see Flandrin et al., 2004
    po = 2.01 + 0.2*(H-0.5) + 0.12*(H-0.5)^2;
    
    nObs = size(data, 1);
    nImf = size(imf, 2);
    dataMean = mean(data);
    
    means = zeros(nImf, 1);
    energy = zeros(nImf, 1);
    noise = zeros(nImf, 1);
    confidence = zeros(nImf, 1);
    period =  zeros(nImf, 1);
    rzcn = zeros(nImf-1, 1);
    numzerlast = 0;
    
    % calculation of the all needed criteria for trend-filtering and data-smoothing
    for i=1:nImf
        [period(i, 1), ~, ~, indzer, numzercur] = period_zero_cross(imf(:, i));
        
        % if IMF has no zero crossing or the complete approaches is used and the current IMF is last then force set numzercur = 0
        if(isempty(indzer) || (i == nImf && (strcmp(method, 'ceemdan') || strcmp(method, 'iceemdan'))))
            numzercur = 0;
        end
        
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
        numzerlast = numzercur;
     end
    
    % trend-filtering and data-smoothing
    % see statistical significance criteria in Flandrin et al., 2004, energy and RZCN criteria in Moghtader et al., 2011 and low-frequency criteria in Afanasyev et al., 2015
    trendidx = 0;
    filteredidx = zeros(1, nImf);
    filtered = zeros(nObs, 1);
    
    for i=nImf:-1:2
        % for complete approaches force include last IMF into trend
        if((i == nImf && (strcmp(method, 'ceemdan') || strcmp(method, 'iceemdan'))) ...
                || (i >= (nImf/2+1) && energy(i, 1) >= confidence(i, 1) && energy(i, 1) > energy(i-1, 1) && (rzcn(i-1, 1) <= rzcntl || rzcn(i-1, 1) >= rzcntr)))
            trendidx = i;
        end
        
        if(energy(i, 1) >= confidence(i, 1))
            filteredidx(1, i) = i;
            filtered = filtered + imf(:, i);
        end
    end
    
    % indexes of statistically significant IMFs
    filteredidx = nonzeros(filteredidx)';
    
    % trend-cyclical component, filtered and detrended data
    filtered = filtered + residue;
    if(strcmp(method, 'emd') && trendidx == 0)
        trend = residue;
    else
        trend = sum(imf(:, trendidx:nImf), 2) + residue;
    end
    ddata = data - trend;
    
    % align min values of the raw and detrended data
    % ddata = ddata + min(data) - min(ddata);
    
    if (plots)
        fontSize = 16;
        fontName = 'Helvetica';
        lineWidth = 1;
    
        %figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        figure;
        for i=1:nImf
            subplot(nImf, 1, i,  'FontName', fontName, 'FontSize', fontSize, 'Box', 'on');
            plot(imf(:, i), 'LineWidth', lineWidth);
            
            if (i == 1)
                title('IMFs');
            end
            
            ylabel(['M', num2str(i)]);
            %ylim([1.01*min(imf(:, i)) 1.01*max(imf(:, i))]);
            xlim([1 nObs]);
            
            if (i < nImf)
                set(gca, 'xticklabel', '');
            end
        end
        
        figure;
        subplot(2, 3, 1, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
        hold on;
        scatter(1:nImf, means, 'filled');
        xlim([0 nImf+1]);
        ylabel('Mean (standardized)');
        grid on;
        hold off;
        
        subplot(2, 3, 2, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
        hold on;
        plot(log2(noise), '-k+', 'LineWidth', lineWidth);
        plot(log2(confidence), '--rx', 'LineWidth', lineWidth);
        scatter(2:nImf, log2(energy(2:nImf)), 'filled');
        xlim([0 nImf+1]);
        legend('"Noise only" model', [num2str((1-alpha)*100), '% confidence interval'], 'Estimated energy');
        ylabel('log_2(Energy)');
        hold off;
        
        subplot(2, 3, 3, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
        hold on;
        scatter(2:nImf, log2(rzcn), 'filled');
        line([0 nImf+1], [log2(rzcntl) log2(rzcntl)], 'LineStyle', '--', 'Color', 'r');
        line([0 nImf+1], [log2(rzcntr) log2(rzcntr)], 'LineStyle', '--', 'Color', 'r');
        xlim([0 nImf+1]);
        ylabel('log_2(RZCN)');
        legend('Ratio of the zero-crossing numbers', [num2str((1-alpha)*100), '% confidence interval']);
        hold off;
        
        subplot(2, 3, 4:6, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
        hold on;
        plot(data, 'LineWidth', lineWidth, 'Color', [0.7 0.7 0.7], 'LineWidth', lineWidth/2);
        plot(trend, 'LineWidth', lineWidth, 'Color', 'r', 'LineWidth', lineWidth);
        xlim([1 nObs]);
        legend('Data', 'Trend');
        hold off;
    end
end