function [ filtered, trend, idx, h, noise, confidence, energy ] = trend_emd_statistical( imfs, plots, alpha )
%TREND_EMD_STATISTICAL Summary of this function goes here
%   Detailed explanation goes here
%   Flandrin, P., Goncalves, P., Rilling, G., 2004. Detrending and denoising with empirical mode decomposition. Proceedings of Eusipco, Wien (Austria), 1581-1584.
    
     if(nargin == 0 || ~ismatrix(imfs))
        error('Input data must be non empty matrix');
    end
    
    if (nargin < 2)
        plots = 0;
    end
    if (nargin < 3)
        alpha = 0.05;
    end
    
    % if time is columns then transpose matrix (I assume that the number of observations is greater than number of IMFs always)
    if (size(imfs, 2) > size(imfs, 1))
        imfs = imfs.';
    end
    
    alphaErr = 'Use only alpha equal to 0.1, 0.08, 0.05, 0.03 or 0.01';
    
    % calculation of corrected empirical Hurst exponent (R/S analysis based) for first IMF
    He = round(hurst(imfs(:,1)), 2);
    % mapping of He to H with available coefficients for the confidence interval computation: 0.45 <= He <= 0.55 mapped to H = 0.5 (random process)
    % TODO: to make a request from authors Flandrin et al., 2004 the coefficients for H = 0.1, 0.2, ..., 0.9, 1.0
    if(He < 0.45)
        H = 0.2;
    elseif(He > 0.55)
        H = 0.8;
    else
        H = 0.5;
    end
    
    % see Flandrin et al., 2004
    % a and b are available only for 1% and 5% significance levels, but for the compability with RZCN I have grouped it with some other values
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
    
    nObs = size(imfs, 1);
    nImf = size(imfs, 2);
    
    % see Flandrin et al., 2004
    po = 2.01 + 0.2*(H-0.5) + 0.12*(H-0.5)^2;
    
    filtered = zeros(nObs, 1);
    trend = zeros(nObs, 1);
    noise = zeros(1, nImf);
    confidence = zeros(1, nImf);
    energy = zeros(1, nImf);
    idx = zeros(1, nImf);
    
    % determine trend indexes using statistical criteria
    for i=1:nImf
        if(i == 1)
            noise(1) = sum(imfs(:, 1).^2, 1)/nObs;
        else
            noise(i) = (noise(1) / beta) * po^(-2*(1-H)*i);
        end
        confidence(i) = 2^(log2(noise(i)) + 2^(a*i + b));
        energy(i) = sum(imfs(:, i).^2, 1)/nObs;
        
        if(energy(i) >= confidence(i))
            idx(i) = i;
            trend = trend + imfs(:, i);
        else
            filtered = filtered + imfs(:, i);
        end
    end
    
    if (plots)
        fontSize = 16;
        fontName = 'Helvetica';
        lineWidth = 1.5;
        markerSize = 8;
        
        h = figure;
        subplot(1, 1, 1, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
        hold on;
        plot(log2(noise), '-k+', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
        plot(log2(confidence), '--rx', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
        scatter(2:nImf, log2(energy(2:nImf)), 'filled');
        xlim([0 nImf+1]);
        legend('"Noise only" model', [num2str((1-alpha)*100), '% confidence interval'], 'Estimated energy', 'Location', 'northwest');
        ylabel('log_2(Energy)');
        hold off;
    end

end
