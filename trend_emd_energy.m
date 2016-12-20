function [ trend, filtered, idx, h, energy ] = trend_emd_energy( imfs, plots )
%TREND_EMD_ENERGY Summary of this function goes here
%   Detailed explanation goes here
%   Moghtader, A., Borgnat, P., Flandrin, P., 2011. Trend filtering: empirical mode decomposition versus l1 and Hodrick-Prescott. Advances in Adaptive Data Analysis 3 (1 and 2), 41–61.

     if(nargin == 0 || ~ismatrix(imfs))
        error('Input data must be non empty matrix');
    end
    
    if (nargin < 2)
        plots = 0;
    end
    
    % if time is columns then transpose matrix (I assume that the number of observations is greater than number of IMFs always)
    if (size(imfs, 2) > size(imfs, 1))
        imfs = imfs.';
    end
    
    nObs = size(imfs, 1);
    nImf = size(imfs, 2);
    
    energy = zeros(1, nImf);
    idx = zeros(1, nImf);
    
    % determine trend indexes using energy criteria
    for i=1:nImf
        energy(i) = sum(imfs(:, i).^2, 1)/nObs;
        
        if(i > 1 && energy(i) >= energy(i-1))
            idx(i) = i;
        end
    end
    
    idx_min = min(nonzeros(idx));
    
    trend = sum(imfs(:, idx_min:nImf), 2);
    filtered = sum(imfs(:, 1:idx_min-1), 2);
    
    if (plots)
        fontSize = 16;
        fontName = 'Helvetica';
        
        h = figure;
        subplot(1, 1, 1, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
        hold on;
        scatter(1:nImf, log2(energy), 'filled');
        xlim([0 nImf+1]);
        ylabel('log_2(Energy)');
        hold off;
    end

end
