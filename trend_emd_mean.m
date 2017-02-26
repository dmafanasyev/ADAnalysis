function [ trend, filtered, idx, means, h] = trend_emd_mean( imfs, plots, threshold)
%TREND_EMD_MEAN Summary of this function goes here
%   Detailed explanation goes here
%   Flandrin, P., Goncalves, P., Rilling, G., 2004. Detrending and denoising with empirical mode decomposition. Proceedings of Eusipco, Wien (Austria), 1581-1584.

    if(nargin == 0 || ~ismatrix(imfs))
        error('Input data must be non empty matrix');
    end
    
    if (nargin < 2)
        plots = 0;
    end
    if (nargin < 3)
        threshold = 0.05;
    end
    
    % if time is columns then transpose matrix (this is assumed that the number of observations is greater than number of IMFs always)
    if (size(imfs, 2) > size(imfs, 1))
        imfs = imfs.';
    end
    
    [nObs, nImf] = size(imfs);
    
    data = sum(imfs, 2);
    means = zeros(nImf-1, 1);
    idx = zeros(1, nImf);
    
    % determine trend indexes using the rule of thumb: "standardized mean of partial reconstructed signal departs significantly from zero"
    for i = 1:nImf-1
        means(i,1) = mean(sum(imfs(:,1:i), 2));
    end
    means = means./sum(means,1);
    
    for i = 1:nImf-1
        if(means(i,1) >= threshold)
            idx(i,1) = i;
        end
    end
    idx(nImf,1) = nImf;
    
    idx_min = min(nonzeros(idx));
    
    trend = sum(imfs(:, idx_min:nImf), 2);
    filtered = sum(imfs(:, 1:idx_min-1), 2);
    
    if (plots)
        fontSize = 16;
        fontName = 'Helvetica';
        lineWidth = 1;
        
        h = figure;
        subplot(1, 2, 1, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
        hold on;
        scatter(1:nImf-1, means, 'filled');
        xlim([0 nImf]);
        ylabel('Mean (standardized)');
        title('Mean of partial reconstructed signal');
        grid on;
        hold off;
        
        subplot(1, 2, 2, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
        hold on;
        plot(data, 'LineWidth', lineWidth/2, 'Color', [0.7 0.7 0.7]);
        plot(trend, 'LineWidth', lineWidth, 'Color', 'r');
        xlim([1 nObs]);
        legend('Data', 'Trend');
        title('Original data and trend');
        hold off;
    end

end
