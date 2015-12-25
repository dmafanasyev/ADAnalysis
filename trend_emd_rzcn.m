function [ filtered, trend, idx, h, rzcn ] = trend_emd_rzcn( imfs, plots, alpha )
%TREND_EMD_ENERGY Summary of this function goes here
%   Detailed explanation goes here
%   Moghtader, A., Borgnat, P., Flandrin, P., 2011. Trend filtering: empirical mode decomposition versus l1 and Hodrick-Prescott. Advances in Adaptive Data Analysis 3 (1 and 2), 41–61.

     if(nargin == 0 || ~ismatrix(imfs))
        error('Input data must be non empty matrix');
    end
    
    if (nargin < 2)
        plots = 0;
    end
    if (nargin < 3)
        alpha = 0.05;
    end
    
    alphaErr = 'Use only alpha equal to 0.1, 0.08, 0.05, 0.03 or 0.01';
    
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
    else
        error(alphaErr);
    end
    
    % if time is columns then transpose matrix (I assume that the number of observations is greater than number of IMFs always)
    if (size(imfs, 2) > size(imfs, 1))
        imfs = imfs.';
    end
    
    nImf = size(imfs, 2);
    
    period = zeros(1, nImf);
    rzcn = zeros(1, nImf);
    idx = zeros(1, nImf);
    
    % determine trend indexes using statistical criteria
    for i=1:nImf
        [period(i, 1), ~, ~, indzer, numzercur] = period_zero_cross(imfs(:, i));
        
        % if current "IMF" is a last (residue) or number of zero crossings less than 2 then force set numzercur = 0
        if(i == nImf || size(indzer, 2) < 2)
            numzercur = 0;
        end
        
        if(i > 1)
            rzcn(i) = numzerlast / numzercur;
            if (isinf(rzcn(i)) || isnan(rzcn(i)))
                rzcn(i) = rzcntr;
            end
            
            if(rzcn(i) <= rzcntl || rzcn(i) >= rzcntr)
                idx(i) = i;
            end
        end
        
        numzerlast = numzercur;
    end
    
    idx_min = min(nonzeros(idx));
    
    trend = sum(imfs(:, idx_min:nImf), 2);
    filtered = sum(imfs(:, 1:idx_min-1), 2);
    
    h = [];
    if (plots)
        fontSize = 16;
        fontName = 'Helvetica';
        
        h = figure;
        subplot(1, 1, 1, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
        hold on;
        scatter(1:nImf, log2(rzcn), 'filled');
        line([0 nImf+1], [log2(rzcntl) log2(rzcntl)], 'LineStyle', '--', 'Color', 'r');
        line([0 nImf+1], [log2(rzcntr) log2(rzcntr)], 'LineStyle', '--', 'Color', 'r');
        xlim([0 nImf+1]);
        ylabel('log_2(RZCN)');
        legend('Ratio of the zero-crossing numbers', [num2str((1-alpha)*100), '% confidence interval']);
        hold off;
    end

end
