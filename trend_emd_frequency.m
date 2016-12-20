function [ trend, filtered, idx ] = trend_emd_frequency( imfs )
%TREND_EMD_FREQUENCY Summary of this function goes here
%   Detailed explanation goes here
%   Afanasyev, D., Fedorova, E., Popov, V., 2015. Fine structure of the price-demand relationship in the electricity market: multi-scale correlation analysis. Energy Economics 51, 215-226.

    if(nargin == 0 || ~ismatrix(imfs))
        error('Input data must be non empty matrix');
    end
    
    % if time is columns then transpose matrix (I assume that the number of observations is greater than number of IMFs always)
    if (size(imfs, 2) > size(imfs, 1))
        imfs = imfs.';
    end
    
    nImf = size(imfs, 2);
    
    idx = zeros(1, nImf);
    
    % determine trend indexes using low-frequency criteria
    for i=1:nImf
        if(i >= nImf/2+1)
            idx(i) = i;
        end
    end
    
    idx_min = min(nonzeros(idx));
    
    trend = sum(imfs(:, idx_min:nImf), 2);
    filtered = sum(imfs(:, 1:idx_min-1), 2);
end
