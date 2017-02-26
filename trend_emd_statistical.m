function [ trend, filtered, idx, energy, confidence ] = trend_emd_statistical( imfs, plots, prints, alpha )
%TREND_EMD_STATISTICAL Detrending and denoising of data using statistical criteria
%   Flandrin, P., Goncalves, P., Rilling, G., 2004. Detrending and denoising with empirical mode decomposition. Proceedings of Eusipco, Wien (Austria), 1581-1584.
    
    if(nargin == 0 || ~ismatrix(imfs))
        error('Input data must be non empty matrix');
    end
    
    if (nargin < 2)
        plots = 0;
    end
    if (nargin < 3)
        prints = 0;
    end
    if (nargin < 4)
        alpha = 0.05;
    end
    
    % determine trend indexes using statistical criteria
    [idx, energy, confidence] = imfs_sign_test(imfs, plots, prints, alpha);
    
    trend = sum(imfs(:,idx),2);
    filtered = sum(imfs,2) - trend;
    
end
