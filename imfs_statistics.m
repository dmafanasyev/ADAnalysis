function [ stat ] = imfs_statistics( data, imfs, doprint, alpha)
    %IMFS_STATISTICS Compute descriptive statistics of IMFs
    %   Input:
    %       data            - source data vector
    %       imfs            - IMFs matrix
    %       doprint        - flag for LaTeX table print, default is 0 (no print)
    %       alpha           - confidence level of statistical tests
    %   
    %   Output:
    %       stat     - descriptive statistics of IMFs:
    %                       mean period,
    %                       variation,
    %                       variation as % of source data vector variation,
    %                       Pearson correlation with source data vector,
    %                       p-value of Pearson correlation,
    %                       log2 of IMF energy,
    %                       log2 of fGn energy confidence level (critical value for given alpha),
    %                       indicator of fGn test passing.
    %   
    %   Copyright (c) 2015 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0  2016.01.05: initial version
    %       1.1  2016.07.31: added the test of H0 about statistical differences of IMFs from fractional Gaussian noise (fGn)
    %       
    
    if(nargin < 2 || isempty(data) || isempty(imfs))
        error('Both source data vector and IMFs matrix must be not empty');
    end
    if(nargin < 3)
        doprint = 0;
    end
    if(nargin < 4)
        alpha = 0.05;
    end
    
    nImfs = size(imfs, 2);
    
    stat = zeros(nImfs, 8);
    
    for i =1:nImfs
        stat(i,1) = round(period_zero_cross(imfs(:, i)), 0);
        [corr, corrp] = corrcoef(imfs(:,i), data);
        stat(i,4) = round(corr(1,2), 2);
        stat(i,5) = round(corrp(1,2), 2);
    end
    
    stat(:,2) = var(imfs, 0, 1)';
    stat(:,3) = round(100*(stat(:,2)./var(data)), 1);
    
    [~, stat(:,6), stat(:,7)] = imfs_sign_test(imfs, 0, 0, alpha);
    stat(:,6:7) = log2(stat(:,6:7));
    stat(:,8) = (stat(:,6) > stat(:,7)).*1;
    
    if(doprint)
        rowTitles = cell(1,nImfs);
        for i = 1:nImfs
            rowTitles{1,i} = num2str(i);
        end
        colTitles = {'IMF \\#', 'Mean period', 'Variation', 'Variation, %%', 'Pearson correlation', 'p-value', 'Energy (log\\_2)', ['Critical value (log\\_2, ', num2str(alpha*100), '\\%%-level)'], 'fGn test passed'};
        tblTitle = 'IMFs statistics';
        formatSpec = {'%.1f', '', '%.1f', '%.2f', '%.2f', '%.2f', '%.2f', '%i'};
        
        print_latex_table([stat], colTitles, rowTitles, formatSpec, tblTitle);
    end
end
