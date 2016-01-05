function [ stat ] = imfs_statistics( data, imfs, doPrint, formatSpec)
    %IMFS_STATISTICS Compute descriptive statistics of IMFs
    %   Input:
    %       data            - source data vector
    %       imfs            - IMFs matrix
    %       doPrint        - flag for LaTeX table print, default is 0 (no print)
    %       formatSpec  - specification of number output format
    %   
    %   Output:
    %       stat     - descriptive statistics of IMFs:
    %                       mean period,
    %                       variation,
    %                       variation as % of source data vector variation,
    %                       Pearson correlation with source data vector,
    %                       p-value of Pearson correlation.
    %   
    %   Copyright (c) 2015 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0  2016.01.05: initial version
    %       
    
    if(nargin < 2 || isempty(data) || isempty(imfs))
        error('Both source data vector and IMFs matrix must be not empty');
    end
    
    if(nargin < 3)
        doPrint = 0;
    end
    
    if(nargin < 4)
        formatSpec = '%f';
    end
    
    nImfs = size(imfs, 2);
    
    stat = zeros(nImfs, 5);
    
    for i =1:nImfs
        stat(i,1) = period_zero_cross(imfs(:, i));
        [corr, corrp] = corrcoef(imfs(:,i), data);
        stat(i,4) = corr(1,2);
        stat(i,5) = corrp(1,2);
    end
    
    stat(:,2) = var(imfs, 0, 1)';
    stat(:,3) = 100*(stat(:,2)./var(data));
    
    if(doPrint)
        rowTitles = cell(1,nImfs);
        for i = 1:nImfs
            rowTitles{1,i} = num2str(i);
        end
        colTitles = {'Mean period', 'Variation', 'Variation (% from total)', 'Pearson correlation', 'Pearson correlation p-value'};
        tblTitle = 'IMFs statistics';
        
        print_latex_table(stat, colTitles, rowTitles, formatSpec, tblTitle);
    end
end
