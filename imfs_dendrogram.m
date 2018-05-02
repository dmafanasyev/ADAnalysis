function [ h ] = imfs_dendrogram( imfSet, method, metric, threshold, titles, lables )
    %IMFS_DENDROGRAM Plot dendrogram of IMFs.
    %   
    %   Input:
    %       imfSet - the 1 x K cell array of IMFs (K - the number of time-series),
    %                each cell imfSet{k} is the N x I matrix,
    %                N - the number of observations, I - the number of IMFs
    %       method - the method to measure the distance between clusters,
    %                see documentation for "linkage" Matlab function
    %       metric - the method to measure the distance between objects in clusters,
    %                see documentation for "pdist" Matlab function
    %       threshold - the value of threshold for unique colors in the dendrogram plot,
    %                   see documentation for "dendrogram" Matlab function, default is 0 (only one color)
    %       titles - the 1 x K cell array of time-series titles that used as subplot titles, default is "Y_k"
    %       lables - the 1 x K cell array of IMFs lables for each time-series from imfSet,
    %                each cell lables{k} is the 1 x I cell array of lables, default is "i"
    %
    %   Output:
    %       h - the graphical object handler of the IMFs dendrogram
    %
    %   Copyright (c) 2016-2018 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0  2016.01.05: initial version
    %       1.1  2018.05.02: documentation and default titles added
    %

    setSize = size(imfSet,2);
    
    if(nargin < 4)
        threshold = 0;
    end
    
    if(nargin < 5)
        for i = 1:setSize
            titles{i} = ['Y_', num2str(i)];
        end
    end
    
    if(nargin < 6)
        for i = 1:setSize
            imfNum = size(imfSet{i},2);
            for j = 1:imfNum
                lables{i}{j} = num2str(j);
            end
        end
    end
    
    if(setSize == 1)
        useWide = 0;
        axRowN = 1;
        axColN = axRowN;
    elseif(setSize == 2)
        useWide = 0;
        axRowN = 1;
        axColN = setSize;
    elseif(setSize == 3)
        useWide = 1;
        axRowN = 2;
        axColN = 2;
    elseif(rem(setSize, 2) == 0)
        useWide = 0;
        axRowN = setSize/2;
        axColN = axRowN;
    else
        useWide = 1;
        axRowN = (setSize+1)/2;
        axColN = axRowN - 1;
    end
    
    fontSize = 8;
    fontName = 'Helvetica';
    lineWidth = 0.5;

    h = figure;
    for i = 1:setSize
        dist = pdist(imfSet{i}', metric);
        tree = linkage(dist, method);
        if(useWide)
            if(i==1)
                ax = subplot(axRowN, axColN, 1:axColN);
            else
                ax = subplot(axRowN, axColN, i+axColN-1);
            end
        else
            ax = subplot(axRowN, axColN, i);
        end
        
        hd = dendrogram(tree,  'ColorThreshold', threshold, 'Labels', lables{i}, 'Reorder', optimalleaforder(tree, dist));
        
%         if((ischar(threshold) && strcmp(threshold, 'default')) || (~ischar(threshold) && threshold > 0))
%             hd = dendrogram(tree,  'ColorThreshold', threshold, 'Labels', lables, 'Reorder', optimalleaforder(tree, dist));
%         else
%             hd = dendrogram(tree,  'Labels', lables, 'Reorder', optimalleaforder(tree, dist));
%         end
        set(hd, 'LineWidth', lineWidth);
        if(threshold == 0)
            set(hd, 'Color', [0, 0.447, 0.741]);
        end
        set(ax, 'FontSize', fontSize);
        set(ax, 'FontName', fontName);
        set(ax, 'Box', 'on');
        set(ax, 'YGrid', 'on');
        title(ax, titles{i});
    end

end