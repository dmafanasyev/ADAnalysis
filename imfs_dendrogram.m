function [ h ] = imfs_dendrogram( imfSet, method, metric, threshold, titles, lables )
%IMFS_DENDROGRAM Summary of this function goes here
%   Detailed explanation goes here

    setSize = size(imfSet,2);
    
    if(nargin < 4)
        threshold = 0;
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
    
    fontSize = 16;
    fontName = 'Helvetica';

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
        set(hd, 'LineWidth', 1.5);
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