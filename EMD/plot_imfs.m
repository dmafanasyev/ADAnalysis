function [ h ] = plot_imfs( data, imfs, hasResidue )
    %PLOT_IMFS Plot original time-series (signal) and its IMFs
    %   
    %   Input:
    %       data - source time-series (signal)
    %       imfs - intrinsic mode functions of time-series (signal)
    %       hasResidue - boolean flag, if 1 (default) then last column of imfs is residue
    %   Output:
    %       h - axis handler
    %
    %   Copyright (c) 2015 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0     2015.09.21: initial version
    %       1.1     2016.03.13: the last column of imfs interpreted as residue (see hasResidue flag)

    if(nargin < 3)
        hasResidue = 1;
    end

    nObs = size(imfs, 1);
    nImf = size(imfs, 2);
    odd = rem(nImf, 2);
    
    nCols = 2;
    nRows = round(nImf/nCols)+1+1;

    fontName = 'Helvetica';
    fontSize = 8;
    lineWidth = 0.5;

    %h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    h = figure;
    
    % plot original signal
    subplot(nRows, nCols, 1:2*nCols, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on');
    hold on;
    plot(data, 'LineWidth', lineWidth);
    ylabel('x_t');
    xlim([1 nObs]);
    %set(gca, 'xticklabel', '');
    hold off;
    
    % plot IMFs and residue (optionaly)
    for i=1:nImf
        if(i < nRows-1)
            axisInd = 2*i-1;
        else
            axisInd = 2*i-(nImf+odd);
        end
         axisInd = axisInd+4;
    
        subplot(nRows, nCols, axisInd, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on');
        hold on;
        plot(imfs(:, i), 'LineWidth', lineWidth);
        xlim([1 nObs]);
        
        if(i == nImf && hasResidue)
            ylabel('R');
        else
            ylabel(['IMF_{', num2str(i), '}']);
        end        

        if (i ~= nRows-2 && i ~=nImf)
            set(gca, 'xticklabel', '');
        end
        hold off;
    end
end