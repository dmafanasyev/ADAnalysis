function [ h ] = plot_imfs( data, imfs, residue )
    %PLOT_IMFS Plot original signal (time-series), its IMFs and the residue (optionaly)
    %   
    %   Input:
    %       data - source time-series or signal
    %       imfs - intrinsic mode functions of data
    %   Output:
    %       h - axis handler
    %
    %   Copyright (c) 2015 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0     2015.09.21: initial version

    plotResidue = 0;

    if(nargin >= 3 && ~isempty(residue))
        plotResidue = 1;
        imfs = [imfs, residue];
    end;

    nObs = size(imfs, 1);
    nImf = size(imfs, 2);
    odd = rem(nImf, 2);
    
    nCols = 2;
    nRows = round(nImf/nCols)+1;

    fontName = 'Helvetica';
    fontSize = 16;

    %h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    h = figure;
    
    % plot original signal
    subplot(nRows, nCols, 1:2, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on');
    hold on;
    plot(data, 'LineWidth', 1.2);
    ylabel('Signal');
    xlim([1 nObs]);
    hold off;
    
    % plot IMFs and residue (optionaly)
    for i=1:nImf;
        if(i < nRows)
            axisInd = 2*i-1;
        else
            axisInd = 2*i-(nImf+odd);
        end
         axisInd = axisInd+2;
    
        subplot(nRows, nCols, axisInd, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on');
        hold on;
        plot(imfs(:, i), 'LineWidth', 1.2);
        xlim([1 nObs]);
        
        if(i == nImf && plotResidue)
            ylabel('Residue');
        else
            ylabel(['IMF_{', num2str(i), '}']);
        end        

        if (i ~= nRows-1 && i ~=nImf)
            set(gca, 'xticklabel', '');
        end
        hold off;
    end
end