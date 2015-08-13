function [ h ] = plot_imfs( imfs, residue )
    %PLOT_IMFS Plot IMFs
    %   Input:
    %       imfs - IMFs of decomposition
    %   Output:
    %       h - axis handler

    plotResidue = 0;
    
    if(nargin >= 2 && ~isempty(residue))
        plotResidue = 1;
    end;
    
    nObs = size(imfs, 1);
    nImf = size(imfs, 2);


figure;
for i=1:nImf;
    subplot(nImf, 1, i);
    
    if (i == 1)
        title('IMFs');
    end
    
    hold on;
        grid on;
        plot(imfs(:, i));
        ylabel(['IMF_', num2str(i)]);
    
    if (i < nImf);
        
        set(gca, 'xticklabel', '');
        hold off;
    else
        hold on;
        grid on;
        plot(imfs(:, i));
        ylabel('Trend');
        set(gca,'XTick', xTickDates);
        datetick('x', xTickDateFormat, 'keepticks');
        hold off;
    end
    
    xlim([1 nObs]);
end