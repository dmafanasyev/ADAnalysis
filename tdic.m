function [ tdc, tlc, tic, instant] = tdic( imfs1, imfs2, bootstrap, plots, prints, subranges, periodn, alpha)
    %TDIC Calculate time-dependent intrinsic correlation.
    %   
    %   For the multi-scale correlation analysis purpose are calculated next types of the IMFs pairs correlations:
    %   (1) simple liner Pearson correlation,
    %   (2) time-dependent local correlation with fixed sliding window (TDLC), mean (optionaly - bootstraped median) of TDLC (see Papadimitriou et al., 2006),
    %   (3) time-dependent intrinsic correlation (TDIC), mean (optionaly - bootstraped median) of TDIC and instantaneous characteristics, using Hilbert transformation (see Chen et al., 2010; Afanasyev et al., 2015).
    %   
    %   NOTE: Unlike Chen et al. (2010) here used Hilbert transformation for
    %   instantaneous periods calculation instead general zero crossing
    %   method. TDIC calculated only for one instantaneous period number (see periodn input parameter),
    %   not for range from minimumu to maximum period (= half of time-series length).
    %   
    %   Input:
    %       imfs1, imfs2 - IMFs matrixes
    %       bootstrap - boolean flag (0 or 1) for bootstraping of TDIC median, default is 0 - simple mean
    %       plots - boolean flag (0 or 1) for graphics plot, default is 0
    %       prints - boolean flag (0 or 1) for print result in LaTeX table, default is 0
    %       subranges - the set of time subranges for correlation averaging, first and last index per row: [1 100; 101 110; 111 150]
    %       periodn - number of instantaneous periods for including into sliding window, default is 1
    %       alpha - significance level, default is 0.05
    %
    %   Output:
    %       tdc - time-dependent intrinsic correlation (TDIC):
    %           tdc.corr - correlation
    %           tdc.prob - probability (significance level) of correlation
    %           tdc.lo - lower boundary of correlation
    %           tdc.up - upper boundary of correlation
    %           tdc.mean - mean level of time-dependent correlation (based on simple mean or bootstraped median), structure is the same as tdc
    %       tlc - time-dependent local correlation with fixed sliding window (TDLC) equal to the max mean period of IMFs pair:
    %           structure is the same as variable tdc
    %       tic - liner Pearson correlation (time-independent), strucutre is same the same as tdc.mean
    %       instant - instantaneous characteristics:
    %           instant.imfs<N> - instantaneous characteristics of IMFs<N> (see method hilbert_transform()):
    %               instant.imfs<N>.h - analytic signal
    %               instant.imfs<N>.ampl- instantaneous amplitude
    %               instant.imfs<N>.phase - instantaneous phase
    %               instant.imfs<N>.omega - instantaneous frequency
    %               instant.imfs<N>.period - instantaneous period
    %	
    %   This function is partial Octave compatible (except bootstrap support).
    %   It requires next packages to be installed (they will be loaded automaticaly in the properly order):
    %       nan (http://octave.sourceforge.net/nan/index.html)
    %       statistics (http://octave.sourceforge.net/statistics/)
    %       signal (http://octave.sourceforge.net/signal/)
    %	
    %   References:
    %       Chen, N., Wu, Z., Huang, N., 2010. The time-dependent intrinsic correlation based on the empirical mode decomposition. Advances in Adaptive Data Analysis 2 (2), 223-265.
    %       Afanasyev, D., Fedorova, E., Popov, V., 2015. Fine structure of the price-demand relationship in the electricity market: multi-scale correlation analysis. Energy Economics 51, 215-226.
    %       Papadimitriou, S., Sun, J., Yu, P., 2006. Local correlation tracking in time series. ICDM, 456-465.
    %
    %   Copyright (c) 2014-2015 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0  2014.04.05: initial version
    %       1.1  2015.04.13: added time-dependent local correlation, refactor internal routines
    %       1.11 2015.12.25: added compatibility with Octave environment
    %       1.2 2016.12.11: the code moved from long-titled time_dependent_intrinsic_corr.m routine
    %       1.3 2017.01.25: added averaging by subperiods, print of result to LaTeX table
    %   
    
    if (nargin < 2)
        error('Two IMFs matrix must be specified');
    end
    if(size(imfs1, 1) ~= size(imfs2, 1))
        error('The length of both IMFs matrix must be the equal');
    end;
    
    if(size(imfs1, 2) ~= size(imfs2, 2))
        error('The number of columns at IMFs matrix must be the equal');
    end;
    
    if (nargin < 3)
        bootstrap = 0;
    end
    if (nargin < 4)
        plots = 0;
    end
    if (nargin < 5)
        prints = 0;
    end
    if (nargin < 6)
        subranges = [];
    end
    if (nargin < 7)
        periodn = 1;
    end
    if (nargin < 8)
        alpha = 0.05;
    end
    
    if(periodn <= 0 || fix(periodn) ~= periodn)
        periodn = 1;
        warning('Number of instantaneous periods for including into sliding window must be positive integer number. Set it to default value 1.');
    end
    
    if(is_run_octave())
        % unload required packages first to provide the properly loading order then
        pkg unload nan;
        pkg unload statistics;
        pkg unload signal;
        
        pkg load nan;
        pkg load statistics;
        pkg load signal;
        
        if(bootstrap == 1)
          warning('Bootstrap is not supported under Octave environment now. The simple mean mode will be used.');
          bootstrap = 0;
        end
    end
    
    [nObs, nCols] = size(imfs1);
    
    [instant.imfs1.h, instant.imfs1.ampl, instant.imfs1.phase, instant.imfs1.omega, instant.imfs1.period] = hilbert_transform(imfs1, 1);
    [instant.imfs2.h, instant.imfs2.ampl, instant.imfs2.phase, instant.imfs2.omega, instant.imfs2.period] = hilbert_transform(imfs2, 1);
    
    tic.corr = nan(nCols, 1);
    tic.prob = nan(nCols, 1);
    tic.lo = nan(nCols, 1);
    tic.up = nan(nCols, 1);
    
    tlc.corr = nan(nObs, nCols);
    tlc.prob = nan(nObs, nCols);
    tlc.lo = nan(nObs, nCols);
    tlc.up = nan(nObs, nCols);
    
    tlc.mean.corr = nan(nCols, 1);
    tlc.mean.prob = nan(nCols, 1);
    tlc.mean.lo = nan(nCols, 1);
    tlc.mean.up = nan(nCols, 1);
    
    tdc.corr = nan(nObs, nCols);
    tdc.prob = nan(nObs, nCols);
    tdc.lo = nan(nObs, nCols);
    tdc.up = nan(nObs, nCols);
    
    tdc.mean.corr = nan(nCols, 1);
    tdc.mean.prob = nan(nCols, 1);
    tdc.mean.lo = nan(nCols, 1);
    tdc.mean.up = nan(nCols, 1);
    
    period_mean = nan(nCols, 2);
    
    nSub = size(subranges, 1);
    if(nSub > 0)
        tdc.sub = cell(nSub,1);
        for s = 1:nSub
            tdc.sub{s}.corr = nan(subranges(s,2)-subranges(s,1)+1, nCols);
            tdc.sub{s}.mean.corr = nan(nCols, 1);
            tdc.sub{s}.mean.prob = nan(nCols, 1);
            tdc.sub{s}.mean.lo = nan(nCols, 1);
            tdc.sub{s}.mean.up = nan(nCols, 1);
        end
    end
    
    for j = 1:nCols
      % time independent correlation
      tic = pearson_full_corr_int(tic, imfs1, imfs2, j, alpha);
      
      period_mean(j, 1) = period_zero_cross(imfs1(:,j));
      period_mean(j, 2) = period_zero_cross(imfs2(:,j));
      halfWinSizeMean = periodn*round(max(period_mean(j, 1), period_mean(j, 2))/2);
      %halfWinSizeMean = periodn*round(max(period_zero_cross(imfs1(:,j)), period_zero_cross(imfs2(:,j)))/2);
      
      % time dependent local and intrinsic correlation
      for t = 1:nObs
          tlc = pearson_window_corr_int(tlc, imfs1, imfs2, t, j, nObs, halfWinSizeMean, alpha);
          
          halfWinSize = periodn*round(max(instant.imfs1.period(t, j), instant.imfs2.period(t, j))/2);
          tdc = pearson_window_corr_int(tdc, imfs1, imfs2, t, j, nObs, halfWinSize, alpha);
      end
      
      tlc = mean_corr_int(tlc, j, bootstrap, alpha);
      tdc = mean_corr_int(tdc, j, bootstrap, alpha);
      
      % averaging for subranges
      for s = 1:nSub
          tdc.sub{s}.corr(:,j) = tdc.corr(subranges(s,1):subranges(s,2), j);
          tdc.sub{s} = mean_corr_int(tdc.sub{s}, j, bootstrap, alpha);
      end
    end
    
    if(plots>0)
        %plot_corr_dynamic_int(tlc, nCols,  nObs, 'Time-dependent local correlations');
        plot_corr_dynamic_int(tdc, nCols,  nObs, 'Time-dependent intrinsic correlations');
        %plot_corr_mean_int( tic,tlc, tdc, nCols);
    end
    
    if(prints>0)
        print_result_int(tic, tlc, tdc, period_mean, nCols);
    end
end

%%%%% Internal routines %%%%%
function [ r, p, rlo, rup ] = pearson_corr_int( x, y, alpha )
    % Calculate Pearson correlation
    if(is_run_octave())
        [r, p, rlo, rup] = corrcoef([x y], 'Mode', 'Pearson', 'alpha', alpha);
    else
        [r, p, rlo, rup] = corrcoef(x, y, 'alpha', alpha);
    end
    
    r = r(1,2);
    p = p(1,2);
    rlo = rlo(1,2);
    rup = rup(1,2);
end

function [ corr_obj ] = pearson_full_corr_int( corr_obj, data1, data2, j, alpha )
    % Calculate Pearson correlation for full sample
    [corr_obj.corr(j), corr_obj.prob(j), corr_obj.lo(j), corr_obj.up(j)] = pearson_corr_int(data1(:, j), data2(:, j), alpha);

    if(isnan(corr_obj.lo(j)))
        corr_obj.lo(j) = corr_obj.corr(j);
    end
    if(isnan(corr_obj.up(j)))
      corr_obj.up(j) = corr_obj.corr(j);
    end
end

function [ corr_obj ] = pearson_window_corr_int( corr_obj, data1, data2, t, j, nObs, halfWinSize, alpha )
    % Calculate Pearson correlation for given window
    ts = t - halfWinSize;
    if(ts < 1)
        ts = 1;
    end

    te = t + halfWinSize;
    if(te > nObs)
        te = nObs;
    end
	
    [corr_obj.corr(t, j), corr_obj.prob(t, j), corr_obj.lo(t, j), corr_obj.up(t, j)] = pearson_corr_int(data1(ts:te, j), data2(ts:te, j), alpha);

    if(isnan(corr_obj.lo(t, j)))
        corr_obj.lo(t, j) = corr_obj.corr(t, j);
    end
    if(isnan(corr_obj.up(t, j)))
		  corr_obj.up(t, j) = corr_obj.corr(t, j);
    end
end

function [ corr_obj ] = mean_corr_int( corr_obj, j, bootstrap, alpha )
    % Calculate mean/median value of time-dependent Pearson correlation
    if(bootstrap > 0)
      % bootstrap median, CI and p-value of t-test
      bootnum = 10000;
      bootstat = bootstrp(bootnum, @median, corr_obj.corr(:, j));
      ci = bootci(bootnum, {@median, corr_obj.corr(:, j)}, 'alpha', alpha);

      corr_obj.mean.corr(j, 1) = mean(bootstat);
      corr_obj.mean.lo(j, 1) = ci(1,1);
      corr_obj.mean.up(j, 1) = ci(2,1);
      [~, corr_obj.mean.prob(j, 1)] = ttest(bootstat, 0, 'alpha', alpha);
    else
      % simple mean, CI and p-value of t-test
      corr_obj.mean.corr(j, 1) = mean(corr_obj.corr(:, j));
      [~, corr_obj.mean.prob(1,j), ci] = ttest(corr_obj.corr(:, j), 0, 'alpha', alpha);
      corr_obj.mean.lo(j, 1) = ci(1);
      corr_obj.mean.up(j, 1) = ci(2);
    end
end

function [] = plot_corr_dynamic_int( corr_obj, nCols,  nObs, title_str )
    % Plot graph of time-dependent Pearson correlation
    figure;
    for i=1:nCols
      subplot(nCols, 1, i, 'FontSize', 8);
      plot(corr_obj.corr(:, i));
      
      if (i == 1)
          %title(title_str);
      end
      
      xlim([1 nObs]);
      ylim([-1 1]);
      
      if (i < nCols)
          set(gca, 'xticklabel', '');
      end
    end
end

function [] = plot_corr_mean_int( tic,tlc, tdc, nCols)
    % Plot mean values of 3 types of time-dependent Pearson correlation
    x = 1:nCols;

    figure;
    subplot(3,1,1, 'FontSize', 8);
    hold on;
    grid on;
    box on;
    errorbar(x, tic.corr, tic.corr-tic.lo, tic.up-tic.corr, 'LineStyle', 'none', 'Marker', 'x', 'LineWidth', 0.5);
    xlim([0 nCols+1]);
    ylim([-1 1]);
    title('Time-independent correlations');
    hold off;

    subplot(3,1,2, 'FontSize', 8);
    hold on;
    grid on;
    errorbar(x, tlc.mean.corr, tlc.mean.corr-tlc.mean.lo, tlc.mean.up-tlc.mean.corr, 'LineStyle', 'none', 'Marker', 'x', 'LineWidth', 0.5);
    xlim([0 nCols+1]);
    ylim([-1 1]);
    title('Mean of the time-dependent local correlations');
    hold off;

    subplot(3,1,3, 'FontSize', 8);
    hold on;
    grid on;
    errorbar(x, tdc.mean.corr, tdc.mean.corr-tdc.mean.lo, tdc.mean.up-tdc.mean.corr, 'LineStyle', 'none', 'Marker', 'x', 'LineWidth', 0.5);
    xlim([0 nCols+1]);
    ylim([-1 1]);
    title('Mean of the time-dependent intrinsic correlations');
    hold off;
end

function [] = print_result_int(tic, tlc, tdc, period_mean, nRows)
    nSub = size(tdc.sub, 1);
    
    tbl = cell(nRows, 6 + nSub);
    rFormatSpec = '%.2f';
    
    for i = 1:nRows
        tbl{i,1} = int2str(i);
        tbl{i,2} = int2str(period_mean(i,1));
        tbl{i,3} = int2str(period_mean(i,2));
        tbl{i,4} = [num2str(tic.corr(i,1), rFormatSpec), '\\textsuperscript{', pvalue_to_asterisks(tic.prob(i,1)), '}'];
        tbl{i,5} = [num2str(tlc.mean.corr(i,1), rFormatSpec), '\\textsuperscript{', pvalue_to_asterisks(tlc.mean.prob(i,1)), '}'];
        tbl{i,6} = [num2str(tdc.mean.corr(i,1), rFormatSpec), '\\textsuperscript{', pvalue_to_asterisks(tdc.mean.prob(i,1)), '}'];
        for s = 1:nSub
            tbl{i,6+s} = [num2str(tdc.sub{s,1}.mean.corr(i,1), rFormatSpec), '\\textsuperscript{', pvalue_to_asterisks(tdc.sub{s,1}.mean.prob(i,1)), '}'];
        end
    end
    
    colTitles = {'$i$', '$\\overline{T}_{1}$', '$\\overline{T}_{2}$', '$r$', '$\\overline{r}$', '$\\overline{\\rho}$'};
    
    for s = 1:nSub
        colTitles{1,6+s} = ['$\\overline{\\rho}^{sub ', num2str(s), '}$'];
    end
    
    tblTitle = 'Multi-scale correlation analysis';
    tblNotes = ['$\\overline{T}$ -- the mean period of IMF, $r$ -- coefficient of linear correlation (Pearson), ', ...
                    '$\\overline{r}$ -- the mean of time-dependent local correlation, $\\overline{\\rho}$ -- the mean of time-dependent intrinsic correlation. ', ...
                    'Significance levels: \\textsuperscript{***} -- 1\\\%, \\textsuperscript{**} -- 5\\\%, \\textsuperscript{*} -- 10\\\%.'];
    
    print_latex_table(tbl, colTitles, {}, '', tblTitle, tblNotes);
end