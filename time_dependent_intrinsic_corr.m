function [ tdc, tlc, tic, instant] = time_dependent_intrinsic_corr( data1, data2, bootstrap, plots, periodn, alpha )
    %TIME_DEPENDENT_INTRINSIC_CORR Calculate time-dependent intrinsic correlation
    %   
    %   For the multi-scale correlation analysis purpose are calculated next types of the IMFs pairs correlations:
    %   (1) simple liner Pearson correlation,
    %   (2) time-dependent local correlation with fixed sliding window (TDLC), mean (optionaly - bootstraped median) of TDLC (see Papadimitriou et al., 2006),
    %   (3) time-dependent intrinsic correlation (TDIC), mean (optionaly - bootstraped median) of TDIC and instantaneous characteristics, using Hilbert transformation (see Chen et al., 2010; Afanasyev et al., 2014).
    %   
    %   NOTE: Unlike Chen et al. (2010) here used Hilbert transformation for
    %   instantaneous periods calculation instead general zero crossing
    %   method. TDIC calculated only for minimum instantaneous period, not
    %   for range from minimumu to maximum period (= half of time-series
    %   length).
    %   
    %   Input:
    %       data1, data2 - IMFs matrixes
    %       bootstrap - boolean flag (0 or 1) for bootstraping of TDIC median, default is 0 - simple mean
    %       plots - boolean flag (0 or 1) for graphics plot, default is 0
    %       periodn - number of instantaneous periods for including to sliding window, default is 1
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
    %           instant.d<N> - instantaneous characteristics of data<N> (see method hilbert_transform()):
    %               instant.d<N>.h - analytic signal
    %               instant.d<N>.ampl- instantaneous amplitude
    %               instant.d<N>.phase - instantaneous phase
    %               instant.d<N>.omega - instantaneous frequency
    %               instant.d<N>.period - instantaneous period
    %
    %   References:
    %       Chen, N., Wu, Z., Huang, N., 2010. The time-dependent intrinsic correlation based on the empirical mode decomposition. Advances in Adaptive Data Analysis 2 (2), 223–265.
    %       Afanasyev, D., Fedorova, E., Popov, V., 2015. Fine structure of the price-demand relationship in the electricity market: multi-scale correlation analysis. Energy Economics 51, 215-226.
    %       Papadimitriou, S., Sun, J., Yu, P., 2006. Local correlation tracking in time series. ICDM, 456–465.
    %
    %   Copyright (c) 2014-2015 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0 2014.04.05: initial version
    %       1.1 2015.04.13: added time-dependent local correlation, refactor internal routines
    %   
    
    if (nargin < 2)
        error('Two IMFs matrix must be specified');
    end
    if(size(data1, 1) ~= size(data2, 1))
        error('The length of both IMFs matrix must be the equal');
    end;
    
    if(size(data1, 2) ~= size(data2, 2))
        error('The number of columns at IMFs matrix must be the equal');
    end;
    
    if (nargin < 3)
        bootstrap = 0;
    end
    if (nargin < 4)
        plots = 0;
    end
    if (nargin < 5)
        periodn = 1;
    end
    if (nargin < 6)
        alpha = 0.05;
    end
    
    [nObs, nCols] = size(data1);
    
    [instant.d1.h, instant.d1.ampl, instant.d1.phase, instant.d1.omega, instant.d1.period] = hilbert_transform(data1, 1);
    [instant.d2.h, instant.d2.ampl, instant.d2.phase, instant.d2.omega, instant.d2.period] = hilbert_transform(data2, 1);
    
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
    
    for j = 1:nCols
        % time independent correlation
        [r, p, rlo, rup] = corrcoef(data1(:,j), data2(:,j), 'alpha', alpha);
        tic.corr(j, 1) = r(1,2);
        tic.prob(j, 1) = p(1,2);
        tic.lo(j, 1) = rlo(1,2);
        tic.up(j, 1) = rup(1,2);
        
        halfWinSizeMean = periodn*round(max(period_zero_cross(data1(:,j)), period_zero_cross(data2(:,j)))/2);
        
        % time dependent intrinsic correlation
        for t = 1:nObs
            tlc = pearson_window_corr_int(tlc, data1, data2, t, j, nObs, halfWinSizeMean, alpha);
            
            halfWinSize = periodn*round(max(instant.d1.period(t, j), instant.d2.period(t, j))/2);
            tdc = pearson_window_corr_int(tdc, data1, data2, t, j, nObs, halfWinSize, alpha);
        end
        
        tlc = mean_corr_int(tlc, j, bootstrap, alpha);
        tdc = mean_corr_int(tdc, j, bootstrap, alpha);
    end
    
    if(plots>0)
        plot_corr_dynamic_int(tlc, nCols,  nObs, 'Time-dependent local correlations');
        plot_corr_dynamic_int(tdc, nCols,  nObs, 'Time-dependent intrinsic correlations');
        plot_corr_mean_int( tic,tlc, tdc, nCols);
    end
end

%%%%% Internal routines %%%%%
function [ corr_obj ] = pearson_window_corr_int( corr_obj, data1, data2, t, j, nObs, halfWinSize, alpha )
    % Calculate Pearson correlation for given window
    ts = t - halfWinSize;
    if(ts<1)
        ts = 1;
    end

    te = t + halfWinSize;
    if(te>nObs)
        te = nObs;
    end

    [r, p, rlo, rup] = corrcoef(data1(ts:te, j), data2(ts:te, j), 'alpha', alpha);

    corr_obj.corr(t, j) = r(1,2);
    corr_obj.prob(t, j) = p(1,2);
    if(isnan(rlo(1,2)))
        corr_obj.lo(t, j) = r(1,2);
    else
        corr_obj.lo(t, j) = rlo(1,2);
    end
    if(isnan(rup(1,2)))
        corr_obj.up(t, j) = r(1,2);
    else
        corr_obj.up(t, j) = rup(1,2);
    end
end

function [ corr_obj ] = mean_corr_int( corr_obj, j, bootstrap,alpha )
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
        [~, tdr.mean.prob(1,j), ci] = ttest(corr_obj.corr(:, j), 0, 'alpha', alpha);
        corr_obj.mean.lo(j, 1) = ci(1);
        corr_obj.mean.up(j, 1) = ci(2);
    end
end

function [] = plot_corr_dynamic_int( corr_obj, nCols,  nObs, title_str )
    % Plot graph of time-dependent Pearson correlation
    figure;
    for i=1:nCols;
        subplot(nCols, 1, i);
        plot(corr_obj.corr(:, i));
        
        if (i == 1)
            title(title_str);
        end
        
        xlim([1 nObs]);
        ylim([-1 1]);
        
        if (i < nCols);
            set(gca, 'xticklabel', '');
        end
    end
end

function [] = plot_corr_mean_int( tic,tlc, tdc, nCols)
    % Plot mean values of 3 types of time-dependent Pearson correlation
    x = 1:nCols;

    figure;
    subplot(3,1,1);
    hold on;
    grid on;
    box on;
    errorbar(x, tic.corr, tic.corr-tic.lo, tic.up-tic.corr, 'LineStyle', 'none', 'Marker', 'x', 'LineWidth', 2);
    xlim([0 nCols+1]);
    ylim([-1 1]);
    title('Time-independent correlations');
    hold off;

    subplot(3,1,2);
    hold on;
    grid on;
    errorbar(x, tlc.mean.corr, tlc.mean.corr-tlc.mean.lo, tlc.mean.up-tlc.mean.corr, 'LineStyle', 'none', 'Marker', 'x', 'LineWidth', 2);
    xlim([0 nCols+1]);
    ylim([-1 1]);
    title('Mean of the time-dependent local correlations');
    hold off;

    subplot(3,1,3);
    hold on;
    grid on;
    errorbar(x, tdc.mean.corr, tdc.mean.corr-tdc.mean.lo, tdc.mean.up-tdc.mean.corr, 'LineStyle', 'none', 'Marker', 'x', 'LineWidth', 2);
    xlim([0 nCols+1]);
    ylim([-1 1]);
    title('Mean of the time-dependent intrinsic correlations');
    hold off;
end