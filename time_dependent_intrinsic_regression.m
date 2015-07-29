function [ tdr, tlr, instant] = time_dependent_intrinsic_regression( y, x, yImf, xImf, bootstrap, plots, periodn, alpha )
    %TIME_DEPEND_CORR Calculate liner Pearson correlation, time-dependent intrinsic correlation (TDIC), median (optionaly - bootstraping median) of TDIC and instantaneous characteristics (using Hilbert transformation).
    %   
    %   Input:
    %       data1 & data2 - IMFs matrixes
    %       bootstrap - boolean flag (0 or 1) for bootstraping of TDIC median, default is 0
    %       plots - boolean flag (0 or 1) for graphics plot, default is 0
    %       periodn - number of instantaneous periods for including in sliding window, default is 1
    %       alpha - confidence level, default is 0.05 (5%)
    %
    %   Output:
    %       tdc - time-dependent correlation:
    %           tdc.corr - correlation
    %           tdc.prob - probability (significance level) of correlation
    %           tdc.lo - lower boundary of correlation
    %           tdc.up - upper boundary of correlation
    %           tdc.mean - mean level of time-dependent correlation, structure is the same as tdc
    %       tic - liner Pearson correlation (time-independent), strucutre is the same as tdc.mean
    %       instant - instantaneous characteristics:
    %           instant.d<N> - instantaneous characteristics of data<N> (see method hilbert_transform()):
    %               instant.d<N>.h - analytic signal
    %               instant.d<N>.ampl- instantaneous amplitude
    %               instant.d<N>.phase - instantaneous phase
    %               instant.d<N>.omega - instantaneous frequency
    %               instant.d<N>.period - instantaneous period
    %
    %   Reference(s):
    %   [1] Chen, N., Wu, Z., Huang, N., 2010. The time-dependent intrinsic correlation based on the empirical mode decomposition. Advances in Adaptive Data Analysis 2 (2), 223–265.
    %   [2] Afanasyev, D., Fedorova, E., Popov, V., 2014. Fine structure of the price-demand relationship in the electricity market: multi-scale
    %   correlation analysis. Preprint submitted to Energy Economics, Sep 2014.
    %
    %   Copyright (c) 2014-2015 by Dmitriy O. Afanasyev
    %   Versions:
    %   v0.1 2014.11.25: initial version
    %   v0.2 2015.06.28: 
    
    if (nargin < 4)
        error('Dependent vector Y, matrix of independent value X and IMFs for both must be specified');
    end
    if (nargin < 5)
        bootstrap = 0;
    end
    if (nargin < 6)
        plots = 0;
    end
    if (nargin < 7)
        periodn = 1;
    end
    if (nargin < 8)
        alpha = 0.05;
    end
    
    [nObs, nImfs] = size(yImf);
    nFactors = size(xImf, 2);
    
    % check number of observations and IMFs at Y and X
    for k = 1:nFactors
        if(nObs ~= size(xImf{1, k}, 1))
            error(['The number of observations at Y IMFs matrix and all of the X IMFs matrix must be the equal ', num2str(nObs)]);
        end
        
        if(nImfs ~= size(xImf{1, k}, 2))
            error(['The number of modes  of Y IMFs matrix and all of the X IMFs matrix must be the equal ', num2str(nImfs)]);
        end
    end
    
    if(bootstrap > 0)
        bootnum = 5000;
    end
    
    [instant.y.h, instant.y.ampl, instant.y.phase, instant.y.omega, instant.y.period] = hilbert_transform(yImf, 1);
    
    for k = 1:nFactors
        [instant.x{1,k}.h, instant.x{1,k}.ampl, instant.x{1,k}.phase, instant.x{1,k}.omega, instant.x{1,k}.period] = hilbert_transform(xImf{1,k}, 1);
    end
    
%     tdr.corr = nan(nObs, nCols);
%     tdr.prob = nan(nObs, nCols);
%     tdr.lo = nan(nObs, nCols);
%     tdr.up = nan(nObs, nCols);
%     
%     tdr.mean.corr = nan(nCols, 1);
%     tdr.mean.prob = nan(nCols, 1);
%     tdr.mean.lo = nan(nCols, 1);
%     tdr.mean.up = nan(nCols, 1);
    
    periodX = nan(1,nFactors);
    periodMeanX = nan(1,nFactors);

    for j = 1:nImfs
        
        %for k = 1:nFactors
        %    periodMeanX(1,k) = period_zero_cross(xImf{1,k}(:,j));
        %end
        %halfWinSizeMean = round(periodn*max([period_zero_cross(yImf(:,j)), periodMeanX])/2);
        halfWinSizeMean = round(periodn*period_zero_cross(yImf(:,j))/2);
        tlr.win(1,j) = halfWinSizeMean+1;
        
        % time dependent liner regression
        for t = 1:nObs
            
            %for k = 1:nFactors
            %    periodX(1,k) = instant.x{1,k}.period(t,j);
            %end
            
            %halfWinSizeMin = round(periodn*max([instant.y.period(t, j), periodX])/2);
            halfWinSizeMin = round(periodn*instant.y.period(t, j)/2);
            
            %halfWinSize = min(halfWinSize, halfWinSizeMean);
            
%             ts = t - halfWinSize;
%             if(ts<1)
%                 ts = 1;
%             end
%             
%             te = t + halfWinSize;
%             if(te>nObs)
%                 te = nObs;
%             end
%             
%             nx_win = te-ts+1;
            nx_win = 0;
            halfWinSize = 0;
            while(nx_win <= nFactors+1)
                halfWinSize = halfWinSize + halfWinSizeMin;
                ts = t - halfWinSize;
                if(ts<1)
                    ts = 1;
                end
                
                te = t + halfWinSize;
                if(te>nObs)
                    te = nObs;
                end
                
                nx_win = te-ts+1;
            end
            
            x_win = nan(nx_win,nFactors);
            x_curr = nan(1,nFactors);
            for k = 1:nFactors
                x_win(1:nx_win,k) = xImf{1,k}(ts:te,j);
                x_curr(1,k) = xImf{1,k}(t,j);
            end
            
            y_win = yImf(ts:te,j);
            
            tdr.win(t,j) = nx_win;
            
            [tdr.b{j}(t,1:nFactors+1), stats] = robustfit(x_win, y_win, 'ols');
            tdr.se{j}(t,1:nFactors+1) = stats.se';
            tdr.ts{j}(t,1:nFactors+1) = stats.t';
            tdr.tsp{j}(t,1:nFactors+1) = stats.p';
            b_ci = abs(tinv(alpha/2, nx_win-nFactors-1)) * tdr.se{j}(1,1:end);
            tdr.blo{j}(t,1:nFactors+1) = tdr.b{j}(1,1:end) - b_ci;
            tdr.bup{j}(t,1:nFactors+1) = tdr.b{j}(1,1:end) + b_ci;
            tdr.R2(t,j) =  1 - (sum(stats.resid.^2, 1) / sum((y_win - mean(y_win)).^2, 1));
            tdr.F(t,j) = (tdr.R2(t,j)/nFactors) / ((1 - tdr.R2(t,j))/(nx_win - nFactors - 1));
            tdr.Fp(t,j) = 1 - fcdf(tdr.F(t,j), nFactors, nx_win-nFactors-1);
            %[~, tdr.LBQp{t,j}] = lbqtest(stats.resid, 'lags', 1:round(log(nx_win)));% m ~ log(T), see Tsay, R. S. Analysis of Financial Time Series. 2nd Ed. Hoboken, NJ: John Wiley & Sons, Inc., 2005.
            %[~, tdr.ARCHp(t,j)] = archtest(stats.resid);
            tdr.fitImf(t,j) = tdr.b{j}(t,:)*[1; x_curr'];
        end
        
        if(bootstrap > 0)
            % bootstrap median, CI and p-value of t-test for beta coefficients and R2
            for k = 1:nFactors+1
                bootstat = bootstrp(bootnum, @median, tdr.b{j}(:,k));
                tdr.mean.b{j}(k,1)  = mean(bootstat);
                tdr.mean.bse{j}(k,1)  = std(bootstat);
                %[~, tdr.mean.bp{j}(k,1)] = ttest(bootstat, 0, 'alpha', alpha);
                %tdr.mean.bts{j}(k,1) = tdr.mean.b{j}(k,1)/(tdr.mean.bse{j}(k,1)/sqrt(bootnum));
                tdr.mean.bts{j}(k,1) = tdr.mean.b{j}(k,1)/(tdr.mean.bse{j}(k,1));
                tdr.mean.bp{j}(k,1) = 1 - tcdf(abs(tdr.mean.bts{j}(k,1)), bootnum-1);
                
                ci = bootci(bootnum, {@median, tdr.b{j}(:,k)}, 'alpha', alpha);
                tdr.mean.blo{j}(k,1) = ci(1,1);
                tdr.mean.bup{j}(k,1) = ci(2,1);
            end
            
            bootstat = bootstrp(bootnum, @median, tdr.R2(:,j));
            tdr.mean.R2(1,j) = mean(bootstat);
        else
            % simple mean, CI and p-value of t-test for beta coefficients and R2
            for k = 1:nFactors+1
                tdr.mean.b{j}(k,1)  = median(tdr.b{j}(:,k));
                [~, tdr.mean.bp{j}(k,1), ci] = ttest(tdr.b{j}(:,k), 0, 'alpha', alpha);
                tdr.mean.blo{j}(k,1) = ci(1);
                tdr.mean.bup{j}(k,1) = ci(2);
            end
            
            tdr.mean.R2(1,j) = mean(tdr.R2(:,j));
            [~, tdr.mean.R2p(1,j), ci] = ttest(tdr.R2(:,j), 0, 'alpha', alpha);
            tdr.mean.R2lo(1,j) = ci(1);
            tdr.mean.R2up(1,j) = ci(2);
        end
    end
    
    tdr.summary.fit = sum(tdr.fitImf, 2);
    tdr.summary.R2 = 1 - sum((y(:)-tdr.summary.fit(:)).^2)/sum((y(:)-mean(y(:))).^2);
    [tdr.summary.MAE, tdr.summary.MAPE, tdr.summary.MSE, tdr.summary.RMSE] = mean_errors(y, tdr.summary.fit);
    
%     if(plots>0)
%         figure;
%         for k=1:nImfs;
%             subplot(nImfs, 1, k);
%             plot(tdr.corr(:, k));
%             
%             if (k == 1)
%                 title('Time dependent intrinsic correlations');
%             end
%             
%             xlim([1 nObs]);
%             ylim([-1 1]);
%             
%             if (k < nImfs);
%                 set(gca, 'xticklabel', '');
%             end
%         end
%         
%         x = 1:nImfs;
%         
%         figure;
%         subplot(2,1,1);
%         hold on;
%         grid on;
%         box on;
%         errorbar(x, tic.corr, tic.corr-tic.lo, tic.up-tic.corr, 'LineStyle', 'none', 'Marker', 'x', 'LineWidth', 2);
%         xlim([0 nImfs+1]);
%         ylim([-1 1]);
%         title('Time independent correlations');
%         hold off;
%         
%         subplot(2,1,2);
%         hold on;
%         grid on;
%         errorbar(x, tdr.mean.corr, tdr.mean.corr-tdr.mean.lo, tdr.mean.up-tdr.mean.corr, 'LineStyle', 'none', 'Marker', 'x', 'LineWidth', 2);
%         xlim([0 nImfs+1]);
%         ylim([-1 1]);
%         title('Median of the time dependent correlations');
%         hold off;
%         
%     end
end
