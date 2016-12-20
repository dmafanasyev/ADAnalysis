function [ tdr, tir, instant] = tdir( y, x, yImf, xImf, bootstrap, plots, factorsName, periodn, alpha )
    %TDIR Calculate time-dependent intrinsic regression (TDIR), bootstraping median of the coefficients estimations and instantaneous characteristics (using Hilbert transformation).
    %   
    %   TODO: write docs
    %
    %   Copyright (c) 2016 by Dmitriy O. Afanasyev
    %   Versions:
    %   v0.1 2014.11.25: initial version
    %   v0.2 2015.06.28: model summary on origin data level, visualization and results print (latex) 
    %   v0.3 2016.08.03: indicative significance function for coefficient avereging; JB, LBQ and ARCH tests
    %   v0.3 2016.09.10: HAC covariance estimation (Newey-West form)
    %
    
    if (nargin < 4)
        error('Dependent vector Y, matrix of independent value X and IMFs for both must be specified');
    end
    
    [nObs, nImfs] = size(yImf);
    nFactors = size(xImf, 2);
    
    if (nargin < 5)
        bootstrap = 0;
    end
    if (nargin < 6)
        plots = 0;
    end
    if (nargin < 7)
        factorsName = cell(1,nFactors);
        for k = 1:nFactors
            factorsName{k} = num2str(k);
        end
    end
    if (nargin < 8)
        periodn = 1;
    end
    if (nargin < 9)
        alpha = 0.05;
    end
    
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
        if(bootstrap == 1)
            bootnum = 10000;
        else
            bootnum = bootstrap;
        end
    end
    
    % switch off some unnecessary warnings
    warning('off', 'stats:jbtest:PTooBig');
    warning('off', 'stats:jbtest:PTooSmall');
    
    [instant.y.h, instant.y.ampl, instant.y.phase, instant.y.omega, instant.y.period] = hilbert_transform(yImf, 1);
    
    for k = 1:nFactors
        [instant.x{1,k}.h, instant.x{1,k}.ampl, instant.x{1,k}.phase, instant.x{1,k}.omega, instant.x{1,k}.period] = hilbert_transform(xImf{1,k}, 1);
    end
    
    periodX = nan(1,nFactors);

    for j = 1:nImfs
    %for j = nImfs:-1:1
        
        %intercept = (j == nImfs);
        intercept = (1==1);
        
        % Time-independent regression
        x_j = nan(nObs, nFactors);
        for k = 1:nFactors
            % if some columns is zero (the factor is excluded from regression) then switch off warrning message about the rank deficient
            if(xImf{1,k}(:,j) == zeros(nObs,1))
                warning('off', 'stats:LinearModel:RankDefDesignMat');
            end
            x_j(:,k) = xImf{1,k}(:,j);
        end
        
        mdl = fitlm(x_j, yImf(:,j), 'Intercept', intercept);
        tir.b{j}(1,1:nFactors+intercept) = mdl.Coefficients.Estimate;
        try
            covEst = hac(mdl, 'type', 'HAC', 'weights', 'BT', 'bandwidth', floor(4*(nObs/100)^(2/9)) + 1, 'display', 'off');% Newey-West (Bartlett kernel)
        catch exception
            covEst = mdl.CoefficientCovariance;
        end
        tir.bse{j}(1,1:nFactors+intercept) = sqrt(diag(covEst))';
        tir.bts{j}(1,1:nFactors+intercept) = tir.b{j}(1,1:nFactors+intercept)./tir.bse{j}(1,1:nFactors+intercept);
        tir.bp{j}(1,1:nFactors+intercept) = 2*(1-tcdf(abs(tir.bts{j}(1,1:nFactors+intercept)), mdl.DFE));
        ci = abs(tinv(alpha/2, mdl.DFE)) * tir.bse{j}(1,1:nFactors+intercept);
        tir.blo{j}(1,1:nFactors+intercept) = tir.b{j}(1,1:nFactors+intercept) - ci;
        tir.bup{j}(1,1:nFactors+intercept) = tir.b{j}(1,1:nFactors+intercept) + ci;
        tir.var(1,j) = var(mdl.Residuals.Raw);
        tir.R2(1,j) = mdl.Rsquared.Ordinary;
        tir.F(1,j) = (tir.R2(1,j)/nFactors) / ((1 - tir.R2(1,j))/(nObs - nFactors - 1));
        tir.Fp(1,j) = 1 - fcdf(tir.F(1,j), nFactors, nObs-nFactors-1);
        [~, tir.JBp(1,j), tir.JB(1,j)] = jbtest(mdl.Residuals.Standardized); % H0: normal distribution with unknown mean and variance
        lbqLags = min(12, nObs-1);
        [~, tir.LBQp(1,j), tir.LBQ(1,j)] = lbqtest(mdl.Residuals.Standardized, 'Lags', lbqLags, 'dof', max(1, lbqLags - nFactors));% H0: about no jointly autocorrelation
        [~, tir.ARCHp(1,j), tir.ARCH(1,j)] = archtest(mdl.Residuals.Standardized, 'Lags', 1);% H0: no conditional heteroscedasticity (ARCH effect)
        
        % Time-dependent regression
        for t = 1:nObs
            
            for k = 1:nFactors
                if(xImf{1,k}(:,j) == zeros(nObs,1))
                    periodX(1,k) = 0;
                    warning('off', 'stats:LinearModel:RankDefDesignMat');
                else
                    periodX(1,k) = instant.x{1,k}.period(t,j);
                end
            end
            
            halfWinSizeMin = round(periodn*max([instant.y.period(t, j), periodX])/2);
            %halfWinSizeMin = round(periodn*instant.y.period(t, j)/2);
            
            nx_win = 0;
            halfWinSize = 0;
            while(nx_win <= nFactors+intercept)
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
            
            mdl = fitlm(x_win, y_win, 'Intercept', intercept);
            warning('on', 'stats:LinearModel:RankDefDesignMat');
            
%             tdr.b{j}(t,1:nFactors+intercept) = mdl.Coefficients.Estimate;
%             tdr.bse{j}(t,1:nFactors+intercept) = mdl.Coefficients.SE;
%             tdr.bts{j}(t,1:nFactors+intercept) = mdl.Coefficients.tStat;
%             tdr.bp{j}(t,1:nFactors+intercept) = mdl.Coefficients.pValue;
%             ci = coefCI(mdl, alpha);
%             tdr.blo{j}(t,1:nFactors+intercept) = ci(1:nFactors+intercept,1);
%             tdr.bup{j}(t,1:nFactors+intercept) = ci(1:nFactors+intercept,2);
            
            tdr.b{j}(t,1:nFactors+intercept) = mdl.Coefficients.Estimate';
            try
                covEst = hac(mdl, 'type', 'HAC', 'weights', 'BT', 'bandwidth', floor(4*(nx_win/100)^(2/9)) + 1, 'display', 'off');% Newey-West (Bartlett kernel)
            catch exception
                covEst = mdl.CoefficientCovariance;
            end
            tdr.bse{j}(t,1:nFactors+intercept) = sqrt(diag(covEst))';
            tdr.bts{j}(t,1:nFactors+intercept) = tdr.b{j}(t,1:nFactors+intercept)./tdr.bse{j}(t,1:nFactors+intercept);
            tdr.bp{j}(t,1:nFactors+intercept) = 2*(1-tcdf(abs(tdr.bts{j}(t,1:nFactors+intercept)), mdl.DFE));
            ci = abs(tinv(alpha/2, mdl.DFE)) * tdr.bse{j}(t,1:nFactors+intercept);
            tdr.blo{j}(t,1:nFactors+intercept) = tdr.b{j}(t,1:nFactors+intercept) - ci;
            tdr.bup{j}(t,1:nFactors+intercept) = tdr.b{j}(t,1:nFactors+intercept) + ci;
            tdr.var(t,j) = var(mdl.Residuals.Raw);
            tdr.R2(t,j) = mdl.Rsquared.Ordinary;
            tdr.F(t,j) = (tdr.R2(t,j)/nFactors) / ((1 - tdr.R2(t,j))/(nx_win - nFactors - 1));
            tdr.Fp(t,j) = 1 - fcdf(tdr.F(t,j), nFactors, nx_win-nFactors-1);
            [~, tdr.JBp(t,j), tdr.JB(t,j)] = jbtest(mdl.Residuals.Standardized); % H0: normal distribution with unknown mean and variance
            lbqLags = min(12, nx_win-1);
            [~, tdr.LBQp(t,j), tdr.LBQ(t,j)] = lbqtest(mdl.Residuals.Standardized, 'Lags', lbqLags, 'dof', max(1, lbqLags - nFactors));% H0: no jointly autocorrelation
            [~, tdr.ARCHp(t,j), tdr.ARCH(t,j)] = archtest(mdl.Residuals.Standardized, 'Lags', 1);% H0: no conditional heteroscedasticity (ARCH effect)
            
%             [tdr.b{j}(t,1:nFactors+1), stats] = robustfit(x_win, y_win, 'ols', 1, 'on');
%             tdr.bse{j}(t,1:nFactors+1) = stats.se';
%             tdr.bts{j}(t,1:nFactors+1) = stats.t';
%             tdr.bp{j}(t,1:nFactors+1) = stats.p';
%             b_ci = abs(tinv(alpha/2, nx_win-nFactors-1)) * tdr.bse{j}(t,:);
%             tdr.blo{j}(t,1:nFactors+1) = tdr.b{j}(t,:) - b_ci;
%             tdr.bup{j}(t,1:nFactors+1) = tdr.b{j}(t,:) + b_ci;
%             tdr.R2(t,j) =  max(0, 1 - (sum(stats.resid.^2, 1) / sum((y_win - mean(y_win)).^2, 1)));
%             tdr.F(t,j) = (tdr.R2(t,j)/nFactors) / ((1 - tdr.R2(t,j))/(nx_win - nFactors - 1));
%             tdr.Fp(t,j) = 1 - fcdf(tdr.F(t,j), nFactors, nx_win-nFactors-1);
            
            % compute fitted values (only for model that significantly differ from intercept model)
            if(round(tdr.Fp(t,j), 2) <= alpha)
                if(intercept)
                    tdr.fitImf(t,j) = tdr.b{j}(t,:)*[1; x_curr'];
                else
                    tdr.fitImf(t,j) = tdr.b{j}(t,:)*x_curr';
                end
            else
                tdr.fitImf(t,j) = 0;
            end
        end
        
        % set all non significant coefficients estimations to zero
        for k = 1:nFactors+intercept
            sign = (round(tdr.bp{j}(:,k), 2) <= alpha);
            tdr.b{j}(:,k) = sign.*tdr.b{j}(:,k);
        end
        
        if(bootstrap > 0)
            % bootstrap median, CI and p-value of t-test for beta coefficients
            for k = 1:nFactors+intercept
                bootstat = bootstrp(bootnum, @median, tdr.b{j}(:,k));
                tdr.mean.b{j}(k,1)  = mean(bootstat);
                
                bootstat = bootstrp(bootnum, @median, tdr.bse{j}(:,k));
                tdr.mean.bse{j}(k,1)  = mean(bootstat);
                
                tdr.mean.bts{j}(k,1) = tdr.mean.b{j}(k,1)/(tdr.mean.bse{j}(k,1));
                tdr.mean.bp{j}(k,1) = 1 - tcdf(abs(tdr.mean.bts{j}(k,1)), bootnum-1);
            end
        else
            % simple mean, CI and p-value of t-test for beta coefficients
            for k = 1:nFactors+intercept
                tdr.mean.b{j}(k,1)  = mean(tdr.b{j}(:,k));
                tdr.mean.bse{j}(k,1)  = mean(tdr.bse{j}(:,k));
                tdr.mean.bts{j}(k,1) = tdr.mean.b{j}(k,1)/(tdr.mean.bse{j}(k,1));
                tdr.mean.bp{j}(k,1) = 1 - tcdf(abs(tdr.mean.bts{j}(k,1)), bootnum-1);
            end
        end
    end
    
    % summary for origin data level
    tdr.summary.fit = sum(tdr.fitImf, 2);
    tdr.summary.R2 = 1 - sum((y(:)-tdr.summary.fit(:)).^2)/sum((y(:)-mean(y(:))).^2);
    [tdr.summary.MAE, tdr.summary.MAPE, tdr.summary.MSE, tdr.summary.RMSE] = mean_errors(y, tdr.summary.fit);
    
    % estimation of the source data liner regression coefficients based on
    % TDIR coefficients
%     tdr.rec.b = nan(1,nFactors+1);
%     x_j = nan(nObs, nFactors+1);
%     
%     for j =1:nImfs
%         intercept = (j == nImfs);
%         
%         if(intercept)
%             x_j(:,1) = ones(nObs,1);
%         else
%             x_j(:,1) = zeros(nObs,1);
%         end
%         
%         for k = 1:nFactors
%             x_j(:,k+1) = xImf{1,k}(:,j);
%         end
%         
%         tdr.rec.b  = tdr.rec.b + tdr.b{j}.*(x_j./x);
%     end
    
    % plot results
    if(plots>0)
         fontName = 'Helvetica';
         fontSize = 16;
         
         % plot model parameters dynamic (without intercept)
         axisInd = 0;
         figure;
         for j = 1:nImfs
             
             %intercept = (j == nImfs);
             intercept = (1==1);
             
             for k = 1:nFactors
                 
                 signDummy = (round(tdr.bp{j}(:,k+intercept), 2) <= alpha);
                 
                 axisInd = axisInd + 1;
                 
                 subplot(nImfs, nFactors, axisInd, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on');
                 hold on;
                 box on;
                 grid on;
                 plot(tdr.b{j}(:,k+intercept).*signDummy);
                 scatter((1:nObs)'.*(~signDummy), tdr.b{j}(:,k+intercept).*(~signDummy), 2, 'r', '.');
                 
                 plot(tdr.blo{j}(:,k+intercept), 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--');
                 plot(tdr.bup{j}(:,k+intercept), 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--');
                 
                 %line([1 n], [tdr.mean.b{imfInd}(factorInd,1) tdr.mean.b{imfInd}(factorInd,1)], 'Color', 'r');
                 %line([1 nObs], [tdr.mean.blo{imfInd}(factorInd,1) tdr.mean.blo{imfInd}(factorInd,1)], 'Color', 'r', 'LineStyle', '-.');
                 %line([1 nObs], [tdr.mean.bup{imfInd}(factorInd,1) tdr.mean.bup{imfInd}(factorInd,1)], 'Color', 'r', 'LineStyle', '-.');
                 xlim([1 nObs]);
                 if(j==1)
                     title(['\beta_{', factorsName{k}, '}']);
                 end
                 if(k==1)
                     ylabel(['CIMF_{', num2str(j),'}']);
                 end
                 
                 hold off;
             end
         end
         
         %plot source vs fitted data
         figure;
         subplot(1, 1, 1,  'FontName', fontName, 'FontSize', fontSize, 'Box', 'on');
         hold on;
         plot((1:nObs)', y, (1:nObs)', tdr.summary.fit);
         xlim([1 nObs]);
         %title('Source and fitted data');
         legend('Source data', ['Fitted data (R^2 = ', num2str(tdr.summary.R2, '%.2f'), ', MAPE = ', num2str(tdr.summary.MAPE, '%.1f'), '%)']);
         hold off;
         
         % print table with model parameters
         bColTitles = cell(1,nFactors+1);
         bRowTitles = cell(1,nImfs);
         bTable = cell(nImfs,nFactors+1);
         formatSpec = '%.3f';
         for j = 1:nImfs
             
             %intercept = (j == nImfs);
             intercept = (1 == 1);
             
             for k = 1:nFactors+1
                 if(intercept)
                     bTable{j, k} = [num2str(tdr.mean.b{j}(k,1), formatSpec), ...
                         '\\textsuperscript{', pvalue_to_asterisks(tdr.mean.bp{j}(k,1)), '}', ...
                         ' (', num2str(tdr.mean.bse{j}(k,1), formatSpec), ')'];
                 else
                     if(k == 1)
                         bTable{j, k} = '--';
                     else
                         bTable{j, k} = [num2str(tdr.mean.b{j}(k-1,1), formatSpec), ...
                             '\\textsuperscript{', pvalue_to_asterisks(tdr.mean.bp{j}(k-1,1)), '}', ...
                             ' (', num2str(tdr.mean.bse{j}(k-1,1), formatSpec), ')'];
                     end
                     
                 end
                 
                 if(j == 1)
                     if(k == 1)
                         factorName = '0';
                     else
                         factorName = factorsName{k-1};
                     end
                     bColTitles{k} = ['$\\beta_{', factorName, '}$'];
                 end
             end
             
             %     bTable{imfInd, factorInd+1} = num2str(tdr.mean.R2(imfInd)*100, '%.1f');
             %     bColTitles{factorInd+1} = ['$R^2$'];
             
             bRowTitles{j} = num2str(j);
         end
         print_latex_table( bTable, bColTitles, bRowTitles);
    end
     
    % switch on some warnings
    warning('on', 'stats:jbtest:PTooBig');
    warning('on', 'stats:jbtest:PTooSmall');
end
