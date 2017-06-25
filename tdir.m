function [ tdr, tir, instant] = tdir(yImf, xImf, bootstrap, showResults, factorsName, periodn, alpha, vifMax )
    %TDIR Calculate time-dependent intrinsic regression (TDIR), bootstraping median of the coefficients estimations and instantaneous characteristics (using Hilbert transformation).
    %   
    %   TODO: docs
    %
    %   Copyright (c) 2014-2017 by Dmitriy O. Afanasyev
    %   Versions:
    %   v0.1 2014.11.25: initial version
    %   v0.2 2015.06.28: model summary on original data level, visualization and results print (LaTeX)
    %   v0.3 2016.08.03: indicative significance function for coefficient avereging; JB, LBQ and ARCH tests for residuals
    %   v0.4 2016.09.10: HAC covariance estimation (Newey-West form)
    %   v0.5 2017.06.04: multicollinearity (VIF) and endogeneity analysis, AD and t tests for residuals, nMAE and nRMSE
    %   v0.6 2017.06.09: one step-ahead forecast
    %
    
    if (nargin < 2)
        error('IMFs for dependent and independent values must be specified');
    end
    
    [nObs, nImfs] = size(yImf);
    nFactors = size(xImf, 2);
    
    if (nargin < 3)
        bootstrap = 0;
    end
    if (nargin < 4)
        showResults = 0;
    end
    if (nargin < 5)
        factorsName = cell(1,nFactors);
        for k = 1:nFactors
            factorsName{k} = num2str(k);
        end
    end
    if (nargin < 6)
        periodn = 1;
    end
    if (nargin < 7)
        alpha = 0.05;
    end
    if (nargin < 8)
        vifMax = 10;
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
    
    y = sum(yImf,2);
    
    [instant.y.h, instant.y.ampl, instant.y.phase, instant.y.omega, instant.y.period] = hilbert_transform(yImf, 1);
    
    for k = 1:nFactors
        [instant.x{1,k}.h, instant.x{1,k}.ampl, instant.x{1,k}.phase, instant.x{1,k}.omega, instant.x{1,k}.period] = hilbert_transform(xImf{1,k}, 1);
    end
    
    periodX = nan(1,nFactors);

    for j = 1:nImfs
        %intercept = (j == nImfs);
        intercept = (1==1);
        
        % Time-independent regression
        x_j = nan(nObs, nFactors);
        for k = 1:nFactors
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
        tir.R2(1,j) = mdl.Rsquared.Ordinary;
        [~,  ~, ~, ~, ~, tir.nMAE(1,j), tir.nRMSE(1,j)] = mean_errors(yImf(:,j), mdl.Fitted);
        tir.F(1,j) = (tir.R2(1,j)/nFactors) / ((1 - tir.R2(1,j))/(nObs - nFactors - 1));
        tir.Fp(1,j) = 1 - fcdf(tir.F(1,j), nFactors, nObs-nFactors-1);
        tir.resMean(1,j) = mean(mdl.Residuals.Raw);
        tir.resVar(1,j) = sum(mdl.Residuals.Raw.^2,1)/(nObs - nFactors - 1);
        [~, tir.resADp(1,j), tir.resAD(1,j)] = adtest(mdl.Residuals.Standardized, 'Distribution', makedist('normal','mu',0,'sigma',1), 'Asymptotic', (nObs > 120)); % H0: normal distribution with mean equal to 0 and variance 1
        [~, tir.resTp(1,j)] = ttest(mdl.Residuals.Standardized, 0); % H0: normal distribution with mean equal to 0 and unknown variance
        [~, tir.resJBp(1,j), tir.resJB(1,j)] = jbtest(mdl.Residuals.Standardized, alpha); % H0: normal distribution with unknown mean and variance
        lbqLags = min(12, nObs-1);
        [~, tir.resLBQp(1,j), tir.resLBQ(1,j)] = lbqtest(mdl.Residuals.Standardized, 'Lags', lbqLags, 'dof', max(1, lbqLags - nFactors));% H0: no jointly autocorrelation
        if(intercept)
            [tir.resDWp(1,j), tir.resDW(1,j)] = dwtest(mdl.Residuals.Standardized, [ones(nObs,1) x_j]);% H0: residuals are uncorrelated;
        else
            [tir.resDWp(1,j), tir.resDW(1,j)] = dwtest(mdl.Residuals.Standardized, x_j);% H0: residuals are uncorrelated;
        end
        [~, tir.ARCHp(1,j), tir.ARCH(1,j)] = archtest(mdl.Residuals.Standardized, 'Lags', 1);% H0: no conditional heteroscedasticity (ARCH effect)
        [r_tmp, p_tmp] = corrcoef([mdl.Residuals.Raw x_j]);
        tir.resXR{j}(1,1:nFactors) = r_tmp(1,2:end);
        tir.resXRp{j}(1,1:nFactors) = p_tmp(1,2:end);
        tir.VIF{j}(1,1:nFactors) = diag(inv(corrcoef(x_j)))';
        
        % Time-dependent regression
        for t = 1:nObs
            for k = 1:nFactors
                periodX(1,k) = instant.x{1,k}.period(t,j);
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
            x_next = nan(1,nFactors);
            for k = 1:nFactors
                x_win(1:nx_win,k) = xImf{1,k}(ts:te,j);
                x_curr(1,k) = xImf{1,k}(t,j);
                if(t ~= nObs)
                    x_next(1,k) = xImf{1,k}(t+1,j);
                end
            end
            
            y_win = yImf(ts:te,j);
            
            mdl = fitlm(x_win, y_win, 'Intercept', intercept);
            
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
            tdr.R2(t,j) = mdl.Rsquared.Ordinary;
            [~,  ~, ~, ~, ~, tdr.nMAE(t,j), tdr.nRMSE(t,j)] = mean_errors(y_win, mdl.Fitted);
            tdr.F(t,j) = (tdr.R2(t,j)/nFactors) / ((1 - tdr.R2(t,j))/(nx_win - nFactors - 1));
            tdr.Fp(t,j) = 1 - fcdf(tdr.F(t,j), nFactors, nx_win-nFactors-1);
            tdr.resMean(t,j) = mean(mdl.Residuals.Raw);
            tdr.resVar(t,j) = sum(mdl.Residuals.Raw.^2,1)/(nx_win - nFactors - 1);
            [~, tdr.resADp(t,j), tdr.resAD(t,j)] = adtest(mdl.Residuals.Standardized, 'Distribution', makedist('normal','mu',0,'sigma',1), 'Asymptotic', (nx_win > 120)); % H0: normal distribution with mean equal to 0 and variance 1
            [~, tdr.resTp(t,j)] = ttest(mdl.Residuals.Standardized, 0); % H0: normal distribution with mean equal to 0 and unknown variance
            [~, tdr.resJBp(t,j), tdr.resJB(t,j)] = jbtest(mdl.Residuals.Standardized, alpha); % H0: normal distribution with unknown mean and variance
            lbqLags = min(12, nx_win-1);
            [~, tdr.resLBQp(t,j), tdr.resLBQ(t,j)] = lbqtest(mdl.Residuals.Standardized, 'Lags', lbqLags, 'dof', max(1, lbqLags - nFactors));% H0: no jointly autocorrelation
            if(intercept)
                [tdr.resDWp(t,j), tdr.resDW(t,j)] = dwtest(mdl.Residuals.Standardized, [ones(nx_win,1) x_win]);% H0: residuals are uncorrelated;
            else
                [tdr.resDWp(t,j), tdr.resDW(t,j)] = dwtest(mdl.Residuals.Standardized, x_win);% H0: residuals are uncorrelated;
            end
            [~, tdr.resARCHp(t,j), tdr.resARCH(t,j)] = archtest(mdl.Residuals.Standardized, 'Lags', 1);% H0: no conditional heteroscedasticity (ARCH effect)
            [r_tmp, p_tmp] = corrcoef([mdl.Residuals.Raw x_win]);
            tdr.resXR{j}(t,1:nFactors) = r_tmp(1,2:end);
            tdr.resXRp{j}(t,1:nFactors) = p_tmp(1,2:end);
            tdr.VIF{j}(t,1:nFactors) = diag(inv(corrcoef(x_win)))';
            tdr.win(t,j) = nx_win;
            
            % compute fitted and one step-ahead forecasted values (for model that significantly differ from intercept model, otherwise use mean of independent value)
            if(round(tdr.Fp(t,j), 2) <= alpha)
                if(intercept)
                    tdr.fitted(t,j) = tdr.b{j}(t,:)*[1; x_curr'];
                    if(t ~= nObs)
                        tdr.forecasted(t,j) = tdr.b{j}(t,:)*[1; x_next'];
                    end
                else
                    tdr.fitted(t,j) = tdr.b{j}(t,:)*x_curr';
                end
            else
                tdr.fitted(t,j) = mean(y_win);
                if(t ~= nObs)
                    tdr.forecasted(t,j) = yImf(t,j);
                end
            end
            
%             if(t ~= nObs)
%                 if(intercept)
%                     tdr.forecasted(t,j) = tdr.b{j}(t,:)*[1; x_next'];
%                 else
%                     tdr.forecasted(t,j) = tdr.b{j}(t,:)*x_next';
%                 end
%             end
        end
        
        % set all non significant coefficients estimations to zero
        for k = 1:nFactors+intercept
            sign = (round(tdr.bp{j}(:,k), 2) <= alpha);
            tdr.b{j}(:,k) = sign.*tdr.b{j}(:,k);
        end
        
        % averaged coefficients and SEs
        if(bootstrap > 0)
            % bootstraped median
            for k = 1:nFactors+intercept
                bootstat = bootstrp(bootnum, @median, tdr.b{j}(:,k));
                tdr.mean.b{j}(k,1)  = mean(bootstat);
                
                bootstat = bootstrp(bootnum, @median, tdr.bse{j}(:,k));
                tdr.mean.bse{j}(k,1)  = mean(bootstat);
                
                tdr.mean.bts{j}(k,1) = tdr.mean.b{j}(k,1)/(tdr.mean.bse{j}(k,1));
                tdr.mean.bp{j}(k,1) = 1 - tcdf(abs(tdr.mean.bts{j}(k,1)), bootnum-1);
            end
        else
            % simple median
            for k = 1:nFactors+intercept
                tdr.mean.b{j}(k,1)  = median(tdr.b{j}(:,k));
                tdr.mean.bse{j}(k,1)  = median(tdr.bse{j}(:,k));
                tdr.mean.bts{j}(k,1) = tdr.mean.b{j}(k,1)/(tdr.mean.bse{j}(k,1));
                tdr.mean.bp{j}(k,1) = 1 - tcdf(abs(tdr.mean.bts{j}(k,1)), bootnum-1);
            end
        end
    end
    
    % averaged indicators
    tdr.mean.R2 = median(tdr.R2);
    tdr.mean.nMAE = median(tdr.nMAE);
    tdr.mean.nRMSE = median(tdr.nRMSE);
    tdr.mean.Fp = median(tdr.Fp);
    tdr.mean.FpN = round(100*sum((tdr.Fp <= alpha),1)/nObs);
    tdr.mean.resMean = median(tdr.resMean);
    tdr.mean.resADp = median(tdr.resADp);
    tdr.mean.resADpN = round(100*sum((tdr.resADp >= alpha),1)/nObs);
    tdr.mean.resTp = median(tdr.resTp);
    tdr.mean.resTpN = round(100*sum((tdr.resTp >= alpha),1)/nObs);
    tdr.mean.resJBp = median(tdr.resJBp);
    tdr.mean.resJBpN = round(100*sum((tdr.resJBp >= alpha),1)/nObs);
    
    for j = 1:nImfs
        tdr.mean.resXRp(j,:) = median(tdr.resXRp{j});
        tdr.mean.resXRpN(j,:) = round(100*sum((tdr.resXRp{j} >= alpha),1)/nObs);
        tdr.mean.VIF(j,:) = median(tdr.VIF{j});
        tdr.mean.VIFN(j,:) = round(100*sum((tdr.VIF{j} < vifMax),1)/nObs);
    end
    
    % summary for original data level
    tdr.original.fitted = sum(tdr.fitted, 2);
    [tdr.original.MAE, tdr.original.MAPE, tdr.original.MSE, tdr.original.RMSE, tdr.original.WAPE, tdr.original.nMAE, tdr.original.nRMSE] = mean_errors(y, tdr.original.fitted);
    
    tdr.original.forecasted = sum(tdr.forecasted, 2);
    [tdr.original.fMAE, tdr.original.fMAPE, tdr.original.fMSE, tdr.original.fRMSE, tdr.original.fWAPE, tdr.original.fnMAE, tdr.original.fnRMSE] = mean_errors(y(2:end,1), tdr.original.forecasted);
    
    % show results
    if(showResults>0)
         fontName = 'Helvetica';
         fontSize = 16;
         
         % plot the model parameters time-series (without intercept)
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
         
         %plot actual vs fitted data
         figure;
         subplot(1, 1, 1,  'FontName', fontName, 'FontSize', fontSize, 'Box', 'on');
         hold on;
         plot((1:nObs)', y, (1:nObs)', tdr.original.fitted);
         xlim([1 nObs]);
         legend('Actual data', 'Fitted data');
         hold off;
         
         % print table with averaged model parameters
         bTable = cell(nImfs,nFactors+1);
         bColTitles = cell(1,nFactors+1);
         bRowTitles = cell(1,nImfs);
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
             bRowTitles{j} = num2str(j);
         end
         print_latex_table(bTable, bColTitles, bRowTitles, '', 'Averaged model parameters.');
         
         % print table with average models quality measures and indicators of the residuals diagnostic
        dTable = cell(nImfs, 6);
        dRowTitles = cell(1,nImfs);
        for j = 1:nImfs
            dTable{j,1} = num2str(round(tdr.mean.nMAE(j)), '%i');
            dTable{j,2} = num2str(round(100*tdr.mean.R2(j)), '%i');
            dTable{j,3} = [num2str(tdr.mean.Fp(j), '%.2f'), ' (', num2str(tdr.mean.FpN(j), '%i'), ')'];
            dTable{j,4} = [num2str(tdr.mean.resADp(j), '%.2f'), ' (', num2str(tdr.mean.resADpN(j), '%i'), ')'];
            dTable{j,5} = [num2str(tdr.mean.resTp(j), '%.2f'), ' (', num2str(tdr.mean.resTpN(j), '%i'), ')'];
            dTable{j,6} = [num2str(tdr.mean.resJBp(j), '%.2f'), ' (', num2str(tdr.mean.resJBpN(j), '%i'), ')'];
            dRowTitles{j} = num2str(j);
        end
        print_latex_table(dTable, ...
            {'$\\overline{nMAE}$', '$\\overline{R^2}$', '$\\overline{p}_F (N_{\\overline{p}_F})$', '$\\overline{p}_{AD} (N_{\\overline{p}_{AD}})$', ...
            '$\\overline{p}_t (N_{\\overline{p}_t})$', '$\\overline{p}_{JB} (N_{\\overline{p}_{JB}})$'}, ...
            dRowTitles, '', 'Averaged the models quality measures and indicators of the residuals diagnostic.');
        
        % print table with multicollinearity and endogeneity analysis results
        mrTable = cell(nImfs, 2*nFactors);
        mrColTitles = cell(1,2*nFactors);
        mrRowTitles = cell(1,nImfs);
        for j = 1:nImfs
            for k =1:nFactors
                if(j == 1)
                    mrColTitles{1,k} = ['$\\overline{VIF}_{', factorsName{k}, '} (N_{\\overline{VIF}_{', factorsName{k}, '}})$'];
                    mrColTitles{1,k+nFactors} = ['$\\overline{p}_{', factorsName{k}, 'r} (N_{\\overline{p}_{', factorsName{k}, 'r}})$'];
                end
                mrTable{j,k} = [num2str(tdr.mean.VIF(j,k), '%.2f'), ' (', num2str(tdr.mean.VIFN(j,k), '%i'), ')'];
                mrTable{j,k+nFactors} = [num2str(tdr.mean.resXRp(j,k), '%.2f'), ' (', num2str(tdr.mean.resXRpN(j,k), '%i'), ')'];
            end
            mrRowTitles{j} = num2str(j);
        end
        print_latex_table(mrTable, mrColTitles, mrRowTitles, '', 'Analysis of multicollinearity and endogeneity.', '');
    end
     
    % switch on some warnings
    warning('on', 'stats:jbtest:PTooBig');
    warning('on', 'stats:jbtest:PTooSmall');
end
