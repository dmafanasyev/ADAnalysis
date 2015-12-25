function [modes, its]=ceemdan_par(x,Nstd,NR,MaxIter)
    % Compute Complete Ensemble Empirical Mode Decomposition with Adaptive Noise using Parallel Computing Toolbox.
    %
    % WARNING: for this code works it is necessary to include in the same
    % directoy the file emd.m developed by Rilling and Flandrin.
    % This file is available at http://perso.ens-lyon.fr/patrick.flandrin/emd.html.
    % We use the default stopping criterion.
    % We use the last modification: 3.2007
    %
    %----------------------------------------------------------------------
    %   INPUTs
    %   x: signal to decompose
    %   Nstd: noise standard deviation
    %   NR: number of realizations
    %   MaxIter: maximum number of sifting iterations allowed.
    %
    %  OUTPUTs
    %  modes: contain the obtained modes in a matrix with the rows being the modes
    %   its: contain the sifting iterations needed for each mode for each realization (one row for each realization)
    % -------------------------------------------------------------------------
    %  Syntax
    %
    %  modes=ceemdan(x,Nstd,NR,MaxIter)
    %  [modes its]=ceemdan(x,Nstd,NR,MaxIter)
    %
    %--------------------------------------------------------------------------
    % This algorithm was presented at ICASSP 2011, Prague, Czech Republic
    % Plese, if you use this code in your work, please cite the paper where the
    % algorithm was first presented.
    % If you use this code, please cite:
    %
    % M.E.TORRES, M.A. COLOMINAS, G. SCHLOTTHAUER, P. FLANDRIN,
    %  "A complete Ensemble Empirical Mode decomposition with adaptive noise,"
    %  IEEE Int. Conf. on Acoust., Speech and Signal Proc. ICASSP-11, pp. 4144-4147, Prague (CZ)
    %
    % Source version 1.0 (Torres ME, Colominas MA, Schlotthauer G, Flandrin P.) was modified by Dmitriy O. Afanasyev for perform parallel computations
    % on multicore computers and computer clusters. See parfor-loop (http://www.mathworks.com/help/distcomp/parfor.html).
    %
    % -------------------------------------------------------------------------
    % Date: June 06,2011
    % Authors:  Torres ME, Colominas MA, Schlotthauer G, Flandrin P.
    % This version was run on Matlab 7.10.0 (R2010a)
    % For problems with the code, please contact the authors:
    % To:  macolominas(AT)bioingenieria.edu.ar
    % CC:  metorres(AT)santafe-conicet.gov.ar
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    % Date: April 22,2014
    % Authors:  Dmitriy O. Afanasyev, dmafanasyev(AT)gmail.com
    % Version 1.1
    % Change Notes: Optimized for parallel computations. This version was run on Matlab R2013b
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    % Date: September 01,2015
    % Authors:  Dmitriy O. Afanasyev, dmafanasyev(AT)gmail.com
    % Version 1.2
    % Change Notes: Some minor optimization for preallocation of variables.
    % -------------------------------------------------------------------------
    
    x=x(:)';
    desvio_x=std(x);
    x=x/desvio_x;
    
    aux=zeros(size(x));
    iter=zeros(NR,round(log2(length(x))+5));
    
    white_noise = cell(1,NR);
    for i=1:NR
        white_noise{i}=randn(size(x));%creates the noise realizations
    end;
    
    modes_white_noise = cell(1,NR);
    for i=1:NR
        modes_white_noise{i}=emd(white_noise{i});%calculates the modes of white gaussian noise
%         modes_white_noise{i} = emdc([], white_noise{i});
    end;
    
    parfor i=1:NR %calculates the first mode
        temp=x+Nstd*white_noise{i};
        [temp, ~, it]=emd(temp,'MAXMODES',1,'MAXITERATIONS',MaxIter);
        %[temp, it] = emd_first_mode(temp, MaxIter);
        temp=temp(1,:);
        aux=aux+temp/NR;
        iter(i,1)=it;
    end;
    
    modes=aux; %saves the first mode
    k=2;
    aux=zeros(size(x));
    acum=sum(modes,1);
    
    while nnz(diff(sign(diff(x-acum))))>2 %calculates the rest of the modes
%         disp(num2str(nnz(diff(sign(diff(x-acum))))));
        parfor i=1:NR
            tamanio=size(modes_white_noise{i});
            if tamanio(1) >= k
                noise=modes_white_noise{i}(k-1,:);
                noise=noise/std(noise);
                noise=Nstd*noise;
                
%                 [temp, it] = emd_first_mode(x-acum+std(x-acum)*noise, MaxIter);
%                 if(isempty(it))
%                     it = 0;
%                     temp=x-acum;
%                 else
%                     temp=temp(1,:);
%                 end
                
                try
                    [temp, ~, it]=emd(x-acum+std(x-acum)*noise,'MAXMODES',1,'MAXITERATIONS',MaxIter);
                    temp=temp(1,:);
                catch
                    it=0;
                    temp=x-acum;
                end;
            else
%                 [temp, it] = emd_first_mode(x-acum, MaxIter);
%                 temp=temp(1,:);
                
                [temp, ~, it]=emd(x-acum,'MAXMODES',1,'MAXITERATIONS',MaxIter);
                temp=temp(1,:);
            end;
            aux=aux+temp/NR;
            iter(i, k) = it;
        end;
        modes=[modes;aux];
        aux=zeros(size(x));
        acum=zeros(size(x));
        acum=sum(modes,1);
        k=k+1;
    end;
    modes=[modes;(x-acum)];
    [a, ~]=size(modes);
    iter=iter(:,1:a);
    modes=modes*desvio_x;
    its=iter;
end

function [mode, it] = emd_first_mode(x, MaxIter)
    if(isinf(MaxIter) || isnan(MaxIter))
        [mode, it] = emdc([], x, [], 1);
    else
        [mode, it] = emdc_fix([], x, MaxIter, 1);
    end
end