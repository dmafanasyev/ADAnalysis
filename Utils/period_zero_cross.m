function [ period, indmin, indmax, indzer, numzer ] = period_zero_cross( data )
    %PERIOD_ZERO_CROSS Calculate period estimation using zero-crossing method.
    %
    %   Input:
    %       data - time-series to processing (row is time)
    %       
    %   Output:
    %       period - period estimation
    %       indmin - indexes of minimum values
    %       indmax - indexes of maximum values
    %       indzer - indexes of zero-crossings
    %       numzer - full number of zero-crossings (real + one per two extremums without zero-crossing)
    %
    %   Copyright (c) 2014-2015 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0 2014.04.12: initial version
    %       1.1 2015.08.04: added output params for indexes of min, max & zero
    %       1.2 2015.09.01: added output param 'numzer' and use it for period calculation
    %                       in the case when time-series or its some cycles has no zero-crossings,
    %                       calculate the "virtual" zero-crossings - one per two extremums without zero-crossing
    %                               
    
        
    [indmin, indmax, indzer] = extr(data');
    
    indext = sort([indmin indmax]);
    numzerext = 0;
    extNum = size(indext, 2);
    for j = 1:extNum-1
        if(nnz(indzer > indext(j) & indzer < indext(j+1)) == 0)
            numzerext = numzerext + 1;
        end
    end
    numzer = size(indzer,2) + numzerext;
    
    period = 4 * (size(data,1) / (size(indmin, 2) + size(indmax, 2) + numzer));
end

